"""
Multiscale, overlapping-neighborhood quantile aggregation of tile features.

For each tile, a neighborhood of R nearest tiles (by x,y) is built with a
KDTree. Neighborhoods overlap by construction (kNN), and we sweep a whole
spectrum of R values without re-querying the tree: kNN neighbor lists are
nested (the R=5 nearest tiles are a prefix of the R=20 nearest tiles), so a
single query at max(R) lets every smaller R just take a slice.

Per-neighborhood aggregation (feature quantiles across the R member tiles)
is done with a single vectorized np.quantile call over a gathered
(n_tiles, R, n_features) array, instead of pandas groupby.quantile over a
duplicated long-format frame. This avoids both the row-duplication memory
cost and the per-group/per-quantile overhead of groupby.

Output: one row per tile (the neighborhood's center), columns are a
MultiIndex (R, feature, q) for every requested scale, concatenated together.
"""

import os
import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
from tqdm import tqdm
import math
from matplotlib import pyplot as plt
import scanpy as sc
import warnings
import json
import tifffile
from matplotlib.path import Path

def contourMask(contour, xy):
    contour = np.asarray(contour).T
    xy = np.asarray(xy)
    path = Path(contour)
    return path.contains_points(xy, radius=1e-9)

def loadFilterFlat(fname, target_mpp=0.5):
    dft = pd.read_csv(fname, index_col=0)
    dff = dft.iloc[:, 7:]
    xy = dft.iloc[:, [6, 5]].values
    with open(os.path.dirname(os.path.dirname(fname)) + '/info.json', 'r') as f:
        info = json.loads(f.read())
    image = info['image']
    with tifffile.TiffFile(image) as f:
        fshape = f.pages[0].shape
    
    scale = float(info['mpp']) / target_mpp
    dims = fshape[1], fshape[0]
    roifile = info['roifile']
    with open(roifile.replace('/dev-komp/data-200/', '/roi/data_batch_0/'), 'r') as f:
        contour = json.loads(f.read())
    contour = np.array([contour['0']['points'], contour['1']['points']])
    contour = (scale * contour * np.array(dims)[:, None]).astype(int)
    cmin = contour.min(axis=1)
    contour[0] -= cmin[0]
    contour[1] -= cmin[1]
    m = contourMask(contour, xy)
    return m

def hex_grid(xy_min, xy_max, R, f=1.5):
    x_min, y_min = xy_min
    x_max, y_max = xy_max
    s = f * R                # center-to-center spacing
    dy = s * math.sqrt(3) / 2  # row-to-row vertical spacing

    points = []
    row = 0
    y = y_min
    while y <= y_max:
        x_offset = (s / 2) if (row % 2) else 0.0
        x = x_min + x_offset
        while x <= x_max:
            points.append((x, y))
            x += s
        y += dy
        row += 1
    return points

def multiscale_neighborhood_quantiles(
    fname,
    delta=224, # spacing in pixels between tile centers
    f=1.75, # hex grid spacing factor
    Rs=[5, 10, 20, 35, 45, 70, 100],
    qs=np.linspace(0.05, 0.95, 10, endpoint=True),
    xy_cols=(6, 5),
    feat_start_col=7,
    workers=-1,
):
    """
    Parameters
    ----------
    fname : str
        CSV path; first column is the index (tile id).
    Rs : sequence of int
        R spectrum; ranges for each step are
        unioned and deduped, e.g. {5,10,15,20,25,30}.
    qs : array-like
        Quantile levels in [0, 1].
    xy_cols : (int, int)
        Positional (iloc) columns for (x, y).
    feat_start_col : int
        Feature columns start at this positional index and run to the end.
    workers : int
        Parallel workers for the KDTree query (-1 = all cores).

    Returns
    -------
    pd.DataFrame indexed like the input, columns MultiIndex (R, feature, q).
    """
    dft = pd.read_csv(fname, index_col=0)
    dff = dft.iloc[:, feat_start_col:]
    print(f"Loaded {len(dft)} tiles, {dff.shape[1]} features")
    
    if True:
        mask = loadFilterFlat(fname)
        xy = dft.iloc[:, list(xy_cols)].loc[mask].to_numpy(dtype=float)
        dft = dft.loc[mask]
        dff = dff.loc[mask]
        print(f"Filtered to {len(dft)} tiles inside the ROI")
    else:
        xy = dft.iloc[:, list(xy_cols)].to_numpy(dtype=float)

    idx = dft.index
    cols = dff.columns
    qs = np.asarray(qs, dtype=float)
    n = xy.shape[0]
    xy_min, xy_max = xy.min(axis=0), xy.max(axis=0)
    # print(xy_min, xy_max)
    df = []
    df_coords = []

    tree = cKDTree(xy)
    for R in tqdm(Rs):
        if R==Rs[-1]:
            # Put one point in the image center, to ensure that the largest R is always represented
            grid = np.array([[0.5*(xy_min[0]+xy_max[0]), 0.5*(xy_min[1]+xy_max[1])]])
            R = 300
        else:
            # Make a hexagonal grid with center-to-center spacing of f*R using xy_min, xy_max
            grid = np.array(hex_grid(xy_min-0.5*R*delta, xy_max, R*delta, f=f))

        # Query the tree for radius R around each grid point
        nn_idx = tree.query_ball_point(grid, r=R*delta, workers=workers)
        counts = np.array([len(idxs) for idxs in nn_idx])
        wh = counts > 0
        grid = grid[wh]
        nn_idx = nn_idx[wh]
        counts = counts[wh]
        
        if len(grid) == 0:
            continue

        # Debug: plot the grid and the neighbors
        if False:
            for i, idxs in enumerate(nn_idx):
                if len(idxs) > 0:
                    plt.scatter(xy[idxs, 0], xy[idxs, 1], s=1)
            ax = plt.gca()
            ax.set_aspect('equal', adjustable='box')
            plt.scatter(grid.T[0], grid.T[1], color='red')
            plt.show()
            print("Grid points:", len(grid))

        # Aggregate the features of the neighbors for each grid point
        c = 0
        ind = []
        mlength = np.median([len(idxs) for idxs in nn_idx])
        # print(f"Median neighborhood size: {mlength:.1f} tiles")
        flat_index = np.concatenate(nn_idx)
        for i, idxs in enumerate(nn_idx):
            coords = xy[idxs]
            ind += [c] * len(idxs)
            c += 1
        xy_temp = dft.iloc[:, list(xy_cols)].iloc[flat_index]
        xy_temp.index = ind
        df_xy = xy_temp.groupby(level=0).mean()

        dff_temp = dff.iloc[flat_index]
        dff_temp.index = ind
        df_agg = dff_temp.groupby(level=0).quantile(qs).unstack().reorder_levels([1, 0], axis=1).sort_index(axis=1)
        df_agg.index =  f'R{R}-' + df_agg.index.astype(str)
        
        df.append(df_agg)
        df_coords.append(df_xy)

    if True:
        df0 = dff.groupby(level=0).quantile(qs).unstack().reorder_levels([1, 0], axis=1).sort_index(axis=1)
        df0 = df0.reset_index(drop=True)
        df0.index = f'R0-' + df0.index.astype(str)
        df0.index.name = None
        df.append(df0)

        df_coords0 = dft.iloc[:, list(xy_cols)].groupby(level=0).mean()
        df_coords0.index = df0.index
        df_coords.append(df_coords0)

    return pd.concat(df), pd.concat(df_coords)

def doCluster(out):
    dfs = out[0].copy()
    xy = out[1].copy()
    dfs.columns = dfs.columns.get_level_values(1) + '_' + dfs.columns.get_level_values(0).round(2).astype(str)
    adsub = sc.AnnData(dfs)
    adsub.obs[['x', 'y']] = xy.values
    adsub.obsm['spatial'] = adsub.obs[['x', 'y']].values
    adsub.obs['R'] = adsub.obs.index.str.split('-').str[0].str[1:].astype(float)
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", category=RuntimeWarning)
        sc.pp.highly_variable_genes(adsub, flavor='seurat', n_top_genes=500)
    sc.pp.scale(adsub)
    sc.pp.pca(adsub, n_comps=30, zero_center=False, mask_var='highly_variable')
    sc.pp.neighbors(adsub, use_rep='X_pca')
    sc.tl.umap(adsub)
    sc.tl.leiden(adsub, key_added='cluster', resolution=0.25)
    return adsub

def mapOne(dataPath, ibatch, sample, F):
    out = multiscale_neighborhood_quantiles(f"{dataPath}/results-batch-{ibatch}/{sample}/features/false-{F}-ctranspath_features.tsv.gz",
                                            Rs=[1, 2, 3, 4, 5, 6, 7, 9, 10, 20, 35, 45, 70, 100],
                                            qs=np.linspace(0.05, 0.95, 10), f=1.75)
    ad = doCluster(out)
    # sc.pl.spatial(ad, color=['x', 'y'], spot_size=500, cmap='jet')
    # sc.pl.umap(ad, color=['R', 'x', 'y'], cmap='jet')

    if True:
        thumb = tifffile.imread(f"{dataPath}/results-batch-{ibatch}/{sample}/thumbnail.jpeg")
        fig, axs = plt.subplots(1, 2, figsize=(10,5), width_ratios=[1,2])
        
        axs[0].imshow(thumb)
        axs[0].axis('off')
    
        x, y, R = ad.obsm['X_umap'].T[0], ad.obsm['X_umap'].T[1], ad.obs['R']
        seq = 5, 200
        wh = (R >= seq[1])
        axs[1].scatter(x[wh], y[wh], s=R[wh], edgecolor='none', c='orange')
        wh = (R >= seq[0]) & (R < seq[1])
        axs[1].scatter(x[wh], y[wh], s=R[wh]+1, edgecolor='none', c='k')
        wh = (R > 0) & (R < seq[0])
        axs[1].scatter(x[wh], y[wh], s=R[wh]+1, edgecolor='none', c='k')
        wh = (R == 0)
        axs[1].scatter(x[wh], y[wh], s=R[wh]+1, edgecolor='none', c='blue')
        axs[1].set_aspect('equal', adjustable='box')
        axs[1].axis('off')
        plt.show()
    return ad
