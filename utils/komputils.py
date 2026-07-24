import os
import json
import numpy as np
import pandas as pd
import tifffile
from tqdm import tqdm
from multiprocessing import get_context, Pool
from matplotlib.path import Path
import matplotlib.pyplot as plt
import scanpy as sc
import h5py
import anndata
from anndata._io.specs import _REGISTRY, IOSpec
@_REGISTRY.register_read(h5py.Dataset, IOSpec("null", "0.1.0"))
def read_null(elem, _reader):
    return None

qs = np.linspace(0.05, 0.95, 10, endpoint=True)

def contourMask(contour, xy):
    contour = np.asarray(contour).T
    xy = np.asarray(xy)
    path = Path(contour)
    return path.contains_points(xy, radius=1e-9)

def loadSamplesParallel(image_paths, save_path, F=1, CPU=16, target_mpp=0.5):
    def _q(args):
        s, fname = args
        dft = pd.read_csv(fname, index_col=0)
        dff = dft.iloc[:, 7:]
        xy = dft.iloc[:, [6, 5]].values
        with open(image_paths[s].replace('thumbnail.tiff', 'info.json'), 'r') as f:
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
        se = dff.loc[m].quantile(qs).stack()
        se.name = s
        return s, se

    if os.path.isfile(save_path):
        return pd.read_pickle(save_path)

    jobs = [(s, image_paths[s].replace('thumbnail.tiff', f'features/false-{F}-ctranspath_features.tsv.gz'))
            for s in sorted(image_paths)[:] if s in ad.obs.index]
    with get_context("fork").Pool(CPU) as p:
        ses = dict(tqdm(p.imap(_q, jobs, chunksize=10), total=len(jobs)))
    df_mfeat = pd.DataFrame(ses).T
    df_mfeat.to_pickle(save_path)

    return df_mfeat

def loadSubsetClusters(image_paths, ids, save_path):
    if os.path.isfile(save_path):
        return pd.read_pickle(save_path)
    dfs = []
    for s in tqdm(ids):
        p = image_paths[s]
        base = os.path.dirname(p)
        cname = f'{base}/clusters.csv.gz'
        se = pd.read_csv(cname, index_col=0).iloc[:, 0]
        fname = f'{base}/features/false-1-ctranspath_features.tsv.gz'
        dff = pd.read_csv(fname, index_col=0).iloc[:, 7:]
        assert (dff.index==se.index).all()
        dff.index = se.values
        dff = dff.groupby(level=0).quantile(qs).unstack().reorder_levels([1, 0], axis=1).T.sort_index().T
        dff.index = s + '.cls' + dff.index.astype(str)
        dfs.append(dff)
    dfs = pd.concat(dfs)
    dfs.to_pickle(save_path)
    return dfs

def reCluster(df_mfeat_sub):
    dfs = df_mfeat_sub.copy()
    ind = dfs.columns.copy()
    dfs.columns = dfs.columns.get_level_values(1) + '_' + dfs.columns.get_level_values(0).round(2).astype(str)
    adsub = sc.AnnData(dfs)
    dfs = None
    sc.pp.highly_variable_genes(adsub, flavor='seurat', n_top_genes=500)
    sc.pp.scale(adsub)
    sc.pp.pca(adsub, n_comps=30, zero_center=False, mask_var='highly_variable')
    sc.pp.neighbors(adsub, use_rep='X_pca')
    sc.tl.umap(adsub)
    adsub.obs['sample'] = adsub.obs.index.str.split('.cls').str[0]
    adsub.obs['slide'] = adsub.obs.index.str.split('.oid').str[0]
    df_sub = pd.DataFrame(adsub.obsm['X_umap'], index=adsub.obs.index).rename({0: 'x', 1: 'y'}, axis=1)
    return df_sub, adsub.obs

def subsetUMAP(df_umap, coords):
    x1, x2, y1, y2 = coords
    
    x = df_umap['x']
    y = df_umap['y']

    wh = (x>x1) & (x<x2) & (y>y1) & (y<y2)
    ids = wh[wh].index
    
    plt.scatter(x, y, c='grey', s=2)
    plt.scatter(x[wh], y[wh], c='k', s=1)
    ax = plt.gca()
    ax.set_aspect('equal')
    ax.axis('off')

    ax.text(x1, y2 + 0.05 * (y2 - y1), f'{len(ids)} samples')
    
    plt.plot([x1, x2, x2, x1, x1], [y1, y1, y2, y2, y1], c='k')
    plt.show()

    return ids
