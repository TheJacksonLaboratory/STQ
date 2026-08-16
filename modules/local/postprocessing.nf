
process DIMRED_CLUSTER_MORPH {

    tag "$sample_id"
    label 'python_low_process'
    errorStrategy  { task.attempt <= 2  ? 'retry' : 'finish' }
    publishDir "${params.outdir}/${sample_id}/figures", pattern: 'figures/*/*.png', saveAs: { filename -> "${filename.split("/")[filename.split("/").length - 1]}" }, mode: 'copy', overwrite: params.overwrite_files_on_publish
    publishDir "${params.outdir}/${sample_id}", pattern: 'clusters.csv.gz', mode: 'copy', overwrite: params.overwrite_files_on_publish
    publishDir "${params.outdir}/${sample_id}", pattern: 'umap.csv.gz', mode: 'copy', overwrite: params.overwrite_files_on_publish
    memory { 1.GB + (Float.valueOf(size) / 1000.0).round(2) * params.memory_scale_factor * 3.GB }
    
    input:
    tuple val(sample_id), path(grid_csv), path(grid_json), path(thumb), path(segmentation_csv), path(features_h5ad), val(expansion_factor), val(suffix)
    
    output:
    tuple val(sample_id), file("figures/*/*.png"), emit: figures
    tuple val(sample_id), file("clusters.csv.gz"), emit: clusters
    tuple val(sample_id), file("umap.csv.gz"), emit: umap

    script:
    """
    #!/usr/bin/env python
    
    import os
    os.environ["NUMBA_CACHE_DIR"] = "./tmp"

    import sys
    import json
    import numpy as np
    import pandas as pd
    import scanpy as sc
    import matplotlib.pyplot as plt
    
    plt.rcParams['figure.dpi'] = ${params.plot_dpi}
    
    def loadImFeatures(dpath):
        df_temp = pd.read_csv(dpath, index_col=[0,1], sep=',').xs(1, level='in_tissue')
        df_temp.insert(0, 'original_barcode', df_temp.index.values)
        ad = sc.AnnData(X=df_temp.loc[:, df_temp.columns.str.contains('feat')],
                        obs=df_temp.loc[:, ~df_temp.columns.str.contains('feat')])
        return ad
    
    def loadAdImage():
        thumbnail = plt.imread("${thumb}")
        with open("${grid_json}", 'r') as f:        
            d = json.load(f)
        grid = pd.read_csv("${grid_csv}", index_col=0, header=None)
        image = {'library_id': {'images': {'lowres': thumbnail},
                                    'metadata': {'chemistry_description': None, 'software_version': None},
                                    'scalefactors': {'tissue_lowres_scalef': thumbnail.shape[0]/d['y'],
                                                     'spot_diameter_fullres': d['spot_diameter_fullres']}}}, grid.index.values, grid[[5, 4]].values
        return image

    # Load data    
    ad = sc.read_h5ad("${features_h5ad}")
    df_temp = pd.read_csv("${segmentation_csv}", index_col=0, header=0).reindex(ad.obs.index)
    df_temp.index.name = 'id'
    ad.obs = pd.concat([ad.obs, df_temp], axis=1)
    
    # Load image
    image = loadAdImage()
    ad.uns['spatial'] = image[0]

    ad.obsm['spatial'] = pd.DataFrame(index=image[1], data=image[2]).reindex(ad.obs['original_barcode']).values
    
    # Morphometrics spatial plots
    cols1 = [None] + df_temp.columns[df_temp.columns.isin(['average_perimeter_length', 'average_area', 'average_eccentricity',
                                                  'average_orientation', 'average_cell_type_prob'])].values.tolist()
    cols2 = [None] + df_temp.columns[df_temp.columns.str.contains('count')].values.tolist()
    
    c, r = np.ptp(ad.obs['array_row']), np.ptp(ad.obs['array_col'])
    f = 5
    if r > c:
        figsize = f, f * c/r
    else:
        figsize = f * r/c, f
    
    if not os.path.exists('figures/show/'):
        os.makedirs('figures/show/')
    if not os.path.exists('figures/umap/'):
        os.makedirs('figures/umap/')

    spot_size = ad.uns['spatial']['library_id']['scalefactors']['spot_diameter_fullres']
    print(spot_size)
    
    spot_size *= ${params.grid_spot_horizontal_spacing} / ${params.grid_spot_diamter}
    print(spot_size)
    
    plt.rcParams["figure.figsize"] = figsize
    sc.pl.spatial(ad, img_key='lowres', color=cols1, spot_size=spot_size, cmap='rainbow', ncols=3, show=False, save='/spatial_plot_morphometric.png');
    sc.pl.spatial(ad, img_key='lowres', color=cols2, spot_size=spot_size, cmap='rainbow', ncols=3, show=False, save='/spatial_plot_classification.png');
    
    print(ad.obs)
    print(ad)
    
    # Dimensionality reduction
    sc.pp.highly_variable_genes(ad, flavor='seurat', n_top_genes=500)
    sc.pp.scale(ad)
    sc.pp.pca(ad, n_comps=min(30, ad.shape[0]-1), zero_center=False, use_highly_variable=True)
    sc.pp.neighbors(ad, use_rep='X_pca')
    sc.tl.umap(ad)
    
    # Clustering
    res = 0.5
    sc.tl.leiden(ad, key_added='cluster', resolution=res)
    plt.rcParams["figure.figsize"] = figsize
    sc.pl.spatial(ad, img_key='lowres', color=[None, 'cluster'], spot_size=spot_size, show=False, save='/cluster.png');
    
    # UMAP plots
    plt.rcParams["figure.figsize"] = (3,3)
    sc.pl.umap(ad, color=['cluster'], s=None, show=False, save='/umap_plot_cluster.png');
    sc.pl.umap(ad, color=['cluster'] + cols1, s=None, ncols=3, show=False, save='/umap_plot_morphometric.png');
    sc.pl.umap(ad, color=['cluster'] + cols2, s=None, ncols=3, show=False, save='/umap_plot_classification.png');

    # Save cluster and UMAP coordinates
    ad.obs['cluster'].to_csv('clusters.csv.gz')
    pd.DataFrame(ad.obsm['X_umap'], index=ad.obs.index, columns=['UMAP1', 'UMAP2']).to_csv('umap.csv.gz')
    """
}


process MAKE_SAMPLER_VECTOR {

    tag "$sample_id"
    label 'python_low_process'
    errorStrategy  { task.attempt <= 2  ? 'retry' : 'finish' }
    publishDir "${params.outdir}/${sample_id}", pattern: 'sampler-vector-*.pkl', mode: 'copy', overwrite: params.overwrite_files_on_publish
    memory { 1.GB + (Float.valueOf(size) / 1000.0).round(2) * params.memory_scale_factor * 3.GB }
    
    input:
    tuple val(sample_id), val(meta), path(srgrid), path(features_csv), val(expansion_factor), val(suffix)
    
    output:
    tuple val(sample_id), file("sampler-vector-*.pkl")

    script:
    """
    #!/usr/bin/env python
    
    import os
    os.environ["NUMBA_CACHE_DIR"] = "./tmp"

    import json
    import numpy as np
    import pandas as pd
    import tifffile
    from matplotlib.path import Path

    def contourMask(contour, xy):
        contour = np.asarray(contour).T
        xy = np.asarray(xy)
        path = Path(contour)
        return path.contains_points(xy, radius=1e-9)

    qs = np.linspace(0.05, 0.95, 10, endpoint=True)

    dft = pd.read_csv("${features_csv}", index_col=0)
    dff = dft.iloc[:, 7:]
    xy = dft.iloc[:, [6, 5]].values
    with tifffile.TiffFile('${meta.image}') as f:
        fshape = f.pages[0].shape

    scale = float(${meta.mpp}) / float(${params.target_mpp})
    dims = fshape[1], fshape[0]
    with open('${meta.roifile}', 'r') as f:
        contour = json.loads(f.read())
    contour = np.array([contour['0']['points'], contour['1']['points']])
    contour = (scale * contour * np.array(dims)[:, None]).astype(int)
    cmin = contour.min(axis=1)
    contour[0] -= cmin[0]
    contour[1] -= cmin[1]
    m = contourMask(contour, xy)
    se = dff.loc[m].quantile(qs).stack()

    F = "${expansion_factor}"
    model = "${suffix}"
    se.to_pickle(f'sampler-vector-{F}-{model}.pkl')
    """
}


process MAKE_CLUSTER_VECTORS {

    tag "$sample_id"
    label 'python_low_process'
    errorStrategy  { task.attempt <= 2  ? 'retry' : 'finish' }
    publishDir "${params.outdir}/${sample_id}", pattern: 'sampler-clusters-*.pkl', mode: 'copy', overwrite: params.overwrite_files_on_publish
    memory { 1.GB + (Float.valueOf(size) / 1000.0).round(2) * params.memory_scale_factor * 3.GB }
    
    input:
    tuple val(sample_id), path(thumb), path(segmentation_csv), path(features_h5ad), val(expansion_factor), val(suffix), path(clusters)
    
    output:
    tuple val(sample_id), file("sampler-clusters-*.pkl")

    script:
    """
    #!/usr/bin/env python
    
    import os
    os.environ["NUMBA_CACHE_DIR"] = "./tmp"

    import numpy as np
    import pandas as pd
    import scanpy as sc
    qs = np.linspace(0.05, 0.95, 10, endpoint=True)

    dff = sc.read_h5ad("${features_h5ad}").to_df()
    se = pd.read_csv("${clusters}", index_col=0).iloc[:, 0]
    assert (dff.index==se.index).all()
    dff.index = se.values
    dff = dff.groupby(level=0).quantile(qs).unstack().reorder_levels([1, 0], axis=1).T.sort_index().T
    dff.index = "${sample_id}" + '.cls' + dff.index.astype(str)

    se_counts = pd.read_csv("${segmentation_csv}", index_col=0)['total_count']
    se_counts = se_counts.reindex(se.index).fillna(0).astype(int)
    se_counts.index = se.values
    se_counts = se_counts.groupby(level=0).mean()
    se_counts.index = "${sample_id}" + '.cls' + se_counts.index.astype(str)
    dff = dff.loc[se_counts > float(${params.sampler_clusters_min_cell_count})]

    F = "${expansion_factor}"
    model = "${suffix}"
    dff.to_pickle(f'sampler-clusters-{F}-{model}.pkl')
    """
}


process DIMRED_CLUSTER {

    tag "$sample_id"
    label 'python_low_process'
    errorStrategy  { task.attempt <= 2  ? 'retry' : 'finish' }
    publishDir "${params.outdir}/${sample_id}/figures", pattern: 'figures/*/*.png', saveAs: { filename -> "${filename.split("/")[filename.split("/").length - 1]}" }, mode: 'copy', overwrite: params.overwrite_files_on_publish
    publishDir "${params.outdir}/${sample_id}", pattern: 'clusters.csv.gz', mode: 'copy', overwrite: params.overwrite_files_on_publish
    publishDir "${params.outdir}/${sample_id}", pattern: 'umap.csv.gz', mode: 'copy', overwrite: params.overwrite_files_on_publish
    memory { 1.GB + (Float.valueOf(size) / 1000.0).round(2) * params.memory_scale_factor * 3.GB }
    
    input:
    tuple val(sample_id), path(grid_csv), path(grid_json), path(thumb), path(features_h5ad), val(expansion_factor), val(suffix)
    
    output:
    tuple val(sample_id), file("figures/*/*.png")
    tuple val(sample_id), file("clusters.csv.gz")
    tuple val(sample_id), file("umap.csv.gz")

    script:
    """
    #!/usr/bin/env python
    
    import os
    os.environ["NUMBA_CACHE_DIR"] = "./tmp"

    import sys
    import json
    import numpy as np
    import pandas as pd
    import scanpy as sc
    import matplotlib.pyplot as plt
    
    plt.rcParams['figure.dpi'] = ${params.plot_dpi}
    
    def loadImFeatures(dpath):
        df_temp = pd.read_csv(dpath, index_col=[0,1], sep=',').xs(1, level='in_tissue')
        df_temp.insert(0, 'original_barcode', df_temp.index.values)
        ad = sc.AnnData(X=df_temp.loc[:, df_temp.columns.str.contains('feat')],
                        obs=df_temp.loc[:, ~df_temp.columns.str.contains('feat')])
        return ad
    
    def loadAdImage():
        thumbnail = plt.imread("${thumb}")
        with open("${grid_json}", 'r') as f:        
            d = json.load(f)
        grid = pd.read_csv("${grid_csv}", index_col=0, header=None)
        image = {'library_id': {'images': {'lowres': thumbnail},
                                    'metadata': {'chemistry_description': None, 'software_version': None},
                                    'scalefactors': {'tissue_lowres_scalef': thumbnail.shape[0]/d['y'],
                                                     'spot_diameter_fullres': d['spot_diameter_fullres']}}}, grid.index.values, grid[[5, 4]].values
        return image

    # Load data    
    ad = sc.read_h5ad("${features_h5ad}")
    
    # Load image
    image = loadAdImage()
    ad.uns['spatial'] = image[0]

    ad.obsm['spatial'] = pd.DataFrame(index=image[1], data=image[2]).reindex(ad.obs['original_barcode']).values
    
    c, r = np.ptp(ad.obs['array_row']), np.ptp(ad.obs['array_col'])
    f = 5
    if r > c:
        figsize = f, f * c/r
    else:
        figsize = f * r/c, f
    
    if not os.path.exists('figures/show/'):
        os.makedirs('figures/show/')
    if not os.path.exists('figures/umap/'):
        os.makedirs('figures/umap/')
    
    plt.rcParams["figure.figsize"] = figsize
    
    print(ad.obs)
    print(ad)
    
    spot_size = ad.uns['spatial']['library_id']['scalefactors']['spot_diameter_fullres']
    print(spot_size)
    
    spot_size *= ${params.grid_spot_horizontal_spacing} / ${params.grid_spot_diamter}
    print(spot_size)

    # Dimensionality reduction
    sc.pp.highly_variable_genes(ad, flavor='seurat', n_top_genes=500)
    sc.pp.scale(ad)
    sc.pp.pca(ad, n_comps=min(30, ad.shape[0]-1), zero_center=False, use_highly_variable=True)
    sc.pp.neighbors(ad, use_rep='X_pca')
    sc.tl.umap(ad)
    
    # Clustering
    res = 0.5
    sc.tl.leiden(ad, key_added='cluster', resolution=res)
    plt.rcParams["figure.figsize"] = figsize
    sc.pl.spatial(ad, img_key='lowres', color=[None, 'cluster'], spot_size=spot_size, show=False, save='/cluster.png');
    
    # UMAP plots
    plt.rcParams["figure.figsize"] = (3,3)
    sc.pl.umap(ad, color=['cluster'], s=None, show=False, save='/umap_plot_cluster.png');

    # Save cluster and UMAP coordinates
    ad.obs['cluster'].to_csv('clusters.csv.gz')
    pd.DataFrame(ad.obsm['X_umap'], index=ad.obs.index, columns=['UMAP1', 'UMAP2']).to_csv('umap.csv.gz')
    """
}

