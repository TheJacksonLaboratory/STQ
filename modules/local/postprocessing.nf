
process DIMRED_CLUSTER {

    tag "$sample_id"
    label 'python_low_process'
    maxRetries 0
    errorStrategy  { task.attempt <= maxRetries  ? 'retry' : 'finish' }
    publishDir "${params.outdir}/${sample_id}/figures", pattern: 'figures/*/*.png', saveAs: { filename -> "${filename.split("/")[filename.split("/").length - 1]}" }, mode: 'copy', overwrite: true
    memory { 1.GB + (Float.valueOf(size) / 1000.0).round(2) * params.memory_scale_factor * 3.GB }
    
    input:
    tuple val(sample_id), path(grid_csv), path(grid_json), path(thumb), path(segmentation_csv), path(features_h5ad), val(expansion_factor), val(suffix)
    
    output:
    tuple val(sample_id), file("figures/*/*.png")

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
    
    plt.rcParams['figure.dpi'] = 150
    
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
    
    plt.rcParams["figure.figsize"] = figsize
    sc.pl.spatial(ad, img_key='lowres', color=cols1, spot_size=None, cmap='rainbow', ncols=3, show=False, save='/spatial_plot_morphometric.png');
    sc.pl.spatial(ad, img_key='lowres', color=cols2, spot_size=None, cmap='rainbow', ncols=3, show=False, save='/spatial_plot_classification.png');
    
    print(ad.obs)
    print(ad)
    
    # Dimensionality reduction
    sc.pp.highly_variable_genes(ad, flavor='seurat', n_top_genes=500)
    sc.pp.scale(ad)
    sc.pp.pca(ad, n_comps=30, zero_center=False, use_highly_variable=True)
    sc.pp.neighbors(ad, use_rep='X_pca')
    sc.tl.umap(ad)
    
    # Clustering
    res = 0.5
    sc.tl.leiden(ad, key_added='cluster', resolution=res)
    plt.rcParams["figure.figsize"] = figsize
    sc.pl.spatial(ad, img_key='lowres', color=[None, 'cluster'], spot_size=None, show=False, save='/cluster.png');
    
    # UMAP plots
    s = 20
    plt.rcParams["figure.figsize"] = (3,3)
    sc.pl.umap(ad, color=['cluster'], s=s, show=False, save='/umap_plot_cluster.png');
    sc.pl.umap(ad, color=['cluster'] + cols1, s=s, ncols=3, show=False, save='/umap_plot_morphometric.png');
    sc.pl.umap(ad, color=['cluster'] + cols2, s=s, ncols=3, show=False, save='/umap_plot_classification.png');
    """
}
