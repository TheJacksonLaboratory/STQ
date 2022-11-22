import os
import gzip
import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix
from scipy.io import mmwrite 
import scanpy as sc

def gz(fname):

    '''Compress file with gzip and remove source
    '''

    with open(fname) as f_in:
        with gzip.open(fname + '.gz', 'wt') as f_out:
            f_out.writelines(f_in)

    os.remove(fname)

    return

def read_sc_from_mtx(outsPath):
    
    sc_adata = sc.read_mtx(outsPath +'matrix.mtx.gz').T
    
    df_var = pd.read_csv(outsPath + 'features.tsv.gz', header=None, sep='\t', index_col=0)
    df_var.index.name = None
    sc_adata.var = df_var
    
    df_obs = pd.read_csv(outsPath + 'barcodes.tsv.gz', header=None).set_index(0)
    df_obs.index.name = None
    sc_adata.obs = df_obs
    
    print(sc_adata.shape)
    
    return sc_adata

def read_mtx_combine_and_write_mtx(adata1, adata2, saveDataDir=''):
    
    if not os.path.exists(saveDataDir):
        os.makedirs(saveDataDir)

    df = pd.concat([adata1.to_df(), adata2.to_df()], axis=1).fillna(0).astype(int)

    obs = pd.Series(df.index)
    obs.to_csv(saveDataDir + '/barcodes.tsv.gz', sep='\t', index=False, header=False)

    var = pd.concat([adata1.var, adata2.var]).loc[df.columns].reset_index()
    var.to_csv(saveDataDir + '/features.tsv.gz', sep='\t', index=False, header=False)

    fname = saveDataDir + '/matrix.mtx'
    mmwrite(fname, csr_matrix(df.values.T))
    gz(fname)
    
    return
