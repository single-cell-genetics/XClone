"""base functions for XClone RDR module: smoothing.
"""
# Author: Rongting Huang
# Date: 2021/11/24
# update: 2022/07/26

import numpy as np
import scanpy as sc
import scipy.sparse as sp
from .smoothing import WMA_smooth, KNN_smooth

import gc

## Smoothing strategy

def RDR_smoothing_base(Xdata,
                       clip = True,
                       outlayer = "RDR_smooth",
                       cell_anno_key = "cell_type",
                       ref_celltype = "unclassified",
                       WMA_window_size = 50,
                       chrom_key = "chr_arm",
                       KNN_sm = True,
                       KNN_connect_use = "connectivities",
                       verbose = False
                       ):
    """
    For smoothing visualization and CNV states guide.
    """
    ## preprocess
    ## normalization [follow scanpy pipeline]
    # Xdata_norm.X = np.log(Xdata_norm.X/Xdata_norm.X.sum(1) * 10000 + 1) 
    # Notes: this only supported by python<3.7 if issparse
    # sc.pp.normalize_total from the Scanpy library supports both dense and sparse matrices. 
    
    Xdata_norm = Xdata.copy()

    #_is_ref = Xdata_norm.obs[cell_anno_key] == ref_celltype
    # modified for multiple ref_celltype
    if isinstance(ref_celltype, list):
        _is_ref = Xdata_norm.obs[cell_anno_key].isin(ref_celltype)
    else:
        _is_ref = Xdata_norm.obs[cell_anno_key] == ref_celltype

    # Normalize and log transform using scanpy functions
    sc.pp.normalize_total(Xdata_norm, target_sum=10000)
    sc.pp.log1p(Xdata_norm)

    # Check if Xdata_norm.X is a sparse matrix
    if sp.issparse(Xdata_norm.X):
        # # Perform operations in a way that supports sparse matrices
        # sums = np.array(Xdata_norm.X.sum(axis=1)).flatten()
        # Xdata_norm.X = Xdata_norm.X.multiply(1 / sums[:, None])
        # Xdata_norm.X = Xdata_norm.X.multiply(10000)
        # # Create a sparse matrix of ones with the same shape
        # ones_sparse = sp.csr_matrix(np.ones(Xdata_norm.X.shape))
        # # Add the sparse matrix of ones to the original sparse matrix
        # Xdata_norm.X = Xdata_norm.X + ones_sparse
        # Xdata_norm.X.data = np.log(Xdata_norm.X.data)
        
        Xdata_norm.var['ref_mean'] = np.array(Xdata_norm[_is_ref, :].X.mean(axis=0)).flatten()
        ref_data_dense = Xdata_norm[_is_ref, :].X.toarray()
        Xdata_norm.var['ref_var'] = ref_data_dense.var(axis=0)
    else:
        # Xdata_norm.X = np.log(Xdata_norm.X / Xdata_norm.X.sum(axis=1) * 10000 + 1)
        Xdata_norm.var['ref_mean'] = Xdata_norm[_is_ref, :].X.mean(0)
        Xdata_norm.var['ref_var'] = Xdata_norm[_is_ref, :].X.var(0)
      

    adata_tmp = Xdata_norm.copy()
    adata_tmp.X = adata_tmp.X - Xdata_norm.var['ref_mean'].values.reshape(1, -1)
    if clip:
        adata_tmp.X = np.clip(adata_tmp.X, -3, 3)
    adata_tmp.X = np.array(adata_tmp.X)

    ## WMA smoothing
    adata_tmp = WMA_smooth(adata_tmp, layer=None, out_layer = "WMA_smoothed", chrom_key=chrom_key, 
        gene_coordinate_key = 'start', method='pyramidinal', window_size = WMA_window_size, verbose=verbose)
    
    Xdata.layers["WMA_smoothed"] = adata_tmp.layers["WMA_smoothed"].copy()
    
    ## KNN smoothing
    if KNN_sm:
        adata_tmp = KNN_smooth(adata_tmp, run_KNN = False, KNN_Xlayer = None, KNN_connect_use = KNN_connect_use,
               layer = "WMA_smoothed", out_layer = outlayer)

        Xdata.layers[outlayer] = adata_tmp.layers[outlayer].copy()
    
    del Xdata_norm
    del adata_tmp
    gc.collect()

    return Xdata