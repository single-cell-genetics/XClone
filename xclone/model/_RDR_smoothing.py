"""base functions for XClone RDR module: smoothing.
"""
# Author: Rongting Huang
# Date: 2021/11/24
# update: 2022/07/26

import numpy as np
from .smoothing import WMA_smooth, KNN_smooth

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
    Xdata_norm = Xdata.copy()
    Xdata_norm.X = np.log(Xdata_norm.X/Xdata_norm.X.sum(1) * 10000 + 1)

    _is_ref = Xdata_norm.obs[cell_anno_key] == ref_celltype
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

    return Xdata