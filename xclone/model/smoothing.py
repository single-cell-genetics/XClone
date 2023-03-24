"""Base functions for XClone smoothing for RDR and BAF.
"""
# Author: Yuanhua Huang, Rongting Huang
# Date: 2022/03/26
# update: 2023/03/21

import numpy as np
from scipy.sparse import csc_matrix

def make_WMA_connectivity(adata, chrom_key='chr_arm', 
                          gene_coordinate_key = 'start', 
                          method='pyramidinal', window_size=101, 
                          verbose=False):
    """
    Generate a gene-by-gene connectivity matrix for 
    weighted moving average smoothing
    """
    chr_arms = np.unique(adata.var[chrom_key])
    
    full_weights = np.arange(1, window_size+1).reshape(1, -1)
    
    # Simple moving average
    if method == 'simple':
        full_weights = full_weights / full_weights
    
    # Exponential moving average
    elif method == 'exponential':
        full_weights = 0.8**(full_weights[::-1])
    
    col_list = np.array([])
    row_list = np.array([])
    dat_list = np.array([])
    for _arm in chr_arms:
        _g_idx = np.where(adata.var[chrom_key].values == _arm)[0]
        _pos = adata.var[gene_coordinate_key][_g_idx].values
        
        _g_orders = np.argsort(_pos)
        gene_idx = _g_idx[_g_orders]
        
        for ig in range(len(gene_idx)):
            _idx_use = np.arange(ig + 1 - window_size, ig + 1)
            _idx_use = _idx_use[_idx_use >= 0]
            
            _weights = full_weights[:, -(len(_idx_use)):]
            _weights = _weights / np.sum(_weights)
            
            # forward
            col_list = np.append(col_list, np.repeat(gene_idx[ig], len(_idx_use)))
            row_list = np.append(row_list, gene_idx[_idx_use])
            dat_list = np.append(dat_list, _weights / 2)
            
            # backward
            col_list = np.append(col_list, np.repeat(gene_idx[::-1][ig], len(_idx_use)))
            row_list = np.append(row_list, gene_idx[::-1][_idx_use])
            dat_list = np.append(dat_list, _weights / 2)
            
        if verbose:
            print(_arm, len(gene_idx))
        
    WMA_connectivity = csc_matrix((dat_list, (row_list, col_list)), 
                                  shape=(adata.shape[1], adata.shape[1]))
    
    return WMA_connectivity


def WMA_smooth(adata, layer=None, out_layer='smoothed', chrom_key='chr_arm', 
               gene_coordinate_key = 'start', method='pyramidinal', 
               window_size=101, connect_key='WMA_connect', run_WMA=None, 
               verbose=True):
    """Perform WMA smoothing
    """
    if connect_key not in adata.varp or run_WMA == True:
        _WMA_mat = make_WMA_connectivity(
            adata, chrom_key=chrom_key, 
            gene_coordinate_key=gene_coordinate_key, 
            method=method, window_size=window_size, verbose=verbose
        )
        adata.varp[connect_key] = _WMA_mat
        print('make WMA connectivities matrix, saved in varp[%s].' %(connect_key))
    else:
        print('%s exists for direct use.' %(connect_key))
        
    if layer is None:
        X_smoothed = adata.X + 0
    else:
        X_smoothed = adata.layers[layer] + 0
        
    X_smoothed = X_smoothed @ adata.varp[connect_key]
    adata_out = adata.copy()
    adata_out.layers[out_layer] = X_smoothed
    
    return adata_out


from .base_utils import normalize
from ._BAF_process import extra_preprocess_BAF

def KNN_smooth(Xdata, run_KNN = False, KNN_Xlayer = None, KNN_connect_use = "connectivities",
               layer = "BAF", out_layer='smoothed'):
    """
    Smoothing info between cells across K-nearest neibourgh

    run_KNN : True, need specify KNN_Xlayer for running KNN.
    """
    # if 'connectivities' not in Xdata.obsp:
        # print('[XClone]Warning: No KNN connectivities available, pls generate')
    if run_KNN == True:
        update_Xdata = extra_preprocess_BAF(Xdata, KNN_Xlayer, run_KNN = True, copy = True)
        print("[XClone] generate KNN connectivities from KNN_Xlayer: ", KNN_Xlayer)
    else:
        update_Xdata = Xdata.copy()

    # normalize connectivities 
    # connectivities = normalize(update_Xdata.obsp[KNN_connect_use].A).copy()
    connectivities = update_Xdata.obsp[KNN_connect_use].copy()
    update_Xdata.layers[out_layer] = connectivities @ update_Xdata.layers[layer]
    return update_Xdata
