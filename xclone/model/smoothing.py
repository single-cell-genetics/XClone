"""Base functions for XClone smoothing for RDR and BAF.
"""
# Author: Yuanhua Huang, Rongting Huang
# Date: 2022/03/26
# update: 2022/07/26

import numpy as np

def WMA_smooth(adata, layer=None, out_layer='smoothed', chrom_key='chr_arm', 
    gene_coordinate_key = 'start', method='pyramidinal', window_size=101, 
    verbose=True):
    
    chr_arms = np.unique(adata.var[chrom_key])
    
    full_weights = np.arange(1, window_size+1).reshape(1, -1)
    
    # Simple moving average
    if method == 'simple':
        full_weights = full_weights / full_weights
    
    # Exponential moving average
    elif method == 'exponential':
        full_weights = 0.8**(full_weights[::-1])
    
    if layer is None:
        X_smoothed = adata.X + 0
    else:
        X_smoothed = adata.layers[layer] + 0
    
    for _arm in chr_arms:
        _g_idx = np.where(adata.var[chrom_key].values == _arm)[0]
        _pos = adata.var[gene_coordinate_key][_g_idx].values
        
        _g_orders = np.argsort(_pos)
        _X_forward = X_smoothed[:, _g_idx[_g_orders]]
        _X_backward = _X_forward[:, ::-1]
        _X_forward_WMA = _X_forward * 0
        _X_backward_WMA = _X_backward * 0
        
        for ig in range(_X_forward.shape[1]):
            _idx_use = np.arange(ig + 1 - window_size, ig + 1)
            _idx_use = _idx_use[_idx_use >= 0]
            
            _weights = full_weights[:, -(len(_idx_use)):]
            _weights = _weights / np.sum(_weights)
            
            _X_forward_WMA[:, ig] = np.sum(
                _X_forward[:, _idx_use] * _weights, axis=1
            )
            
            _X_backward_WMA[:, ig] = np.sum(
                _X_backward[:, _idx_use] * _weights, axis=1
            )
        
        _X_WMA = (_X_forward_WMA + _X_backward_WMA[:, ::-1]) / 2
        X_smoothed[:, _g_idx[_g_orders]] = _X_WMA
        
        if verbose:
            print(_arm, len(_g_idx), round(_X_forward.mean(), 4), 
                  round(_X_WMA.mean(), 4))
                
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


## Smoothening by weighted moving average (similar to inferCNV)
## yuanhua
def WMA_chrom(adata, chrom_key='chr_arm', gene_coordinate_key = 'start', 
              method='pyramidinal', window_size=101, verbose=True):
    
    chr_arms = np.unique(adata.var[chrom_key])
    
    full_weights = np.arange(1, window_size+1).reshape(1, -1)
    
    # Simple moving average
    if method == 'simple':
        full_weights = full_weights / full_weights
    
    # Exponential moving average
    elif method == 'exponential':
        full_weights = 0.8**(full_weights[::-1])
    
    X_raw = adata.X + 0
    
    for _arm in chr_arms:
        _g_idx = np.where(adata.var[chrom_key].values == _arm)[0]
        _pos = adata.var[gene_coordinate_key][_g_idx].values
        
        _g_orders = np.argsort(_pos)
        _X_forward = X_raw[:, _g_idx[_g_orders]]
        _X_backward = _X_forward[:, ::-1]
        _X_forward_WMA = _X_forward * 0
        _X_backward_WMA = _X_backward * 0
        
        for ig in range(_X_forward.shape[1]):
            _idx_use = np.arange(ig + 1 - window_size, ig + 1)
            _idx_use = _idx_use[_idx_use >= 0]
            
            _weights = full_weights[:, -(len(_idx_use)):]
            _weights = _weights / np.sum(_weights)
            
            _X_forward_WMA[:, ig] = np.sum(
                _X_forward[:, _idx_use] * _weights, axis=1
            )
            
            _X_backward_WMA[:, ig] = np.sum(
                _X_backward[:, _idx_use] * _weights, axis=1
            )
        
        _X_WMA = (_X_forward_WMA + _X_backward_WMA[:, ::-1]) / 2
        X_raw[:, _g_idx[_g_orders]] = _X_WMA
        
        if verbose:
            print(_arm, len(_g_idx), round(_X_forward.mean(), 4), 
                  round(_X_WMA.mean(), 4))
                
    adata_rv = adata.copy()
    adata_rv.X = X_raw
        
    return adata_rv