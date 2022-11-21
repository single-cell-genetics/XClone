"""Base functions for XClone data preprocessing
"""
# Author: Rongting Huang
# Date: 2021/05/04
# update: 2021/07/22 for xclonedata format

# import os
import sys
import warnings

from anndata import AnnData

import logging
# logging.debug('debug message')
# logging.info('info message')
# logging.warn('warn message')
# logging.error('error message')
# logging.critical('critical message') 


import numpy as np
import scipy as sp
import pandas as pd
from scipy.sparse import issparse

## Part I: Common used Function

def data_check(Xdata, Xlayer):
    """
    check negative number or Nan vaule in raw data
    for RDR module and BAF module.
    """
    if sp.sparse.issparse(Xdata.layers[Xlayer]):
        X_mtx = Xdata.layers[Xlayer].copy().A
    else:
        X_mtx = Xdata.layers[Xlayer].copy()
    
    # detect nan Value
    nan_count = np.isnan(X_mtx).sum()
    if nan_count > 0:
        print("[XClone data check in layer:]", Xlayer)
        raise ValueError("[XClone data] NaN value in counts! Pls check!")
    
    # detect negative Value
    if np.any(X_mtx < 0):
        print("[XClone data check in layer:]", Xlayer)
        raise ValueError("[XClone data] Negative value in counts! Pls check!")
    
    return None

def valid_cell(Xdata, cell_anno_key = "cell_type", filter_na_anno = True, verbose = True):
    """
    Function:
    Filter the cells without annotation.
    """
    if filter_na_anno:
        valid_cells = Xdata.obs[cell_anno_key] == Xdata.obs[cell_anno_key]
        update_Xdata = Xdata[valid_cells,:].copy()

        if verbose:
            cells_num = Xdata.shape[0]
            print("Keep valid cells: Filter out %d cells / %d total cells, remain %d valid cells with annotation" 
               %(cells_num - valid_cells.sum(), cells_num, valid_cells.sum()))
    else:
        update_Xdata = Xdata.copy()
        update_Xdata.obs[cell_anno_key].fillna("unannotated", inplace = True)
    return update_Xdata

## RDR used Function
def gene_filter(Xdata, cell_detection_rate = 0.05, verbose = True):
    """
    Function:
    ---------
    # 1) detect nan value in Xdata;
    2) do gene filtering base on cell_detection_rate; 
    default: keep genes have expression across at least 5% cells;

    Parameters:
    -----------
    Xdata:
    cell_detection_rate: float, default: 0.05
    verbose: bool

    Return:
    ------
    update_Xdata:
    Flag_: `array like` bool
    
    Examples
    --------
    update_Xdata, Flag_ = gene_filter(Xdata)
    """
    if sp.sparse.issparse(Xdata.X):
        X_mtx = Xdata.copy().X.A
    else:
        X_mtx = Xdata.copy().X
    
    # filter genes based on detection rate
    No_detection_rate_ = (X_mtx == 0).astype(int).sum(axis=0) / X_mtx.shape[0]

    Flag_ = (1- No_detection_rate_) > cell_detection_rate
    
    update_Xdata = Xdata[:, Flag_].copy()

    if verbose:
        gene_num = X_mtx.shape[1]
        print("[XClone-RDR preprocessing] Filter out %d genes / %d total genes, remain %d genes" 
               %(gene_num - Flag_.sum(), gene_num, Flag_.sum()))
    # return update_Xdata, Flag_
    return update_Xdata

## BAF used Function


###===================================
## todo:
# BAF data DP coverage cutoff 20
# actually just can be used in BAF merge_Xdata
####==================================
def DP_coverage_check(Xdata, Xlayer, threshold = 20):
    """
    can be applied on BAF_adata or merge_Xdata.
    """
    if issparse(Xdata.layers[Xlayer]):
        flag_ = Xdata.layers[Xlayer].A.sum(axis=0) > threshold
    else:
        flag_ = Xdata.layers[Xlayer].sum(axis=0) > threshold

    update_Xdata = Xdata.copy()[:, flag_]
    return update_Xdata
    

## Part I: Simple preprocessing functions
### filter data


def filter_nulldata(X,X_name="X", axis=1, output_format="np.arr"):
    """
    function: filter the all zero/nan rows or columns;
    X:input data, can be np.array or np.matrix or sp.sparse_matrix.
    default axis: 1,filter the rows without information; 0, filter the cols
    here, default filter the rows (features); if for cols(cellls), axis=0
    default data name: X_name="X", can use RDR,AD,DP, etc.
    return filtered data (and index): output format can be "np.arr","np.mtx",
    "sp.sparse_mtx"(default:csr sparse matirx)
    """
    if type(X) == np.matrix or issparse(X):
        X=X.A
    # if isinstance(X, AnnData):
    #     X = X.X
    #     ## If there is only one layer in the AnnData
    print(X_name,X.shape)
    X_mask = np.all(np.isnan(X) | np.equal(X, 0), axis=axis)
    X_filter_idx = ~X_mask
    if axis==1:
        X_filter = X[X_filter_idx]
    elif axis==0:
        X_filter = X[:,X_filter_idx]
    print(X_name," filter 0 and nan:", X_filter.shape)
    if output_format=="np.mtx":
        X_filter = np.mat(X_filter)
    if output_format=="sp.sparse_mtx":
        X_filter = sp.sparse.csr_matrix(X_filter)
    if output_format==None:
        return X_filter_idx
    return X_filter_idx, X_filter

def filter_2nulldata(X,Y,X_name="X", Y_name="Y",axis=1, output_format="np.arr"):
    """
    function: filter the all zero/nan rows or columns for 2 dataset;
    X,Y:input data, can be np.array or np.matrix or sp.sparse_matrix.
    default axis: 1, filter the rows without information; 0, filter the cols
    default data name: X_name="X",Y_name="Y", can use AD,DP, etc.
    return filtered data:output format can be "np.arr","np.mtx",
    "sp.sparse_mtx"(default:csr sparse matirx)
    -----
    Notes:deprecated in xclonedata format
    """
    if type(X) == np.matrix or issparse(X):
        X=X.A
    if type(Y) == np.matrix or issparse(Y):
        Y=Y.A
    print(X_name, X.shape)
    print(Y_name, Y.shape)
    # X_mask = np.all(np.isnan(X) | np.equal(X, 0), axis=axis)
    Y_mask = np.all(np.isnan(Y) | np.equal(Y, 0), axis=axis) # Y is DP in XClone
    # XY_mask = X_mask|Y_mask
    XY_mask = Y_mask
    XY_filter_idx = ~XY_mask
    if axis==1:
        X_filter = X[XY_filter_idx]
        Y_filter = Y[XY_filter_idx]
    elif axis==0:
        X_filter = X[:, XY_filter_idx]
        Y_filter = Y[:, XY_filter_idx]
    print(X_name," filter 0 and nan:", X_filter.shape)
    print(Y_name," filter 0 and nan:", Y_filter.shape)
    ## output format
    if output_format=="np.mtx":
        X_filter = np.mat(X_filter)
        Y_filter = np.mat(Y_filter)
    if output_format=="sp.sparse_mtx":
        X_filter = sp.sparse.csr_matrix(X_filter)
        Y_filter = sp.sparse.csr_matrix(Y_filter)
    print(type(X_filter), type(Y_filter))
    warnings.warn('[XClone-preprocessing]filter_2nulldata-Maybe DeprecationWarning', DeprecationWarning)
    return XY_filter_idx, X_filter, Y_filter

def tidy_Xdata(
    Xdata: AnnData,
    drop_obs: bool = False,
    drop_features: bool = True,
    inplace: bool = True):
    """
    Function:
    filter the all zero/nan rows or columns;
    -----
    params:
    data: xclonedata formatm AnnData
          The (annotated) data matrix of shape `n_obs` × `n_vars`.
          Rows correspond to cells/samples and columns to genes/blocks/features.
    -----
    Examples:
    

    """
    data_mode = Xdata.uns["data_mode"]
    if data_mode == "BAF":
        # AD = Xdata.layers["AD"]
        DP = Xdata.layers["DP"]
        if drop_obs:
            logging.info('info message: [XClone-tidy_Xdata-BAF] check obs')
            # filter_idx, X_filter, Y_filter = filter_2nulldata(AD,DP,X_name="AD", Y_name="DP",axis=1, output_format="np.arr")
            filter_idx0 = filter_nulldata(DP, X_name="BAF", axis=1, output_format=None)
        else:
            filter_idx0 = np.ones(DP.shape[0], dtype= bool)
        if drop_features:
            logging.info('info message: [XClone-tidy_Xdata-BAF] check features')
            filter_idx1 = filter_nulldata(DP, X_name="BAF", axis=0, output_format=None)
        else:
            filter_idx1 = np.ones(DP.shape[1], dtype= bool)
        Xdata.uns["log"]
    elif data_mode == "RDR":
        RDR = Xdata.X
        if drop_obs:
            logging.info('info message: [XClone-tidy_Xdata-RDR] check obs')
            filter_idx0 = filter_nulldata(RDR, X_name="RDR", axis=1, output_format=None)
        else:
            filter_idx0 = np.ones(RDR.shape[0], dtype= bool)
        if drop_features:
            logging.info('info message: [XClone-tidy_Xdata-RDR] check features')
            filter_idx1 = filter_nulldata(RDR, X_name="RDR", axis=0, output_format=None)
        else:
            filter_idx1 = np.ones(RDR.shape[1], dtype= bool)
    else:
        logging.error('error message: [XClone-tidy_Xdata] Not xclonedata format')
        sys.exit(1)
    log_msg = "remove %d obs and %d features,"  % ((~filter_idx0).sum(), (~filter_idx0).sum())
    log_msg += "tidydata (%d, %d)" % (filter_idx0.sum(), filter_idx1.sum())
    Xdata.uns["log"]["tidy_Xdata"] = log_msg
    return Xdata[filter_idx0, filter_idx1]

# import scanpy as sc
import numpy as np
def filter_pre(
    Xdata: AnnData,
    inplace: bool = False,
    **kwargs):
    """
    Function: add some annotation for obs/features
    ----
    Examples:
    import xclone
    BAF_adata = xclone.pp.load_TNBC1_BAF()
    RDR_adata = xclone.pp.load_TNBC1_RDR()
    BAF_adata = xclone.pp.tidy_Xdata(BAF_adata)
    RDR_adata = xclone.pp.tidy_Xdata(RDR_adata)
    xclone.pp.filter_pre(BAF_adata)
    xclone.pp.filter_pre(RDR_adata)
    """
    data_mode = Xdata.uns["data_mode"]
    if data_mode == "BAF":
        # AD = Xdata.layers["AD"]
        # DP = Xdata.layers["DP"]
        for layer_ in Xdata.layers:
            Xdata.obs[layer_+'_counts_per_obs'] = Xdata.layers[layer_].sum(axis=1)
            Xdata.obs[layer_+'_features_per_obs'] = np.count_nonzero(Xdata.layers[layer_].A, axis=1)
            Xdata.var[layer_+'_counts_per_feature'] = Xdata.layers[layer_].A.sum(axis=0)
            Xdata.var[layer_+'_obs_per_feature'] = np.count_nonzero(Xdata.layers[layer_].A, axis=0)
            log_msg = 'add annotation for per obs/features'
    elif data_mode == "RDR":
        Xdata.obs['RDR'+'_counts_per_obs'] = Xdata.X.sum(axis=1)
        Xdata.obs['RDR'+'_features_per_obs'] = np.count_nonzero(Xdata.X.A, axis=1)
        Xdata.var['RDR'+'_counts_per_feature'] = Xdata.X.A.sum(axis=0)
        Xdata.var['RDR'+'_obs_per_feature'] = np.count_nonzero(Xdata.X.A, axis=0)
        log_msg = 'add annotation for per obs/features'
    else:
        logging.error('error message: [XClone-filter_obs] Not xclonedata format')
        sys.exit(1)
    Xdata.uns["log"]["filter_pre"] = log_msg
    return Xdata

def filter_obs():
    ## to do list
    pass

def filter_features(
    Xdata: AnnData,
    inplace: bool = True,
    **kwargs):
    ## todo list
    data_mode = Xdata.uns["data_mode"]
    if data_mode == "BAF":
        for layer_ in Xdata.layers:
            Xdata.obs[layer_+'_counts_per_features'] = Xdata.layers[layer_].sum(axis=0)
    elif data_mode == "RDR":
        pass
        # Xdata = sc.pp.filter_genes(Xdata, inplace=inplace, **kwargs)
    else:
        logging.error('error message: [XClone-filter_features] Not xclonedata format')
        sys.exit(1)
    return Xdata

## II chromosome based functions--- subset 

def sub_chr(X, region_index, X_name="X",chr_list=-1,output_format="np.arr"):
    """
    Notes: To do list - not used in this preprocessing version 2021-07-23
    function: extract the subset of chromsomes;
    X:input data, can be np.array or np.matrix or sp.sparse_matrix.
    region_index: the corresponding data of X contains the chr info,1,2,3..;
    X_name: default data name "X", can use RDR,AD,DP, etc.
    chr_list: list of the the extracted chrs, support for [1,2,...,22];need update for X,Y;
    output_format: can be "np.arr","np.mtx","sp.sparse_mtx"(default:csr sparse matirx)
    return the extracted subset for testing/analysis
    """
    if type(X) == np.matrix or issparse(X):
        X=X.A
    for i in range(len(chr_list)):
        chr_idx = region_index[:,0] ==chr_list[i]
        print("chr", chr_list[i], sum(chr_idx))
        if i ==0:
            X_use = X[chr_idx,:]
        else:
            X_use = np.vstack((X_use, X[chr_idx, :]))
    if output_format=="np.mtx":
        X_use = np.mat(X_use)
    if output_format=="sp.sparse_mtx":
        X_use = sp.sparse.csr_matrix(X_use)
    print(X_name + "subset:",X_use.shape,"format:"+output_format)
    return X_use


def select_chr_region(X, region_index, X_name="X", mode="genome", select_list=-1, output_format="np.arr"):
    """
    ## base function for sub_features
    ## for blocks_based
    function: extract the subset of chromsomes; for hg19/38 genes/blocks;
    X:input data, can be np.array or np.matrix or sp.sparse_matrix.
    region_index: the corresponding data of X contains the chr info: 1,2,3.. and arm info: p, q;
    X_name: default data name "X", can use RDR,AD,DP, etc.
    select_list:chr list of the the extracted chrs, support for ['1','2',...,'22','X','Y'];have been updated for X,Y;
    ['1p','1q','2p',...,'22p','22q','Xp','Xq']----maybe need a better solution for select chr region
    gene_lst: update version in select_feature() function!
    output_format: can be "np.arr","np.mtx","sp.sparse_mtx"(default:csr sparse matirx)
    return the extracted subset for testing/analysis
    """
    if type(X) == np.matrix or issparse(X):
        X=X.A
    ## regions mode1
    if mode == "genome":
        print("regions mode1: genome...")
        if output_format=="np.mtx":
            X = np.mat(X)
        if output_format=="sp.sparse_mtx":
            X = sp.sparse.csr_matrix(X)
        
        select_idx = pd.Series(np.ones(X.shape[0], dtype= bool))
        print(X_name + " whole genome:", X.shape,"format:" + output_format)
        return select_idx, X
    ## regions mode2
    if mode == "select_chr":
        print("regions mode2: select_chr...")
        region_index["select_index"] = region_index["chr"]
    ## regions mode3
    if mode == "select_chr_arm":
        print("regions mode3: select_chr_arm...")
        region_index["select_index"] = region_index["chr"] + region_index["arm"]
    
    if (mode == "select_chr") or (mode == "select_chr_arm"):
        print("regions mode2 or mode3: select_chr or chr_arm...")
        for i in range(len(select_list)):
            select_flag = region_index["select_index"] == select_list[i]
            print("chr", select_list[i], sum(select_flag))
            if i ==0:
                X_use = X[select_flag,:]
                select_idx = select_flag
            else:
                X_use = np.vstack((X_use, X[select_flag, :]))
                select_idx = select_idx|select_flag
    
    ## regions mode4
    if mode == "select_chr_region":
        print("regions mode4: select_chr_region...")
        for i in range(len(select_list)):
            tmp_chr = select_list[i].split(":")[0]
            tmp_region = select_list[i].split(":")[1]
            tmp_start_loc = int(tmp_region.split("-")[0])
            tmp_end_loc = int(tmp_region.split("-")[1])
            chr_flag = region_index["chr"] == tmp_chr
            start_flag =  region_index["start"] >= tmp_start_loc
            end_flag = region_index["stop"] <= tmp_end_loc
            select_flag =  chr_flag & start_flag & end_flag
            print("chr", select_list[i], sum(select_flag))
            if i ==0:
                X_use = X[select_flag,:]
                select_idx = select_flag
            else:
                X_use = np.vstack((X_use, X[select_flag, :]))
                select_idx = select_idx|select_flag
    if output_format==None:
        return select_idx
    if output_format=="np.mtx":
        X_use = np.mat(X_use)
    if output_format=="sp.sparse_mtx":
        X_use = sp.sparse.csr_matrix(X_use)
    print(X_name + " subset:", X_use.shape, "format:" + output_format)
    # return X_use
    return select_idx, X_use

def select_features(X, feature_index=None, X_name="X", regions_mode="genome", chr_list=-1,
                   include_state=None, exclude_state=None, 
                   include_genes_lst=None, exclude_genes_lst=None, output_format="np.arr"):
    """
    ## base function for sub_features
    ## for genes_based
    function: extract the subset of chromsomes; v2 fit for hg38_genes and hg19_genes; not for blocks now;
    X:input data, can be np.array or np.matrix or sp.sparse_matrix.
    feature_index: the corresponding data of X contains the chr info: 1,2,3.. and arm info: p, q;
    X_name: default data name "X", can use RDR,AD,DP, etc.
    regions_mode: can be "hg19_blocks", "hg19_genes", "hg38_blocks", "hg38_genes", default: "genome";
    chr_list: chr list of the the extracted chrs, support for ['1','2',...,'22','X','Y'], have been updated for X,Y;
    ['1p','1q','2p',...,'22p','22q','Xp','Xq']----maybe need a better solution for select specific chr region(wait for update);
    gene_lst: include house keeping genes and exclude cell cycle genes;
    output_format: can be "np.arr","np.mtx","sp.sparse_mtx" (default:csr sparse matirx)
    return the extracted subset and the index(cells and features used) for testing/analysis

    ## usage
    RDR_idx, RDR_sub = select_feature(RDR, regions, "RDR", regions_mode, chr_lst, 
                            include_state=include_genes, exclude_state=exclude_genes,
                            include_genes_lst=include_genes_lst, exclude_genes_lst=exclude_genes_lst)
    RDR_idx, RDR_sub = select_feature(RDR, regions, "RDR", regions_mode, chr_lst, 
                            include_state=True, exclude_state=Flase,
                            include_genes_lst=include_genes_lst)
    """
    print("load select_feature function...")
    if type(X) == np.matrix or issparse(X):
        X=X.A
    ## select genes/items
    if include_state==True:
        if exclude_state==True:
            select_lst = set(include_genes_lst) - set(exclude_genes_lst)
            print("include %d genes, exclude %d genes, select %d  genes" 
                    %(len(include_genes_lst),len(exclude_genes_lst), len(select_lst)))
        else:
            select_lst = include_genes_lst
            print("include %d  genes" %(len(select_lst)))
    else:
        if exclude_state==True:
            select_lst = exclude_genes_lst
            print("exclude %d  genes" %(len(select_lst))) # print("exclude {:d}  genes".format(len(select_lst)))
        else:
            select_lst = feature_index["GeneName"].tolist()
            print("filtering no genes")  
    ## regions mode1
    if regions_mode == "genome":
        print("regions mode1: genome...")
        print(X_name + " whole genome:", X.shape, "format:" + str(output_format))

        select_idx =  feature_index["GeneName"].isin(select_lst)
        select_count = select_idx.sum() # select_idx-logistic value
        print("select count: %d" %(select_count))
        X_use = X[select_idx, :]
        print(X_name + " whole genome:", X_use.shape)
        ## output format
        if output_format== None:
            return select_idx
        if output_format=="np.mtx":
            X_use = np.mat(X_use)
        if output_format=="sp.sparse_mtx":
            X_use = sp.sparse.csr_matrix(X_use)
        return select_idx, X_use
    
    ## regions mode2
    if regions_mode == "select_chr":
        print("regions mode2: select_chr...")
        print(X_name + " before select_chr:", X.shape, "format:" + str(output_format))
        feature_index["select_index"] = feature_index["chr"]
    
    ## regions mode3
    if regions_mode == "select_chr_arm":
        print("regions mode2: select_chr_arm...")
        print(X_name + " before select_chr_arm:", X.shape, "format:" + str(output_format))
        feature_index["select_index"] = feature_index["chr"] + feature_index["arm"]
    
    if (regions_mode == "select_chr") or (regions_mode == "select_chr_arm"):
        for i in range(len(chr_list)):
            chr_flag = feature_index["select_index"] == chr_list[i]
            gene_flag = feature_index["GeneName"].isin(select_lst)
            select_flag = chr_flag&gene_flag
            print("chr", chr_list[i], sum(select_flag))
            if i ==0:
                X_use = X[select_flag,:]
                select_idx = select_flag
            else:
                X_use = np.vstack((X_use, X[select_flag, :]))
                select_idx = select_idx|select_flag
    ## regions mode4
    if regions_mode == "select_chr_region":
        for i in range(len(chr_list)):
            chr_flag = feature_index["select_index"] == chr_list[i]
            tmp_chr = chr_list[i].split(":")[0]
            tmp_region = chr_list[i].split(":")[1]
            tmp_start_loc = int(tmp_region.split("-")[0])
            tmp_end_loc = int(tmp_region.split("-")[1])
            chr_flag = feature_index["chr"] == tmp_chr
            start_flag =  feature_index["start"] >= tmp_start_loc
            end_flag = feature_index["stop"] <= tmp_end_loc
            gene_flag = feature_index["GeneName"].isin(select_lst)
            select_flag =  chr_flag & start_flag & end_flag & gene_flag
            print("chr", chr_list[i], sum(select_flag))
            if i ==0:
                X_use = X[select_flag,:]
                select_idx = select_flag
            else:
                X_use = np.vstack((X_use, X[select_flag, :]))
                select_idx = select_idx|select_flag
    ## output format
    if output_format==None:
        return select_idx
    if output_format=="np.mtx":
        X_use = np.mat(X_use)
    if output_format=="sp.sparse_mtx":
        X_use = sp.sparse.csr_matrix(X_use)
    print(X_name + " subset:", X_use.shape, "format:" + output_format)
    return select_idx, X_use

def sub_features(
    Xdata: AnnData,
    regions_mode,
    chr_lst,
    inplace: bool = True,
    **kwargs):
    """
    Function:
    
    -----
    params:
    data: xclonedata formatm AnnData
          The (annotated) data matrix of shape `n_obs` × `n_vars`.
          Rows correspond to cells/samples and columns to genes/blocks/features.
    regions_mode:
    chr_lst:
    ------
    optional params:
    include_state, default=None
    exclude_state, default=None
    include_genes_lst, default=None
    exclude_genes_lst, default=None
    -----
    Examples:
    hk_genes = xclone.pp.load_hk_genes()
    cc_genes = xclone.pp.load_cc_genes()
    test1 = xclone.pp.sub_features(RDR_adata,regions_mode="select_chr",chr_lst=['1','2'])
    test2 = xclone.pp.sub_features(RDR_adata,regions_mode="select_chr",chr_lst=['1','2'], 
                                   include_state=True, include_genes_lst=hk_genes["GeneName"])
    test3 = xclone.pp.sub_features(RDR_adata,regions_mode="select_chr",chr_lst=['1','2'], 
                                   include_state=True, include_genes_lst=hk_genes["GeneName"], 
                                   exclude_state=True, exclude_genes_lst=cc_genes["GeneName"])
    """
    data_mode = Xdata.uns["data_mode"]
    genome_mode = Xdata.uns["genome_mode"]
    if data_mode == "BAF":
        if "genes" in genome_mode:
            select_idx= select_features(Xdata.X.T, Xdata.var, "BAF", regions_mode=regions_mode, chr_list=chr_lst, output_format=None, **kwargs)
        elif "blocks" in genome_mode:
            select_idx = select_chr_region(Xdata.X.T, Xdata.var, X_name="BAF", mode=regions_mode, select_list=chr_lst, output_format=None)
        else:
            logging.error('error message: [XClone-sub_feature] Not right genome_mode info in xclonedata')
            sys.exit(1)
        # TO DO LIST UPDATE select_feature, and can remove the first params
        # select_idx = select_feature(Xdata.X.T, Xdata.var, "BAF", regions_mode=regions_mode, chr_list=chr_lst,include_state=None, exclude_state=None, 
                #    include_genes_lst=None, exclude_genes_lst=None, output_format=None)
    elif data_mode == "RDR":
        if "genes" in genome_mode:
            select_idx= select_features(Xdata.X.T, Xdata.var, "RDR", regions_mode=regions_mode, chr_list=chr_lst, output_format=None, **kwargs)
        elif "blocks" in genome_mode:
            select_idx = select_chr_region(Xdata.X.T, Xdata.var, X_name="RDR", mode=regions_mode, select_list=chr_lst, output_format=None)
        else:
            logging.error('error message: [XClone-sub_feature] Not right genome_mode info in xclonedata')
            sys.exit(1)
    else:
        logging.error('error message: [XClone-sub_feature] Not xclonedata format')
        sys.exit(1)
    log_msg = '[sub_features] ' + str(chr_lst)
    Xdata.uns["log"]["sub_feature"] = log_msg
    return Xdata[:,select_idx]

def sub_CellCluster(mtx_data, mtx_barcodes_file,node_barcodes_lst):
    """
    ## can be deprecated in preprocessing anndata
    Function: get the cells subset for specific node/cells in specific barcodeslst
    
    mtx_data: like AD/DP/RDR
    mtx_barcodes_file: also the output from xcltk
    node_barcodes_lst:the cells need to be extrcted
    
    Example:
    mtx_barcodes_file = baf_dat_dir + "cellSNP.samples.tsv"
    AD_node = sub_CellCluster(AD, mtx_barcodes_file,barcodes_8805_lst)
    DP_node = sub_CellCluster(DP, mtx_barcodes_file,barcodes_8805_lst)
    
    """
    mtx_barcodes = pd.read_table(mtx_barcodes_file,header = None)
    mtx_barcodes.columns = ["cell"]
    node_flag = mtx_barcodes["cell"].isin(node_barcodes_lst)
    node_idx = mtx_barcodes[node_flag].index
    sub_cluster_mtx_T = mtx_data.T[node_idx]
    sub_cluster_mtx = sub_cluster_mtx_T.T
    num_cells = sub_cluster_mtx.shape[1]
    num_blocks = sub_cluster_mtx.shape[0]
    print("sub_cells:",num_cells,";","blocks/genes:",num_blocks)
    return sub_cluster_mtx

def sub_cells(
    Xdata: AnnData,
    cellbarcodes_lst,
    exclude: bool = False):
    """
    Function:
    
    -----
    params:
    data: xclonedata formatm AnnData
          The (annotated) data matrix of shape `n_obs` × `n_vars`.
          Rows correspond to cells/samples and columns to genes/blocks/features.
    cellbarcode_lst:
    exclude: default False, means subsetting the cells in the cellbarcodes_lst;
             if True, exclude the cells
    -----
    optional params:

    -----
    Examples:
    clone1_barcodes_lst = ['AAACCTGCACCTTGTC-1', 'AAACGGGAGTCCTCCT-1']
    BAF_adata1 = sub_cells(BAF_adata, clone1_barcodes_lst)
    BAF_adata1 = sub_cells(BAF_adata, clone1_barcodes_lst,exclude=True)
    """
    mtx_barcodes = Xdata.obs.index
    select_idx = mtx_barcodes.isin(cellbarcodes_lst)
    if exclude == True:
        select_idx = ~select_idx
    return Xdata[select_idx,:]