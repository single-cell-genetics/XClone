"""Utility functions for XClone results analysis
"""
# Author: Rongting Huang
# Date: 2021/03/15
# update: 2021/07/22

## data preprocessing
import numpy as np
import scipy as sp
from scipy.sparse import issparse

import warnings

def filter_data(X,X_name="X", axis=1, output_format="np.arr"):
    """
    function: filter the all zero/nan rows or columns;
    X:input data, can be np.array or np.matrix or sp.sparse_matrix.
    default axis: 1,filter the rows without information; 0, filter the cols
    default data name:X_name="X", can use RDR,AD,DP, etc.
    return filtered data:output format can be "np.arr","np.mtx",
    "sp.sparse_mtx"(default:csr sparse matirx)
    """
    if type(X) == np.matrix or issparse(X):
        X=X.A
    print(X_name,X.shape)
    X_mask = np.all(np.isnan(X) | np.equal(X, 0), axis=axis)
    if axis==1:
        X_filter = X[~X_mask]
    elif axis==0:
        X_filter = X[:,~X_mask]
    print(X_name," filter 0 and nan:", X_filter.shape)
    if output_format=="np.mtx":
        X_filter = np.mat(X_filter)
    if output_format=="sp.sparse_mtx":
        X_filter = sp.sparse.csr_matrix(X_filter)
    warnings.warn('[XClone-analysis]-DeprecationWarning', DeprecationWarning)
    return X_filter

def filter_2data(X,Y,X_name="X", Y_name="Y",axis=1, output_format="np.arr"):
    """
    function: filter the all zero/nan rows or columns for 2 dataset;
    X,Y:input data, can be np.array or np.matrix or sp.sparse_matrix.
    default axis: 1, filter the rows without information; 0, filter the cols
    default data name:X_name="X",Y_name="Y", can use AD,DP, etc.
    return filtered data:output format can be "np.arr","np.mtx",
    "sp.sparse_mtx"(default:csr sparse matirx)
    """
    if type(X) == np.matrix or issparse(X):
        X=X.A
    if type(Y) == np.matrix or issparse(Y):
        Y=Y.A
    print(X_name,X.shape)
    print(Y_name,Y.shape)
    # X_mask = np.all(np.isnan(X) | np.equal(X, 0), axis=axis)
    Y_mask = np.all(np.isnan(Y) | np.equal(Y, 0), axis=axis)
    # XY_mask = X_mask|Y_mask
    XY_mask = Y_mask
    if axis==1:
        X_filter = X[~XY_mask]
        Y_filter = Y[~XY_mask]
    elif axis==0:
        X_filter = X[:,~XY_mask]
        Y_filter = Y[:,~XY_mask]
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
    warnings.warn('[XClone-analysis]-DeprecationWarning', DeprecationWarning)
    return X_filter, Y_filter

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
    warnings.warn('[XClone-analysis]-DeprecationWarning', DeprecationWarning)
    return X_filter_idx, X_filter

def filter_2nulldata(X,Y,X_name="X", Y_name="Y",axis=1, output_format="np.arr"):
    """
    function: filter the all zero/nan rows or columns for 2 dataset;
    X,Y:input data, can be np.array or np.matrix or sp.sparse_matrix.
    default axis: 1, filter the rows without information; 0, filter the cols
    default data name: X_name="X",Y_name="Y", can use AD,DP, etc.
    return filtered data:output format can be "np.arr","np.mtx",
    "sp.sparse_mtx"(default:csr sparse matirx)
    """
    if type(X) == np.matrix or issparse(X):
        X=X.A
    if type(Y) == np.matrix or issparse(Y):
        Y=Y.A
    print(X_name, X.shape)
    print(Y_name, Y.shape)
    # X_mask = np.all(np.isnan(X) | np.equal(X, 0), axis=axis)
    Y_mask = np.all(np.isnan(Y) | np.equal(Y, 0), axis=axis) # Y is DP here
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
    warnings.warn('[XClone-analysis]-DeprecationWarning', DeprecationWarning)
    return XY_filter_idx, X_filter, Y_filter

def sub_chr(X, region_index, X_name="X",chr_list=-1,output_format="np.arr"):
    """
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

## block based
def select_chr_region(X, region_index, X_name="X", mode="genome", select_list=-1, output_format="np.arr"):
    """
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

    if output_format=="np.mtx":
        X_use = np.mat(X_use)
    if output_format=="sp.sparse_mtx":
        X_use = sp.sparse.csr_matrix(X_use)
    print(X_name + " subset:", X_use.shape, "format:" + output_format)
    # return X_use
    return select_idx, X_use


## gene_based
def select_feature(X, feature_index=None, X_name="X", regions_mode="genome", chr_list=-1,
                   include_state=None, exclude_state=None, 
                   include_genes_lst=None, exclude_genes_lst=None, output_format="np.arr"):
    """
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
    ## select genes
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
        print(X_name + " whole genome:", X.shape, "format:" + output_format)

        select_idx =  feature_index["GeneName"].isin(select_lst)
        select_count = select_idx.sum() # select_idx-logistic value
        print("select count: %d" %(select_count))
        X_use = X[select_idx, :]
        print(X_name + " whole genome:", X_use.shape)
        ## output format
        if output_format=="np.mtx":
            X_use = np.mat(X_use)
        if output_format=="sp.sparse_mtx":
            X_use = sp.sparse.csr_matrix(X_use)
        return select_idx, X_use
    
    ## regions mode2
    if regions_mode == "select_chr":
        print("regions mode2: select_chr...")
        print(X_name + " before select_chr:", X.shape, "format:" + output_format)
        feature_index["select_index"] = feature_index["chr"]
    
    ## regions mode3
    if regions_mode == "select_chr_arm":
        print("regions mode2: select_chr_arm...")
        print(X_name + " before select_chr_arm:", X.shape, "format:" + output_format)
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
    if output_format=="np.mtx":
        X_use = np.mat(X_use)
    if output_format=="sp.sparse_mtx":
        X_use = sp.sparse.csr_matrix(X_use)
    print(X_name + " subset:", X_use.shape, "format:" + output_format)
    return select_idx, X_use



import pandas as pd
import matplotlib.pyplot as plt

def get_node_barcodeslst(dloupe_node_file, node_id):
    """
    Function: get the barcodes list for specific node (cell cluster)
    
    dloupe_node_file: path of the file, defalut for dloupe output csv file,
    which contains node_id, barcodes, num_cells, num_noisy and blabla...
    node_id: int, the ID of the node need to be extrcted.
    
    example:
    mkn45_10x_dloupe = "/storage/yhhuang/research/mito/mkn45/fulldepth/mode2/
    mkn45_5k_dloupe-group_10438-region_1_1-5120000_Y_56320001-59373566-heatmap_copy number.csv"
    barcodes_10050_lst = get_node_barcodeslst(mkn45_10x_dloupe, 10050)
    """
    data_10x_cluster = pd.read_table(dloupe_node_file,sep=",")
    node_data = data_10x_cluster[data_10x_cluster["node_id"] == node_id]
    node_barcodes = node_data["barcodes"].str.split(";")
    node_barcodes_lst = list(node_barcodes)[0]
    print("node:", node_id,"num_cells:",len(node_barcodes_lst))
    warnings.warn('[XClone-analysis]-DeprecationWarning', DeprecationWarning)
    return node_barcodes_lst


def get_non_node_barcodeslst(dloupe_node_file, node_id):
    """
    Function: get the barcodes list for excluded node (cell cluster)
    specific for CellRanger dloupe output csv file.
    
    dloupe_node_file: path of the file, defalut for dloupe output csv file,
    which contains node_id, barcodes, num_cells, num_noisy and blabla...
    node_id: int, the ID of the excluded node need to be extrcted.
    
    example:
    mkn45_10x_dloupe = "/storage/yhhuang/research/mito/mkn45/fulldepth/mode2/
    mkn45_5k_dloupe-group_10438-region_1_1-5120000_Y_56320001-59373566-heatmap_copy number.csv"
    barcodes_non_node_lst = get_non_node_barcodeslst(mkn45_10x_dloupe, 10050)
    """
    data_10x_cluster = pd.read_table(dloupe_node_file,sep=",")
    Flag_ = data_10x_cluster["node_id"] != node_id
    node_data = data_10x_cluster[Flag_]
    node_barcodes = node_data["barcodes"].str.split(";")
    node_barcodes_lst = []
    for idx_ in node_data.index:
        for i in range(len(node_barcodes[idx_])):
            node_barcodes_lst.append(node_barcodes[idx_][i])     
    print("non_node:", node_id,"num_cells:",len(node_barcodes_lst))
    warnings.warn('[XClone-analysis]-DeprecationWarning', DeprecationWarning)
    return node_barcodes_lst

def processing_dlouple_file(dloupe_node_file, save_status = False, save_file = "xxxx_scDNA_cellranger_node.txt"):
    """
    Function: get the barcodes list from the cellranger scDNA dloupe heatmap csv file
    for each node (cell cluster)
    
    dloupe_node_file: path of the file, defalut for dloupe output csv file,
    which contains node_id, barcodes, num_cells, num_noisy and blabla...
    save_status: default False for saving the result for future processing
    
    return:
    barcodes and corresponding node_id, order is the same with the heatmap in loupe scDNA Browser.
    
    example:
    mkn45_10x_dloupe = "/storage/yhhuang/research/mito/mkn45/fulldepth/mode2/
    mkn45_5k_dloupe-group_10438-region_1_1-5120000_Y_56320001-59373566-heatmap_copy number.csv"
    cellranger_resutls = processing_dlouple_file(mkn45_10x_dloupe, 10050)
    """
    data_10x = pd.read_table(dloupe_node_file, sep=",")
    data_10x["node_barcodes_lst"] = data_10x["barcodes"].str.split(";")
    node_barcodes_lst = []
    node_lst = []
    for idx_ in data_10x.index:
        node_barcodes_lst_tmp = data_10x["node_barcodes_lst"][idx_]
        node_tmp = data_10x["node_id"][idx_]
        for i in range(len(node_barcodes_lst_tmp)):
            node_barcodes_lst.append(node_barcodes_lst_tmp[i])
            node_lst.append(node_tmp)
    cellranger_data = {'barcodes':node_barcodes_lst,
                     'node_id': node_lst}
    cellranger_results = pd.DataFrame(cellranger_data)
    if save_status:
        cellranger_results.to_csv(save_file, sep="\t", index=False)
    warnings.warn('[XClone-analysis]-DeprecationWarning', DeprecationWarning)
    return cellranger_results

def get_barcodeslst(cellranger_results_file, node_id_lst, include=True, return_index=False):
    """
    Function:
    cellranger_results_file: processed and saved cellbarcodes and corresponding node_id
    
    usage:
    node_barcodes_lst = get_barcodeslst(cellranger_results_file,[2710,2701],include=False)
    """
    cellranger_results = pd.read_csv(cellranger_results_file, sep="\t")
    if include:
        flag_ = cellranger_results["node_id"].isin(node_id_lst)
    else:
        flag_ = ~(cellranger_results["node_id"].isin(node_id_lst))
    barcodes_lst = cellranger_results_test["barcodes"][flag_].tolist()
    print("node_id_lst: ", node_id_lst)
    print("include_node_status: ", include)
    print("get barcodeslst cell numbers: ", len(barcodes_lst))
    barcodes_lst_index = flag_.index
    if return_index:
        return barcodes_lst_index, barcodes_lst
    warnings.warn('[XClone-analysis]-DeprecationWarning', DeprecationWarning)
    return barcodes_lst

def annotate_node(cellranger_results, node_cut_id):
    """
    Function: annotate the cellranger clusters with the node_cut_id
    cellranger_results: the dataframe with the cellbarcodes and correspoding node_id;
    node_cut_id: the node_id lst for cutting the clusters, [ ) format
    return:
    new dataframe with annotated cell clusters in the new column "cluster"
    
    usage:
    cellranger_cut_results = annotate_node(cellranger_results, [2536])
    cellranger_cut_results = annotate_node(cellranger_results, [2536,2543])
    
    """
    flag_ = cellranger_results["node_id"].drop_duplicates().isin(node_cut_id)
    idx_ = cellranger_results["node_id"].drop_duplicates().isin(node_cut_id).index
    cut_idx_ = idx_[flag_]
    cluster_num = len(cut_idx_)
    cellranger_results["cluster"] = cluster_num
    print(cluster_num)
    for i in reversed(range(cluster_num)):
        tmp_cut = cut_idx_[i]
        cellranger_results["cluster"][cellranger_results.index < tmp_cut ] = i
    warnings.warn('[XClone-analysis]-DeprecationWarning', DeprecationWarning)
    return cellranger_results


def sub_CellCluster(mtx_data, mtx_barcodes_file,node_barcodes_lst):
    """
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
    warnings.warn('[XClone-analysis]-DeprecationWarning', DeprecationWarning)
    return sub_cluster_mtx

def get_region_position(regions,chr):
    """
    Function: get the start position for chromosome in the genome coordinate
    regions:chr\t block_start \t block_end for each line
    chr: int, specific chr num
    Example:
    baf_dat_dir = '/storage/yhhuang/users/rthuang/processed_data/xcltk/xianjie-cpos/mkn45_020821/mkn45-500k/mkn45-500k-csp-post/phase-snp-even/'
    regions = np.genfromtxt(baf_dat_dir + "blocks.50kb.tsv", delimiter="\t")
    get_region_position(regions,4)
    """
    flag_regions = list(regions[:,0] == chr)
    first_index = flag_regions.index(True)
    ## TODO THERE SHOULD BE ANOTHER FUNCTION TO GET THE FIRST INDEX DIRECTLY #rt
    print("Start position for chr",chr,": ",first_index)
    warnings.warn('[XClone-analysis]-DeprecationWarning', DeprecationWarning)
    return first_index


## comparison analysis-confusion matrix and visualization
# import pandas as pd
# import numpy as np
def compare_clusters(model_proba, mtx_barcodes_file, node_barcodes_lst):
    """
    Function:
    
    Example:
    #### data
    filename = "mkn45-500k-BAF-chr2/_model_baf_sub2_6.model"
    _model_baf_sub2_6 = joblib.load(filename)
    
    baf_dat_dir = '/storage/yhhuang/users/rthuang/processed_data/xcltk/xianjie-cpos/mkn45_020821/mkn45-500k/mkn45-500k-csp-post/phase-snp-even/'
    mtx_barcodes_file = baf_dat_dir + "cellSNP.samples.tsv"
    
    mkn45_10x_dloupe = "/storage/yhhuang/research/mito/mkn45/fulldepth/mode2/mkn45_5k_dloupe-group_10438-region_1_1-5120000_Y_56320001-59373566-heatmap_copy number.csv"
    barcodes_8805_lst = get_node_barcodeslst(mkn45_10x_dloupe, 8805)
    
    ## usage
    compare_clusters(_model_baf_sub2_6.proba_, mtx_barcodes_file, barcodes_8805_lst)
    
    """
    cluster_index_lst = np.argmax(model_proba, axis=1)
    uni_cluster = np.unique(cluster_index_lst)
    mtx_barcodes_lst = pd.read_table(mtx_barcodes_file,header = None)
#     print(mtx_barcodes_lst)
    for uni_id in uni_cluster:
        print("cluster",uni_id, ":", np.sum(cluster_index_lst == uni_id))
        uni_cluster_barcodes = get_cluster_barcodes_lst(uni_id, cluster_index_lst, mtx_barcodes_lst)
#         print(uni_cluster_barcodes)
        inter_set = (set(uni_cluster_barcodes)&set(node_barcodes_lst))
#         print(inter_set)
        print("inter_set", len(inter_set))

def get_cluster_barcodes_lst(cluster_id, cluster_index_lst, mtx_barcodes_lst):
    """
    Function:
    
    Example:
    uni_cluster_barcodes = get_cluster_barcodes_lst(1, cluster_index_lst, mtx_barcodes_lst)
    """
    if len(cluster_index_lst) != len(mtx_barcodes_lst):
        print("error! No matching in cluster_index and mtx_barcodes")
    flag_ = cluster_index_lst ==  cluster_id
    barcodes_ = mtx_barcodes_lst[flag_]
    return list(barcodes_[0])


def identified_node_lst(mtx_barcodes_lst, barcodes_node_lst, node_ID="node_10050"):
    """
    Function: generate list for confusion matrix building
    mtx_barcodes_lst: ref list with the same order as the matrix used
    barcodes_node_lst: barcodes for specified node, like 10050
    
    Example:
    baf_dat_dir = '/storage/yhhuang/users/rthuang/processed_data/xcltk/xianjie-cpos/mkn45_021021/mkn45-fulldepth/mkn45-fulldepth-csp-post/phase-snp-even/'
    mtx_barcodes_file = baf_dat_dir + "cellSNP.samples.tsv"
    mtx_barcodes_lst = pd.read_table(mtx_barcodes_file,header = None)

    mkn45_10x_dloupe = "/storage/yhhuang/research/mito/mkn45/fulldepth/mode2/mkn45_5k_dloupe-group_10438-region_1_1-5120000_Y_56320001-59373566-heatmap_copy number.csv"
    barcodes_10050_lst = get_node_barcodeslst(mkn45_10x_dloupe, 10050)
    node_lst = identified_node_lst(mtx_barcodes_lst, barcodes_10050_lst, node_ID="node_10050")
    """
    bool_node_lst = mtx_barcodes_lst.isin(barcodes_node_lst)
    print("node_cells_num",bool_node_lst.sum())
    node_lst_ = bool_node_lst.applymap(lambda x: node_ID if x == True else "non_" + node_ID)[0].values
    return node_lst_

### from vireoSNP
### https://github.com/single-cell-genetics/vireo/blob/master/vireoSNP/utils/base_utils.py#L3
def get_confusion(ids1, ids2):
    """Get confusion matrix
    
    Parameters
    ----------
    ids1: numpy.array or list
        id list in the first annotation
    ids2: numpy.array or list
        id list in the second annotation
        
    Return
    ------
    (confuse_mat, ids1_uniq, ids2_uniq)
    confuse_mat[i, j]: 
        number of samples have ids1 == ids1_uniq[i]
        and ids2 == id2_uniq[j]
    """
    if type(ids1) == list: ids1 = np.array(ids1)
    if type(ids2) == list: ids2 = np.array(ids2)
    
    ids1_uniq = np.unique(ids1)
    ids2_uniq = np.unique(ids2)
    
    confuse_mat = np.zeros((len(ids1_uniq), len(ids2_uniq)), dtype=int)
    for i, _id1 in enumerate(ids1_uniq):
        for j, _id2 in enumerate(ids2_uniq):
            confuse_mat[i, j] = np.sum((ids1 == _id1) * (ids2 == _id2))
            
    return confuse_mat, ids1_uniq, ids2_uniq

def get_confuse_mat_df(confuse_mat,index_names=None,clolums_names=None):
    """
    Function: construct df for confuse_mat for heatmap plot
    ========================
    Example:
    confuse_mat_df = get_confuse_mat_df(confuse_mat,index_names=["cluster0", "cluster1"],clolums_names=["node_10050", "non node_10050"])
    """
    confuse_mat_df = pd.DataFrame(confuse_mat)
    if index_names:
        confuse_mat_df.index = index_names
    if clolums_names:
        confuse_mat_df.columns = clolums_names
    return confuse_mat_df
    

# def confuse_heatmap(confuse_mat_df, cmap = "Blues", version1 = True, version2 =True, save_pdf = None):
#     """
#     Function:
#     confusion matrix heatmap
    
#     confuse_mat_df: count- confusion matrix in pd.DataFrame
    
#     Example：
#     confuse_heatmap(confuse_mat_df, version1 = False, save_pdf = True)
#     """
#     import matplotlib.pyplot as plt
#     import seaborn as sns
#     ## setting
#     font = {'family' : 'DejaVu Sans',
#             'size'   : 15}
#     plt.rc('font', **font)
#     if version1:
#         sns.heatmap(confuse_mat_df, annot=True, fmt="d", cmap=cmap)
#     if version2:
#         ## data normalization
#         plot_df = (confuse_mat_df.T/confuse_mat_df.sum(axis=1).values).T
#         ## plot percentage
#         plt.figure(figsize=[8,6])
#         ax = sns.heatmap(plot_df, cmap=cmap) # annot=True,
#         ## annot original count
#         height, width = np.array(confuse_mat_df).shape
#         text_colors = ['black', 'white']
#         for x in range(width):
#             for y in range(height):
#                 ax.annotate(str(np.array(confuse_mat_df)[y][x]), xy=(x+0.5, y+0.5), 
#                             horizontalalignment='center',
#                             verticalalignment='center', color=text_colors[int(np.array(plot_df)[y][x] > 0.5)], fontsize = 15)
#         plt.xlabel('CellRanger nodes')
#         plt.ylabel('CNV clones')
#         # plt.yticks(rotation=45)
#         # plt.xticks(rotation=45)
#         # plt.xticks(range(3), confusion_matrix.columns) #, rotation=315)
#         # plt.yticks(range(len(norm_conf)), set(clone_id))
#         plt.tight_layout()
#         plt.title('Confusion matrix for comparing CNV clones with node_10050', fontsize = 18)
#         if save_pdf:
#             plt.savefig('conf_mat_CNV_2clones-node10050_update.pdf', dpi=200, bbox_inches='tight')

## confusion matrix heatmap
def confuse_heatmap(confuse_mat_df, cmap = "Blues", version1 = True, version2 =True, plot_xlabel = None,
                    plot_ylabel = None, plt_title = None, save_file_name = None):
    """
    Function:
    confusion matrix heatmap
    confuse_mat_df: count- confusion matrix in pd.DataFrame
    
    Example：            
    confuse_mat, ids1_uniq, ids2_uniq = get_confusion(CopyKAT_lst, expression_lst)
    confuse_mat_df = get_confuse_mat_df(confuse_mat,
                                        index_names=["cloneA", "cloneB","Normal"],
                                        clolums_names=["cloneA", "cloneB","Normal"])
    confuse_heatmap(confuse_mat_df,  
                    plot_xlabel = "CopyKAT",
                    plot_ylabel = "Ground Truth", 
                    plt_title = "Concordance in subclone identification", 
                    save_file_name = 'CopyKAT_vs_Groudtruth1.pdf')
    """
    import matplotlib.pyplot as plt
    import seaborn as sns
    ## setting
    font = {'family' : 'DejaVu Sans',
            'size'   : 15}
    plt.rc('font', **font)
    if version1:
        plt.figure(figsize=[8,6])
        sns.heatmap(confuse_mat_df, annot=True, fmt="d", cmap=cmap)
    if version2:
        ## data normalization
        plot_df = (confuse_mat_df.T/confuse_mat_df.sum(axis=1).values).T
        ## plot percentage
        plt.figure(figsize=[8,6])
        ax = sns.heatmap(plot_df, cmap=cmap) # annot=True,
        ## annot original count
        height, width = np.array(confuse_mat_df).shape
        text_colors = ['black', 'white']
        for x in range(width):
            for y in range(height):
                ax.annotate(str(np.array(confuse_mat_df)[y][x]), xy=(x+0.5, y+0.5), 
                            horizontalalignment='center',
                            verticalalignment='center', color=text_colors[int(np.array(plot_df)[y][x] > 0.5)], fontsize = 15)
        if plot_xlabel is None:
            pass
        else:
            plt.xlabel(plot_xlabel)
        if plot_ylabel is None:
            pass
        else:
            plt.ylabel(plot_ylabel)
        # plt.yticks(rotation=45)
        # plt.xticks(rotation=45)
        # plt.xticks(range(3), confusion_matrix.columns) #, rotation=315)
        # plt.yticks(range(len(norm_conf)), set(clone_id))
        plt.tight_layout()
        if plt.title is None:
            pass
        else:
            plt.title(plt_title, fontsize = 18) #plt.title('Concordance in subclone identification', fontsize = 18)
        if save_file_name is None:
            pass
        else:
            plt.savefig(save_file_name, dpi=300, bbox_inches='tight') 
             #plt.savefig('CopyKAT_vs_Groudtruth.pdf', dpi=300, bbox_inches='tight')
              


## visualization version1
def BAF_plot(data,xlabel,ylabel,title,display_state="show"):
    """
    Function: plot BAF-heatmap
    display_state: default:show  alternative: "save",saved as jpg
    Example:
    BAF_plot(data=AD.T/DP.T,xlabel="blocks in genome",ylabel="cells",title="node 8805",display_state="show")
    """
    fig = plt.figure(figsize=(6, 5))
    im = plt.imshow(data, interpolation='none', 
                aspect='auto', cmap='coolwarm') #seismic, bwr
    plt.colorbar()
    plt.tight_layout()
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)
    if display_state == "show":
        plt.show()
    if display_state == "save":
        plt.savefig(title + ".jpg")

from scipy.sparse import issparse

def plot_BAF_scatter(AD,DP,xlabel,ylabel,title,display_state="show"):
    """
    Function: plot BAF for subcluster-scatter plot
    
    Example1:
    plot_BAF_scatter(AD_sub4_node,DP_sub4_node,xlabel="chr4",ylabel="BAF",title="BAF_node8805_chr4")
    Example2:
    plot_BAF_scatter(AD_node,DP_node,xlabel="genome",ylabel="BAF",title="BAF_node8805_genome")
    ### TODO-get the accurate position and then display it(maybe try to use the pd.Series's index)
    
    """
    ## preprocessing-merge all cells and filter nan
    print("AD shape:", AD.shape)
    print("DP shape:", DP.shape)
    if type(AD) == np.matrix or issparse(AD):
        AD = AD.A
    if type(DP) == np.matrix or issparse(DP):
        DP = DP.A
    print("AD type:", type(AD))
    print("DP type:", type(DP))
    ### merge all cells to get BAF
    AD_sum = AD.sum(axis=1)
    DP_sum = DP.sum(axis=1)
    ### Filter nan value
    AD_FLAG = ~(AD_sum == 0)
    DP_FLAG = ~(DP_sum == 0)
    FLAG_ = AD_FLAG & DP_FLAG
    AD_sum_filter = AD_sum[FLAG_]
    DP_sum_filter = DP_sum[FLAG_]
    print("filter AD",len(AD_sum_filter))
    print("filter DP",len(DP_sum_filter))
    
    x_len = len(AD_sum_filter)
    y_data = AD_sum_filter/DP_sum_filter
    plt.figure(figsize=(24, 6))
    plt.scatter(range(x_len), y_data, alpha=0.5)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)


def compare_plot_BAF(AD,DP, model_proba, mtx_barcodes_file,regions_file, chr_lst,xlabel1="genome",xlabel2="chr4", plot_title="BAF_cluster_chr4(6-clusters)"):
    """
    Example:
    compare_plot_BAF(AD,DP, model_proba, mtx_barcodes_file,regions_file, [4])
    
    """
    cluster_index_lst = np.argmax(model_proba, axis=1)
    uni_cluster = np.unique(cluster_index_lst)
    mtx_barcodes_lst = pd.read_table(mtx_barcodes_file,header = None)
    regions = np.genfromtxt(regions_file, delimiter="\t")
    for uni_id in uni_cluster:
        uni_cluster_barcodes = get_cluster_barcodes_lst(uni_id, cluster_index_lst, mtx_barcodes_lst)
        AD_cluster_tmp = sub_CellCluster(AD, mtx_barcodes_file, uni_cluster_barcodes)
        DP_cluster_tmp = sub_CellCluster(DP, mtx_barcodes_file, uni_cluster_barcodes)
        AD_sub_cluster = sub_chr(AD_cluster_tmp, regions,"AD", chr_lst)
        DP_sub_cluster = sub_chr(DP_cluster_tmp, regions,"DP", chr_lst)
        plot_BAF_scatter(AD_cluster_tmp, DP_cluster_tmp, xlabel=xlabel1, ylabel="BAF", title= plot_title)
        plot_BAF_scatter(AD_sub_cluster, DP_sub_cluster, xlabel=xlabel2, ylabel="BAF", title= plot_title)



## visualization version2
def get_BAF_plot_data(regions_file, AD, DP,filter_nan=True):
    """
    Function: get the fitted format for plot BAF scatter
    regions_file:contain the chr and also the blocks info for each chr
    AD: AD information or all cells/subcluster cells, mtx
    DP: DP information or all cells/subcluster cells, mtx
    
    Example:
    
    """
    regions_df = pd.read_table(regions_file, header = None)
    regions_df.columns = ["chrom","chr_block_start","chr_block_stop"]
#     regions_df["chrom"] = regions_df["chrom"].astype('category')
    regions_df["genome_pos"] = regions_df.index
    regions_df["pos_in_chr"] = regions_df["chr_block_stop"]/50000
    regions_df["pos_in_chr"] = regions_df["pos_in_chr"].astype(np.int64)
    ## add AD and DP
    ### merge AD and DP for cells to get the bulk level
    AD_sum = AD.sum(axis=1)
    DP_sum = DP.sum(axis=1)
    regions_df["AD"] = pd.Series(AD_sum.A[:,0])
    regions_df["DP"] = pd.Series(DP_sum.A[:,0])
    ## filter tha nan
    if filter_nan:
#         AD_flag = ~(regions_df['AD'] == 0)
        DP_flag = ~(regions_df['DP'] == 0)
#         FLAG_ = AD_flag & DP_flag
        FLAG_ = DP_flag
        regions_df = regions_df[FLAG_]
    ### get the BAF_value
    regions_df.loc[:,("BAF")] = regions_df.loc[:,("AD")] / regions_df.loc[:,("DP")]
    ## return the dataframe for plot
    return regions_df


def BAF_scatter_plot2(BAF_df,plot_title, chr_lst=['2','4'], sub_chr_plot=False, filenotes = "filename"):
    import limix_plot as lp
    ## todo- can change code in limix_plot@Rongtingting
#     df_plot = BAF_df[["chrom","genome_pos", "BAF"]]
#     df_plot = df_plot.rename(columns={"genome_pos": "pos"})
    df_plot = BAF_df[["chrom","pos_in_chr", "BAF"]]
    df_plot = df_plot.rename(columns={"pos_in_chr": "pos"})
    if sub_chr_plot:
        df_plot = df_plot[df_plot["chrom"].isin(chr_lst)]
        plt = lp.get_pyplot()
        manhattan(df_plot,colora="#75E6DA", colorb="#EF7C8E", pts_kws = {"markersize": 8, "marker":'.'})  
        plt.axhline(0.5, color='#167D7F')
        plt.title(plot_title, fontsize=18)
        plt.savefig(filenotes + "sub_chr.pdf")
#         plt.xlabel(xlabel)
#         plt.ylabel(ylabel)
    else:
        manhattan(df_plot,colora="#F79489", colorb="#2E8BC0")
        plt = lp.get_pyplot()
        plt.axhline(0.5, color='red')
        plt.title(plot_title, fontsize=18)
        plt.savefig(filenotes + "genome.pdf")


def compare_BAF_plot2(AD,DP, model_proba, mtx_barcodes_file,regions_file, chr_lst, plot_title="BAF_cluster_chr2+4(2-clusters)", cluster_note=["Clone1 (xx cells)","Clone2 (xx cells)"]):
    """
    Function:
    integration of BAF scatter plot for different clusters for specific regions and genome
    Example:
    compare_BAF_plot2(AD,DP, model_proba, mtx_barcodes_file,regions_file, ['2','4'],)
    
    """
    cluster_index_lst = np.argmax(model_proba, axis=1)
    uni_cluster = np.unique(cluster_index_lst)
    mtx_barcodes_lst = pd.read_table(mtx_barcodes_file,header = None)
    regions = np.genfromtxt(regions_file, delimiter="\t")
    cluster_cnt = 0
    for uni_id in uni_cluster:
        uni_cluster_barcodes = get_cluster_barcodes_lst(uni_id, cluster_index_lst, mtx_barcodes_lst)
        AD_cluster_tmp = sub_CellCluster(AD, mtx_barcodes_file, uni_cluster_barcodes)
        DP_cluster_tmp = sub_CellCluster(DP, mtx_barcodes_file, uni_cluster_barcodes)
        ## prepare data in dataframe for plot
        BAF_df = get_BAF_plot_data(regions_file, AD_cluster_tmp, DP_cluster_tmp)
        BAF_scatter_plot2(BAF_df, plot_title= plot_title + cluster_note[cluster_cnt], filenotes = str(uni_id))
        BAF_scatter_plot2(BAF_df, plot_title= plot_title + cluster_note[cluster_cnt], chr_lst=chr_lst, sub_chr_plot=True, filenotes = str(uni_id))
        cluster_cnt +=1

        
def node_BAF_plot2(AD,DP,mtx_barcodes_file,node, regions_file, chr_lst, plot_title1="BAF_cluster_chr2+4(node-xxx)", plot_title2="BAF_cluster_chr2+4(non_node-xxx)"):
    mkn45_10x_dloupe = "/storage/yhhuang/research/mito/mkn45/fulldepth/mode2/mkn45_5k_dloupe-group_10438-region_1_1-5120000_Y_56320001-59373566-heatmap_copy number.csv"
#     barcodes_10050_lst = get_node_barcodeslst(mkn45_10x_dloupe, 10050)
    barcodes_node_lst = get_node_barcodeslst(mkn45_10x_dloupe, node)
    barcodes_non_node_lst = get_non_node_barcodeslst(mkn45_10x_dloupe, node)

    AD_node = sub_CellCluster(AD, mtx_barcodes_file,barcodes_node_lst)
    DP_node = sub_CellCluster(DP, mtx_barcodes_file,barcodes_node_lst)
    BAF_df = get_BAF_plot_data(regions_file, AD_node, DP_node)
    BAF_scatter_plot2(BAF_df, plot_title= plot_title1, filenotes = "node10050")
    BAF_scatter_plot2(BAF_df, plot_title= plot_title1, chr_lst=chr_lst, sub_chr_plot=True, filenotes = "node10050")
    ## 
    AD_non_node = sub_CellCluster(AD, mtx_barcodes_file,barcodes_non_node_lst)
    DP_non_node = sub_CellCluster(DP, mtx_barcodes_file,barcodes_non_node_lst)
    BAF_df_ = get_BAF_plot_data(regions_file, AD_non_node, DP_non_node)
    BAF_scatter_plot2(BAF_df_, plot_title= plot_title2, filenotes = "non_node")
    BAF_scatter_plot2(BAF_df_, plot_title= plot_title2, chr_lst=chr_lst, sub_chr_plot=True, filenotes = "non_node")
    
     
def get_specific_pos():
    pass

def BAF_scatter_plot_for_specific_region():
    pass
    
## manhattan plot for BAF
### update from https://github.com/Rongtingting/limix-plot
### [orogin: https://github.com/limix/limix-plot]
###======================================== manhattan plot for BAF[start]
def manhattan(data, colora="#5689AC", colorb="#21334F", 
              anno_pv_max=None, pts_kws=None, ax=None):
    """
    Produce a manhattan plot.
    Parameters
    ----------
    data : DataFrame, dict
        DataFrame containing the chromosome, base-pair positions, and
        BAF-values.
    colora : matplotlib color
        Points color of the first group.
    colorb : matplotlib color
        Points color of the second group.
    anno_pv_max : float
        Threshold of maximum p value for annotating data['id'] in the plot.
    pts_kws : dict, optional
        Keyword arguments forwarded to the matplotlib function used for
        plotting the points.
    ax : matplotlib Axes, optional
        The target handle for this figure. If ``None``, the current axes is
        set.
    Example
    -------
    .. plot::
        >>> import limix_plot as lp
        >>> from numpy import log10
        >>>
        >>> df = lp.load_dataset('gwas')
        >>> df = df.rename(columns={"chr": "chrom"})
        >>> df = df.rename(columns={"pv": "BAF"})
        >>> print(df.head())
            chrom     pos       BAF
        234    10  224239  0.00887
        239    10  229681  0.00848
        253    10  240788  0.00721
        258    10  246933  0.00568
        266    10  255222  0.00593
        >>> lp.manhattan(df)
        >>> plt = lp.get_pyplot()
        >>> _ = plt.axhline(-log10(1e-7), color='red')
        >>> _ = plt.ylim(2, plt.ylim()[1])
    """
    from numpy import log10, unique, where
    from xarray import DataArray
    import pandas as pd

#     plt = get_pyplot()
    plt.figure(figsize=(24, 6))

    if isinstance(data, pd.DataFrame):
        data = DataArray(
            data=data["BAF"],
            dims=["candidate"],
            coords={k: ("candidate", data[k]) for k in data.columns},
        )
    else:
        data = DataArray(data=data)

    if len(data) == 0:
        raise ValueError("DataFrame is empty.")

    if pts_kws is None:
        pts_kws = dict()

    ax = plt.gca() if ax is None else ax

    data["chrom"] = data["chrom"].astype(str)
    data["pos"] = data["pos"].astype(int)
    chr_order = _chr_precedence(data)
    data["order"] = ("candidate", [chr_order[i] for i in data["chrom"].values])

    data = data.sortby(["order", "pos"])

    data = _abs_pos(data)

    if "markersize" not in pts_kws:
        pts_kws["markersize"] = 2
    if "marker" not in pts_kws:
        pts_kws["marker"] = "."
    if "linestyle" not in pts_kws:
        pts_kws["linestyle"] = ""

    colors = {0: colora, 1: colorb}

    for i, c in enumerate(unique(data["order"])):
        ok = data["order"] == c
        pts_kws["color"] = colors[i % 2]
        x = data.loc[ok]["abs_pos"]
        y = data.loc[ok].values
        ax.plot(x, y, **pts_kws)
        
        if anno_pv_max is not None:
            _idx = where(y > anno_pv_max)[0]
            for _ii in _idx:
                if 'id' in data.coords:
                    _txt = data['id'].loc[ok].loc[_ii].values
                else:
                    _txt = (str(data['chrom'].loc[ok].loc[_ii].values) + "_" + 
                            str(data['pos'].loc[ok].loc[_ii].values))
                ax.annotate(_txt, (x[_ii], y[_ii]), ha='center')

    ax.set_xlim(data["abs_pos"].min(), data["abs_pos"].max())
#     ax.set_ylim(0, ax.get_ylim()[1])
    ax.set_ylim(-0.01, 1.01)

    ax.set_ylabel("BAF", fontsize=15)
    ax.set_xlabel("chromosome", fontsize=15)

    u = unique(data["chrom"].values)
    chrom_labels = sorted(u, key=lambda x: chr_order[x])
    _set_ticks(ax, _chrom_bounds(data), chrom_labels)


def _plot_points(ax, data, alpha, null_style, alt_style):
    from numpy import log10

    null = data.loc[data.values >= alpha, :]
    alt = data.loc[data.values < alpha, :]

    ax.plot(null["abs_pos"], null.values, ".", ms=7, **null_style)
    ax.plot(alt["abs_pos"], alt.values, ".", ms=7, **alt_style)


def _set_ticks(ax, chrom_bounds, chrom_labels):
    from numpy import asarray, mean

    n = len(chrom_bounds)
    xticks = asarray([mean(chrom_bounds[i]) for i in range(n)])
    ax.set_xticks(xticks)
    ax.set_xticklabels(chrom_labels, fontsize=15)
    
    ax.tick_params(axis='both', which='major', labelsize=15)
#     ax.tick_params(axis='both', which='minor', labelsize=8)


def _abs_pos(data):
    from numpy import cumsum, flipud, unique

    order = unique(data["order"].values)
    chrom_ends = [data["pos"].values[data["order"].values == c].max() for c in order]

    offset = flipud(cumsum(chrom_ends)[:-1])

    data["abs_pos"] = data["pos"].copy()

    order = list(reversed(order))
    for i, oi in enumerate(offset):
        ix = data["order"] == order[i]
        data["abs_pos"].values[ix] = data.loc[ix]["abs_pos"] + oi

    return data


def _chrom_bounds(data):
    from numpy import unique

    order = unique(data["order"])
    v = []
    for c in order:
        vals = data["abs_pos"][data["order"] == c]
        v += [(vals.min(), vals.max())]
    return v


def _isint(i):
    try:
        int(i)
    except ValueError:
        return False
    else:
        return True


def _chr_precedence(data):
    from numpy import unique

    uchr = unique(data["chrom"].values)
    nchr = [int(i) for i in uchr if _isint(i)]
    if len(nchr) > 0:
        offset = max(nchr)
    else:
        offset = -1
    precedence = {str(i): i for i in nchr}
    schr = sorted([i for i in uchr if not _isint(i)])
    for i, s in enumerate(schr):
        precedence[s] = offset + i + 1
    return precedence

###======================================== manhattan plot for BAF[end]