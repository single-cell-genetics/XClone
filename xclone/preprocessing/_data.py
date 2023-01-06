"""Base functions for XClone data preprocessing
"""
# Author: Rongting Huang
# Date: 2021/07/15
# update: 2022/07/13

import os
import time
import numpy as np
import pandas as pd
import scipy as sp
from scipy import io
# from scipy.io import mmread
from scipy.sparse import vstack
from anndata import AnnData

from ._anno_data import load_anno


## Part I: mtx preprocessing

def get_Xmtx(X, genome_mode):
    """
    Function: prepare data format for XClone
    ------
    params:
    X: csr_mtx /csr_mtx path
    genome_mode: hg38_genes/hg38_blocks/hg19_genes/hg19_blocks
    default: hg38_genes
    ------
    usage:
    from xclone.model.preprocessing_utils import get_Xmtx
    dat_dir = '/storage/yhhuang/users/rthuang/processed_data/xcltk/
    xianjie-cpos/kat_022621/TNBC1-csp-post/phase-snp-even/'
    X = dat_dir + 'kat-csp-post.50kb.block.AD.mtx'
    Xmtx = get_Xmtx(X, "hg38_blocks")
    """
    # X can be file path/or loaded sparse matrix
    if sp.sparse.issparse(X):
        X_data = X   
    elif os.path.exists(X):
        X_data = sp.io.mmread(X).tocsr()
    ## use chr1-22+XY, 
    ## can be updated if xcltk output change[Note1]
    if genome_mode=="hg19_genes":
        X_data = X_data[0:32696,:]
    if genome_mode=="hg38_genes":
        X_data = X_data[0:33472,:]
    return X_data.T

def resort_mtx_bychr(mtx_file, features_file, assign_sort_index=None, out_file = None, keep_all = True):
    """
    Function:
    sort the mtx as the assignned chromosome index.
    For thr output gene based mtx is not in order.
    
    out_file: "resort_matrix.mtx"
    
    Example1:
    import xclone
    mtx_file = rdr_data_dir  + 'matrix.mtx'
    features_file = rdr_data_dir + "features.tsv"
    out_RDR_file = rdr_data_dir + "resort_matrix.mtx"
    resort_mtx = xclone.pp.resort_mtx_bychr(mtx_file, features_file, out_file = out_RDR_file, keep_all=True)
    
    Example2:
    --------
    AD_file = baf_data_dir  + 'xcltk.AD.mtx'
    DP_file = baf_data_dir  + 'xcltk.DP.mtx'
    features_file = baf_data_dir  +  "xcltk.region.tsv"
    out_AD_file =  baf_data_dir  + "resort_AD.mtx"
    out_DP_file = baf_data_dir  +  "resort_DP.mtx"

    resort_AD = xclone.pp.resort_mtx_bychr(AD_file, features_file, out_file = out_AD_file, keep_all=True)
    resort_DP = xclone.pp.resort_mtx_bychr(DP_file, features_file, out_file = out_DP_file, keep_all=True)
    """
    sort_index=['1', '2','3', '4', '5', '6', '7', '8', '9','10', '11', '12', '13', '14', '15', '16', '17', '18', '19','20', '21', '22', 'X', 'Y']
    
    mtx_init  = sp.io.mmread(mtx_file).tocsr()
    print("mtx init:", mtx_init.shape)
    regions = pd.read_table(features_file,header = None)
    try:
        regions.columns = ["chr_info", "start", "end"]
        chr_init = regions["chr_info"].str.split("chr").str[1]
    except:
        regions.columns = ["chr_info", "start", "end", "all_info"]
        chr_init = regions["chr_info"]
    else:
        pass
    # chr_init = regions["chr_info"].str.split("chr").str[1]
    regions["chr_init"] = chr_init
    print("original order: \n", regions.drop_duplicates("chr_init")["chr_init"])
    # drop_duplicates default keep_first
    
    if assign_sort_index == None:
        assign_sort_index = sort_index
        print("default chr index for ordering: 1-22,XY")
    else:
        print(assign_sort_index)
    
    cnt = 0
    for i in assign_sort_index:
        print(i)
        tmp_mtx = mtx_init[chr_init == i]
        if cnt==0:
            resort_mtx = tmp_mtx
        else:
            resort_mtx = vstack([resort_mtx,tmp_mtx])
        cnt+=1
    if cnt == len(assign_sort_index):
        print("sorted %d chromosomes" %cnt)
    else:
        print("Wrong match with assign_sort_index")
    
    print("sorted part:", resort_mtx.shape)
    
    if keep_all:
        flag = chr_init.isin(assign_sort_index)
        tmp_mtx = mtx_init[~flag]
        resort_mtx = vstack([resort_mtx,tmp_mtx])

    if out_file == None:
        print("please assign a file path for output")
    else:
        sp.io.mmwrite(out_file, resort_mtx)
    
    print("output:", resort_mtx.shape)
    return resort_mtx


def process_barcodes(barcodes_infile, barcodes_outfile):
    """
    special for BCH869 dataset-smart-seq
    Tools for merging the bamfiles keep the .sort
    """
    barcodes = pd.read_table(barcodes_infile, header = None, index_col=0)
    
    barcodes.index = barcodes.index.str.split(".").str[0]
    barcodes.to_csv(barcodes_outfile, sep='\t', index=True, header=False)

def check_RDR_BAF_anno(BAF_barcodes_file, RDR_barcodes_file):
    """
    """
    BAF_barcodes = pd.read_table(BAF_barcodes_file, header = None, index_col=0)
    RDR_barcodes = pd.read_table(RDR_barcodes_file, header = None, index_col=0)
    
    BAF_cell_num = BAF_barcodes.shape[0]
    RDR_cell_num = RDR_barcodes.shape[0]
    
    if BAF_cell_num == RDR_cell_num:
        if (BAF_barcodes.index == RDR_barcodes.index).sum() == BAF_cell_num:
            print("[XClone cell anno checking]: RDR and BAF cell anno in the same order.")
            success_flag = True
        else:
            print("[XClone cell anno checking]: RDR and BAF cell anno not matched. pls check!")
            print("""[XClone data hint]: Pls resort RDR mtx in the same order with BAF barcodes 
                by FUNC 'xclone.pp.resort_mtx_bycell'.""")
            success_flag = False
    else:
        print("[XClone cell anno checking]: RDR and BAF cells number not matched! Pls check!")
        print("""[XClone data hint]: Pls resort RDR mtx in the same order with BAF barcodes 
                by FUNC 'xclone.pp.resort_mtx_bycell'.""")
        success_flag = False
    return success_flag


def resort_mtx_bycell(BAF_barcodes_file, RDR_barcodes_file, RDR_mtx_file, out_mtx_file):
    """
    resort RDR mtx by BAF_barcodes order.
    """
    BAF_barcodes = pd.read_table(BAF_barcodes_file, header = None, index_col=0)
    RDR_barcodes = pd.read_table(RDR_barcodes_file, header = None, index_col=0)

    RDR_mtx_init  = sp.io.mmread(RDR_mtx_file).tocsr()
    print("mtx init:", RDR_mtx_init.shape)

    RDR_mtx_df = pd.DataFrame(RDR_mtx_init.A.T)
    RDR_mtx_df["cell_barcodes"] = RDR_barcodes.index
    RDR_mtx_df.set_index('cell_barcodes', inplace=True)

    update_mtx = BAF_barcodes.merge(RDR_mtx_df, left_index=True, right_index=True, how = "left")
    
    out_mtx = update_mtx.to_numpy().T.tocsr()

    if out_mtx_file == None:
        print("[XClone warning] Please assign a file path for output.")
    else:
        sp.io.mmwrite(out_mtx_file, out_mtx)
    print("[XClone warning] Need update RDR barcodes file after mtx resorting.")
    
    return out_mtx

## Part II: XClone data format-AnnData
def xclonedata(X, data_mode, 
               mtx_barcodes_file, 
               regions_anno_file = None, 
               genome_mode = "hg38_genes", 
               data_notes = None):
    """
    Function: 
    ------
    prepare data format for XClone
    
    Params:
    -------
    X: csr_mtx/csr_mtx path, can be list
    data_mode: 'BAF' OR 'RDR'
    mtx_barcodes_file: barcodes_file path
    genome_mode: hg38_genes/hg38_blocks/hg19_genes/hg19_blocks
    default: hg38_genes
    
    Example:
    --------
    * Load xclonedata
    >>> import xclone
    >>> dat_dir = '/storage/yhhuang/users/rthuang/processed_data/xcltk/
    xianjie-cpos/kat_022621/TNBC1-csp-post/phase-snp-even/'
    >>> AD_file = dat_dir + 'kat-csp-post.50kb.block.AD.mtx'
    >>> DP_file = dat_dir + 'kat-csp-post.50kb.block.DP.mtx'
    >>> mtx_barcodes_file = dat_dir + "cellSNP.samples.tsv"
    >>> BAF_adata = xclone.pp.xclonedata([AD_file, DP_file], 'BAF', mtx_barcodes_file, 
    "hg38_blocks", "TNBC1 scRNA-seq data in copyKAT")
    >>> RDR_adata = xclonedata(RDR_file, 'RDR', mtx_barcodes_file, "hg38_genes")
    """
    ## data loading
    ### obs anno-mtx cell barcoedes
    if isinstance(mtx_barcodes_file, pd.DataFrame):
        cell_anno = mtx_barcodes_file
    elif os.path.exists(mtx_barcodes_file):
        cell_anno = pd.read_table(mtx_barcodes_file, header = None, index_col=0)
    cell_anno.index.name = None
    ### var anno
    if regions_anno_file is None:
        regions_anno = load_anno(genome_mode)
    else:
        regions_anno = pd.read_table(regions_anno_file, header = None, index_col=0)

    ## initialize the data in AnnData format
    if data_mode == 'BAF':
        AD = get_Xmtx(X[0], genome_mode)
        DP = get_Xmtx(X[1], genome_mode)
        X_adata = AnnData(AD, obs=cell_anno, var=regions_anno) # dtype='int32'
        X_adata.layers["AD"] = AD
        X_adata.layers["DP"] = DP  
    elif data_mode =='RDR':
        RDR = get_Xmtx(X, genome_mode)
        X_adata = AnnData(RDR, obs=cell_anno, var=regions_anno) # dtype='int32'
        X_adata.layers["raw_expr"] = RDR
    
    ## unstructed anno
    X_adata.uns["log"] = dict([('init_data', str(X_adata.shape))])
    X_adata.uns["log"]["data_mode"] = data_mode
    X_adata.uns["log"]["data_notes"] = data_notes
    X_adata.uns["log"]["genome_mode"] = genome_mode
    X_adata.uns["data_mode"] = data_mode
    X_adata.uns["genome_mode"] = genome_mode
    if data_notes is None:
        data_notes = time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())
    else:
        data_notes = data_notes + ": " + time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())
    X_adata.uns["data_notes"] = data_notes
    return X_adata

def extra_anno(Xdata, anno_file, 
               barcodes_key = "cell",
               cell_anno_key = "cell_type", 
               sep = ",", 
               copy=True):
    """
    Func:
    add extra annotation for Xdata.obs
    Params:
    --------
    anno_file: first column is cell barcodes, default col_name is "cell"; other
            cols contains different type of annotations, like celltypes or clones.
    cell_anno_key: default "cell_type", also can be list ["cell_type", "Clone ID"]

    Example:
    >>> import xclone
    >>> anno_dir = "/storage/yhhuang/research/xomics/GEX5/GX109-T1c/"
    >>> anno_file = anno_dir + "GX109-T1c_yh.csv"
    >>> RDR_adata = xclone.pp.extra_anno(RDR_adata, anno_file, barcodes_key = "cell", 
    cell_anno_key = "cell_type", sep =",") 

    """
    Xdata = Xdata.copy() if copy else Xdata
    # annotation data loading
    anno_data = pd.read_csv(anno_file, sep = sep)
    # merge annotaion
    anno_data_update = Xdata.obs.merge(anno_data.set_index(barcodes_key), left_index=True, right_index=True, how = "left")
    
    if isinstance(cell_anno_key, list):
        for key_ in cell_anno_key:
            anno_data_update[key_] = anno_data_update[key_].astype('category')
    else:
        anno_data_update[cell_anno_key] = anno_data_update[cell_anno_key].astype('category')
    
    Xdata.obs = anno_data_update
    return Xdata if copy else None

## Part III: check data validity
from ._preprocessing import data_check, valid_cell, gene_filter

def check_RDR(Xdata, cell_anno_key = "cell_type", 
              cell_detection_rate = 0.05, filter_na_anno = True, verbose = True):
    """
    After data loading, check data quality.
    """
    data_check(Xdata, Xlayer = "raw_expr")
    print("[XClone data preprocessing] check RDR raw dataset value: success")
    
    update_Xdata = valid_cell(Xdata, cell_anno_key, filter_na_anno, verbose = verbose)
    print("[XClone data preprocessing] check RDR cell annotation: success")
    
    update_Xdata = gene_filter(update_Xdata, cell_detection_rate, verbose = verbose)
    print("[XClone data preprocessing] detect RDR genes: done")

    return update_Xdata

def check_BAF(Xdata, cell_anno_key = "cell_type", filter_na_anno = True, verbose = True):
    """
    After data loading, check data quality.
    """
    data_check(Xdata, Xlayer = "AD")
    data_check(Xdata, Xlayer = "DP")
    print("[XClone data preprocessing] check BAF raw dataset value: success")
    
    update_Xdata = valid_cell(Xdata, cell_anno_key, filter_na_anno, verbose = verbose)
    print("[XClone data preprocessing] check BAF cell annotation: success")
    return update_Xdata

def check_RDR_BAF_cellorder(RDR_Xdata, BAF_Xdata):
    """
    * need to make sure the BAF and RDR are in the same cell order.

    Function:
    --------
    preprocessing step.
    check the annotation for both RDR and BAF data.
    actually should be the same as the anno_file and `cell_index` order are the same.
    get the same cell shape for RDR and BAF, and in the same order.
    """
    RDR_cellnum = RDR_Xdata.obs.shape[0]
    BAF_cellnum = BAF_Xdata.obs.shape[0]
    if RDR_cellnum == BAF_cellnum:
        check_ordernum = (RDR_Xdata.obs.index == BAF_Xdata.obs.index).sum()
        if RDR_cellnum == check_ordernum:
            print("[XClone data checking]: RDR and BAF in same cell order.")
            success_flag = True
        else:
            print("[XClone data checking]: RDR and BAF in different cell order! Pls check!")
            success_flag = False
    return success_flag


## Part IV: check data before combination analysis.
def check_data_combine(RDR_Xdata, BAF_Xdata):
    """
    Function:
    after-processing step.

    check cell order
    check var
    check probability states, etc.

    Example:
    
    """
    ## check cell order
    success_flag = check_RDR_BAF_cellorder(RDR_Xdata, BAF_Xdata)
    ## check var
    if success_flag:
        ## check probability
        if "posterior_mtx" in RDR_Xdata.layers.keys():
            pass
        else:
            raise ValueError("[XClone data warning] No layer 'posterior_mtx' exists in RDR module.")
        if "posterior_mtx" in BAF_Xdata.layers.keys():
            pass
        else:
            raise ValueError("[XClone data warning] No layer 'posterior_mtx' exists in BAF module.")
    else:
        raise ValueError("[XClone data warning]-pls check cell order before combination step.")

def exclude_XY_adata(Xdata):
    """
    exclude XY in chromosome analysis.
    """
    flag_ = ~(Xdata.var["chr"].isin(["X", "Y"]))
    print("[XClone] exclude chr X&Y in the analysis.")
    return Xdata[:, flag_].copy()

def check_RDR_BAF_chrmapping(RDR_Xdata, BAF_Xdata):
    """
    mainly for chr mapping check.
    """
    rdr_chr = RDR_Xdata.var["chr"].drop_duplicates().reset_index(drop = True)
    baf_chr = BAF_Xdata.var["chr"].drop_duplicates().reset_index(drop = True)
    
    success_flag = rdr_chr.equals(baf_chr)
    return success_flag