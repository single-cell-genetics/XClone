"""Base functions for XClone RDR processing
"""

# Author: Rongting Huang
# Date: 2021-09-20
# update: 2022-07-12

import numpy as np
import pandas as pd
import anndata as ad
import scipy as sp

import datetime
import multiprocessing

from scipy import stats
from scipy.special import logsumexp

from ..model.base_utils import normalize
from ..plot._data import reorder_data_by_cellanno

# About different filtering functions in data processing.
# 1. filtering genes based on cell detection rate-FUNC `gene_filter``
# 2. filtering genes based on ref cells-for raW ratio based following analysis- FUNC `Xdata_RDR_preprocess`
# 3. filtering nan value in any mtx for genes.(like nan in emm_prob) to solve num issue

## Part-I data preprocessing part in RDR module
# YH comment: combine the multiple output adata objects into one
def Xdata_RDR_preprocess(Xdata, 
                         filter_ref_ave = 0,
                         cell_anno_key = "cell_type", 
                         ref_celltype = None,
                         var_key = "ref_avg",
                         obs_key = "counts_ratio",
                         mode = "ALL"):
    """
    Xdata_RDR_preprocess
    actually can replace V1_rdr_Xdata, V2_rdr_Xdata 20211102
    
    Function:
    ---------
    1) Filter genes based on ref_ave value; 
       at least filter 0 ref_ave values to avoid numeric issues;
    2) Transform original matrix into celltypebased matrix 
    and cellbased matrix with reference in the first row;
    Output the anndata format for visualization.
    3) calculate 'counts ratio'
    
    Params:
    --------
    Xdata: anndata.
    filter_ref_ave: int. if filter_ref_ave is None, skip the filtering step.
    ref_celltype: char. assigned reference celltype.
    cell_anno_key: char. colname in Xdata.var for cell annotation.
    mode: "FILTER", "ALL". default: "ALL" for whole preprocessing pipeline.
          if mode is FILTER, then return the filtered anndata;
          if mode is default, then return anndata and other for following analysis.
    
    Return:
    ------
    update_Xdata: filter the ref values and update the Xdata;
    ref_obs_bulk_ad: first row is ref_bulk, 
                    celltype based anndata, bulk [ref + obs] anndata;
    ref_obs_ad: first row is ref_mean, then add the obs cell anndata;
    ref_obs_ad1: first row is ref_mean, another version keep the original ref cells info.
    rr_ad_celltype: log raw_ratio 
    rr_ad_cell: log raw_ratio

    Examples
    --------
    >>> update_Xdata = xclone.model.Xdata_RDR_preprocess(test_adata, filter_ref_ave = 2, 
        ref_celltype = "unclassified", mode = "FILTER")
    >>> update_Xdata = xclone.model.Xdata_RDR_preprocess(test_adata, filter_ref_ave = None, 
        ref_celltype = "unclassified")
    >>> update_Xdata, ref_obs_bulk_ad, ref_obs_ad, ref_obs_ad1, rr_ad_celltype, rr_ad_cell = 
        xclone.model.Xdata_preprocess(test_adata, filter_ref_ave = 0, ref_celltype ="unclassified", 
        cell_anno_key = "cell_type", order_cell = True)
    """
    # check sparse matrix
    if sp.sparse.issparse(Xdata.X):
        pass
    else:
        raise ValueError("Xdata should be sparse matrix! Pls check!")
    
    # if mode == "Counts_ratio":
    #     ## cell counts_ratio to be compared with learned libratio
    #     Xdata.obs[obs_key] = Xdata.X.A.sum(axis=1) / Xdata.var['ref_avg'].sum()
    #     return Xdata
    
    # reference data preprocessing
    if ref_celltype is None:
        raise ValueError("ref_celltype should be assigned! Pls check!")

    ref_flag = Xdata.obs[cell_anno_key] == ref_celltype
    ref_Xdata = Xdata[ref_flag,:]

    Xdata.var['ref_avg'] = ref_Xdata.X.A.mean(axis=0)

    
    ## filter genes based on ref_average counts
    ### if filter_ref_ave is None, skip the filtering step*
    ### if mode == "FILTER", just do the filter step*
    gene_num = Xdata.var.shape[0]
    if filter_ref_ave is None:
        gene_flag = np.ones((gene_num), dtype=bool)
    else:
        gene_flag = Xdata.var['ref_avg'] > filter_ref_ave
    
    # mode="FILTER"
    update_Xdata = Xdata[:, gene_flag].copy()
    ## cell counts_ratio to be compared with learned libratio
    update_Xdata.obs[obs_key] = update_Xdata.X.A.sum(axis=1) / update_Xdata.var['ref_avg'].sum()
    
    if mode == "FILTER":
        
        print("[XClone-RDR preprocessing] Filter out %d genes / %d total genes, remain %d genes" 
               %(gene_num - gene_flag.sum(), gene_num, gene_flag.sum()))
        return update_Xdata
    
    # mode="ALL"
    ## process both cellbased data and celltype based data.
    celltype_ad = rr_ad_celltype_processing(update_Xdata, ref_celltype, cell_anno_key)
    is_ref = celltype_ad.obs[cell_anno_key] == ref_celltype

    celltype_ad.var["ref_avg"] = celltype_ad[is_ref,:].X.toarray()[0]

    rr_cell_ad = rr_ad_cell_processing(update_Xdata, ref_celltype, cell_anno_key)
    update_Xdata.layers["raw_ratio"] = rr_cell_ad.X

    return update_Xdata, celltype_ad

### Sub part-For visualization -RDR raw ratio

def rr_ad_celltype_processing(Xdata, ref_celltype, cell_anno_key):
    """
    For pseudo bulk. celltype-based count.
    version:0.0.2
    """
    ref_flag = Xdata.obs[cell_anno_key] == ref_celltype
    ref_Xdata = Xdata[ref_flag,:]
    obs_Xdata = Xdata.copy()

    ref_bulk = ref_Xdata.X.A.sum(axis=0)
    ref_normalization_term = ref_bulk.sum()
    ref_norm = ref_bulk/ref_normalization_term

    # 01-celltype based anndata
    ## prepare the obs_bulk
    groups = obs_Xdata.obs.groupby(cell_anno_key).indices ## dict
    cell_type_order = []
    
    cnt = 0
    for tmp_celltype, inds in groups.items():
        cell_type_order.append(tmp_celltype)
        if cnt==0:
            obs_all_bulk = obs_Xdata[inds,:].X.A.sum(axis=0)
        else:
            obs_all_bulk = np.vstack((obs_all_bulk, obs_Xdata[inds,:].X.A.sum(axis=0)))
        cnt+=1
    
    celltype_obs = pd.DataFrame({cell_anno_key:cell_type_order})
    obs_all_bulk_ad = ad.AnnData(obs_all_bulk, var=obs_Xdata.var.copy(), obs = celltype_obs)

    ## calculate ratio-bulk
    ## return rr_ad_celltype for log raw ratio visualization and choosing normal chrs
    cell_lib = obs_all_bulk.sum(axis=1, keepdims=True) # equivalent normalization term like ref_bulk
    # raw_ratio = np.log2((obs_all_bulk+0.1) / (cell_lib*ref_norm)) ## add small value 0.1 for visualization 
    raw_ratio = np.log2((obs_all_bulk + 1e-8) / (cell_lib*ref_norm)) ## add small value 1e-6 for visualization

    obs_all_bulk_ad.layers["raw_ratio"] = raw_ratio
    return obs_all_bulk_ad

def rr_ad_cell_processing(Xdata, ref_celltype, cell_anno_key):
    """
    version:0.0.2
    """
    ref_flag = Xdata.obs[cell_anno_key] == ref_celltype
    ref_Xdata = Xdata[ref_flag,:]

    obs_Xdata = Xdata
    
    ref_bulk = ref_Xdata.X.A.sum(axis=0)
    ref_normalization_term = ref_bulk.sum()
    ref_norm = ref_bulk/ref_normalization_term
    
    # 02-cell based anndata
    ## calculate ratio-cell
    cell_lib = obs_Xdata.X.A.sum(axis=1, keepdims=True)
    # raw_ratio = np.log((obs_Xdata.X.A + 0.1) / (cell_lib*ref_norm)) ## add small value 0.1 for visualization
    raw_ratio = np.log((obs_Xdata.X.A + 1e-8) / (cell_lib*ref_norm)) ## add small value 1e-6 for visualization

    rr_cell_ad = ad.AnnData(raw_ratio, var=obs_Xdata.var.copy(), obs = obs_Xdata.obs.copy())
    
    print("output anndata is not sparse matrix.")
    return rr_cell_ad


## Part-II libsize ratio estimation part in RDR module
## see _RDR_libratio

# strategy notes(pipeline here):
# 1) preprocess the Xdata here for following analysis.
# 2) select 4 normal chrs based on the celltype bulk Xdata.
# 3) calculate libratio using GLM for celltype/cell(use the extended selected chrs).
# 4) get emm_prob from RDR VB model and put back to HMM framework.


## Part-III gene-specific Dispersions estimation part in RDR module
## see _RDR_dispersion

# strategy notes
# 1. view celltype information to find which celltype have the most cells
# 2. select the celltype with most cell nums
# 3. and estimate gene-specific dispersions based on the selected celltype(do not perform 
# this across wholedataset to avoid counting in the biological difference) 


### subpart-todo
# actually, inf or nan in learned libratio or dispersion; 
# need do some filtering before the following analysis.

# and need select rational estimated libsize/common total counts
# need select rational estimated gene-specific dispersions


## Part-IV Estimate CNV states' ratio based via different gene expression group

## see _RDR_CNVratio

## todo CNV estimation iteration


## Part-V reference biological ratio prediction based on NMF

## see _RDR_NMF

## Part-VI merging each part in RDR module|pipeline construction|HMM

# two main strategies
# 1) calculate NB prob as emm_prob and integrate into HMM model
# 2) intergrate RDR CNV_prob (as emm_prob) into HMM model


# pipeline in RDR module

# 0.preprocessing
# 1.libsize ratio (done)
# 2.gene-specific overdispersion (done)
# 3.emm_prob prepare
# 4.HMM
# 5.estimate CNV ratio based on the step3-4 iterations, cal logliklihood.
# 6.final visualization

def gene_length_scale(Xdata, avg_key = "ref_avg", scale_avg_key = "ref_avg_scale"):
    """
    For smart-seq counts: take gene length into account.
    """

    Xdata.var["gene_length"] = Xdata.var["stop"] - Xdata.var["start"]

    Xdata.var["length_ratio"] = Xdata.var["gene_length"]/Xdata.var["gene_length"].median()

    Xdata.var[scale_avg_key] = Xdata.var[avg_key] * Xdata.var["length_ratio"]
    
    print("[Xclone preprocessing for smart-seq counts]: ", scale_avg_key ,"added.")


    # ## cell counts_ratio to be compared with learned libratio
    # Xdata.obs['counts_ratio'] = Xdata.X.A.sum(axis=1) / Xdata.var['ref_avg'].sum()

    return Xdata

