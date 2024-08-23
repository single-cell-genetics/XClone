"""Base functions for XClone spatial transcriptomics CNV analysis.
"""
import numpy as np

def CalculateSampleCNVProb(Xdata, brk = "chr_arm", Xlayer = "prob1_merge"):
    """
    For combine data;
    For rdr data;
    For BAF data
    """
    brk_item = Xdata.var[brk].drop_duplicates(keep="first")
    obs_key_lst = []
    
    for brk_ in brk_item:
    
        tmp_region_flag = Xdata.var[brk] == brk_
        tmp_p_cnv = Xdata[:, tmp_region_flag].layers[Xlayer][:,:,[0,1,3]].sum(axis = 2).mean(axis = 1)
        obs_key = "p_cnv_" + brk_
        print(obs_key)
        obs_key_lst.append(obs_key)
        Xdata.obs[obs_key] = tmp_p_cnv
        Xdata.uns["p_cnv_lst"] = obs_key_lst
    
    return Xdata
 
def CalculateSampleBAF(Xdata, brk = "chr_arm"):
    """
    For combine data;
    For rdr data;
    For BAF data
    """
    brk_item = Xdata.var[brk].drop_duplicates(keep="first")
    BAF_lst = []
    BAF_smooth_lst = []
    
    for brk_ in brk_item:
    
        tmp_region_flag = Xdata.var[brk] == brk_
        Xlayer = "fill_BAF_phased"
        tmp_BAF = Xdata[:, tmp_region_flag].layers[Xlayer].mean(axis=1)
        obs_key = "BAF_" + brk_
        BAF_lst.append(obs_key)
        Xdata.obs[obs_key] = tmp_BAF
        Xdata.uns["BAF_lst"] = BAF_lst
        
        
        Xlayer = "BAF_phased_KNN_WMA"
        tmp_smooth_BAF = Xdata[:, tmp_region_flag].layers[Xlayer].mean(axis=1)
        obs_key = "BAF_smooth_" + brk_
        BAF_smooth_lst.append(obs_key)
        Xdata.obs[obs_key] = tmp_smooth_BAF
        Xdata.uns["BAF_smooth_lst"] = BAF_smooth_lst
    
    return Xdata

def CalculateSampleRDR(Xdata, brk = "chr_arm"):
    """
    For rdr data;
    """
    brk_item = Xdata.var[brk].drop_duplicates(keep="first")
    RDR_smooth_lst = []
    RDR_smooth_log_lst = []
    
    for brk_ in brk_item:
    
        tmp_region_flag = Xdata.var[brk] == brk_
        Xlayer = "RDR_smooth"
        tmp_smooth_RDR = Xdata[:, tmp_region_flag].layers[Xlayer].mean(axis=1)
        obs_key = "RDR_smooth_" + brk_ + " (log)"
        RDR_smooth_log_lst.append(obs_key)
        Xdata.obs[obs_key] = tmp_smooth_RDR
        Xdata.uns["RDRsmooth_log_lst"] = RDR_smooth_log_lst
        
        
        tmp_smooth_RDR = np.exp(Xdata[:, tmp_region_flag].layers[Xlayer]).mean(axis=1)
        obs_key = "RDR_smooth_" + brk_
        RDR_smooth_lst.append(obs_key)
        Xdata.obs[obs_key] = tmp_smooth_RDR
        Xdata.uns["RDR_smooth_lst"] = RDR_smooth_lst
    
    return Xdata


def spatial_analysis(visium_adata, xclone_adata, baf_adata, rdr_adata=None):
    """_summary_
    Sometimes only include combined cnv analysis and baf analysis.

    Args:
        visium_adata (_type_): _description_
        xclone_adata (_type_): _description_
        baf_adata (_type_): _description_
        rdr_adata (_type_): _description_
    """

    adata_sp = visium_adata[visium_adata.obs.index.isin(xclone_adata.obs.index),:]
    # cnv
    adata_sp.obs[xclone_adata.uns["p_cnv_lst"]] = xclone_adata.obs[xclone_adata.uns["p_cnv_lst"]]
    adata_sp.uns["p_cnv_lst"] = xclone_adata.uns["p_cnv_lst"]
    # baf
    adata_sp.obs[baf_adata.uns["BAF_lst"]] = baf_adata.obs[baf_adata.uns["BAF_lst"]]
    adata_sp.uns["BAF_lst"] = baf_adata.uns["BAF_lst"]
    adata_sp.obs[baf_adata.uns["BAF_smooth_lst"]] = baf_adata.obs[baf_adata.uns["BAF_smooth_lst"]]
    adata_sp.uns["BAF_smooth_lst"] = baf_adata.uns["BAF_smooth_lst"]

    # rdr
    if rdr_adata is not None:
        adata_sp.obs[rdr_adata.uns["RDR_smooth_lst"]] = rdr_adata.obs[rdr_adata.uns["RDR_smooth_lst"]]
        adata_sp.uns["RDR_smooth_lst"] = rdr_adata.uns["RDR_smooth_lst"]

        adata_sp.obs[rdr_adata.uns["RDRsmooth_log_lst"]] = rdr_adata.obs[rdr_adata.uns["RDRsmooth_log_lst"]]
        adata_sp.uns["RDRsmooth_log_lst"] = rdr_adata.uns["RDRsmooth_log_lst"]

    return adata_sp