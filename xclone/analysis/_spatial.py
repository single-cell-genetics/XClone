"""Base functions for XClone spatial transcriptomics CNV analysis.
"""

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
        obs_key = "BAF_smooth" + brk_
        BAF_smooth_lst.append(obs_key)
        Xdata.obs[obs_key] = tmp_smooth_BAF
        Xdata.uns["BAF_smooth_lst"] = BAF_smooth_lst
    
    return Xdata