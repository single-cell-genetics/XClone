"""Base functions for XClone BAF processing
"""

# Author: Rongting Huang
# Date: 2021/10/18
# update: 2022/04/08

import numpy as np
import pandas as pd
import anndata as ad
import multiprocessing
import datetime

import gc

from .phasing import Local_Phasing, Global_Phasing
# from .analysis_utils import select_chr_region

## Functions for data processing
def select_chr_region(Xdata, region_loc, region_key="chr"):
    """
    Func:
    select region from whole genome.
    """
    FLAG = Xdata.var[region_key] == region_loc
    return Xdata[:,FLAG]

## Functions for phasing
### local phasing
def process_bin(idx, AD, DP):
    """
    Func:
    bin: default for 100 genes
    """ 
    ## Z: allele flipping probability
    ad_sum, ad_sum1, dp_sum, Z, thetas, _logLik_new = Local_Phasing(AD, DP)
    is_flip = np.argmax(Z, axis =1)
    
    ## important *
    BD = DP-AD
    AD_phased = AD.toarray() * np.expand_dims(1-is_flip, axis=1) + BD.toarray() * np.expand_dims(is_flip, axis=1)

    RV_bin = {}
    RV_bin["bin_idx"] = idx
    RV_bin['bin_idx_lst'] = np.repeat(idx, AD_phased.shape[0])
    RV_bin['Z'] = Z
    RV_bin['flip'] = is_flip
    RV_bin["AD_phased"] = AD_phased
    RV_bin['ad_bin_softcnt'] = ad_sum[0, :]
    RV_bin['ad_bin'] = ad_sum1[0, :]
    RV_bin['dp_bin'] = dp_sum[0, :]
    RV_bin['theta_bin'] = np.array(thetas)[:, 0]
    RV_bin['logLik'] = _logLik_new
    return RV_bin

def process_region(AD_region, DP_region, phasing_len = 100, nproc=1, last_bin_min = 30):
    """
    Func:
    region: default chr_based. 
    """
    n_bins = int(AD_region.shape[0] / phasing_len)
    last_bin_len = AD_region.shape[0] % phasing_len

    # merge last bin if too short
    merge_last_bin = last_bin_len < last_bin_min and last_bin_len > 0
    

    ## do phasing
    if nproc > 1:
        result = []
        pool = multiprocessing.Pool(processes=nproc)

        for ib in range(n_bins):
            if merge_last_bin and ib == n_bins - 1:
                idx = range(ib * phasing_len, AD_region.shape[0])
            else:
                idx = range(ib * phasing_len, (ib + 1) * phasing_len)
            result.append(pool.apply_async(process_bin, (ib, AD_region[idx, :], DP_region[idx, :]), callback=None))

        if not merge_last_bin and last_bin_len > 0:
            idx = range(-last_bin_len, 0)
            result.append(pool.apply_async(process_bin, (n_bins, AD_region[idx, :], DP_region[idx, :]), callback=None))

        pool.close()
        pool.join()
        result = [res.get() for res in result]

    else:
        result = []
        for ib in range(n_bins):
            if merge_last_bin and ib == n_bins - 1:
                idx = range(ib * phasing_len, AD_region.shape[0])
            else:
                idx = range(ib * phasing_len, (ib + 1) * phasing_len)
            RV_bin = process_bin(ib, AD_region[idx, :], DP_region[idx, :])
            result.append(RV_bin)

        if not merge_last_bin and last_bin_len > 0:
            idx = range(-last_bin_len, 0)
            RV_bin = process_bin(n_bins, AD_region[idx, :], DP_region[idx, :])
            result.append(RV_bin)
      
    
    ## resolve result
    for i, RV_bin in zip(range(len(result)), result):
        if i == 0:
            AD_phased = RV_bin["AD_phased"]
            ad_bin_softcnt = RV_bin["ad_bin_softcnt"]
            ad_bin = RV_bin["ad_bin"]
            dp_bin = RV_bin["dp_bin"]
            theta_bin = RV_bin["theta_bin"]
            bin_idx = RV_bin["bin_idx"]
            bin_idx_lst = RV_bin["bin_idx_lst"]
            allele_flip_local = RV_bin["flip"]
        else:
            AD_phased = np.vstack((AD_phased, RV_bin["AD_phased"]))
            ad_bin_softcnt = np.vstack((ad_bin_softcnt, RV_bin["ad_bin_softcnt"]))
            ad_bin = np.vstack((ad_bin, RV_bin["ad_bin"]))
            dp_bin = np.vstack((dp_bin, RV_bin["dp_bin"]))
            theta_bin = np.vstack((theta_bin, RV_bin["theta_bin"]))
            bin_idx = np.append(bin_idx, RV_bin["bin_idx"])
            bin_idx_lst = np.append(bin_idx_lst, RV_bin["bin_idx_lst"])
            allele_flip_local = np.append(allele_flip_local, RV_bin["flip"])

    ## resolve results for global phasing input
    RV_region = {
        "AD_phased": AD_phased,
        "bin_idx": bin_idx,
        "ad_bin_softcnt": ad_bin_softcnt,
        "ad_bin": ad_bin,
        "dp_bin": dp_bin,
        "theta_bin": theta_bin,
        "allele_flip_local": allele_flip_local,
        "bin_idx_lst": bin_idx_lst
    }
    # return AD_phased
    return RV_region


def BAF_Local_phasing(Xdata, chr_lst = None, 
                      region_key = "chr", 
                      phasing_len = 100, 
                      region_nproc=1, 
                      bin_nproc=1, 
                      feature_mode = "GENE", # feature_mode ="BLOCK"
                      var_add = None,
                      verbose = False):
    """
    Func:
    phasing_len: default for 100 genes.

    """
    from scipy import sparse
    start_t = datetime.datetime.now()

    if chr_lst is None:
        chr_lst = Xdata.var[region_key].drop_duplicates(keep="first")

    if region_nproc > 1:
        result = []
        pool = multiprocessing.Pool(processes = region_nproc)
        for chr_ in chr_lst:
            tmp_Xdata = select_chr_region(Xdata, chr_, region_key)
            AD_region = tmp_Xdata.layers["AD"].T
            DP_region = tmp_Xdata.layers["DP"].T
            result.append(pool.apply_async(process_region,(AD_region, DP_region, phasing_len, bin_nproc), callback = None))
        pool.close()
        pool.join()
        result = [res.get() for res in result]
    else:
        result = []
        for chr_ in chr_lst:
            tmp_Xdata = select_chr_region(Xdata, chr_, region_key)
            AD_region = tmp_Xdata.layers["AD"].T
            DP_region = tmp_Xdata.layers["DP"].T
            RV_region = process_region(AD_region, DP_region, phasing_len = phasing_len, nproc=bin_nproc)
            result.append(RV_region)

    ## process the data from the list to a dataset
    for i , RV_region in zip(range(len(result)), result):
        if i == 0:
            AD_phased = RV_region["AD_phased"]

            theta_bin = RV_region["theta_bin"]
            bin_idx = RV_region["bin_idx"]
            bin_idx_lst = RV_region["bin_idx_lst"]
            ad_bin_softcnt = RV_region["ad_bin_softcnt"]
            ad_bin = RV_region["ad_bin"]
            dp_bin = RV_region["dp_bin"]
            allele_flip_local = RV_region["allele_flip_local"]

        else:
            AD_phased = np.vstack((AD_phased, RV_region["AD_phased"]))
            
            theta_bin = np.vstack((theta_bin, RV_region["theta_bin"]))
            bin_idx = np.append(bin_idx, RV_region["bin_idx"])
            bin_idx_lst = np.append(bin_idx_lst, RV_region["bin_idx_lst"])
            ad_bin_softcnt = np.vstack((ad_bin_softcnt, RV_region["ad_bin_softcnt"]))
            ad_bin = np.vstack((ad_bin, RV_region["ad_bin"]))
            dp_bin = np.vstack((dp_bin, RV_region["dp_bin"]))
            allele_flip_local = np.append(allele_flip_local, RV_region["allele_flip_local"])
         
    end_t = datetime.datetime.now()
    elapsed_sec = (end_t - start_t).total_seconds()
    
    print("[XClone-Local_phasing] time_used: " + "{:.2f}".format(elapsed_sec) + "seconds")

    update_Xdata = Xdata.copy()

    update_Xdata.layers["AD_phased"] = AD_phased.T
    update_Xdata.var["bin_idx"] = bin_idx_lst
    update_Xdata.var["allele_flip_local"] = allele_flip_local
    try:
        update_Xdata.var["allele_flip_local"].replace({0: False, 1: True}, inplace = True)
    except Exception as e:
        print("[XClone Warning]", e)
    else:
        print("[XClone hint] get allele flip status from local phasing.")
    
    update_Xdata.obsm["theta_bin"] = theta_bin.T
    update_Xdata = process_bin_id(update_Xdata, region_key, verbose)

    ##check theta_bin reuslts first and there are nan value
    ## save for visualization
    ad_bin_softcnt = sparse.csr_matrix(ad_bin_softcnt)
    ad_bin = sparse.csr_matrix(ad_bin)
    dp_bin = sparse.csr_matrix(dp_bin)

    ## generate merge_Xdata var
    merge_var = update_Xdata.var.copy()
    merge_var.drop_duplicates("bin_idx_cum", keep="first", inplace=True)
    # print(merge_var)
    merge_var.rename(columns={'stop': 'gene1_stop'}, inplace=True)

    merge_var2 = update_Xdata.var.copy()
    merge_var2.drop_duplicates("bin_idx_cum", keep="last", inplace=True)
    # print(merge_var2)
    merge_var["stop"] = merge_var2["stop"].copy().tolist()
    merge_var["bin_stop_arm"] = merge_var2["arm"].copy().tolist()
    merge_var["bin_stop_chr_arm"] = merge_var2["chr_arm"].copy().tolist()
    merge_var["bin_stop_band"] = merge_var2["band"].copy().tolist()

    
    merge_var3 =  update_Xdata.var.copy()
    def get_bin_genes(merge_var, group_key = "bin_idx_cum"):
        GeneName_lst = []
        GeneID_lst = []
        GeneName_dict = {}
        GeneID_dict = {}
        dict_ = merge_var.groupby(group_key).groups
        for key, idx_ in dict_.items():
            tmp_genename_lst = merge_var.loc[idx_]['GeneName'].copy().tolist()
            tmp_geneid_lst = merge_var.loc[idx_]['GeneID'].copy().tolist()
            GeneName_lst.append(tmp_genename_lst)
            GeneID_lst.append(tmp_geneid_lst)
            GeneName_dict[key] = tmp_genename_lst
            GeneID_dict[key] = tmp_geneid_lst
        
        return GeneName_lst, GeneID_lst, GeneName_dict, GeneID_dict
    if feature_mode == "GENE":
        GeneName_lst, GeneID_lst, GeneName_dict, GeneID_dict = get_bin_genes(merge_var3, group_key = "bin_idx_cum")
        merge_var["GeneName_lst"] = [",".join(x) for x in GeneName_lst]
        merge_var["GeneID_lst"] = [",".join(x) for x in GeneID_lst]
        
        var_keep_lst = ["chr", "start", "stop", "arm", "chr_arm", "band", "gene1_stop", 
        "bin_stop_arm","bin_stop_chr_arm", "bin_stop_band", "bin_idx", "bin_idx_cum",
        "GeneName_lst", "GeneID_lst"]
        if var_add is not None:
            if var_add in Xdata.var.columns:
                var_keep_lst.append(var_add)
                # print("add var: ", var_add)

        merge_var = merge_var[var_keep_lst]
        merge_var["bin_genes_cnt"] = merge_var["GeneName_lst"].str.len()
    else:
        var_keep_lst = ["chr", "start", "stop", "arm", "chr_arm", "band", 
        "bin_stop_arm","bin_stop_chr_arm", "bin_stop_band", "bin_idx", "bin_idx_cum"]
        if var_add is not None:
            if var_add in Xdata.var.columns:
                var_keep_lst.append(var_add)
                # print("add var: ", var_add)

        merge_var = merge_var[var_keep_lst]


    merge_Xdata = ad.AnnData(ad_bin.T, var = merge_var, obs = update_Xdata.obs.copy())
    
    # ## in case that the list saved as string in anndata
    # merge_Xdata.uns["GeneName_lst"] = GeneName_dict
    # merge_Xdata.uns["GeneID_lst"] = GeneID_dict

    ## soft phasing
    merge_Xdata.layers["ad_bin_softcnt"] = ad_bin_softcnt.T
    ## hard phasing
    merge_Xdata.layers["ad_bin"] = ad_bin.T
    merge_Xdata.layers["dp_bin"] = dp_bin.T

    merge_Xdata.uns["local_phasing_key"] = region_key
    merge_Xdata.uns["local_phasing_len"] = phasing_len

    del Xdata
    gc.collect()
    
    return update_Xdata, merge_Xdata
    

def process_bin_id(Xdata, region_key = "chr", verbose = False):
    """
    Func used in Func `BAF_Local_phasing`.
    """
    update_Xdata = Xdata.copy()
    ## whole genome bin idx
    region_bin_counts = np.array(update_Xdata.var.drop_duplicates(region_key, keep="last")["bin_idx"]) + 1
    cumbin_counts = np.cumsum(region_bin_counts)[:-1]
    if verbose:
        print(cumbin_counts)
    update_Xdata.var["bin_idx_cum"] = update_Xdata.var["bin_idx"].copy()
    region_brk_item = update_Xdata.var[region_key].drop_duplicates(keep="first")
    
    brk1 = region_brk_item[0]
    flag_ = update_Xdata.var[region_key] == brk1
    bin_idx_cum = update_Xdata[:, flag_].var["bin_idx"]
    
    for cumcount_, brk_ in zip(cumbin_counts, region_brk_item[1:]):
        if verbose:
            print(cumcount_, brk_)
        flag_ = update_Xdata.var[region_key] == brk_
        tmp_Xdata = update_Xdata[:, flag_]
        bin_idx_cum = np.append(bin_idx_cum, tmp_Xdata.var["bin_idx_cum"] + cumcount_)
    update_Xdata.var["bin_idx_cum"] = bin_idx_cum
    
    
    del Xdata
    gc.collect()
    
    return update_Xdata

### global phasing
'''
def BAF_Global_phasing(Xdata, bin_Xdata):
    """
    ##
    need test
    todo
    1) get filp status for bins (genes)
    2) add plots for global phasing part (Xlayer = "BAF_globbal_phased" or global_phased = True)

    """
    start_t = datetime.datetime.now()

    ## Global_Phasing
    is_flips, distances, theta_new = Global_Phasing(Xdata.obsm["theta_bin"].T)

    ## apply global phasing method on original AD_phased
    
    bin_value_counts = Xdata.var["bin_idx_cum"].value_counts()
    bin_item = Xdata.var["bin_idx_cum"].drop_duplicates(keep="first")
    
    ## check is_flips the same length with bin_item
    if is_flips.shape[0] == bin_item.shape[0]:
        pass
    else:
        raise ValueError("[XClone] Pls check bin_item!")

    for i, bin_id in enumerate(bin_item):
        if i == 0:
            gene_flip = np.repeat(is_flips[i], bin_value_counts[bin_id])
        else:
            gene_flip = np.append(gene_flip, np.repeat(is_flips[i], bin_value_counts[bin_id]))

    # get gene_flip for AD_Phased
    AD_phased = Xdata.layers["AD_phased"]
    BD_phased = Xdata.layers["DP"] - Xdata.layers["AD_phased"]
    AD_global_phased = AD_phased + 0
    AD_global_phased[:, gene_flip] = BD_phased[: , gene_flip] + 0

    update_Xdata = Xdata.copy()
    update_Xdata.layers["AD_global_phased"] = AD_global_phased
    update_Xdata.var["allele_flip_global"] = gene_flip
    print("[XClone hint] get allele flip status from global phasing.")
    if {'allele_flip_local', 'allele_flip_global'}.issubset(update_Xdata.var.columns):
        update_Xdata.var["allele_flip"] =  update_Xdata.var["allele_flip_local"] ^ update_Xdata.var["allele_flip_global"]
        print("[XClone hint] get final allele flip status.")

    ## apply global phasing method on AD_phased bins
    ad_bin_softcnt = bin_Xdata.layers["ad_bin_softcnt"]
    ad_bin = bin_Xdata.layers["ad_bin"]
    dp_bin = bin_Xdata.layers["dp_bin"]
    
    bd_bin_softcnt = dp_bin - ad_bin_softcnt
    ad_bin_softcnt_phased = ad_bin_softcnt + 0
    ad_bin_softcnt_phased[: , is_flips] = bd_bin_softcnt[:, is_flips] + 0

    bd_bin = dp_bin - ad_bin
    ad_bin_phased = ad_bin + 0
    ad_bin_phased[: , is_flips] = bd_bin[:, is_flips] + 0

    bin_Xdata.layers["ad_bin_softcnt_phased"] = ad_bin_softcnt_phased
    ## the ad_bin_softcnt_phased is weighted counts(bin counts is weighted)(soft phasing)
    bin_Xdata.layers["ad_bin_phased"] = ad_bin_phased
    ## ad_bin_phased is hard phasing counts.

    end_t = datetime.datetime.now()
    elapsed_sec = (end_t - start_t).total_seconds()
    print("[XClone-Global_phasing] time_used: " + "{:.2f}".format(elapsed_sec) + "seconds")

    del Xdata
    gc.collect()
    
    return update_Xdata, bin_Xdata

'''
    
def BAF_Global_phasing(Xdata, bin_Xdata):
    """
    Global phasing per chromosome.
    For each chromosome, perform global phasing on the corresponding theta_bin,
    then combine is_flips for all chromosomes.
    """

    start_t = datetime.datetime.now()

    # Prepare
    chr_list = bin_Xdata.var['chr'].drop_duplicates(keep="first")
    theta_bin = Xdata.obsm["theta_bin"]
    chr_var = bin_Xdata.var['chr']

    is_flips_all = []
    distances_all = []
    theta_new_all = []

    for chr_ in chr_list:
        chr_mask = chr_var == chr_
        # Get the indices for this chromosome
        chr_indices = np.where(chr_mask)[0]
        # Get the theta_bin for this chromosome (shape: n_bins_chr, n_cells)
        theta_chr = theta_bin[:, chr_indices]
        # Global phasing for this chromosome
        is_flips_chr, distances_chr, theta_new_chr = Global_Phasing(theta_chr.T)
        is_flips_all.append(is_flips_chr)
        # distances_all.append(distances_chr)
        # theta_new_all.append(theta_new_chr)

    # Concatenate results for all chromosomes
    is_flips = np.concatenate(is_flips_all)
    # distances = np.concatenate(distances_all)
    # theta_new = np.concatenate(theta_new_all, axis=0)

    # The rest of the function is unchanged
    bin_value_counts = Xdata.var["bin_idx_cum"].value_counts()
    bin_item = Xdata.var["bin_idx_cum"].drop_duplicates(keep="first")
    
    if is_flips.shape[0] == bin_item.shape[0]:
        pass
    else:
        raise ValueError("[XClone] Pls check bin_item!")

    for i, bin_id in enumerate(bin_item):
        if i == 0:
            gene_flip = np.repeat(is_flips[i], bin_value_counts[bin_id])
        else:
            gene_flip = np.append(gene_flip, np.repeat(is_flips[i], bin_value_counts[bin_id]))

    # get gene_flip for AD_Phased
    AD_phased = Xdata.layers["AD_phased"]
    BD_phased = Xdata.layers["DP"] - Xdata.layers["AD_phased"]
    AD_global_phased = AD_phased + 0
    AD_global_phased[:, gene_flip] = BD_phased[: , gene_flip] + 0

    update_Xdata = Xdata.copy()
    update_Xdata.layers["AD_global_phased"] = AD_global_phased
    update_Xdata.var["allele_flip_global"] = gene_flip
    print("[XClone hint] get allele flip status from global phasing.")
    if {'allele_flip_local', 'allele_flip_global'}.issubset(update_Xdata.var.columns):
        update_Xdata.var["allele_flip"] =  update_Xdata.var["allele_flip_local"] ^ update_Xdata.var["allele_flip_global"]
        print("[XClone hint] get final allele flip status.")

    # apply global phasing method on AD_phased bins
    ad_bin_softcnt = bin_Xdata.layers["ad_bin_softcnt"]
    ad_bin = bin_Xdata.layers["ad_bin"]
    dp_bin = bin_Xdata.layers["dp_bin"]
    
    bd_bin_softcnt = dp_bin - ad_bin_softcnt
    ad_bin_softcnt_phased = ad_bin_softcnt + 0
    ad_bin_softcnt_phased[: , is_flips] = bd_bin_softcnt[:, is_flips] + 0

    bd_bin = dp_bin - ad_bin
    ad_bin_phased = ad_bin + 0
    ad_bin_phased[: , is_flips] = bd_bin[:, is_flips] + 0

    bin_Xdata.layers["ad_bin_softcnt_phased"] = ad_bin_softcnt_phased
    bin_Xdata.layers["ad_bin_phased"] = ad_bin_phased

    end_t = datetime.datetime.now()
    elapsed_sec = (end_t - start_t).total_seconds()
    print("[XClone-Global_phasing] time_used: " + "{:.2f}".format(elapsed_sec) + "seconds")

    del Xdata
    import gc
    gc.collect()
    
    return update_Xdata, bin_Xdata

def BAF_Global_phasing_rev(Xdata, bin_Xdata):
    """
    Global phasing per chromosome (reversed).
    For each chromosome, reverse the bin order before performing global phasing,
    then combine is_flips for all chromosomes.
    """

    start_t = datetime.datetime.now()

    # Prepare
    chr_list = bin_Xdata.var['chr'].drop_duplicates(keep="first")
    theta_bin = Xdata.obsm["theta_bin"]
    chr_var = bin_Xdata.var['chr']

    is_flips_all = []

    for chr_ in chr_list:
        chr_mask = chr_var == chr_
        chr_indices = np.where(chr_mask)[0]
        # Reverse the bin indices for this chromosome
        chr_indices_rev = chr_indices[::-1]
        # Get reversed theta_bin (shape: n_bins_chr, n_cells)
        theta_chr_rev = theta_bin[:, chr_indices_rev]
        # Global phasing on reversed bins
        is_flips_chr_rev, _, _ = Global_Phasing(theta_chr_rev.T)
        # Reorder flips back to original bin order
        is_flips_chr = is_flips_chr_rev[::-1]
        is_flips_all.append(is_flips_chr)

    # Concatenate results for all chromosomes
    is_flips = np.concatenate(is_flips_all)

    # Map flips to genes
    bin_value_counts = Xdata.var["bin_idx_cum"].value_counts()
    bin_item = Xdata.var["bin_idx_cum"].drop_duplicates(keep="first")

    if is_flips.shape[0] != bin_item.shape[0]:
        raise ValueError("[XClone] Pls check bin_item!")

    gene_flip = np.concatenate([
        np.repeat(is_flips[i], bin_value_counts[bin_id])
        for i, bin_id in enumerate(bin_item)
    ])

    # Apply global phasing to AD_phased
    AD_phased = Xdata.layers["AD_phased"]
    BD_phased = Xdata.layers["DP"] - AD_phased
    AD_global_phased = AD_phased.copy()
    AD_global_phased[:, gene_flip] = BD_phased[:, gene_flip]

    update_Xdata = Xdata.copy()
    update_Xdata.layers["AD_global_phased_rev"] = AD_global_phased
    update_Xdata.var["allele_flip_global_rev"] = gene_flip

    print("[XClone hint] get allele flip status from reversed global phasing.")
    if {'allele_flip_local', 'allele_flip_global_rev'}.issubset(update_Xdata.var.columns):
        update_Xdata.var["allele_flip_rev"] = update_Xdata.var["allele_flip_local"] ^ update_Xdata.var["allele_flip_global_rev"]
        print("[XClone hint] get final allele flip status (reversed).")

    # Apply global phasing to bin-level AD
    ad_bin_softcnt = bin_Xdata.layers["ad_bin_softcnt"]
    ad_bin = bin_Xdata.layers["ad_bin"]
    dp_bin = bin_Xdata.layers["dp_bin"]

    bd_bin_softcnt = dp_bin - ad_bin_softcnt
    ad_bin_softcnt_phased = ad_bin_softcnt.copy()
    ad_bin_softcnt_phased[:, is_flips] = bd_bin_softcnt[:, is_flips]

    bd_bin = dp_bin - ad_bin
    ad_bin_phased = ad_bin.copy()
    ad_bin_phased[:, is_flips] = bd_bin[:, is_flips]

    bin_Xdata.layers["ad_bin_softcnt_phased_rev"] = ad_bin_softcnt_phased
    bin_Xdata.layers["ad_bin_phased_rev"] = ad_bin_phased

    end_t = datetime.datetime.now()
    elapsed_sec = (end_t - start_t).total_seconds()
    print("[XClone-Global_phasing_rev] time_used: {:.2f} seconds".format(elapsed_sec))

    del Xdata
    gc.collect()
    return update_Xdata, bin_Xdata



## Visualization
## after global phasing and apply moving average and KNN smoothing 
## mainly for visualization but not for methods

# strategy notes*
# BAF_smoothed1 <- KNN_connectivity @ AF 
# BAF_smoothed2 <- BAF_smoothed1 @ WMA_connectivity



def BAF_fillna(Xdata, Xlayer = "BAF", out_layer = "fill_BAF"):
    """
    """
    Xmtx_used = Xdata.layers[Xlayer].copy()
    Xmtx_used = np.nan_to_num(Xmtx_used, nan=0.5)
    Xdata.layers[out_layer] = Xmtx_used
    return Xdata


def BAF_remove_ref(Xdata, Xlayer = "BAF", out_layer = "reref_BAF", 
                   anno_key = "cell_type", ref_cell = "unclassified"):
    """
    apply logit func (BAF center is 0.5)
    """
    from scipy.special import logit, expit

    Xmtx_used = Xdata.layers[Xlayer].copy()
    Xmtx_used = logit(Xmtx_used)

    is_ref_ = Xdata.obs[anno_key] == ref_cell
    ref_ = Xmtx_used[is_ref_].mean(axis = 0)

    Xdata.layers[out_layer] = expit(Xmtx_used - ref_)

    return Xdata



