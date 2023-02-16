"""Base functions for XClone BAF and RDR combination:
XClone CNV states combination.
# Option 1: merge genes into bin before combination.
# Option 2: extend bins to genes scale before combination.
"""

# Author: Rongting Huang
# Date: 2021/11/25
# update: 2022/09/25

import numpy as np
import time

from scipy.special import logsumexp
from .utils import pairwise1
from .base_utils import normalize
from ..preprocessing._data import check_RDR_BAF_cellorder, check_RDR_BAF_chrmapping
from ..preprocessing._Xdata_manipulation import Xdata_cell_selection

from ..plot.CNV_plot import remove_Xdata_layers


# Option1: 

## get RDR_merge_Xdata
def cells_mapping(Xdata,
                  BAF_merge_Xdata):
    """
    Function:
    ----------
    1) check cells in the same order before merge genes to bin.
    2) init new RDR_merge_Xdata
    """
    cells_num = Xdata.shape[0]
    BAF_cells_num = BAF_merge_Xdata.shape[0]
    if cells_num == BAF_cells_num:
        RDR_merge_Xdata = BAF_merge_Xdata.copy()
    else:
        Xdata.obs["cellbarcodes"] = Xdata.obs.index
        cell_lst_ = Xdata.obs["cellbarcodes"].copy().tolist()

        RDR_merge_Xdata = BAF_merge_Xdata.copy()
        RDR_merge_Xdata.obs["cellbarcodes"] = RDR_merge_Xdata.obs.index
        
        RDR_merge_Xdata = Xdata_cell_selection(RDR_merge_Xdata, 
                                           select = True, 
                                           cell_anno_key = "cellbarcodes",
                                           cell_lst = cell_lst_,
                                           update_uns=False)
    
    # remove layers in RDR_merge_Xdata for memory saving.
    RDR_merge_Xdata = remove_Xdata_layers(RDR_merge_Xdata, copy = True)
    # rename var items inherit from BAF
    RDR_merge_Xdata.var.rename(columns={'GeneName_lst': 'BAF_GeneName_lst', 
    'GeneID_lst': 'BAF_GeneID_lst', 'bin_genes_cnt': 'BAF_bin_genes_cnt'}, inplace=True)
    # del things from BAF
    if "GeneName_lst" in RDR_merge_Xdata.uns.keys():
        del RDR_merge_Xdata.uns["GeneName_lst"]
    if "GeneID_lst" in RDR_merge_Xdata.uns.keys():
        del RDR_merge_Xdata.uns["GeneID_lst"]

    data_notes = time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()) 
    RDR_merge_Xdata.uns["data_notes"] = data_notes

    return RDR_merge_Xdata

# merge RDR module from genes to bins to match with BAF_merge_Xdata
def gene_to_bin_mapping(Xdata,
                        BAF_merge_Xdata,
                        Xlayer = "posterior_mtx_log",
                        merge_layer = "RDR_merge_prob",
                        return_prob = False):
    """
    Function:
    ----------
    merge log probability from genes to bins.
    And the output contains both log prob and prob.

    Parameters:
    -----------
    Xlayer: probability in log
    "posterior_mtx"  results mapping

    Example:
    --------
    >>> import xclone
    >>> prob_log = xclone.model.gene_to_bin_mapping(RDR_adata, 
        BAF_merge_Xdata, return_prob=True)
    """
    print("[XCLone] RDR merge genes to bin")
    ## check cells number
    RDR_merge_Xdata = cells_mapping(Xdata, BAF_merge_Xdata)
    ## check cell order
    success_flag = check_RDR_BAF_cellorder(RDR_merge_Xdata, BAF_merge_Xdata)
    fail_flag = not success_flag
    if fail_flag:
        print("[XClone-combination]gene to bin mapping:")
        raise ValueError("pls check input: cell orders of the two Xdata are not matched.")
    ## check chr mapping
    success_flag = check_RDR_BAF_chrmapping(RDR_merge_Xdata, BAF_merge_Xdata)
    fail_flag = not success_flag
    if fail_flag:
        print("[XClone-combination] gene to bin mapping:")
        raise ValueError("[XClone-combination] pls check input: chr nums of the two Xdata are not matched.")
    
    Xdata.var["gene_index"] = [int(x) for x in Xdata.var.index]
    BAF_merge_Xdata.var["brk_gene_index"] = [int(x) for x in BAF_merge_Xdata.var.index]
    phasing_len = BAF_merge_Xdata.uns["local_phasing_len"]
    last_item = BAF_merge_Xdata.var["brk_gene_index"][-1] + phasing_len
    brk_index = BAF_merge_Xdata.var["brk_gene_index"].copy().tolist() + [last_item]

    phasing_key = BAF_merge_Xdata.uns["local_phasing_key"]
    Xdata.var[phasing_key].drop_duplicates(keep="last").index

    # output Xdata merge results
    merge_results = []
    merge_genes_lst = []
    merge_geneID_lst = []

    merge_genes_dict = {}
    merge_geneID_dict = {}

    i = 0 
    for idx1, idx2 in pairwise1(brk_index):
        flag_ = np.vectorize(lambda x: idx1 <= x < idx2)(Xdata.var["gene_index"])
        if flag_.sum() == 0:
            print("No genes in this bin:", idx1, idx2, ", inherit previous bin.")
            # use previous bins prob-adopt this strategy here| for LOH corrected next step
            # todo update, remove bins to keep another version of outpout if need.| keep in RDR merged anndata.
            # tmp_res = tmp_res
            # np.zeros((cell_nums, 1, states_nums))
            # tmp_res[:,:,1] = 1
            tmp_genes = [""]
            tmp_geneID = [""]

        else:
            # tmp_res = Xdata.layers[Xlayer][:,flag_,:].sum(axis = 1, keepdims= True)
            tmp_res = np.exp(Xdata.layers[Xlayer])[:,flag_,:].mean(axis = 1, keepdims= True)
            # tmp_res = np.median(np.exp(Xdata.layers[Xlayer])[:,flag_,:], axis = 1, keepdims= True)

            tmp_genes = Xdata.var["GeneName"][flag_].copy().tolist()
            tmp_geneID = Xdata.var["GeneID"][flag_].copy().tolist()

        merge_results.append(tmp_res)
        
        merge_genes_lst.append(tmp_genes)
        merge_geneID_lst.append(tmp_geneID)
        
        merge_genes_dict[i] = tmp_genes
        merge_geneID_dict[i] = tmp_geneID

        i+=1
    
    res_log = np.hstack(np.log(merge_results))
    
    # normalised the new merged probabilty
    res_log += -logsumexp(res_log, axis=2, keepdims=True)
    res_ = np.exp(res_log)

    ## output for visualization
    BAF_merge_Xdata.layers[merge_layer] = res_
    RDR_merge_Xdata.layers[merge_layer] = res_

    layer_ = merge_layer + "_log"
    BAF_merge_Xdata.layers[layer_] = res_log
    RDR_merge_Xdata.layers[layer_] = res_log

    
    RDR_merge_Xdata.var["GeneName_lst"] = merge_genes_lst
    RDR_merge_Xdata.var["GeneID_lst"] = merge_geneID_lst
    RDR_merge_Xdata.var["bin_genes_cnt"] = RDR_merge_Xdata.var["GeneName_lst"].str.len()

    RDR_merge_Xdata.uns["GeneName_lst"] = merge_genes_dict
    RDR_merge_Xdata.uns["GeneID_lst"] = merge_geneID_dict
  
    BAF_merge_Xdata.var["RDR_GeneName_lst"] = merge_genes_lst
    BAF_merge_Xdata.var["RDR_GeneID_lst"] = merge_geneID_lst
    BAF_merge_Xdata.uns["RDR_GeneName_lst"] = merge_genes_dict
    BAF_merge_Xdata.uns["RDR_GeneID_lst"] = merge_geneID_dict

    if return_prob == True:
        return res_log
    else:
        return BAF_merge_Xdata, RDR_merge_Xdata

def CNV_prob_combination11(RDR_Xdata,
                           BAF_Xdata,
                           RDR_layer = "RDR_merge_prob",
                           BAF_layer = "posterior_mtx"):
    """
    deprecated.
    use merged to bin level RDR prob combine with BAF prob.

    """
    rdr_prob_ = RDR_Xdata.layers[RDR_layer].copy()
    baf_prob_ = BAF_Xdata.layers[BAF_layer].copy()

    rdr_prob_ = np.expand_dims(rdr_prob_, axis = -1)
    baf_prob_ = np.expand_dims(baf_prob_, axis = -2)

    combine_base_prob = rdr_prob_ * baf_prob_
    RDR_Xdata.layers["combine_base_prob"] = combine_base_prob
    BAF_Xdata.layers["combine_base_prob"] = combine_base_prob

    ## combination adjust
    # RDR_merge_Xdata-corrected for copy loss identification
    loss_corrected_prob = copyloss_corrected(combine_base_prob)

    RDR_Xdata.layers["loss_corrected_prob"] = loss_corrected_prob
    BAF_Xdata.layers["loss_corrected_prob"] = loss_corrected_prob
    return RDR_Xdata, BAF_Xdata

# Option2: get extend BAF module from bins to genes scale to match with RDR_Xdata.

def bin_to_gene_mapping(BAF_merge_Xdata,
                        RDR_Xdata,
                        Xlayer = "posterior_mtx_log",
                        extend_layer = "BAF_extend_prob",
                        return_prob = False):
    """
    Function:
    ---------
    output BAF information in RDR gene resolution.
    extend BAF Xdata from bins to genes.
    """

    print("[XClone] BAF extend bins to genes.")
    ## check cell order
    success_flag = check_RDR_BAF_cellorder(RDR_Xdata, BAF_merge_Xdata)
    fail_flag = not success_flag
    if fail_flag:
        print("[XClone-combination]gene to bin mapping:")
        raise ValueError("[XClone-combination] pls check input: cell orders of the two Xdata are not matched.")
    ## check chr mapping
    success_flag = check_RDR_BAF_chrmapping(RDR_Xdata, BAF_merge_Xdata)
    fail_flag = not success_flag
    if fail_flag:
        print("[XClone-combination]gene to bin mapping:")
        raise ValueError("[XClone-combination] pls check input: chr nums of the two Xdata are not matched.")

    RDR_Xdata.var["gene_index"] = [int(x) for x in RDR_Xdata.var.index]
    BAF_merge_Xdata.var["brk_gene_index"] = [int(x) for x in BAF_merge_Xdata.var.index]
    phasing_len = BAF_merge_Xdata.uns["local_phasing_len"]
    last_item = BAF_merge_Xdata.var["brk_gene_index"][-1] + phasing_len
    brk_index = BAF_merge_Xdata.var["brk_gene_index"].copy().tolist() + [last_item]

    phasing_key = BAF_merge_Xdata.uns["local_phasing_key"]
    RDR_Xdata.var[phasing_key].drop_duplicates(keep="last").index

    # output Xdata extend results
    use_mtx = BAF_merge_Xdata.layers[Xlayer]
    bin_cnt = BAF_merge_Xdata.var.shape[0]
    # extend_results = []
    mtx_ndim= use_mtx.ndim
    
    def extend_base(brk_index, RDR_Xdata, use_mtx, bin_cnt, dimensions = 3):
        """
        """
        extend_results = []
        if dimensions == 3:
            cnt = 0 
            for idx1, idx2 in pairwise1(brk_index):
                flag_ = np.vectorize(lambda x: idx1 <= x < idx2)(RDR_Xdata.var["gene_index"])
                repeat_cnt = flag_.sum()
                if repeat_cnt == 0:
                    print("No genes in this bin:", idx1, idx2, ", skip this bin.")
                else:
                    tmp_extend_mtx = np.expand_dims(use_mtx[:,cnt,:], axis = 1).repeat(repeat_cnt, axis = 1)
                    extend_results.append(tmp_extend_mtx)
                cnt+=1
        if dimensions == 2:
            cnt = 0 
            for idx1, idx2 in pairwise1(brk_index):
                flag_ = np.vectorize(lambda x: idx1 <= x < idx2)(RDR_Xdata.var["gene_index"])
                repeat_cnt = flag_.sum()
                if repeat_cnt == 0:
                    print("No genes in this bin:", idx1, idx2, ", skip this bin.")
                else:
                    tmp_extend_mtx = np.expand_dims(use_mtx[:,cnt], axis = 1).repeat(repeat_cnt, axis = 1)
                    extend_results.append(tmp_extend_mtx)
                cnt+=1

        if cnt != bin_cnt:
            raise ValueError("[XClone extend base]Pls ensure all bins are processing to extend.")
        return extend_results
    
    extend_results = extend_base(brk_index, RDR_Xdata, use_mtx, bin_cnt, mtx_ndim)

    extend_res = np.hstack(extend_results)

    RDR_Xdata.layers[extend_layer] = extend_res
    # combined_Xdata = RDR_Xdata.copy()
    # combined_Xdata.layers[extend_layer] = extend_res
    if return_prob == True:
        return extend_res
    else:
        return RDR_Xdata
        # return combined_Xdata

def copyloss_corrected(Xdata, Xlayer, mode = 1):
    """
    (1) can use bafprob to find out allele-specific copy loss
    (2) adjust copy loss(RDR)-copy neutral(BAF) situation.
    """
    prob_ = Xdata.layers[Xlayer].copy()
    
    ## BAF 3 states
    if prob_.shape[-1] == 3:
        # strategy 1: default used.
        if mode == 1:
            prob_[:,:,0,0] += prob_[:,:,0,1]/3
            prob_[:,:,0,2] += prob_[:,:,0,1]/3
            prob_[:,:,1,1] += prob_[:,:,0,1]/3

            prob_[:,:,0,1] = 0
    
        # strategy 2:
        if mode == 2:
            prob_[:,:,1,1] += prob_[:,:,0,1]
            prob_[:,:,0,1] = 0
    
    ## BAF 5 states
    if prob_.shape[-1] == 5:
        if mode == 1:
            prob_[:,:,0,1] += prob_[:,:,0,2]/3
            prob_[:,:,0,3] += prob_[:,:,0,2]/3
            prob_[:,:,1,2] += prob_[:,:,0,2]/3

            prob_[:,:,0,2] = 0
        if mode == 2:
            prob_[:,:,1,2] += prob_[:,:,0,2]
            prob_[:,:,0,2] = 0
    
    return prob_

def copygain_corrected(Xdata, Xlayer, mode = 1):
    """
    not sure if necessary; can improve in next release.
    Not used in this release ("0.3.4")
    can use bafprob to find out allele-specific copy gain.
    However, due to narrow distance between allele-specific copy gain and allele balance,
    it is hard to do now.
    """
    prob_ = Xdata.layers[Xlayer].copy()
    ## BAF 3 states
    ### default mode = 1
    if prob_.shape[-1] == 3:
        prob_[:,:,2,0] += prob_[:,:,2,1]/3
        prob_[:,:,2,2] += prob_[:,:,2,1]/3
        prob_[:,:,1,1] += prob_[:,:,2,1]/3

        prob_[:,:,2,1] = 0
    
    ## BAF 5 states
    if prob_.shape[-1] == 5:
        if mode == 1:
            prob_[:,:,2,1] += prob_[:,:,2,2]/3
            prob_[:,:,2,3] += prob_[:,:,2,2]/3
            prob_[:,:,1,2] += prob_[:,:,2,2]/3

            prob_[:,:,2,2] = 0
        
        if mode == 2:
            prob_[:,:,1,2] += prob_[:,:,1,1]/2
            prob_[:,:,2,1] += prob_[:,:,1,1]/2
            prob_[:,:,1,1] = 0
            prob_[:,:,1,2] += prob_[:,:,1,3]/2
            prob_[:,:,2,3] += prob_[:,:,1,3]/2
            prob_[:,:,1,3] = 0
        
        if mode == 3:
            prob_[:,:,1,2] += prob_[:,:,1,1]*1/3
            prob_[:,:,2,1] += prob_[:,:,1,1]*2/3
            prob_[:,:,1,1] = 0
            prob_[:,:,1,2] += prob_[:,:,1,3]*1/3
            prob_[:,:,2,3] += prob_[:,:,1,3]*2/3
            prob_[:,:,1,3] = 0
    return prob_

def CNV_prob_combination(Xdata,
                         RDR_layer = "posterior_mtx",
                         BAF_layer = "BAF_extend_post_prob",
                         copyloss_correct = True,
                         copyloss_correct_mode = 1,
                         copygain_correct = True,
                         copygain_correct_mode = 2):
    """
    Combine RDR prob with BAF prob.
    """
    rdr_prob_ = Xdata.layers[RDR_layer].copy()
    baf_prob_ = Xdata.layers[BAF_layer].copy()

    rdr_prob_ = np.expand_dims(rdr_prob_, axis = -1)
    baf_prob_ = np.expand_dims(baf_prob_, axis = -2)

    ## combination base
    combine_base_prob = rdr_prob_ * baf_prob_
    Xdata.layers["combine_base_prob"] = combine_base_prob
    
    ## combination adjust
    ## copy loss corrected-strategy1(default mode 1)
    if copyloss_correct & copygain_correct:
        corrected_prob = copyloss_corrected(Xdata, "combine_base_prob", mode = copyloss_correct_mode)
        Xdata.layers["corrected_prob"] = corrected_prob
        corrected_prob = copygain_corrected(Xdata, "corrected_prob", mode = copygain_correct_mode)
        Xdata.layers["corrected_prob"] = corrected_prob
    elif copyloss_correct:
        corrected_prob = copyloss_corrected(Xdata, "combine_base_prob", mode = copyloss_correct_mode)
        Xdata.layers["corrected_prob"] = corrected_prob
    elif copygain_correct:
        corrected_prob = copygain_corrected(Xdata, "combine_base_prob", mode = copygain_correct_mode)
        Xdata.layers["corrected_prob"] = corrected_prob
    else:
        Xdata.layers["corrected_prob"] = Xdata.layers["combine_base_prob"].copy()

    prob1_merge = CNV_prob_merge(Xdata, "corrected_prob")
    Xdata.layers["prob1_merge"] = prob1_merge

    return Xdata

def CNV_prob_merge(Xdata,
                   Xlayer):
    """
    Merge states prob for evaluation.
    """
    prob_ = Xdata.layers[Xlayer].copy()
    
    ## BAF 3 states
    if prob_.shape[-1] == 3:
        copy_loss = prob_[:,:,0,:].sum(axis = -1)
        loh = prob_[:,:,1,0] + prob_[:,:,1,2]
        copy_neutral = prob_[:,:,1,1] 
        copy_gain = prob_[:,:,2,:].sum(axis = -1)
    
    ## BAF 5 states
    if prob_.shape[-1] == 5:
        copy_loss = prob_[:,:,0,:].sum(axis = -1)
        loh = prob_[:,:,1,0] + prob_[:,:,1,4]
        copy_neutral = prob_[:,:,1,2]
        copy_gain_less = prob_[:,:,1,1] + prob_[:,:,1,3]
        copy_gain = prob_[:,:,2,:].sum(axis = -1) + copy_gain_less

    prob_merge = np.stack([copy_loss, loh, copy_neutral, copy_gain], axis = -1)
    return prob_merge

def CNV_prob_merge_for_plot(Xdata,
                            Xlayer = "corrected_prob"):
    """
    Merge states prob for more detailed evaluation and visualization.
    """
    prob_ = Xdata.layers[Xlayer].copy()
    ## BAF 3 states
    if prob_.shape[-1] == 3:
        copy_loss = prob_[:,:,0,:].sum(axis = -1)
        loh = prob_[:,:,1,0] + prob_[:,:,1,2]
        copy_neutral = prob_[:,:,1,1] 
        copy_gain = prob_[:,:,2,:].sum(axis = -1)

        copy_loss_A = prob_[:,:,0,2]
        copy_loss_B = prob_[:,:,0,0]
    
        loh_A = prob_[:,:,1,2]
        loh_B = prob_[:,:,1,0]
    
    ## BAF 5 states
    if prob_.shape[-1] == 5:
        copy_loss = prob_[:,:,0,:].sum(axis = -1)
        loh = prob_[:,:,1,0] + prob_[:,:,1,4]
        copy_neutral = prob_[:,:,1,2]
        copy_gain_less = prob_[:,:,1,1] + prob_[:,:,1,3]
        copy_gain = prob_[:,:,2,:].sum(axis = -1) + copy_gain_less

        copy_loss_A = prob_[:,:,0,3] + prob_[:,:,0,4]
        copy_loss_B = prob_[:,:,0,0] + prob_[:,:,0,1]
    
        loh_A = prob_[:,:,1,4]
        loh_B = prob_[:,:,1,0]

    plot_prob_merge1 = np.stack([copy_loss, loh, copy_neutral, copy_gain], axis = -1)
    Xdata.layers["plot_prob_merge1"] = plot_prob_merge1
    plot_prob_merge2 = np.stack([copy_loss_A, copy_loss_B, loh, copy_neutral, copy_gain], axis = -1)
    Xdata.layers["plot_prob_merge2"] = plot_prob_merge2
    plot_prob_merge3 = np.stack([copy_loss_A, copy_loss_B, loh_A, loh_B, copy_neutral, copy_gain], axis = -1)
    Xdata.layers["plot_prob_merge3"] = plot_prob_merge3
    plot_prob_merge4 = np.stack([copy_loss, loh_A, loh_B, copy_neutral, copy_gain], axis = -1)
    Xdata.layers["plot_prob_merge4"] = plot_prob_merge4
    
    return Xdata