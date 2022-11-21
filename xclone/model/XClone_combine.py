"""Base functions for XClone BAF and RDR combination
Base functions for XClone CNV states combination.
"""

# Author: Rongting Huang
# Date: 2021/11/25
# update: 2022/09/25

import numpy as np
import time

from scipy.special import logsumexp
from .utils import pairwise1
from .base_utils import normalize
from ..preprocessing._data import check_RDR_BAF_cellorder
from ..preprocessing._Xdata_manipulation import Xdata_cell_selection

from ..plot.CNV_plot import remove_Xdata_layers


# utils
def class2onehot(prob_, states_class_num):
    index_ = np.argmax(prob_, axis=-1)
    one_hot_ = np.eye(states_class_num)[index_]
    return one_hot_

# Part I: RDR_merge_Xdata

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
        

# Part II: BAF_merge_Xdata-corrected for LOH identification
## (based on information provided by RDR)
## LOH identification and evaluation[important]
def LOH_corrected(bafprob_, rdrprob_, mode = 0):
    """
    states table here for rdr:
    copy loss, copy neutral, copy gain
        states_table = np.array([[1,0,0],
                            [0,1,0],
                            [0,0,1]])

    mode = 0: mask copy gain for allele-specific loss visualization.
    mode = 1: mask copy loss for allele-specific gain visualization.
    [notes: mode = 1 not for the datasets used now]
    mode = 2: mask copy gain and loss and keep the corrected loh.

    """
    lohprob_ = bafprob_.copy()
    rdr_states = class2onehot(rdrprob_, 3)
    
    if mode == 0:
        ## comfirm not copy gain from rdr information, mask wrong loh
        rdr_comfirm = (~np.all(rdr_states == np.array([0,0,1]), axis=-1))
    if mode == 1:
        ## comfirm not copy loss from rdr information, mask wrong loh
        rdr_comfirm = (~np.all(rdr_states == np.array([1,0,0]), axis=-1))
    if mode == 2:
        ## keep true LOH information, BAF loss + RDR neutral
        rdr_comfirm = (np.all(rdr_states == np.array([0,1,0]), axis=-1))

    lohprob_[:,:,0] = lohprob_[:,:,0] * rdr_comfirm
    lohprob_[:,:,2] = lohprob_[:,:,2] * rdr_comfirm
    
    lohprob_corrected = normalize(lohprob_)
    return lohprob_corrected

def LOH_merge(lohprob):
    """
    merge 2 type of LOH into 1 for evaluation.
    original order: loh1, copy neutral, loh2
    """
    lohprob_merge = lohprob[:,:,[0,2]].copy()

    lohprob_merge = np.max(lohprob_merge, axis= -1, keepdims = True)
    lohprob_merge = np.insert(lohprob_merge, 1, values = lohprob[:,:,1],  axis=-1)

    lohprob_merge =normalize(lohprob_merge)

    return lohprob_merge

def BAF_LOH_merge(BAF_merge_Xdata, BAF_layer = "posterior_mtx"):
    """
    original order: copyloss1, copy neutral, copyloss2
    """
    BAF_prob_ = BAF_merge_Xdata.layers[BAF_layer].copy()
    loh_merge = LOH_merge(BAF_prob_)
    BAF_merge_Xdata.layers["lohprob_merge"] = loh_merge
    
    return BAF_merge_Xdata

def BAF_LOH_corrected(RDR_merge_Xdata,
                  BAF_merge_Xdata, 
                  RDR_layer = "RDR_merge_prob",
                  BAF_layer = "posterior_mtx"):
    """
    corrected BAF calling results based on RDR results.
    (rdr copy gain-remove loh)
    """
    RDR_prob_ = RDR_merge_Xdata.layers[RDR_layer].copy()
    BAF_prob_ = BAF_merge_Xdata.layers[BAF_layer].copy()

    lohprob_corrected_g = LOH_corrected(BAF_prob_, RDR_prob_, mode = 0)
    lohprob_corrected_g_merge = LOH_merge(lohprob_corrected_g)
    lohprob_corrected_l = LOH_corrected(BAF_prob_, RDR_prob_, mode = 1)
    lohprob_corrected_l_merge = LOH_merge(lohprob_corrected_l)
    
    lohprob_corrected = LOH_corrected(BAF_prob_, RDR_prob_, mode = 2)
    lohprob_corrected_merge = LOH_merge(lohprob_corrected)

    ## lohprob_corrected based on rdr information
    ## lohprob_corrected_merge for LOH evaluation
    BAF_merge_Xdata.layers["lohprob_corrected"] = lohprob_corrected
    BAF_merge_Xdata.layers["lohprob_corrected_merge"] = lohprob_corrected_merge

    BAF_merge_Xdata.layers["lohprob_corrected_g"] = lohprob_corrected_g
    BAF_merge_Xdata.layers["lohprob_corrected_g_merge"] = lohprob_corrected_g_merge

    BAF_merge_Xdata.layers["lohprob_corrected_l"] = lohprob_corrected_l
    BAF_merge_Xdata.layers["lohprob_corrected_l_merge"] = lohprob_corrected_l_merge

    return BAF_merge_Xdata

# Part III: RDR_merge_Xdata-corrected for copy loss identification
#           RDR_adata-corrected for copy loss identification
## extend Xdata from bins to genes

def bin_to_gene_mapping(BAF_merge_Xdata,
                        RDR_Xdata,
                        Xlayer = "posterior_mtx_log",
                        extend_layer = "BAF_extend_prob",
                        return_prob = False):
    """
    Function:
    ---------
    output BAF information in RDR gene resolution.
    """

    print("[XClone] BAF extend bins to genes.")
    ## check cell order
    success_flag = check_RDR_BAF_cellorder(RDR_Xdata, BAF_merge_Xdata)
    fail_flag = not success_flag
    if fail_flag:
        print("[XClone-combination]gene to bin mapping:")
        raise ValueError("pls check input: cell orders of the two Xdata are not matched.")

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
    
    if return_prob == True:
        return extend_res
    else:
        return RDR_Xdata

def BAF_extend(RDR_Xdata, BAF_merge_Xdata, 
                    extend_layer = "posterior_mtx_log",
                    out_layer = "BAF_posterior"):
    """
    extend BAF bins to genes
    1) based on RDR genes (choose this one)
    2) based on BAF original adata

    maybe todo.
    """
    pass

def copygain_corrected(combine_base_prob):
    """
    not sure if necessary;
    can use bafprob to find out allele-specific copy gain
    """
    pass

def copyloss_corrected(combine_base_prob, mode = 1):
    """
    (1) can use bafprob to find out allele-specific copy loss

    (2) adjust copy loss(RDR)-copy neutral(BAF) situation.
    """
    prob_ = combine_base_prob.copy()
    # np.argsort(prob_, axis = -1)
    
    # strategy 1: 
    if mode == 1:
        prob_[:,:,0,0] += prob_[:,:,0,1]/3
        prob_[:,:,0,2] += prob_[:,:,0,1]/3
        prob_[:,:,1,1] += prob_[:,:,0,1]/3

        prob_[:,:,0,1] = 0
    
    # strategy 2:
    if mode == 2:
        prob_[:,:,1,1] += prob_[:,:,0,1]
        prob_[:,:,0,1] = 0
    if mode == 3:
        pass

    return prob_

def CNV_prob_combination11(RDR_Xdata,
                         BAF_Xdata,
                         RDR_layer = "RDR_merge_prob",
                         BAF_layer = "posterior_mtx"):
    """
    todo

    """
    rdr_prob_ = RDR_Xdata.layers[RDR_layer].copy()
    baf_prob_ = BAF_Xdata.layers[BAF_layer].copy()

    rdr_prob_ = np.expand_dims(rdr_prob_, axis = -1)
    baf_prob_ = np.expand_dims(baf_prob_, axis = -2)

    combine_base_prob = rdr_prob_ * baf_prob_
    RDR_Xdata.layers["combine_base_prob"] = combine_base_prob
    BAF_Xdata.layers["combine_base_prob"] = combine_base_prob

    ## combination adjust

    loss_corrected_prob = copyloss_corrected(combine_base_prob)

    RDR_Xdata.layers["loss_corrected_prob"] = loss_corrected_prob
    BAF_Xdata.layers["loss_corrected_prob"] = loss_corrected_prob
    return RDR_Xdata, BAF_Xdata

def CNV_prob_combination(Xdata,
                         RDR_layer = "posterior_mtx",
                         BAF_layer = "BAF_extend_post_prob"):
    """

    """
    rdr_prob_ = Xdata.layers[RDR_layer].copy()
    baf_prob_ = Xdata.layers[BAF_layer].copy()

    rdr_prob_ = np.expand_dims(rdr_prob_, axis = -1)
    baf_prob_ = np.expand_dims(baf_prob_, axis = -2)

    ## combination base
    combine_base_prob = rdr_prob_ * baf_prob_
    Xdata.layers["combine_base_prob"] = combine_base_prob
    ## combination adjust
    ## copy loss corrected
    loss_corrected_prob1 = copyloss_corrected(combine_base_prob, mode = 1)
    Xdata.layers["loss_corrected_prob1"] = loss_corrected_prob1
    prob1_merge = CNV_prob_merge(Xdata, "loss_corrected_prob1")
    Xdata.layers["prob1_merge"] = prob1_merge

    loss_corrected_prob2 = copyloss_corrected(combine_base_prob, mode = 2)
    Xdata.layers["loss_corrected_prob2"] = loss_corrected_prob2
    prob2_merge = CNV_prob_merge(Xdata, "loss_corrected_prob2")
    Xdata.layers["prob2_merge"] = prob2_merge

    return Xdata

def CNV_prob_merge(Xdata,
                   Xlayer):
    """
    merge states prob for evaluation.
    """
    copy_loss = Xdata.layers[Xlayer][:,:,0,:].sum(axis = -1)
    loh = Xdata.layers[Xlayer][:,:,1,0] + Xdata.layers[Xlayer][:,:,1,2]
    copy_neutral = Xdata.layers[Xlayer][:,:,1,1] 
    copy_gain = Xdata.layers[Xlayer][:,:,2,:].sum(axis = -1)

    prob_merge = np.stack([copy_loss, loh, copy_neutral, copy_gain], axis = -1)
    return prob_merge

def CNV_prob_merge_for_plot(Xdata,
                            Xlayer = "loss_corrected_prob1"):
    """
    merge states prob for evaluation.
    """
    copy_loss = Xdata.layers[Xlayer][:,:,0,:].sum(axis = -1)
    loh = Xdata.layers[Xlayer][:,:,1,0] + Xdata.layers[Xlayer][:,:,1,2]
    copy_neutral = Xdata.layers[Xlayer][:,:,1,1] 
    copy_gain = Xdata.layers[Xlayer][:,:,2,:].sum(axis = -1)

    copy_loss_A = Xdata.layers[Xlayer][:,:,0,0]
    copy_loss_B = Xdata.layers[Xlayer][:,:,0,2]
    
    loh_A = Xdata.layers[Xlayer][:,:,1,0]
    loh_B = Xdata.layers[Xlayer][:,:,1,2]
    
    plot_prob_merge1 = np.stack([copy_loss, loh, copy_neutral, copy_gain], axis = -1)
    Xdata.layers["plot_prob_merge1"] = plot_prob_merge1
    plot_prob_merge2 = np.stack([copy_loss_A, copy_loss_B, loh, copy_neutral, copy_gain], axis = -1)
    Xdata.layers["plot_prob_merge2"] = plot_prob_merge2
    plot_prob_merge3 = np.stack([copy_loss_A, copy_loss_B, loh_A, loh_B, copy_neutral, copy_gain], axis = -1)
    Xdata.layers["plot_prob_merge3"] = plot_prob_merge3
    plot_prob_merge4 = np.stack([copy_loss, loh_A, loh_B, copy_neutral, copy_gain], axis = -1)
    Xdata.layers["plot_prob_merge4"] = plot_prob_merge4
    
    return Xdata

# Part IV: Combination update version. - heuristic
def CNV_states_combination(RDR_merge_Xdata,
                           BAF_merge_Xdata, 
                           RDR_layer = "RDR_merge_prob",
                           BAF_layer = "posterior_mtx"):
    """
    AIM:
    output 4 states: copy gain, copy neutral, copy loss, LOH
    output 6 states: copy gain, copy neutral, copy lossA, copy lossB, LOH1, LOH2

    stable transformation query:
    0: copy neutral
    1 2 3: copy loss
    4: LOH1
    5: LOH2
    6 7 8: copy gain
    and adjusted copy loss 3-->1 or 3-->2
    
    output layers:
    "combined_states": 0-8.
    "combined_states_adjust": 0-8 exclude 3.
    """
    #########################################
    # states_table
    # original order: copy loss, LOH1, LOH2, copy neutral, copy gain
    ## copy neutral(1)  copy loss(3)  LOH(2)  copy gain(3)
    stable = np.array([[0,0,0,1,0],
                       [1,1,0,0,0],
                       [1,0,1,0,0],
                       [1,0,0,1,0],
                       [0,1,0,1,0],
                       [0,0,1,1,0],
                       [0,1,0,0,1],
                       [0,0,1,0,1],
                       [0,0,0,1,1]])
    #########################################

    RDR_prob_ = RDR_merge_Xdata.layers[RDR_layer].copy()
    BAF_prob_ = BAF_merge_Xdata.layers[BAF_layer].copy()
    ## copy of BAF_merge_Xdata but remove layers.
    combined_Xdata = remove_Xdata_layers(BAF_merge_Xdata, copy = True)
    combined_Xdata.uns["log"] = dict([('init_data', str(combined_Xdata.shape))])
    combined_Xdata.uns["log"]["data_mode"] = "combine_X"
    
    combined_Xdata.layers["RDR_posterior_prob"] = RDR_prob_
    combined_Xdata.layers["BAF_posterior_prob"] = BAF_prob_

    def prob_insert(prob_, insert_idx1=1, insert_idx2=1):
        prob_ = np.insert(prob_, insert_idx1, values = 0,  axis=-1)
        prob_ = np.insert(prob_, insert_idx2, values = 0,  axis=-1)
        return prob_
    
    RDR_prob_ = prob_insert(RDR_prob_, 1, 1)
    BAF_prob_ = prob_insert(BAF_prob_[:,:,[0,2,1]], 0, 4)
    # order: copy loss, LOH1, LOH2, copy neutral, copy gain
    print("[XClone combination] RDR and BAF states matching.")

    prob1_ = class2onehot(RDR_prob_, 5)
    prob2_ = class2onehot(BAF_prob_, 5)

    result_ = np.logical_or(prob1_, prob2_) *np.ones_like(prob1_)
    mapping_ = np.dot(result_ , stable.T)

    valid_index = (mapping_ == 2).sum(axis=-1)
    states_index = np.argmax(mapping_, axis=-1)
    
    ## combined states
    ## order copy neutral, copy loss, LOH1, LOH2, copy gain.
    # 0-8 matched with 9 states in stable
    combined_states = valid_index * states_index
    
    combined_Xdata.layers["combined_states"] = combined_states
    print("[XClone combination] combined_states generated.")
    
    #-----------------------adjust----------------------
    combined_Xdata = copyloss_adjust(combined_Xdata,
                    in_layer = "combined_states",
                    out_layer1 = "mask_loss",
                    out_layer2 = "combined_states_adjust")
    print("[XClone combination] combined_states corrected.")
    
    return combined_Xdata

def copyloss_adjust(combined_Xdata,
                    in_layer = "combined_states",
                    out_layer1 = "mask_loss",
                    out_layer2 = "combined_states_adjust"):
    """
    Function:
    ---------
    if rdr is copy loss but BAF identified as copy neutral,
    we adjust copy loss to copy lossA or copy lossB by BAF information, 2nd argmax.

    3-->1  or 3-->2 based on the probability.
    
    Parameters:
    -----------
    RDR_prob_: copy loss, copy neutral, copy gain
    BAF_prob_: copy loss-A, copy neutral, copy loss-B
    combined_Xdata: anndata for combined items.
    in_layer: default, "combined_states".
    """
    RDR_prob_ = combined_Xdata.layers["RDR_posterior_prob"]
    BAF_prob_ = combined_Xdata.layers["BAF_posterior_prob"]

    RDR_mask = ~(np.argmax(RDR_prob_, axis = -1) == 0)
    
    BAF_2nd_index = np.argsort(BAF_prob_, axis= -1)[:,:,1].copy()
    BAF_2nd_index[BAF_2nd_index == 0] = -2
    BAF_2nd_index[BAF_2nd_index == 1] = 0
    BAF_2nd_index[BAF_2nd_index == 2] = -1
    
    BAF_2nd_index[RDR_mask] = 0
    
    combined_Xdata.layers[out_layer1] = BAF_2nd_index
    combined_Xdata.layers[out_layer2] = combined_Xdata.layers[in_layer] + BAF_2nd_index

    # print(np.unique(BAF_2nd_index)) should be array([-2, -1,  0])

    ## mask layer
    state_nochange = (combined_Xdata.layers[out_layer1] == 0).sum()
    state_change1 = (combined_Xdata.layers[out_layer1] == -1).sum()
    state_change2 = (combined_Xdata.layers[out_layer1] == -2).sum()

    check_num = state_nochange + state_change1 + state_change2
    data_size = combined_Xdata.shape[0] * combined_Xdata.shape[1]
    if check_num == data_size:
        print("[XClone] copyloss adjust done!")
        states_keep = (combined_Xdata.layers[out_layer2] == 3).sum()
        loss_transfer =  (combined_Xdata.layers[in_layer] == 3).sum() - states_keep
        lossA_add = (combined_Xdata.layers[out_layer2] == 1).sum() - (combined_Xdata.layers[in_layer] == 1).sum()
        lossB_add = (combined_Xdata.layers[out_layer2] == 2).sum() - (combined_Xdata.layers[in_layer] == 2).sum()
        if states_keep == 0:
            print("[XClone]  %d (%4f %%) copy loss states corrected." %(loss_transfer, loss_transfer/data_size * 100))
            if loss_transfer == lossA_add + lossB_add:
                print("[XClone]  %d (%.4f %%) copy loss states corrected to copy lossA; " %(lossA_add, lossA_add/data_size * 100))
                print("[XClone]  %d (%.4f %%) copy loss states corrected to copy lossB. " %(lossB_add, lossB_add/data_size * 100))
    else:
        print("[XClone] copyloss adjust error! Pls check input.")
    return combined_Xdata

# Part IV: Combination visualization.
def states_identification(combined_Xdata,
                          mode = 0,
                          in_layer = "combined_states_adjust", 
                          out_layer = "visualized_states"):
    """
    mainly for visualization.

    copy neutral; copy loss; copy gain; LOH
    
    copy neutral; copy loss; copy gain; LOH-A; LOH-B

    copy neutral; copy lossA; copy loss B; copy gain; LOH
    copy neutral; copy lossA; copy loss B; copy gain; LOH-A; LOH-B

    mode = 0: merge_copyloss = False, merge_loh = False
    mode = 1: merge_copyloss = False, merge_loh = True
    mode = 2: merge_copyloss = True, merge_loh = False
    mode = 3: merge_copyloss = True, merge_loh = True
    """
    combined_states = combined_Xdata.layers[in_layer]

    ## merge copy gain
    # combined_states = np.where(combined_states == 6, 6, combined_states)
    combined_states = np.where(combined_states == 7, 6, combined_states)
    combined_states = np.where(combined_states == 8, 6, combined_states)
    
    ## 0 copy neutral, 1,copy loss A, 2, copy loss B, 4, loh1, 5, loh2, 6, copy gain
    if mode == 0:
        ## merge_copyloss = False, merge_loh = False
        haplotyped = combined_states.copy()
        
        haplotyped = np.where(haplotyped == 4, 3, haplotyped)
        haplotyped = np.where(haplotyped == 5, 4, haplotyped)
        # copy gain
        haplotyped = np.where(haplotyped == 6, 5, haplotyped)
        combined_Xdata.layers[out_layer + str(mode) + "_haplotyped"] = haplotyped

    elif mode == 1:
        ## merge_copyloss = False, merge_loh = True
        ## copy loss allele specific
        haplotyped = combined_states.copy()
        haplotyped = np.where(haplotyped == 4, 3, haplotyped)
        haplotyped = np.where(haplotyped == 5, 3, haplotyped)
        haplotyped = np.where(haplotyped == 6, 4, haplotyped)

        combined_Xdata.layers[out_layer + str(mode) + "_haplotyped"] = haplotyped
    
    elif mode == 2:
        ## merge_copyloss = True, merge_loh = False
        haplotyped = combined_states.copy()
        # copy loss states
        haplotyped = np.where(haplotyped == 2, 1, haplotyped)
        # 2 loh states
        haplotyped = np.where(haplotyped == 4, 2, haplotyped)
        haplotyped = np.where(haplotyped == 5, 3, haplotyped)
        # copy gain states
        haplotyped = np.where(haplotyped == 6, 4, haplotyped)
        combined_Xdata.layers[out_layer + str(mode) + "_haplotyped"] = haplotyped

    elif mode == 3:
        ## merge_copyloss = True, merge_loh = True
        no_haplotyped = combined_states.copy()
        ## merged copy loss information
        no_haplotyped = np.where(no_haplotyped == 2, 1, no_haplotyped)
        ## merged loh information
        no_haplotyped = np.where(no_haplotyped == 4, 2, no_haplotyped)
        no_haplotyped = np.where(no_haplotyped == 5, 2, no_haplotyped)
        ## merged copy gain information
        no_haplotyped = np.where(no_haplotyped == 6, 3, no_haplotyped)
        combined_Xdata.layers[out_layer + str(mode) + "_no_haplotyped"] = no_haplotyped
    
    return combined_Xdata

def visualize_states_process(combined_Xdata, out_layer = "visualized_states"):
    """
    ## for visualization, be processed separately
    """
    combined_Xdata = states_identification(combined_Xdata, in_layer = "combined_states_adjust", 
                                            mode = 0, out_layer = out_layer)
    combined_Xdata = states_identification(combined_Xdata, in_layer = "combined_states_adjust",
                                            mode = 1, out_layer = out_layer)
    combined_Xdata = states_identification(combined_Xdata, in_layer = "combined_states_adjust",
                                            mode = 2, out_layer = out_layer)
    combined_Xdata = states_identification(combined_Xdata, in_layer = "combined_states_adjust",
                                            mode = 3, out_layer = out_layer)
    return combined_Xdata