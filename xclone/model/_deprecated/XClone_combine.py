# gene_to_bin RDR strategy

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

    lohprob_merge = normalize(lohprob_merge)

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


# Part IV: Combination update version. - heuristic

from ..utils import class2onehot
from ...plot.CNV_plot import remove_Xdata_layers
def CNV_states_combination(RDR_merge_Xdata,
                           BAF_merge_Xdata, 
                           RDR_layer = "RDR_merge_prob",
                           BAF_layer = "posterior_mtx"):
    """
    init version but deprecated.
    Quite interesting in the stable construction part.[record]
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


################################################################################
# testing. deprecated.
## todo error in prob processing.
def BAF_prob_processing(Xdata, 
                        haplotype = True, 
                        in_layer = "posterior_mtx",
                        out_layer = "haplotyped_posterior"):
    """
    deprecated
    Xdata: anndata. BAF_Xdata with CNV calling results.
    haplotype: bool.
    in_layer: probability layer, states order: copy loss, neutral,
              copy loss
    out_layer: processed probability.

    """
    ## process the probability, *not the log format.
    prob_ = Xdata.layers[in_layer].copy()
    if haplotype == True:
        ## support 2 copy loss states
        prob_ = np.insert(prob_, 1, values = prob_[:,:,2], axis=2)
        # ## copy loss1
        # prob_log[:,:,0] = prob_log[:,:,0] - prob_log[:,:,2]
        # ## copy loss2
        # prob_log[:,:,1] = prob_log[:,:,1] - prob_log[:,:,2]
        # ## copy neutral
        # prob_log[:,:,2] = prob_log[:,:,2] - prob_log[:,:,2]
        ## copy gain
        # prob_log[:,:,3] = 0
        prob_[:,:,3] = prob_[:,:,2]

    else:
        ## merge copy loss states
        prob_[:,:,0] = prob_[:,:,0] + prob_[:,:,2]
        ## ratio-compare with neutral state
        ## copy loss
        # prob_log[:,:,0] = prob_log[:,:,0] - prob_log[:,:,1]
        # ## copy neutral
        # prob_log[:,:,1] = prob_log[:,:,1] - prob_log[:,:,1]
        # ## copy gain
        # prob_log[:,:,2] = 0

        ## copy gain
        prob_[:,:,2] = prob_[:,:,1]

    ## normalize
    # prob_log += -logsumexp(prob_log, axis=2, keepdims=True)
    prob_ = normalize(prob_)


    Xdata.layers[out_layer] = prob_
    return Xdata

from .base_utils import normalize
def CNV_states_combination_X(RDR_merge_Xdata,
                           BAF_merge_Xdata, 
                           RDR_layer = "RDR_merge_prob",
                           BAF_layer = "posterior",
                           out_layer = "haplotyped",
                           haplotype= False):
    """
    deprecated
    Two format

    Support 3 states: copy loss, copy neutral, copy gain
    Support 4 states: copy loss1, copy loss2, copy loss, copy neutral, copy gain
    """
    RDR_prob_ = RDR_merge_Xdata.layers[RDR_layer].copy()
    BAF_prob_ = BAF_merge_Xdata.layers[BAF_layer].copy()
    
    if haplotype == True:
        ## support 4 states
        RDR_prob_ = np.insert(RDR_prob_, 1, values = RDR_prob_[:,:,0], axis=2)
        # ## normalize
        # RDR_prob_log += -logsumexp(RDR_prob_log, axis=2, keepdims=True)
        
        # prob_log = np.exp(RDR_prob_log) + np.exp(BAF_prob_log)
        prob_ = RDR_prob_ + BAF_prob_
    else:
        ## support 3 states
        # prob_log = np.exp(RDR_prob_log) + np.exp(BAF_prob_log)
        prob_ = RDR_prob_ + BAF_prob_
    
    # prob_log += -logsumexp(prob_log, axis=2, keepdims=True)
    
    # BAF_merge_Xdata.layers[out_layer] = np.exp(prob_log)
    BAF_merge_Xdata.layers[out_layer] = normalize(prob_log)

    return BAF_merge_Xdata

import numpy as np
## deprecated record
def CNV_states_combination(RDR_merge_Xdata,
                           BAF_merge_Xdata, 
                           RDR_layer = "RDR_merge_prob",
                           BAF_layer = "posterior",
                           out_layer = "combined_states",
                           states_table = 1):
    """
    output 4 states: copy gain, copy neutral, copy loss, LOH
    output 5 states: copy gain, copy neutral, copy loss, LOH1, LOH2
    """
    #########################################
    # states_table = 0: ignore some RDR copy gain and copy loss states; deprecated.
    # states_table = 1: more reasonable, default one.
    #########################################
    # states_table1
    ## copy neutral(1) copy loss(2)  LOH(2)  copy gain(2)
    states_table1 = np.array([[0,0,0,1,0],
                            [1,1,0,0,0],
                            [1,0,1,0,0],
                            [0,1,0,1,0],
                            [0,0,1,1,0],
                            [0,1,0,0,1],
                            [0,0,1,0,1]])
    # states_table2
    ## copy neutral(1)  copy loss(3)  LOH(2)  copy gain(3)
    states_table2 = np.array([[0,0,0,1,0],
                            [1,1,0,0,0],
                            [1,0,1,0,0],
                            [1,0,0,1,0],
                            [0,1,0,1,0],
                            [0,0,1,1,0],
                            [0,1,0,0,1],
                            [0,0,1,0,1],
                            [0,0,0,1,1]])
    if states_table == 0:
        stable = states_table1
    elif states_table == 1:
        stable = states_table2

    #########################################

    RDR_prob_ = RDR_merge_Xdata.layers[RDR_layer].copy()
    BAF_prob_ = BAF_merge_Xdata.layers[BAF_layer].copy()

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
    combined_states = valid_index * states_index
    if states_table == 0:
        ## 0 copy neutral, 1,2 copy loss, 3, loh1, 4, loh2, 5, 6 copy gain
        haplotyped = combined_states.copy()
        
        haplotyped = np.where(haplotyped == 2, 1, haplotyped)
        haplotyped = np.where(haplotyped == 3, 2, haplotyped)
        haplotyped = np.where(haplotyped == 4, 3, haplotyped)
        haplotyped = np.where(haplotyped == 5, 4, haplotyped)
        haplotyped = np.where(haplotyped == 6, 4, haplotyped)

        no_haplotyped = combined_states.copy()
        no_haplotyped = np.where(no_haplotyped == 2, 1, no_haplotyped)
        no_haplotyped = np.where(no_haplotyped == 3, 2, no_haplotyped)
        no_haplotyped = np.where(no_haplotyped == 4, 2, no_haplotyped)
        no_haplotyped = np.where(no_haplotyped == 5, 3, no_haplotyped)
        no_haplotyped = np.where(no_haplotyped == 6, 3, no_haplotyped)

    elif states_table == 1:
        ## 0 copy neutral, 1,2,3 copy loss, 4, loh1, 5, loh2, 6, 7, 8 copy gain
        haplotyped = combined_states.copy()
        
        haplotyped = np.where(haplotyped == 2, 1, haplotyped)
        haplotyped = np.where(haplotyped == 3, 1, haplotyped)
        haplotyped = np.where(haplotyped == 4, 2, haplotyped)
        haplotyped = np.where(haplotyped == 5, 3, haplotyped)
        haplotyped = np.where(haplotyped == 6, 4, haplotyped)
        haplotyped = np.where(haplotyped == 7, 4, haplotyped)
        haplotyped = np.where(haplotyped == 8, 4, haplotyped)

        no_haplotyped = combined_states.copy()
        no_haplotyped = np.where(no_haplotyped == 2, 1, no_haplotyped)
        no_haplotyped = np.where(no_haplotyped == 3, 1, no_haplotyped)
        no_haplotyped = np.where(no_haplotyped == 4, 2, no_haplotyped)
        no_haplotyped = np.where(no_haplotyped == 5, 2, no_haplotyped)
        no_haplotyped = np.where(no_haplotyped == 6, 3, no_haplotyped)
        no_haplotyped = np.where(no_haplotyped == 7, 3, no_haplotyped)
        no_haplotyped = np.where(no_haplotyped == 8, 3, no_haplotyped)

    BAF_merge_Xdata.layers[out_layer] = no_haplotyped
    BAF_merge_Xdata.layers[out_layer + "_haplotyped"] = haplotyped
    return BAF_merge_Xdata
