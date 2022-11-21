
################################################################################
# Part II:
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
