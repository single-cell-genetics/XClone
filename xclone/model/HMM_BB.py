"""Base funcs for feature(gene/block) specific HMM model in XClone BAF module"""

# This file contains a few HMM models and its combinations
# For emmission probs - Beta binomial count model

# Author: Rongting Huang
# Date: 2022-03-28
# update: 2022-05-20

import datetime

import numpy as np
import pandas as pd
import scipy as sp

import anndata as ad
import itertools
from scipy.stats import betabinom
from scipy.special import logsumexp

## Part I: Emmission_prob based on beta binomial for HMM in BAF module
# todo concentration-gene-specific

## todo states-specific [Done]
## steps: 
## 1) get ref_BAF
## 2) Generate gene_specific_BAF matrix as `states` input of FUNC `calculate_Xemm_prob_bb`

def get_BAF_ref(Xdata, Xlayer = "fill_BAF_phased", out_anno = "ref_BAF_phased", 
                anno_key = "cell_type", ref_cell = "unclassified", clipping = False):
    """
    get ref_BAF from different BAF Xlayers.
    For example, the `phased_BAF layer` or `smoothened phased BAF layer`.
    
    get ref gene/bin BAF as gene specific copy neutral BAF states.
    previous version: use theoratical 0.5 as copy neutral BAF states.
    """
    # from scipy.special import logit, expit

    Xmtx_used = Xdata.layers[Xlayer].copy()
    # Xmtx_used = logit(Xmtx_used)

    is_ref_ = Xdata.obs[anno_key] == ref_cell
    ref_ = Xmtx_used[is_ref_].mean(axis = 0)
    
    ## do clipping at ref BAF
    if clipping:
        print("do clipping at ref BAF")
        qt_10 = np.quantile(ref_, 0.1)
        qt_90 = np.quantile(ref_, 0.9)
        bd_low = max(qt_10, 0.3)
        bd_high = min(qt_90, 0.7)
        ref_ = np.clip(ref_, bd_low, bd_high)

    # Xdata.layers[out_layer] = expit(Xmtx_used - ref_)
    Xdata.var[out_anno] = ref_

    return Xdata

def gene_specific_BAF(Xdata,
                      theo_normal_states = 0.5,
                      theo_states= np.array([0.01, 0.99, 1/3, 2/3]), 
                      specific_BAF = "ref_BAF_phased", 
                      rescale = False):
    """
    np.array([0.01, 0.99, 0.5, 1/3, 2/3]) ->
    np.array([[0.01, 0.99, specific-BAF, 1/3, 2/3],
              [0.01, 0.99, specific-BAF, 1/3, 2/3]])

    support 3 states: copy loss and copy neutral.
    np.array([0.01, 0.99, 0.5]) ->
    np.array([[0.01, 0.99, specific-BAF],
              [0.01, 0.99, specific-BAF]])

    and support rescale.
    """
    from scipy.special import logit, expit

    specific_states = Xdata.var[specific_BAF]
    gene_num = Xdata.shape[1]
    # base_states = np.tile([0.01, 0.99, 1/3, 2/3], gene_num).reshape(gene_num, -1)
    base_states = np.tile(theo_states, gene_num).reshape(gene_num, -1)
    if len(theo_states) == 4:
        used_specific_states = np.insert(base_states, 2, values = specific_states, axis=1)
        if rescale:
            delta1 = logit(specific_states) - logit(theo_normal_states)
            used_specific_states[:,0] = expit(logit(theo_states[0])+ delta1)
            used_specific_states[:,1] = expit(logit(theo_states[1])+ delta1)
            used_specific_states[:,3] = expit(logit(theo_states[2])+ delta1)
            used_specific_states[:,4] = expit(logit(theo_states[3])+ delta1)
    elif len(theo_states) == 2:
        used_specific_states = np.insert(base_states, 1, values = specific_states, axis=1)
        if rescale:
            delta1 = logit(specific_states) - logit(theo_normal_states)
            used_specific_states[:,0] = expit(logit(theo_states[0])+ delta1)
            used_specific_states[:,2] = expit(logit(theo_states[1])+ delta1)
    return used_specific_states

## todo maybe need remove ref BAF
def specific_BAF(Xdata, specific_BAF_layer = "reref_BAF"):
    """
    deprecated
    can not be used. not normal state.
    """
    specific_states = np.expand_dims(Xdata.layers[specific_BAF_layer], axis=2)
    
    used_specific_states = np.insert(specific_states, 1, 2/3, axis=2)
    used_specific_states = np.insert(used_specific_states, 1, 1/3, axis=2)
    used_specific_states = np.insert(used_specific_states, 0, 0.99, axis=2)
    used_specific_states = np.insert(used_specific_states, 0, 0.01, axis=2)
    return used_specific_states


def generate_bb_logprob(AD, DP,
                        states = None,
                        concentration = 10,
                        gene_specific = False,
                        verbose = False):
    """
    Function: 
    beta binomial distributions[scipy]

    https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.betabinom.html
    https://en.wikipedia.org/wiki/Beta-binomial_distribution
    https://search.r-project.org/CRAN/refmans/aod/html/betabin.html
    https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.nbinom.html

    states here is (theoretical) BAF value for default.
    overdispersion is gene-specific and pre-learned from GLM fitting.

    
    Parameters
    ----------
    AD: np.array.
    DP: np.array.
    states: default np.array([0.01, 0.99, 0.5, 1/3, 2/3]) 
    # example for copy loss, neutral and copy gain
    can apply other value and can apply gene-specific states from FUNC
    gene_specific_BAF

    concentration:
    overdispersion:
    gene_specific:

    Returns
    -------
    emm_prob_log: betabinom.logpmf, should be in `log` format.
    """
    if states is None:
        states = np.array([0.01, 1/3, 0.5, 2/3, 0.99])
        states = states.reshape(1, 1, -1)
    elif states.ndim == 1:
        states = states.reshape(1, 1, -1)
    else:
        print("[XClone] specific Center states used.")

    # mu_ = states * np.ones(AD.shape)[:,:,np.newaxis]
    mu_ = np.expand_dims(np.ones(AD.shape), axis=2) * states
    
    if gene_specific:
        print("[XClone hint] specific concentrations used.")

    # obs_ = np.repeat(AD, states.shape[-1]).reshape(AD.shape[0], -1, states.shape[-1])
    # _bb_total = np.repeat(DP, states.shape[-1]).reshape(DP.shape[0], -1, states.shape[-1])

    obs_ = np.expand_dims(AD, axis=2)
    _bb_total = np.expand_dims(DP, axis=2)

    _bb_prob = mu_
    _bb_phi = 1/(concentration + 1)

    _bb_alpha = _bb_prob*(1/_bb_phi - 1)
    _bb_beta =  _bb_alpha*(1/_bb_prob - 1)
    
    if verbose:
        print("_bb_alpha", _bb_alpha)
        print("_bb_beta", _bb_beta)

    emm_prob_log = betabinom.logpmf(obs_, _bb_total, _bb_alpha, _bb_beta, loc=0)

    return emm_prob_log

def process_BAF_emm_prob(emm_prob_log):
    """
    copy gain for selecting the fitted type.
    ##todo
    ## 1) merge copy gain status? or just processing in HMM
    ## 2) save all status

    ## not in used now.
    """
    emm_prob_log_update = emm_prob_log.copy()
    emm_prob_log_update = emm_prob_log_update[:,:,0:3]

    copy_gain_ = emm_prob_log[:,:,2:].max(axis=-1, keepdims = True)

    emm_prob_log_update[:,:,2:] = copy_gain_

    return emm_prob_log_update

def calculate_Xemm_prob_bb(Xdata,
                           AD_key = "AD",
                           DP_key = "DP",
                           outlayer = "BAF_emm_prob_log",
                           states = None,
                           states_num = 5,
                           gene_specific = False,
                           concentration = 10,
                           verbose = False):
    """
    Function:
    Usage of FUNC `generate_bb_logprob`.


    Parameters:
    ----------
    Xdata: anndata.
    states:

    Return:
    ------

    Example:
    -------

    """
    # time stamp
    start_time_ = datetime.datetime.now()

    ## prepare the input from Xdata for generate_nb_logprob
    if sp.sparse.issparse(Xdata.layers[AD_key]):
        AD = Xdata.layers[AD_key].A
    else:
        AD = Xdata.layers[AD_key]
    
    if sp.sparse.issparse(Xdata.layers[DP_key]):
        DP = Xdata.layers[DP_key].A
    else:
        DP = Xdata.layers[DP_key]
    
    if states is None:
        if states_num == 5:
            states = np.array([0.01, 1/3, 0.5, 2/3, 0.99])
        elif states_num == 3:
            states = np.array([0.01, 0.5, 0.99])

    ## check the shape for the input
    if gene_specific == True:
        if len(Xdata.var["concentration"]) != AD.shape[1]:
            raise ValueError("[XClone]gene-speciific concentration doesn't match with the gene numbers! Pls check!")
        concentration = Xdata.var["concentration"].values.reshape(1, -1, 1)
    elif concentration is None:
        concentration = 10
 
    ## calculate emm_prob_log
    if states.ndim == 1:
        print("states used:", states)
    else:
        head, *tail = states
        print("states used:", head)
        print(".....")
    emm_prob_log = generate_bb_logprob( AD, DP,
                                        states = states,
                                        gene_specific = gene_specific,
                                        concentration = concentration
                                        )
    emm_prob_log = validating_prob(emm_prob_log)
    # normalised loglikelihood (dosen't affect HMM)
    # emm_prob_log += -logsumexp(emm_prob_log, axis=2, keepdims=True)

    if verbose:
        print("states", states)
        print("emm_prob_log", emm_prob_log)
        print("emm_prob_log shape", emm_prob_log.shape)
    
    Xdata.layers[outlayer] = emm_prob_log
    # emm_prob_log_update = process_BAF_emm_prob(emm_prob_log)
    # Xdata.uns["BAF_emm_prob"] = emm_prob_log_update

    # time stamp
    end_time_= datetime.datetime.now()
    time_used = end_time_ - start_time_
    print("cal emm prob time", time_used.seconds, "seconds")

    return Xdata

## Part II: emm_prob_log processing

def validating_prob(pre_prob, replace_value = np.log(1e-200)):
    """
    Mainly for inf situation.
    Filter cells/genes before KNN/WMA/HMM smoothing.
    
    replace_value：ref min, np.log(1e-320)
    -----
    Notes:   
    https://www.cnblogs.com/herbert/p/3402245.html
    float: 3.4E-38～3.4E+38
    double: 1.7E-308～1.7E+308
    """

    emm_prob_log = pre_prob.copy()
    is_valid_ = np.isfinite(emm_prob_log).all()
    
    if is_valid_:
        print("[XClone]: validated probability, all finite.")
        return emm_prob_log
    else:
        print("[XClone]: not validated probability, infinite exists.")

    cond_ = ~np.isfinite(emm_prob_log)

    idx_ = np.where(cond_)
    emm_prob_log[idx_] = replace_value
    print("[XClone]: replacing the infinite value as ", replace_value)
    print("max and min of prob", emm_prob_log.max(), emm_prob_log.min())
    return emm_prob_log