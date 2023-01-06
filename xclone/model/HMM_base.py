"""Base funcs for feature(gene/block) specific HMM model in XClone."""

# This file contains a few HMM models and its combinations
# for different emmission probs
# Author: Rongting Huang
# Date: 30/07/2021
# update: 19/01/2022

import numpy as np
import datetime
from scipy.special import logsumexp
from .base_utils import normalize

## Part I: base functions for HMM smoothing in XClone
## Forward-backward algorithm
def fwd_bkw_prob1(obs_num, states_num, start_prob=None, trans_prob=None, emm_prob=None):
    """Forward–backward algorithm.
    fwd_bkw_prob: version 0.0.1
    try to fix Numerical issue by normalization

    # Example1:
    # --------
    # reference = np.array([10, 30, 50])
    # observations = np.array([8, 80, 50])
    # states = np.array([0.5, 1.0, 1.5])
    # start_prob = np.array([0.15, 0.7, 0.15])
    # trans_prob = np.array([[0.7, 0.15, 0.15],
    #                       [0.15, 0.7, 0.15],
    #                       [0.15, 0.15, 0.7]])
    # emm_prob = generate_nb_emmission_prob(reference, observations, states)
    # states_num = len(states)
    # obs_num = len(observations)
    # fwd_bkw_prob(obs_num, states_num, start_prob, trans_prob, emm_prob)

    Example2:
    reference = np.array([10, 10, 10,10,10, 10,10])
    # reference = np.array([5., 9, 3, 5, 12, 6, 9])
    observations = np.array([5., 9, 3, 5, 12, 6, 9])
    states = np.array([0.5, 1.0])
    start_prob = np.array([0.8, 0.2])
    trans_prob = np.array([[0.7,0.3],[0.2,0.8]])

    emm_prob = generate_nb_emmission_prob(reference, observations,states)
    # emm_prob = generate_poisson_emmission_prob(reference, observations,states)
    states_num = len(states)
    obs_num = len(observations)
    fwd_bkw_prob(obs_num, states_num, start_prob, trans_prob, emm_prob)
    """
    # Forward part of the algorithm
    # obs_num = len(observations)
    # states_num = len(states)
    # fwd_mtx = np.zeros((obs_num, states_num))
    for g in range(obs_num):
        # f_curr = np.zeros((states_num))
        if g == 0:
            prev_f_sum = start_prob
        else:
            prev_f_sum = np.dot(f_prev, trans_prob)        
        
        f_curr = emm_prob[g, :] * prev_f_sum
        f_curr = normalize(f_curr)
        # print("f_curr:", f_curr)

        if g == 0:
            fwd_mtx = f_curr
        else:
            fwd_mtx = np.vstack((fwd_mtx, f_curr))
        
        f_prev = f_curr.copy()

    # Backward part of the algorithm
    bkw = []
    for g in list(reversed(range(obs_num))):
        # b_curr = np.zeros((states_num))
        if g == obs_num-1:
            # base case for backward part
            b_curr = np.ones((states_num))
        else:
            b_curr = np.dot(trans_prob,(emm_prob[g+1,:] * b_prev).T)
        # print("normalize bkw") # actually, no need to normalize backward part
        b_curr = normalize(b_curr) # but add normalize can fix value problem, and not affect results
        bkw.insert(0,b_curr)
        b_prev = b_curr.copy()
    
    for cnt, b_curr in enumerate(bkw):
        if cnt == 0:
            bkw_mtx = b_curr
        else:
            bkw_mtx = np.vstack((bkw_mtx, b_curr))
    
    ## get two matrix, fwd_mtx and bkw_mtx: genes*states matrix
    ## Merging the two parts
    # print("fwd_mtx:", fwd_mtx)
    # print("bkw_mtx:", bkw_mtx)
    posterior_pre_mtx = fwd_mtx * bkw_mtx
    posterior_sum = posterior_pre_mtx.sum(axis=1).reshape(-1,1)
    posterior_mtx = posterior_pre_mtx/posterior_sum
    # print("posterior_pre_mtx:", posterior_pre_mtx)
    # print('posterior_mtx:', posterior_mtx)
    ## return genes*states matrix
    return posterior_mtx, fwd_mtx, bkw_mtx


## Forward-backward algorithm
def fwd_bkw_prob2(start_prob=None, trans_prob=None, emm_prob_log=None):
    """Forward–backward algorithm.

    fwd_bkw_prob: version 0.0.2
    try to fix Numerical issue by storing the log value, stable version

    # Example1:
    # --------
    # reference = np.array([10, 30, 50])
    # observations = np.array([8, 80, 50])
    # states = np.array([0.5, 1.0, 1.5])
    # start_prob = np.array([0.15, 0.7, 0.15])
    # trans_prob = np.array([[0.7, 0.15, 0.15],
    #                       [0.15, 0.7, 0.15],
    #                       [0.15, 0.15, 0.7]])
    # emm_prob_log = generate_nb_emmission_logprob(reference, observations, states)
    # states_num = len(states)
    # obs_num = len(observations)
    # fwd_bkw_prob2(start_prob, trans_prob, emm_prob_log)

    Example2:
    reference = np.array([10, 10, 10,10,10, 10,10])
    # reference = np.array([5., 9, 3, 5, 12, 6, 9])
    observations = np.array([5., 9, 3, 5, 12, 6, 9])
    states = np.array([0.5, 1.0])
    start_prob = np.array([0.8, 0.2])
    trans_prob = np.array([[0.7,0.3],[0.2,0.8]])

    emm_prob_log = generate_nb_emmission_logprob(reference, observations,states)
    # emm_prob_log = generate_poisson_emmission_logprob(reference, observations,states)
    # states_num = len(states)
    # obs_num = len(observations)
    fwd_bkw_prob2(start_prob, trans_prob, emm_prob_log)
    """
    obs_num = emm_prob_log.shape[0]
    states_num = emm_prob_log.shape[1]

    # Forward part of the algorithm
    trans_prob_log = np.log(trans_prob)
    start_prob_log = np.log(start_prob)
    # emm_prob_log = np.log(emm_prob)
    for g in range(obs_num):
        if g == 0:
            prev_f_sum_log = start_prob_log
        else:
            prev_f_sum_log = logsumexp(f_prev_log.reshape(-1,1) + trans_prob_log, axis=0)
        
        f_curr_log = emm_prob_log[g, :] + prev_f_sum_log
        # normalize
        f_curr_log = f_curr_log - logsumexp(f_curr_log)
        # f_curr_log = np.log(normalize(np.exp(loglik_amplify(f_curr_log, axis=-1)))) ## numerical issue still

        if g == 0:
            fwd_mtx_log = f_curr_log
        else:
            fwd_mtx_log = np.vstack((fwd_mtx_log, f_curr_log))
        
        f_prev_log = f_curr_log.copy()

    # Backward part of the algorithm
    bkw_log = []
    for g in list(reversed(range(obs_num))):
        if g == obs_num-1:
            # base case for backward part
            b_curr_log = np.log(np.ones((states_num)))
        else:
            b_curr_log = logsumexp(trans_prob_log + (emm_prob_log[g+1, :] + b_prev_log).reshape(1, -1), axis=1)
        
        # b_curr_log = b_curr_log - logsumexp(b_curr_log) # no need to normalize
        bkw_log.insert(0,b_curr_log)
        b_prev_log = b_curr_log.copy()
        
    for cnt, b_curr_log in enumerate(bkw_log):
        if cnt == 0:
            bkw_mtx_log = b_curr_log
        else:
            bkw_mtx_log = np.vstack((bkw_mtx_log, b_curr_log))
    
    ## get two matrix, fwd_mtx and bkw_mtx: genes*states matrix
    ## Merging the two parts
    # print("fwd_mtx_log:", fwd_mtx_log)
    # print("exp fwd_mtx_log:", np.exp(fwd_mtx_log))
    # print("bkw_mtx_log:", bkw_mtx_log)
    # print("exp bkw_mtx_log:", np.exp(bkw_mtx_log))

    posterior_pre_mtx_log = fwd_mtx_log + bkw_mtx_log
    posterior_mtx_log = posterior_pre_mtx_log - logsumexp(posterior_pre_mtx_log, axis=-1, keepdims=True)
    # print("posterior_mtx_log:", posterior_mtx_log)
    posterior_mtx = np.exp(posterior_mtx_log)
    # print('posterior_mtx:', posterior_mtx)
    ## return genes*states matrix
    return posterior_mtx, posterior_mtx_log, fwd_mtx_log, bkw_mtx_log
 
def fwd_bkw_prob3(start_prob=None, trans_prob=None, emm_prob_log=None):
    """Forward–backward algorithm.

    fwd_bkw_prob: version 0.0.3 [to be tested]
    1)try to fix Numerical issue by storing the log value, stable version
    2)improve efficiency by applying to matrix.
      suitable for matrix calculation, broadcast to be fitted across cells.

    Parameters:
    ----------
    start_prob: np.array, states shape.
    trans_prob: np.array, states*states matrix.
    emm_prob_log: np.array, cells*genes*states shape.

    Return:
    -------
    # return cells*genes*states array
    posterior_mtx: 
    posterior_mtx_log:
    fwd_mtx_log:
    bkw_mtx_log:

    Example:
    -------
    reference = np.array([10, 10, 10,10,10, 10,10])
    # reference = np.array([5., 9, 3, 5, 12, 6, 9])
    observations = np.array([5., 9, 3, 5, 12, 6, 9])
    states = np.array([0.5, 1.0])
    start_prob = np.array([0.8, 0.2])
    trans_prob = np.array([[0.7,0.3],[0.2,0.8]])

    emm_prob_log = generate_nb_emmission_logprob(reference, observations, states)
    fwd_bkw_prob3(start_prob, trans_prob, emm_prob_log)
    """
    sample_num = emm_prob_log.shape[0] ## cells in XClone   
    obs_num = emm_prob_log.shape[-2] ## genes in XClone
    states_num = emm_prob_log.shape[-1]

    # Forward part of the algorithm
    trans_prob_log = np.log(trans_prob)
    start_prob_log = np.log(start_prob)

    for g in range(obs_num):
        if g == 0:
            prev_f_sum_log = start_prob_log
        else:
            prev_f_sum_log = logsumexp(f_prev_log.swapaxes(1,2) + trans_prob_log, axis=1)
            ## swap genes axis and states axis for calculation
        
        f_curr_log = emm_prob_log[:, g, :] + prev_f_sum_log
        # normalize
        f_curr_log = f_curr_log - logsumexp(f_curr_log)
        # keepdims
        f_curr_log = f_curr_log[:, np.newaxis, :]

        if g == 0:
            fwd_mtx_log = f_curr_log
        else:
            fwd_mtx_log = np.concatenate((fwd_mtx_log, f_curr_log), axis=1)
        
        f_prev_log = f_curr_log.copy()

    # Backward part of the algorithm
    bkw_log = []
    for g in list(reversed(range(obs_num))):
        if g == obs_num-1:
            # base case for backward part
            b_curr_log = np.log(np.ones((sample_num, states_num)))
        else:
            b_curr_log = logsumexp(trans_prob_log + emm_prob_log[:, g+1, :][:,np.newaxis, :] + b_prev_log, axis=-1)
        b_curr_log = b_curr_log[:, np.newaxis, :]
        
        bkw_log.insert(0,b_curr_log)
        b_prev_log = b_curr_log.copy()
        
    for cnt, b_curr_log in enumerate(bkw_log):
        if cnt == 0:
            bkw_mtx_log = b_curr_log
        else:
            bkw_mtx_log = np.concatenate((bkw_mtx_log, b_curr_log), axis=1)
            
    ## get two matrix, fwd_mtx and bkw_mtx: cells*genes*states matrix
    ## Merging the two parts
    posterior_pre_mtx_log = fwd_mtx_log + bkw_mtx_log
    posterior_mtx_log = posterior_pre_mtx_log - logsumexp(posterior_pre_mtx_log, axis=-1, keepdims=True)
    posterior_mtx = np.exp(posterior_mtx_log)

    return posterior_mtx, posterior_mtx_log, fwd_mtx_log, bkw_mtx_log

def fwd_bkw_prob_base(start_prob=None, trans_prob=None, emm_prob_log=None, verbose = False):
    emm_ndim = emm_prob_log.ndim
    # dim=2 means one sample/cell
    # dim=3 means many samples/cells
    if emm_ndim == 2:
        if verbose:
            print("[XClone] fwd_bkw_prob_base func, emm_prob dim=2")
        res = fwd_bkw_prob2(start_prob, trans_prob, emm_prob_log)
    elif emm_ndim == 3:
        if verbose:
            print("[XClone] fwd_bkw_prob_base func, emm_prob dim=3")
        res = fwd_bkw_prob3(start_prob, trans_prob, emm_prob_log)
    return res

########################################################################### end of the base HMM functions
## HMM smoothing
## Part II: emm_prob_log processing

def processing_prob_bygene(emm_prob_log, Xdata):
    """
    Based on some Selection criteria.
    mainly for nan situation.
    Filter genes before HMM smoothing.
    """
    idx_1 = np.where(emm_prob_log == emm_prob_log)
    gene_idx = np.unique(idx_1[1])
    emm_prob_log = emm_prob_log[:,gene_idx,:]

    update_Xdata = Xdata[:,gene_idx].copy()
    
    if len(gene_idx) < Xdata.shape[1]:
        print("Gene level: filter nan emm_prob")
    else:
        print("Gene level: no filtering emm_prob")
    return emm_prob_log, update_Xdata

def processing_prob_bycell(emm_prob_log, Xdata):
    """
    Based on some Selection criteria.
    mainly for nan situation.
    Filter cells before HMM smoothing.
    """
    idx_1 = np.where(emm_prob_log == emm_prob_log)
    cell_idx = np.unique(idx_1[0])
    emm_prob_log = emm_prob_log[cell_idx,:,:]

    update_Xdata = Xdata[cell_idx,:].copy()
    
    if len(cell_idx) < Xdata.shape[0]:
        print("Cell level: filter nan emm_prob")
        print("[XClone] warning: filter cells!")## add warning
    else:
        print("Cell level: no filtering emm_prob")
    return emm_prob_log, update_Xdata

## Part III: HMM smoothing in XClone

import multiprocessing

t = 1e-6
def XC_HMM_base(emm_prob_log, 
                start_prob = np.array([0.1, 0.8, 0.1]), 
                trans_prob = np.array([[1-2*t, t, t],[t, 1-2*t, t],[t, t, 1-2*t]]),
                verbose = True):
    """
    Function:

    Parameters:
    ----------
    Default setting for states = np.array([0.5, 1.0, 1.5]).

    Return:
    ------
    
    Example:
    -------
    reference = np.array([10, 10, 10,10,10, 10,10])
    observations = np.array([5., 9, 3, 5, 12, 6, 9])
    """

    ## run fwd_bwd_prob
    res = fwd_bkw_prob_base(start_prob, trans_prob, emm_prob_log, verbose)
    return res

def show_progress(RV=None):
    return RV      

def brk_HMM_base(emm_prob_log, start_prob, trans_prob, obs_Xdata, brk, verbose = True, nproc=1):
    """
    Function:
    version 0.0.2
    (1) improve the calculation efficiency by applying fwd_bkw at multi-cells scale.
    (2) multiprocessing in chr brk items.

    Parameters:
    ----------
    emm_prob_log: np.array.
    start_prob: np.array.
    trans_prob: np.array.

    Return:
    ------

    Example:
    -------

    """
    brk_item = obs_Xdata.var[brk].drop_duplicates(keep="first")
    if verbose:
        print(brk_item)
    
    if nproc == 1:
        # loop for each brk chrs/chr_arms
        cnt = 0
        for brk_ in brk_item:
            tmp_region_flag = obs_Xdata.var[brk] == brk_
            tmp_res = XC_HMM_base(emm_prob_log[:,tmp_region_flag,:], start_prob, trans_prob, verbose)
            if  cnt == 0:
                res = tmp_res[0]
                res_log = tmp_res[1]
            else:
                res = np.concatenate((res, tmp_res[0]), axis = 1) # merge the genes axis
                res_log = np.concatenate((res_log, tmp_res[1]), axis = 1)
            cnt+=1
    
    elif nproc > 1:
        print("[XClone] multiprocessing for each brk item")
        print("nproc:", nproc)
        pool = multiprocessing.Pool(processes=nproc)
        
        result = []
        idx_kept = []
        for brk_ in brk_item:
            tmp_region_flag = obs_Xdata.var[brk] == brk_
            idx_kept.append(brk_)

            result.append(pool.apply_async(
                XC_HMM_base, (emm_prob_log[:,tmp_region_flag,:], start_prob, trans_prob, verbose),
                callback=show_progress
            ))
        
        pool.close()
        pool.join()

        ## processing result
        # result = [res.get() for res in result]
        cnt = 0
        for tmp_res in result:
            if cnt == 0:
                res = tmp_res.get()[0]
                res_log = tmp_res.get()[1]
            else:
                res = np.concatenate((res, tmp_res.get()[0]), axis = 1) # merge the genes axis
                res_log = np.concatenate((res_log, tmp_res.get()[1]), axis = 1)
            cnt+=1

        if verbose:
            print(idx_kept)

    return res, res_log

def XHMM_smoothing(Xdata,
                   brk = "chr_arm",
                   emm_inlayer = None,
                   start_prob = None,
                   trans_prob = None,
                   verbose = True, 
                   nproc = 1,
                   **kwargs):
    """
    Function:
    Apply HMM smoothing.

    Parameters:
    ----------
    Xdata: anndata. ref_obs bulk counts or ref_obs raw counts.
    Xdata here should be the same one with calculate_Xemm_prob.
    emm_prob_log: can be calculated or just input from any methods.(NB/VB)

    **kwargs: 
    # emm_prob_log = cal_Xemm_prob(Xdata, states = None, 
    #           gene_specific = False, overdispersion = None,
    #           ref_normalization_term = None,
    #           obs_normalization_term = None,
    #           verbose = True)

    Return:
    ------
    
    Example:
    --------
    """
    ## time stamp
    start_time_ = datetime.datetime.now()

    # check emm_prob_log
    if emm_inlayer is None:
        raise ValueError("[XClone]emm_prob_log doesn't exist. Pls check!")
    else:
        emm_prob_log = Xdata.layers[emm_inlayer]
        # check emm_prob_log format with Xdata
        cell_dim_flag = emm_prob_log.shape[0] == Xdata.shape[0]
        feature_dim_flag = emm_prob_log.shape[1] == Xdata.shape[1]
        if cell_dim_flag and feature_dim_flag:
            pass
        else:
            raise ValueError("[XClone]emm_prob_log doesn't match with the Xdata dims! Pls check!")

    if start_prob is None:
        start_prob = np.array([0.1, 0.8, 0.1])
    if trans_prob is None:
        trans_prob = np.array([[1-2*t, t, t],[t, 1-2*t, t],[t, t, 1-2*t]])

    ## filter nan emm_prob_log and update the Xdata
    emm_prob_log, update_Xdata = processing_prob_bycell(emm_prob_log, Xdata) # more strict
    emm_prob_log, update_Xdata = processing_prob_bygene(emm_prob_log, update_Xdata)

    update_Xdata.layers["emm_prob_log_noHMM"] = emm_prob_log
    update_Xdata.layers["emm_prob_noHMM"] = np.exp(emm_prob_log)

    ## HMM smoothing for each brk block-XC_HMM_base
    ## here todo  check  both celltype and cellbased!
    res, res_log = brk_HMM_base(emm_prob_log, start_prob, trans_prob, update_Xdata, brk, verbose, nproc) 
    
    ## time stamp
    end_time_= datetime.datetime.now()
    time_used = end_time_ - start_time_
    print("[XClone HMM smoothing] Time used:", time_used.seconds, "seconds")

    update_Xdata.layers["posterior_mtx"] = res
    update_Xdata.layers["posterior_mtx_log"] = res_log
    return update_Xdata

################################################################end of the HMM smoothing functions