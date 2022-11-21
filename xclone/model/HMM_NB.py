"""Base funcs for feature(gene/block) specific HMM model in XClone RDR module"""

# This file contains a few HMM models and its combinations
# For emmission probs - Negative binomial count model

# Author: Rongting Huang
# Date: 2021-12-13
# update: 2022-05-20

import numpy as np
import pandas as pd
import scipy as sp

import anndata as ad
import itertools
import datetime

from scipy.stats import nbinom
from scipy.stats import poisson

from scipy.sparse import csr_matrix

from scipy.special import logsumexp

from .base_utils import normalize, loglik_amplify

## Part I: Emmission_prob based on negative binomial for HMM in RDR module

def generate_nb_logprob(reference, observations,
                        states = np.array([0.5, 1.0, 1.5]),
                        gene_specific= False, overdispersion = None, 
                        ref_normalization_term = None, obs_normalization_term = None):
    """
    Function: 
    negative binomial distributions[scipy]

    update version for numpy array
    emission prob log format - gene specific distribution-params

    ## Normally, the normalization term is for cell library, but sometimes we generate the emm_prob for each chromosome
    ## so need specify the normalization term in the cell scale if the input(reference and observations) are chr specific.

    ## best way is calculate the emm_prob at first, in the cell scale
    
    ## default None means the input ref and obs are genome scale and 
    ## can achieve the normalization term via sum and calculate the normalized value.
    ## if learned lib ratio(obs_normalization_term) exists, then the ref_normalization_term should be 1.

    Parameters
    ----------
    reference:
    observations:
    overdispersion:
    gene_specific:
    ref_normalization_term:
    obs_normalization_term:
    states: default np.array([0.5, 1.0]) # example for copy loss and neutral

    Returns
    -------
    emm_prob_log: nbinom.logpmf, should be in `log` format.
    
    Example
    -------
    reference = np.array([10, 10, 10,10,10, 10,10])
    # reference = np.array([5., 9, 3, 5, 12, 6, 9])
    observations = np.array([[5., 9, 3, 5, 12, 6, 9]])
    observations = np.array([[5., 9, 3, 5, 12, 6, 9],[5., 9, 3, 5, 12, 6, 9]])

    states = np.array([0.5, 1.0])
    emm_prob_log = generate_nb_logprob_update(reference, observations, states = states)

    emm_prob_log_update = generate_nb_logprob_update(reference, observations, states = states) # default overdispersion 0.1
    emm_prob_log_update = generate_nb_logprob_update(reference, observations, overdispersion = 0.2, states = states)

    overdispersion = np.array([0.1, 0.2, 0.3, 0.1, 0.2, 0.3, 0.1])
    emm_prob_log = generate_nb_logprob_update(reference, observations, 
                                              overdispersion, gene_specfic= True, 
                                              ref_normalization_term = 1, obs_normalization_term = sample_libsize,
                                              states = states)
    """
    if ref_normalization_term is None:
        ref_normalization_term = reference.sum()
    if obs_normalization_term is None:
        obs_normalization_term= observations.sum(axis=1)
    else:
        obs_normalization_term = np.array(obs_normalization_term)

    if gene_specific:
        if len(overdispersion) == observations.shape[1]:
            pass
        else:
            raise ValueError("[XClone]-pls check gene-specific overdispersion")
    elif overdispersion is None:
        overdispersion = 0.1
    
    _GEX_ref = reference/ref_normalization_term
    _GEX_ref = states * _GEX_ref[:,np.newaxis]

    ## extend the GEX_ref the same format with observations
    # for i in range(observations.shape[0]):
    #     if i == 0:
    #         _GEX_ref_extend =  _GEX_ref[np.newaxis,:]
    #     else:
    #         _GEX_ref_extend =  np.vstack((_GEX_ref_extend, _GEX_ref[np.newaxis,:]))
    ## accelerated
    _GEX_ref_extend = np.tile(_GEX_ref, (observations.shape[0],1)).reshape((-1,_GEX_ref.shape[0] ,_GEX_ref.shape[1]))
    
    _mu =  _GEX_ref_extend * obs_normalization_term[:,np.newaxis,np.newaxis]

    ## if overdispersion is gene specfic, then use overdispersion[:,np.newaxis]
    if gene_specific:
        _var = _mu + overdispersion[np.newaxis, :, np.newaxis] * _mu**2
    else:
        _var = _mu + overdispersion * _mu**2
    
    _nb_prob =  _mu / _var
    _nb_total = _mu * _nb_prob / (1 - _nb_prob)

    # for i in range(observations.shape[0]):
    #     obs_tmp = np.repeat(observations[i], states.shape[-1]).reshape(observations.shape[1],states.shape[-1])
    #     if i == 0:
    #         obs_ = obs_tmp[np.newaxis,:]
    #     else:
    #         obs_ = np.vstack((obs_, obs_tmp[np.newaxis,:]))
    ## accelerated
    obs_ = np.repeat(observations, states.shape[-1]).reshape(observations.shape[0], -1, states.shape[-1])
    
    emm_prob_log = nbinom.logpmf(obs_, _nb_total, _nb_prob)

    print("generate a emm_prob matrix for %d states, matrix shape is %s" %(states.shape[-1], str(emm_prob_log.shape)))
    return emm_prob_log

def generate_nb_logprob2(reference, observations, ref_pred_obs,
                        states = np.array([0.5, 1.0, 1.5]),
                        gene_specific= False, overdispersion = None, 
                        ref_normalization_term = None, obs_normalization_term = None):
    """
    Function: 
    negative binomial distributions[scipy]

    update version for numpy array
    emission prob log format - gene specific distribution-params
    update version for ref_bio information predicted obs, using NMF. 2022-03-04

    ## Normally, the normalization term is for cell library, but sometimes we generate the emm_prob for each chromosome
    ## so need specify the normalization term in the cell scale if the input(reference and observations) are chr specific.

    ## best way is calculate the emm_prob at first, in the cell scale
    
    ## default None means the input ref and obs are genome scale and 
    ## can achieve the normalization term via sum and calculate the normalized value.
    ## if learned lib ratio(obs_normalization_term) exists, then the ref_normalization_term should be 1.

    Parameters
    ----------
    reference:
    observations:
    overdispersion:
    gene_specific:
    ref_normalization_term:
    obs_normalization_term:
    states: default np.array([0.5, 1.0]) # example for copy loss and neutral

    Returns
    -------
    emm_prob_log: nbinom.logpmf, should be in `log` format.
    
    Example
    -------
    reference = np.array([10, 10, 10,10,10, 10,10])
    # reference = np.array([5., 9, 3, 5, 12, 6, 9])
    observations = np.array([[5., 9, 3, 5, 12, 6, 9]])
    observations = np.array([[5., 9, 3, 5, 12, 6, 9],[5., 9, 3, 5, 12, 6, 9]])

    states = np.array([0.5, 1.0])
    emm_prob_log = generate_nb_logprob_update(reference, observations, states = states)

    emm_prob_log_update = generate_nb_logprob_update(reference, observations, states = states) # default overdispersion 0.1
    emm_prob_log_update = generate_nb_logprob_update(reference, observations, overdispersion = 0.2, states = states)

    overdispersion = np.array([0.1, 0.2, 0.3, 0.1, 0.2, 0.3, 0.1])
    emm_prob_log = generate_nb_logprob_update(reference, observations, 
                                              overdispersion, gene_specfic= True, 
                                              ref_normalization_term = 1, obs_normalization_term = sample_libsize,
                                              states = states)
    """
    if ref_normalization_term is None:
        ref_normalization_term = reference.sum()
    if obs_normalization_term is None:
        obs_normalization_term= observations.sum(axis=1)
    else:
        obs_normalization_term = np.array(obs_normalization_term)

    if gene_specific:
        if len(overdispersion) == observations.shape[1]:
            pass
        else:
            raise ValueError("[XClone]-pls check gene-specific overdispersion")
    elif overdispersion is None:
        overdispersion = 0.1
    
    _GEX_ref = reference/ref_normalization_term
    _GEX_ref = states * _GEX_ref[:,np.newaxis]

    _GEX_ref_extend = np.tile(_GEX_ref, (observations.shape[0],1)).reshape((-1, _GEX_ref.shape[0], _GEX_ref.shape[1]))
    
    # states_num = len(states)
    # _ref_pred_obs_extend = np.tile(ref_pred_obs, (states_num,1)).reshape((-1, ref_pred_obs.shape[0], ref_pred_obs.shape[1])).swapaxes(0,2)
    
    
    # take libratio (obs_normalization_term) and ref bio information ratio (ref_pred_obs) into account
    _mu =  _GEX_ref_extend * ref_pred_obs[:,:, np.newaxis] * obs_normalization_term[:,np.newaxis,np.newaxis]

    ## if overdispersion is gene specfic, then use overdispersion[:,np.newaxis]
    if gene_specific:
        _var = _mu + overdispersion[np.newaxis, :, np.newaxis] * _mu**2
    else:
        _var = _mu + overdispersion * _mu**2
    
    _nb_prob =  _mu / _var
    _nb_total = _mu * _nb_prob / (1 - _nb_prob)

    obs_ = np.repeat(observations, states.shape[-1]).reshape(observations.shape[0], -1, states.shape[-1])
    
    emm_prob_log = nbinom.logpmf(obs_, _nb_total, _nb_prob)

    print("generate a emm_prob matrix for %d states, matrix shape is %s" %(states.shape[-1], str(emm_prob_log.shape)))
    return emm_prob_log


def calculate_Xemm_probTry(Xdata, 
                           dispersion_key = "dispersion", 
                           states=None, 
                           clip_low=None, 
                           KNN_smooth=True, 
                           outlayer = "emm_prob_log",
                           combine_BAF = False):
    """
    Shape: (n_cell, n_gene, n_state)
    """
    if sp.sparse.issparse(Xdata.X):
        Xmtx = Xdata.X.A
    else:
        Xmtx = Xdata.X
    
    ## Need to ensure the shape of states
    # states = states[0, :].reshape(1, 1, -1)
    states = np.expand_dims(states, axis=0)

    obs_ = np.expand_dims(Xmtx, axis=2)

    _mu = np.expand_dims(Xdata.layers['expected'], axis=2) * states
    _overdisp = Xdata.var[dispersion_key].values.reshape(1, -1, 1)
    _var = _mu + _overdisp * _mu**2

    _nb_prob =  _mu / _var
    _nb_total = _mu * _nb_prob / (1 - _nb_prob)

    emm_prob_log = nbinom.logpmf(obs_, _nb_total, _nb_prob)
    
    if combine_BAF:
        emm_prob_log += Xdata.uns["BAF_emm_prob_log"]
    
    # normalised loglikelihood (dosen't affect HMM)
    emm_prob_log += -logsumexp(emm_prob_log, axis=2, keepdims=True)

    # clip
    if clip_low is not None:
        # clip_low = -5 # should be a safe trial (generally not much difference)
        print('proportion below the threshold: ' 
              %(np.mean(emm_prob_log < clip_low)))
        emm_prob_log = np.clip(emm_prob_log, clip_low, 0)

    ## optional: smooth the emm_prob_log with KNN graph
    if KNN_smooth:
        if 'connectivities' not in Xdata.obsp:
            print('Warning: No KNN connectivities available, skipped.')
            return emm_prob_log

        # from scvelo.preprocessing.neighbors import get_connectivities
        # connectivities = get_connectivities(Xdata, mode='connectivities', 
        #                                     recurse_neighbors=False)
        # normalize connectivities
        # important notes: need normalization and use sparse matrix in calculation.
        # connectivities = Xdata.obsp['connectivities']
        connectivities = normalize(Xdata.obsp['connectivities'])
        for k in range(emm_prob_log.shape[2]):
            emm_prob_log[:, :, k] = connectivities @ csr_matrix(emm_prob_log[:, :, k]).toarray()
    
    Xdata.layers[outlayer] = emm_prob_log
    return Xdata


def calculate_Xemm_prob(Xdata,
                        states = None, 
                        gene_specific = False, overdispersion = None,
                        ref_normalization_term = None,
                        obs_normalization_term = None,
                        verbose = True):
    """
    Function:
    Usage of FUNC `generate_nb_logprob`.


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
    if sp.sparse.issparse(Xdata.X):
        X_mtx = Xdata.X.A
    else:
        X_mtx = Xdata.X
    
    ref_ = X_mtx[0, :]
    obs_ = X_mtx[1:, :]
    
    if states is None:
        states = np.array([0.5, 1.0, 1.5])

    ## check the shape for the input
    if gene_specific == True:
        if len(overdispersion) != len(ref_):
            raise ValueError("[XClone]overdispersion doesn't match with the gene numbers! Pls check!")
    if obs_normalization_term is not None:
        if len(obs_normalization_term) != obs_.shape[0]:
            raise ValueError("[XClone]libratio doesn't match with the cell numbers! Pls check!")
 
    ## calculate emm_prob_log
    emm_prob_log = generate_nb_logprob( ref_, obs_,
                                        states = states,
                                        gene_specific = gene_specific, overdispersion = overdispersion, 
                                        ref_normalization_term = ref_normalization_term,
                                        obs_normalization_term = obs_normalization_term)
    if verbose:
        print("reference", ref_)
        print("observations", obs_)
        print("states", states)
        print("emm_prob_log", emm_prob_log)
    
    # time stamp
    end_time_= datetime.datetime.now()
    time_used = end_time_ - start_time_
    print("cal emm prob time", time_used.seconds, "seconds")
    return emm_prob_log


def calculate_Xemm_prob2(Xdata,
                        states = None, 
                        gene_specific = False, overdispersion = None,
                        ref_normalization_term = None,
                        obs_normalization_term = None,
                        verbose = True):
    """
    Function:
    Usage of FUNC `generate_nb_logprob`.


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
    if sp.sparse.issparse(Xdata.X):
        X_mtx = Xdata.X.A
    else:
        X_mtx = Xdata.X
    
    ref_ = X_mtx[0, :]
    obs_ = X_mtx[1:, :]

    obs_pred_ = Xdata.uns["NMF_pred_obs"]
    
    if states is None:
        states = np.array([0.5, 1.0, 1.5])

    ## check the shape for the input
    if gene_specific == True:
        if len(overdispersion) != len(ref_):
            raise ValueError("[XClone]overdispersion doesn't match with the gene numbers! Pls check!")
    if obs_normalization_term is not None:
        if len(obs_normalization_term) != obs_.shape[0]:
            raise ValueError("[XClone]libratio doesn't match with the cell numbers! Pls check!")
 
    ## calculate emm_prob_log
    emm_prob_log = generate_nb_logprob2(ref_, obs_, obs_pred_, 
                                        states = states,
                                        gene_specific = gene_specific, overdispersion = overdispersion, 
                                        ref_normalization_term = ref_normalization_term,
                                        obs_normalization_term = obs_normalization_term)
    if verbose:
        print("reference", ref_)
        print("observations", obs_)
        print("states", states)
        print("emm_prob_log", emm_prob_log)
    
    # time stamp
    end_time_= datetime.datetime.now()
    time_used = end_time_ - start_time_
    print("cal emm prob time", time_used.seconds, "seconds")
    return emm_prob_log

## Next Part : HMM smoothing in RDR(see HMM_base)