"""Base functions for XClone BAF processing BB_HMM framework.
"""

# Author: Rongting Huang
# Date: 2021-05-04
# update: 2022-05-20

##############################################
## pipeline in BAF module
# (1) BAF data read and check
# (2) BAF data preprocessing-phasing
##############################################

import datetime
import numpy as np
import pandas as pd

from .smoothing import WMA_smooth, KNN_smooth
from .base_utils import normalize
from .base_utils import cal_log_lik

from .HMM_BB import gene_specific_BAF
from .HMM_BB import calculate_Xemm_prob_bb
from .HMM_base import XHMM_smoothing
from ._RDR_CNVratio import CNV_visualization

from scipy.special import logsumexp
from scipy.sparse import csr_matrix
from scipy.special import logit, expit

def extrme_count_capping(Xdata,
                         quantile_pct = 0.99,
                         verbose = False):
    """
    remove extreme counts influence.
    """
    ad_Counts = Xdata.layers["ad_bin"].A.sum(axis=0)
    dp_Counts = Xdata.layers["dp_bin"].A.sum(axis=0)
    if verbose:
        ## visualize counts distribution
        import matplotlib.pylab as plt
        plt.plot(ad_Counts)
        plt.plot(dp_Counts)

    cutoff_ = np.quantile(dp_Counts, quantile_pct)
    flag_ = dp_Counts > cutoff_
    ratio_replace = dp_Counts[flag_] / cutoff_

    Xdata.layers["dp_bin_backup"] = Xdata.layers["dp_bin"].copy()
    Xdata.layers["ad_bin_backup"] = Xdata.layers["ad_bin"].copy()
    # mask to, and insert your replacement values:
    Xdata.layers["dp_bin"][:, flag_] = np.ceil(Xdata.layers["dp_bin"][:, flag_]/ratio_replace)
    Xdata.layers["ad_bin"][:, flag_] = np.ceil(Xdata.layers["ad_bin"][:, flag_]/ratio_replace)

    if verbose:
        ## update counts distribution
        ad_Counts = Xdata.layers["ad_bin"].A.sum(axis=0)
        dp_Counts = Xdata.layers["dp_bin"].A.sum(axis=0)
        plt.plot(ad_Counts)
        plt.plot(dp_Counts)

    return Xdata

def concentration_mapping(Xdata, 
                          concentration_lower = 30,
                          concentration_uppper = 100):
    """
    Return bin specific concentration.
    """
    # ad_Counts = Xdata.layers["ad_bin"].A.sum(axis=0)
    dp_Counts = Xdata.layers["dp_bin"].A.sum(axis=0)
    max_cnt = dp_Counts.max()
    min_cnt = dp_Counts.min()
    
    
    x_ = np.array([max_cnt, min_cnt])
    y_ = np.array([concentration_lower, concentration_uppper])
    
    from scipy.stats import linregress
    linear_res = linregress(x_, y_)
    concentration = linear_res.slope * dp_Counts + linear_res.intercept
    Xdata.var["concentration"] = concentration

    return Xdata
    

def BAF_smoothing(Xdata, 
                  inlayer = "BAF_emm_prob_log",
                  outlayer = "BAF_emm_prob_log_KNN",
                  KNN_smooth = True, 
                  KNN_connectivities_key = "connectivities_expr"):
    """
    ## todo check in and out
    do KNN smoothing across cells.
    do WMA smoothing across genes/bins.

    apply on calculated BAF, or
    apply on calculated emm_prob_log.--but seems only use one type of smoothing.
    ## todo improve and integrate all smoothing on BAF or emm_prob
    """
    start_t = datetime.datetime.now()

    emm_prob_log = Xdata.layers[inlayer].copy()
    print("normalize the input emm_prob_log")
    # normalised loglikelihood (dosen't affect HMM)
    emm_prob_log += -logsumexp(emm_prob_log, axis=2, keepdims=True)
    print("normalized emm_prob_log")

    ## optional: smooth the emm_prob_log with KNN graph
    if KNN_smooth:
        if KNN_connectivities_key not in Xdata.obsp:
            print('Warning: No KNN connectivities available, skipped.')
            return Xdata

        connectivities = normalize(Xdata.obsp[KNN_connectivities_key])
        for k in range(emm_prob_log.shape[2]):
            # emm_prob_log[:, :, k] = connectivities @ emm_prob_log[:, :, k]
            emm_prob_log[:, :, k] = connectivities @ csr_matrix(emm_prob_log[:, :, k]).toarray()
  
    print("generate new layer key value:", outlayer)
    Xdata.layers[outlayer] = emm_prob_log

    end_t = datetime.datetime.now()
    elapsed_sec = (end_t - start_t).total_seconds()
    print("[BAF smoothing] time_used: " + "{:.2f}".format(elapsed_sec) + "seconds")

    return Xdata

def fit_BAF_theoretical_values(Xdata, 
                               chr_lst = None, 
                               celltype_lst = None, 
                               cell_anno = "cell_type", 
                               region_anno = "chr_arm",
                               random_seed = None,
                               n_sample_cells = 200,
                               AD_layer = "ad_bin_phased",
                               DP_layer = "dp_bin",
                               iterations = 1,
                               calculate_ref = True,
                               ref_celltype = "unclassified",
                               verbose = True):
    """
    Sample areas from copy loss/copy gain (ground truth) for 
    fitting of state specific BAF theoretical values.
    Based on Binomial analytical prob(mu_).
    """

    if random_seed is not None:
        rvgen = np.random.RandomState(random_seed)
    else:
        rvgen = np.random
    
    if chr_lst is None:
        is_region = Xdata.var[region_anno] == Xdata.var[region_anno]
    else:
        is_region = Xdata.var[region_anno].isin(chr_lst)

    if celltype_lst is None:
        is_cell = Xdata.obs[cell_anno] == Xdata.obs[cell_anno]
    else:
        is_cell = Xdata.obs[cell_anno].isin(celltype_lst)

    sample_pool = Xdata[:, is_region][is_cell, :].copy()
    n_cells = sample_pool.shape[0]
    theoretical_prob = []
    for i in range(iterations):
        sample_idx = rvgen.randint(0, n_cells, n_sample_cells)

        AD_sum = sample_pool[sample_idx, :].layers[AD_layer].sum()
        DP_sum = sample_pool[sample_idx, :].layers[DP_layer].sum()

        theoretical_prob_value = AD_sum / (AD_sum + DP_sum)
        if verbose:
            print("[XClone]fit_BAF_theoretical_value", i, ":",theoretical_prob_value)
        theoretical_prob.append(theoretical_prob_value)
    if calculate_ref:
        is_ref = Xdata.obs[cell_anno] == ref_celltype
        ref_obs = Xdata[is_ref, :][:, is_region].copy()
        ref_BAF = ref_obs.layers["fill_BAF1_phased"].mean()
        # ref_BAF = (ref_obs.layers[AD_layer] / ref_obs.layers[DP_layer]).mean()
    else:
        ref_BAF = 0.5
    if verbose:
        print("theoretical_prob:", theoretical_prob)
        print("ref_BAF:", ref_BAF)
    return theoretical_prob, ref_BAF

def BAF_theoretical_value_optimization(merge_Xdata,
            AD_layer = "ad_bin_phased",
            DP_layer = "dp_bin",
            ref_anno = "ref_BAF1_phased_clipped",
            random_seed1 = None,
            random_seed2 = None,
            init_delta1 = 0.7, 
            init_delta2 = 1.4,
            region_anno1 = "chr_arm",
            region_anno2 = "chr_arm",
            chr_lst1 = None,
            chr_lst2 = None,
            celltype_lst1 = None,
            celltype_lst2 = None,
            max_iter = 20, 
            min_iter = 10,
            epsilon_conv = 1e-2, 
            fitBAF_verbose = True,
            HMM_verbose = False,
            start_prob = None,
            trans_prob = None,
            KNN_smooth = False,
            nproc = 1,
            log_display = True,
            verbose = True,
            **kwargs):
    """
    ## todo-optimization not yet error now
    init_delta1 = 0.7 for copy gain,
    init_delta2 = 1.4 for copy loss.
    
    kwargs: general params in fit_BAF_theoretical_values if not using default.
    cell_anno
    n_sample_cells
    AD_layer
    DP_layer
    """
    start_t = datetime.datetime.now()
    ## initalization
    ### likelihood init
    Logliklihood = np.zeros(max_iter)

    ### params setting and fit_BAF_ratio
    fit_params = {}
    fit_params['Xdata'] = merge_Xdata
    fit_params['AD_layer'] = AD_layer
    fit_params['DP_layer'] = DP_layer
    fit_params['chr_lst'] = chr_lst1
    fit_params['celltype_lst'] = celltype_lst1
    fit_params['region_anno'] = region_anno1
    fit_params['random_seed'] = random_seed1
    fit_params['iterations'] = max_iter-1
    fit_params['verbose'] = fitBAF_verbose
    fit_params.update(**kwargs)
    theoretical_prob1, ref_BAF1 = fit_BAF_theoretical_values(**fit_params)
    delta1 = np.absolute(logit(theoretical_prob1) - logit(ref_BAF1))
    # delta1 = np.absolute(logit(theoretical_prob1) - logit(0.5))
    print("delta1:", delta1)
    
    fit_params['chr_lst'] = chr_lst2
    fit_params['celltype_lst'] = celltype_lst2
    fit_params['region_anno'] = region_anno2
    fit_params['random_seed'] = random_seed2
    theoretical_prob2, ref_BAF2 = fit_BAF_theoretical_values(**fit_params)
    delta2 = np.absolute(logit(theoretical_prob2) - logit(ref_BAF2))
    # delta2 = np.absolute(logit(theoretical_prob2) - logit(0.5))
    print("delta2:", delta2)

    if start_prob is None:
        start_prob = np.array([0.1, 0.1, 0.6, 0.1, 0.1])
    if trans_prob is None:
        t = 1e-6
        trans_prob = np.array([[1-4*t, t, t, t, t],[t, 1-4*t, t, t, t],[t, t, 1-4*t, t, t], 
        [t, t, t, 1-4*t, t], [t, t, t, t, 1-4*t]])

    ## iteration
    for it in range(max_iter):
        if verbose == True:
            print("[XClone] CNV_optimazation iteration: ", it+1)
        if it == 0: 
            # init round
            # based on 0.5 or based on ref
            copy_gain1 = expit(logit(0.5) - init_delta1)
            copy_gain2 = expit(logit(0.5) + init_delta1)
            copy_loss1 = expit(logit(0.5) - init_delta2)
            copy_loss2 = expit(logit(0.5) + init_delta2)
            theo_states_ = np.array([copy_loss1, copy_loss2, copy_gain1, copy_gain2])
        else:
            copy_gain1 = expit(logit(ref_BAF1) - delta1[it-1])
            copy_gain2 = expit(logit(ref_BAF1) + delta1[it-1])
            copy_loss1 = expit(logit(ref_BAF2) - delta2[it-1])
            copy_loss2 = expit(logit(ref_BAF2) + delta2[it-1])
            # copy_gain1 = expit(logit(0.5) - delta1[it-1])
            # copy_gain2 = expit(logit(0.5) + delta1[it-1])
            # copy_loss1 = expit(logit(0.5) - delta2[it-1])
            # copy_loss2 = expit(logit(0.5) + delta2[it-1])
            theo_states_ = np.array([copy_loss1, copy_loss2, copy_gain1, copy_gain2])
        # used_specific_states = gene_specific_BAF(Xdata, theo_states= np.array([0.2, 0.8, 1/3, 2/3]), 
        #                                 specific_BAF = "ref_BAF1_phased_clipped", rescale=True)
        used_specific_states = gene_specific_BAF(merge_Xdata, theo_states= theo_states_, 
                                    specific_BAF = ref_anno, rescale=True)
        merge_Xdata = calculate_Xemm_prob_bb(merge_Xdata, 
                                         AD_key = AD_layer, 
                                         DP_key = DP_layer, 
                                         concentration = 100, 
                                         uns_out = "bin_phased_BAF_specific_center_emm_prob_log", 
                                         states = used_specific_states)
        if KNN_smooth:
            merge_Xdata = BAF_smoothing(merge_Xdata, 
                  uns_in = "bin_phased_BAF_specific_center_emm_prob_log",
                  uns_out = "bin_phased_BAF_specific_center_emm_prob_log_KNN",
                  KNN_smooth = True)
            emm_prob_log = merge_Xdata.uns["bin_phased_BAF_specific_center_emm_prob_log_KNN"].copy()
        else:
            emm_prob_log = merge_Xdata.uns["bin_phased_BAF_specific_center_emm_prob_log"].copy()

        ### posterior_mtx_log
        update_Xdata = XHMM_smoothing(merge_Xdata, emm_prob_log = emm_prob_log, 
                                    start_prob = start_prob, trans_prob = trans_prob, 
                                    nproc = nproc, verbose = HMM_verbose)
        if log_display == True:
            CNV_visualization(update_Xdata, states_weight = np.array([1, 1, 2, 3, 3]))
            CNV_visualization(update_Xdata, states_weight = np.array([1, 1, 2, 3, 3]), weights = False)

        _logLik = cal_log_lik(update_Xdata.uns["emm_prob_log_noHMMsmoothing"], update_Xdata.uns["posterior_mtx_log"])
        Logliklihood[it] = _logLik

        if it > min_iter:
            if Logliklihood[it] < Logliklihood[it - 1]:
                if verbose:
                    print("[XClone] Warning: Lower bound decreases!\n")
            elif it == max_iter - 1:
                if verbose:
                    print("[XClone] Warning: CNV ration optimization did not converge!\n")
                    print("[XClone] Notes: try to increase the max_iter: ", max_iter, "!\n")
            elif Logliklihood[it] - Logliklihood[it - 1] < epsilon_conv:
                break
    Logliklihood = Logliklihood[:it+1]

    if verbose == True:
        print("iteration_end_round: ", it+1)
        print("Logliklihood: ", Logliklihood)
    
    end_t = datetime.datetime.now()
    elapsed_sec = (end_t - start_t).total_seconds()
    print("[XClone fitting BAF theoratical_values] time_used: " + "{:.2f}".format(elapsed_sec) + "seconds")
    return update_Xdata

## todo- not used yet
def BAF_theoretical_value_CNV_optimization(merge_Xdata,
            AD_layer = "ad_bin_phased",
            DP_layer = "dp_bin",
            ref_anno = "ref_BAF1_phased_clipped",
            random_seed1 = None,
            random_seed2 = None,
            init_delta1 = 0.7, 
            init_delta2 = 1.4,
            region_anno1 = "chr_arm",
            region_anno2 = "chr_arm",
            chr_lst1 = None,
            chr_lst2 = None,
            celltype_lst1 = None,
            celltype_lst2 = None,
            max_iter = 20, 
            min_iter = 10,
            epsilon_conv = 1e-2, 
            fitBAF_verbose = True,
            HMM_verbose = False,
            start_prob = None,
            trans_prob = None,
            KNN_smooth = False,
            nproc = 1,
            log_display = True,
            verbose = True,
            **kwargs):
    """
    init_delta1 = 0.7 for copy gain,
    init_delta2 = 1.4 for copy loss.
    
    kwargs: general params in fit_BAF_theoretical_values if not using default.
    cell_anno
    n_sample_cells
    AD_layer
    DP_layer
    """
    start_t = datetime.datetime.now()
    ## initalization
    ### likelihood init
    Logliklihood = np.zeros(max_iter)

    ### params setting and fit_BAF_ratio
    fit_params = {}
    fit_params['Xdata'] = merge_Xdata
    fit_params['AD_layer'] = AD_layer
    fit_params['DP_layer'] = DP_layer
    fit_params['chr_lst'] = chr_lst1
    fit_params['celltype_lst'] = celltype_lst1
    fit_params['region_anno'] = region_anno1
    fit_params['random_seed'] = random_seed1
    fit_params['iterations'] = max_iter-1
    fit_params['verbose'] = fitBAF_verbose
    fit_params.update(**kwargs)
    theoretical_prob1, ref_BAF1 = fit_BAF_theoretical_values(**fit_params)
    delta1 = np.absolute(logit(theoretical_prob1) - logit(ref_BAF1))
    # delta1 = np.absolute(logit(theoretical_prob1) - logit(0.5))
    print("delta1:", delta1)
    
    fit_params['chr_lst'] = chr_lst2
    fit_params['celltype_lst'] = celltype_lst2
    fit_params['region_anno'] = region_anno2
    fit_params['random_seed'] = random_seed2
    theoretical_prob2, ref_BAF2 = fit_BAF_theoretical_values(**fit_params)
    delta2 = np.absolute(logit(theoretical_prob2) - logit(ref_BAF2))
    # delta2 = np.absolute(logit(theoretical_prob2) - logit(0.5))
    print("delta2:", delta2)

    if start_prob is None:
        start_prob = np.array([0.1, 0.1, 0.6, 0.1, 0.1])
    if trans_prob is None:
        t = 1e-6
        trans_prob = np.array([[1-4*t, t, t, t, t],[t, 1-4*t, t, t, t],[t, t, 1-4*t, t, t], 
        [t, t, t, 1-4*t, t], [t, t, t, t, 1-4*t]])

    ## iteration
    for it in range(max_iter):
        if verbose == True:
            print("[XClone] CNV_optimazation iteration: ", it+1)
        if it == 0: 
            # init round
            # based on 0.5 or based on ref
            copy_gain1 = expit(logit(0.5) - init_delta1)
            copy_gain2 = expit(logit(0.5) + init_delta1)
            copy_loss1 = expit(logit(0.5) - init_delta2)
            copy_loss2 = expit(logit(0.5) + init_delta2)
            theo_states_ = np.array([copy_loss1, copy_loss2, copy_gain1, copy_gain2])
        else:
            copy_gain1 = expit(logit(ref_BAF1) - delta1[it-1])
            copy_gain2 = expit(logit(ref_BAF1) + delta1[it-1])
            copy_loss1 = expit(logit(ref_BAF2) - delta2[it-1])
            copy_loss2 = expit(logit(ref_BAF2) + delta2[it-1])
            # copy_gain1 = expit(logit(0.5) - delta1[it-1])
            # copy_gain2 = expit(logit(0.5) + delta1[it-1])
            # copy_loss1 = expit(logit(0.5) - delta2[it-1])
            # copy_loss2 = expit(logit(0.5) + delta2[it-1])
            theo_states_ = np.array([copy_loss1, copy_loss2, copy_gain1, copy_gain2])
        # used_specific_states = gene_specific_BAF(Xdata, theo_states= np.array([0.2, 0.8, 1/3, 2/3]), 
        #                                 specific_BAF = "ref_BAF1_phased_clipped", rescale=True)
        used_specific_states = gene_specific_BAF(merge_Xdata, theo_states= theo_states_, 
                                    specific_BAF = ref_anno, rescale=True)
        merge_Xdata = calculate_Xemm_prob_bb(merge_Xdata, 
                                         AD_key = AD_layer, 
                                         DP_key = DP_layer, 
                                         concentration = 100, 
                                         uns_out = "bin_phased_BAF_specific_center_emm_prob_log", 
                                         states = used_specific_states)
        if KNN_smooth:
            merge_Xdata = BAF_smoothing(merge_Xdata, 
                  uns_in = "bin_phased_BAF_specific_center_emm_prob_log",
                  uns_out = "bin_phased_BAF_specific_center_emm_prob_log_KNN",
                  KNN_smooth = True)
            emm_prob_log = merge_Xdata.uns["bin_phased_BAF_specific_center_emm_prob_log_KNN"].copy()
        else:
            emm_prob_log = merge_Xdata.uns["bin_phased_BAF_specific_center_emm_prob_log"].copy()

        ### posterior_mtx_log
        update_Xdata = XHMM_smoothing(merge_Xdata, emm_prob_log = emm_prob_log, 
                                    start_prob = start_prob, trans_prob = trans_prob, 
                                    nproc = nproc, verbose = HMM_verbose)
        if log_display == True:
            CNV_visualization(update_Xdata, states_weight = np.array([1, 1, 2, 3, 3]))
            CNV_visualization(update_Xdata, states_weight = np.array([1, 1, 2, 3, 3]), weights = False)

        _logLik = cal_log_lik(update_Xdata.uns["emm_prob_log_noHMMsmoothing"], update_Xdata.uns["posterior_mtx_log"])
        Logliklihood[it] = _logLik

        if it > min_iter:
            if Logliklihood[it] < Logliklihood[it - 1]:
                if verbose:
                    print("[XClone] Warning: Lower bound decreases!\n")
            elif it == max_iter - 1:
                if verbose:
                    print("[XClone] Warning: CNV ration optimization did not converge!\n")
                    print("[XClone] Notes: try to increase the max_iter: ", max_iter, "!\n")
            elif Logliklihood[it] - Logliklihood[it - 1] < epsilon_conv:
                break
    Logliklihood = Logliklihood[:it+1]

    if verbose == True:
        print("iteration_end_round: ", it+1)
        print("Logliklihood: ", Logliklihood)
    
    end_t = datetime.datetime.now()
    elapsed_sec = (end_t - start_t).total_seconds()
    print("[XClone fitting BAF theoratical_values] time_used: " + "{:.2f}".format(elapsed_sec) + "seconds")
    return update_Xdata

class CNVoptimizer():
    """
    """
    def __init__(self):
        pass
    def fit():
        pass