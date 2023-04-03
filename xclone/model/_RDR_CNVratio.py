"""Base functions for XClone RDR processing
strategies for fitting the CNV ratio.
"""

# Author: Rongting Huang
# Date: 2022-01-20
# update: 2022-03-15

## RDR
## Part-IV Estimate CNV states' ratio based via different gene expression group [source _RDR]
import datetime
import numpy as np
import scipy as sp

from .base_utils import normalize
from .base_utils import cal_log_lik

from .HMM_NB import calculate_Xemm_prob2
# from .HMM_NB import calculate_Xemm_prob
from .HMM_NB import calculate_Xemm_probTry
from .HMM_base import XHMM_smoothing

from ..plot._data import reorder_data_by_cellanno
from ..plot._base_xanndata import Xheatmap
from ..plot.CNV_plot import CNV_visualization

def gene_exp_group(Xdata, n_group = 5, ref_log = True, verbose=True):
    """
    Function:
    group genes into different clusters based on the ref expression.

    
    Parameters:
    ----------
    Xdata: anndata.
    
    Return:
    ------
    Xdata: anndata.
    
    Example:
    -------

    """
    ref_exp_raw = np.array(Xdata.var["ref_avg"])
    ref_exp_log = np.log(ref_exp_raw)

    if ref_log:
        ref_exp_ = ref_exp_log
    else:
        ref_exp_ = ref_exp_raw
    
    quantile_step = 1/n_group
    quantile_array = np.arange(0, 1, quantile_step)
    
    ## expression break
    expression_brk = []
    expression_brk.append(ref_exp_.min())
    for qt_value in quantile_array[1:]:
        expression_brk.append(np.quantile(ref_exp_, qt_value))
    expression_brk.append(ref_exp_.max())
    expression_brk = np.array(expression_brk)

    if verbose:
        print("expression_brk", expression_brk)
    
    if ref_log:
        Xdata.uns["ref_log_expression_brk"] = expression_brk
    else:
        Xdata.uns["ref_expression_brk"] = expression_brk
    
    ## select genes group based on the expression level
    group_genes = []
    for i in range(n_group):
        flag1 = ref_exp_ >= expression_brk[i]
        flag2 = ref_exp_ < expression_brk[i+1]
        flag_ = flag1 & flag2
        genes_ = Xdata.var["GeneName"][flag_]
        group_genes.append(genes_)
    
    Xdata.uns["group_genes"] = group_genes
    return Xdata

def hard_assign_prob(init_prob, keep_first = False, axis = -1):
    """
    Function:
    Transform the prob to hard assign prob based on the CNV states.(default axis = -1)

    strategy: Find the maximum value from each row set it to 1 and the rest setted as 0.

    Parameters:
    ----------
    init_prob: np.array. cell * gene * CNV state probability matrix.

    keep_first: bool. If there are multiple ones in a row with the same maximum value and 
    you would like to set only the first one as 1, then you set keep_first = True, otherwise
    keep the default False, and will transform all the equal max values to 1.

    axis: default -1.
    
    Example:
    -------
    init_prob = hard_assign_prob(init_prob)
    """
    if keep_first == False:
        out_prob = (init_prob == init_prob.max(axis=axis, keepdims=True)).astype(float)
    else:
        idx_ = init_prob.argmax(axis=axis)
        out_prob = np.zeros_like(init_prob, dtype=float)
        out_prob[np.arange(init_prob.shape[0]), idx_] = 1
    
    return out_prob

def generate_emm_GT(Xdata,
                    n_cells, n_genes, 
                    copy_gain_chr=None, copy_gain_chr_arm=None, 
                    copy_loss_chr=None, copy_loss_chr_arm=None):
    """
    Function:
    generate ground truth of emm_prob for CNV ratio checking.

    Example:
    emm_prob_GT = generate_emm_GT(Xdata,
                    n_cells=2, n_genes=6, 
                    copy_gain_chr=["8","20"], 
                    copy_loss_chr=["18"], copy_loss_chr_arm=["12p"])
    """
    ## init emm_prob with all normal states
    n_states = 3
    emm_prob = np.ones((n_cells, n_genes, n_states))
    emm_prob[:] =  np.array([0.0, 1.0, 0.0])

    if copy_gain_chr is None:
        raise ValueError("[XClone] Pls specify copy_gain_chr!")
    if copy_loss_chr is None:
        raise ValueError("[XClone] Pls specify copy_losss_chr!")
    
    ## init the flag
    ## set ground truth for chrs or chr_arms
    
    gain_flag = Xdata.var["chr"].isin(copy_gain_chr)
    if copy_gain_chr_arm is None:
        pass
    else:
        gain_flag_tmp = Xdata.var["chr_arm"].isin(copy_gain_chr_arm)
        gain_flag = gain_flag|gain_flag_tmp
    
    loss_flag = Xdata.var["chr"].isin(copy_loss_chr)
    if copy_loss_chr_arm is None:
        pass
    else:
        loss_flag_tmp = Xdata.var["chr_arm"].isin(copy_loss_chr_arm)
        loss_flag = loss_flag|loss_flag_tmp
    
    emm_prob[:,gain_flag,:] = np.array([0.0, 0.0, 1.0])
    emm_prob[:,loss_flag,:] = np.array([1.0, 0.0, 0.0])

    return emm_prob


## learn CNV states first.
## 1) strategy1: GMM mixture model, not good for sparse matrix.(log ratio)
## 2) strategy2: learn outliers (from log ratio) as init CNV states

def guide_CNV_chrs(Xdata, 
                   Xlayer = "RDR_smooth", 
                   anno_key = "chr_arm", 
                   remove_XY = True):
    """
    Function:
    --------
    Find most potential chrs with different CNV states in RDR module.
    default find region at the resolution of chr_arm.
    
    Example:
    -------
    >>> chr_lst, anno_key = guide_CNV_chrs(RDR_adata, 
        Xlayer = "RDR_smooth", anno_key = "chr_arm")
    """
    from collections import OrderedDict
    
    mtx_use = Xdata.layers[Xlayer]
    chr_use = Xdata.var[anno_key].unique()

    chr_dict = {}
    for chr_ in chr_use:
        flag_ = Xdata.var[anno_key] == chr_
        avg_ = mtx_use[:, flag_].mean()
        # avg_ = np.median(mtx_use[:, flag_])
        chr_dict[chr_] = avg_
    
    # remove chrX+Y from the chr_dict if exists.
    if remove_XY:
        for chr_items in ['X', 'Y', 'Xp', 'Xq', 'Yp', 'Yq']:
            try:
                chr_dict.pop(chr_items)
            except:
                pass    
    
    sort_chr_dict = OrderedDict(sorted(chr_dict.items(), key=lambda t: t[1]))
    
    dict_cnt = len(sort_chr_dict)
    dict_keys = list(sort_chr_dict.keys())
    
    ## recommend use chr for different CNV states
    loss_ = dict_keys[0]
    neutral_ = dict_keys[round(dict_cnt/2)+1]
    gain_ = dict_keys[dict_cnt-1]
    
    chr_lst = [loss_, neutral_, gain_]
    
    Xdata.uns["chr_dict"] = dict(sort_chr_dict)
    Xdata.uns["guide_CNV_chrs_use_layers"] = Xlayer
    Xdata.uns["guide_CNV_chrs_use_anno_key"] = anno_key
    print("[XClone] RDR CNV states chrs guiding(copy loss, copy neutral, copy gain):", chr_lst)
    
    return chr_lst

def correct_guide_RDR_cnv(guide_cnv_ratio, threshold = 0.15):
    """
    if guide cnv ratio values are similar,(then maybe no cnv states exists here) 
    return default theoratical value

    todo: may change threshold.
    """
    if guide_cnv_ratio[1] - guide_cnv_ratio[0] < threshold:
        guide_cnv_ratio[0] = 0.5
        print("correct RDR CNV guiding copy loss ratio")
    if guide_cnv_ratio[2] - guide_cnv_ratio[1] < threshold:
        guide_cnv_ratio[2] = 1.5
        print("correct RDR CNV guiding copy gain ratio")
    
    ## add clipping for this version
    if guide_cnv_ratio[0] < 0.3:
        guide_cnv_ratio[0] = 0.3
        print("[XClone] warning: correct RDR CNV guiding copy loss ratio")
        print("[XClone] hints: can change guide_qt_lst")
    if guide_cnv_ratio[2] > 3:
        print("[XClone] warning: correct RDR CNV guiding copy gain ratio")
        print("[XClone] hints: can change guide_qt_lst")
        guide_cnv_ratio[2] = 3
    
    return guide_cnv_ratio

def guide_CNV_states(RDR_Xdata, 
                     Xlayer = "RDR_smooth", 
                     anno_key = "chr", 
                     chr_lst = ["18", "1", "8"], 
                     qt_lst = [0.00001, 0.96, 0.99999],
                     states = ["CNV loss: ", "CNV neutral: ", "CNV gain: "],
                     show_boxplot = True):
    """
    get outliers from copy gain/loss region for CNV states init.
    """
    import matplotlib.pylab as plt
    guide_cnv_ratio = []
    for i, tmp_state in zip(range(len(states)), states):
        
        flag_ = RDR_Xdata.var[anno_key] == chr_lst[i]
        data_ = RDR_Xdata[:, flag_].layers[Xlayer]
        value_ = np.exp(np.quantile(data_, qt_lst[i]))
        print(tmp_state, value_)
        guide_cnv_ratio.append(value_)
        if show_boxplot:
            plt.boxplot(np.exp(data_).reshape(-1,1))
            plt.title(tmp_state.split(":")[0])
            plt.show()
    # # copy gain
    # gain_ = np.exp(np.quantile(RDR_Xdata[:, RDR_Xdata.var[anno_key] == chr_lst[0]].layers[Xlayer], qt_lst[0]))
    # # copy loss
    # loss_ = np.exp(np.quantile(RDR_Xdata[:, RDR_Xdata.var[anno_key] == chr_lst[1]].layers[Xlayer], qt_lst[1]))
    # # copy neutral
    # neutral_ = np.exp(np.quantile(RDR_Xdata[:, RDR_Xdata.var[anno_key] == chr_lst[2]].layers[Xlayer], qt_lst[2]))
    
    # print("CNV gain:", gain_)
    # print("CNV loss:", loss_)
    # print("CNV neutral:", neutral_)
    
    guide_cnv_ratio = np.array(guide_cnv_ratio)
    guide_cnv_ratio = correct_guide_RDR_cnv(guide_cnv_ratio)
    print("[XClone] RDR CNV states ratio guiding(copy loss, copy neutral, copy gain):", guide_cnv_ratio)
    return guide_cnv_ratio

from sklearn.mixture import GaussianMixture
def guide_RDR_CNV_ratio(RDR_Xdata, 
                        Xlayer = "RDR_smooth", 
                        anno_key = "chr_arm", 
                        chr_lst = ["18p", "9q", "8q"], 
                        states = ["CNV loss: ", "CNV neutral: ", "CNV gain: "],
                        means_init = None,
                        n_components = 3,
                        max_iter = 50,
                        **kwargs):
    """
    based on GMM model.
    """
    import matplotlib.pylab as plt
    guide_cnv_ratio = []
    
    # for i, tmp_state in zip(range(len(states)), states):
    if True:
        
        # flag_ = RDR_Xdata.var[anno_key] == chr_lst[i]
        flag_ = RDR_Xdata.var[anno_key].isin(chr_lst)

        data_ = RDR_Xdata[:, flag_].layers[Xlayer].reshape(-1,1)
        
        r = np.random.RandomState(seed=1234)
        params = {}
        params["n_components"] = n_components
        params["max_iter"] = max_iter
        params["random_state"] = r
        params["tol"] = 1e-9
        if means_init is not None:
            params["means_init"] = means_init
        # if i == 0:
        #     params["means_init"] = np.array([0.5]).reshape(-1, 1)
        # if i == 1:
        #     params["means_init"] = np.array([1]).reshape(-1, 1)
        # if i == 2:
        #     params["means_init"] = np.array([1.5]).reshape(-1, 1)


        params.update(**kwargs)
        gmm = GaussianMixture(**params).fit(data_)

        

        sort_means = np.exp(np.sort(gmm.means_[:,0]))
        print(sort_means)
        # print(tmp_state, sort_means)
        # if i == 0:
        #     guide_cnv_ratio.append(sort_means[0])

    return None



def fit_CNV_ratio(Xdata, 
                  init_prob_layer = "emm_prob_log", 
                  hard_assign = True, 
                  n_sample_cells = 10,
                  feature_X = None, 
                  avg_key = "ref_avg",
                  depth_key = "library_ratio_capped", 
                  dispersion_set = None,
                  dispersion_key = "dispersion_capped", 
                  random_seed = None, 
                  verbose=True,
                  NB_kwargs={'disp': True, 'skip_hessian': True}):
    """
    Function:
    version 0.0.2

    Based on selected group genes and GLM model to get the estimate CNV ratio
    for different state, default state is copy loss, neutral, copy gain.
    GLM loop is for different gene groups. 

    We estimate CNV state-specific ratio for each gene group, for each gene we
    sample the cells from the cells pool.
 
    Here use weighted GLM.
    https://www.statsmodels.org/devel/examples/notebooks/generated/glm_weights.html#
    https://www.statsmodels.org/dev/examples/notebooks/generated/glm.html
    https://www.statsmodels.org/dev/_modules/statsmodels/discrete/discrete_model.html#NegativeBinomial

    Parameters:
    -----------
    Xdata: anndata.
    init_prob: np.array.
               probability for each gene across cells for all CNV states.
    hard_assign: bool. 
                 preprocessing for init_prob to get the hard assign prob or not.
    n_sample_cells: int.
    dispersion_set: float.
    random_seed: int.      
    verbose : bool
            Whether print out log info
    test_func: bool.
    
    Return:
    ------
    Xdata: anndata.


    Example:
    -------

    """
    import statsmodels.api as sm
    import scipy as sp
    print("[XClone] fit CNV ratio")

    if sp.sparse.issparse(Xdata.X):
        Xmtx = Xdata.X.A
    else:
        Xmtx = Xdata.X
    
    if random_seed is not None:
        rvgen = np.random.RandomState(random_seed)
    else:
        rvgen = np.random
    
    # time stamp
    start_time_ = datetime.datetime.now()
    
    init_prob = normalize(np.exp(Xdata.layers[init_prob_layer].copy()))
    ## for the first round and no region selection
    ## hard_assign should be true
    if hard_assign == True:
        if verbose:
            print("[XClone] fit_CNV_ratio: hard_assign prob")
        init_prob = hard_assign_prob(init_prob)

    # check input
    if Xdata.shape[0] == init_prob.shape[0]:
        if Xdata.shape[1] == init_prob.shape[1]:
            if verbose:
                print("input data is in matched dimension")
        else:
            raise ValueError("[XClone] gene dimension is not matched! Pls check!")
    else:
        raise ValueError("[XClone] cell dimension is not matched! Pls check!")
    
    ## data by different gene_expression groups
    group_genes_lst = Xdata.uns["group_genes"]
    n_gene_group = len(group_genes_lst)
    n_states = init_prob.shape[2]
    n_cells = Xdata.shape[0]

    model_results = []
    CNV_ratio = []
    ref_bio_ratio = []

    ## GLM model
    for i in range(n_gene_group):
        flag_ = Xdata.var["GeneName"].isin(group_genes_lst[i])

        ref_vec = np.array(Xdata.var[avg_key][flag_])
        tmp_mtx = Xmtx[:,flag_]
        tmp_prob = init_prob[:,flag_, :]
        
        libratio = Xdata.obs[depth_key]
        
        if dispersion_set is None:
            # dispersion_ = Xdata.var[dispersion_key][flag_].median()
            # gene-specific dispersion
            dispersion_ = Xdata.var[dispersion_key][flag_]
            print(dispersion_)
        else:
            dispersion_ = dispersion_set
        
        genes_num = tmp_mtx.shape[1]
        
        ## prepare GLM input
        ref_exp_ = np.repeat(ref_vec, n_sample_cells)
        dispersion_used = np.repeat(dispersion_, n_sample_cells)
        
        sample_cell_idx = rvgen.randint(0, n_cells, n_sample_cells*genes_num)

        gene_idx = np.repeat(range(0, genes_num), n_sample_cells)
        gene_exp_ = tmp_mtx[sample_cell_idx, gene_idx]
        freq_prob_ = tmp_prob[sample_cell_idx, gene_idx, :]
        libratio_used = libratio[sample_cell_idx]

        obs_y = gene_exp_.reshape(-1)
        if feature_X is None:
            feature_x = np.ones(len(obs_y))
        else:
            feature_x1 = np.ones(len(obs_y))
            feature_x2 = feature_X[sample_cell_idx, gene_idx].reshape(-1)
            feature_x = np.vstack((feature_x1,feature_x2)).T
        # exposure_x = ref_exp_.reshape(-1) # version 0.0.0 [no lib ratio counted in]
        # should take libratio into account
        exposure_x = np.array(ref_exp_ * libratio_used).reshape(-1) 
       
        for state_idx in range(n_states):
            freq_w = freq_prob_[:, state_idx]

            try:
            # if True:
                W_NB_glm = sm.GLM(obs_y, 
                                  feature_x,
                                  exposure = exposure_x,
                                  family=sm.families.NegativeBinomial(alpha=dispersion_used),
                                  freq_weights = freq_w)
            
                W_NB_res = W_NB_glm.fit(**NB_kwargs)
            except:
                print("[XClone] GLM error:")
                print(obs_y, feature_x, exposure_x, freq_w)
                model_results.append("Error in this gene group")
                tmp_ratio = np.NaN
                CNV_ratio.append(tmp_ratio)
                if feature_X is None:
                    pass
                else:
                    tmp_ratio1 = np.NaN
                    ref_bio_ratio.append(tmp_ratio1)
            else:
                print("[XClone] GLM success:")
                print(obs_y, feature_x,exposure_x, freq_w)
                model_results.append(W_NB_res)
                tmp_ratio = np.exp(W_NB_res.params[0])
                CNV_ratio.append(tmp_ratio)
                if feature_X is None:
                    pass
                else:
                    tmp_ratio1 = np.exp(W_NB_res.params[1])
                    ref_bio_ratio.append(tmp_ratio1)    

    CNV_ratio = np.array(CNV_ratio).reshape(-1, n_states)

    # time stamp
    end_time_ = datetime.datetime.now()
    time_used = end_time_ - start_time_
    print("Time used", time_used.seconds, "seconds")
    
    if feature_X is None:
        return CNV_ratio, model_results
    else:
        ref_bio_ratio = np.array(ref_bio_ratio).reshape(-1, n_states)

        return CNV_ratio, ref_bio_ratio, model_results


def gene_specific_states(Xdata, CNV_ratio):
    """
    Function:
    Parsing_CNV_ratio:parsing the cnv ratio for different gene groups
    generate group genes-specific CNV ratio array for calculating the NB prob.

    Parameters:
    ----------
    Xdata: anndata. uns["gene_group"] provides gene-specific information.
    CNV_ratio: np.array. different CNV states ratio for different gene groups.


    Return:
    ------
    gene-specific CNV states ratio: np.array.


    Example:
    -------
    xclone.model.gene_specific_states(ref_obs_ad1, np.array(CNV_ratio5).reshape(-1, 3))

    """

    ## init
    # n_cells = Xdata.shape[0]
    n_genes = Xdata.shape[1]
    n_states = CNV_ratio.shape[1]

    states_array = np.zeros((n_genes, n_states))

    ## data
    group_genes_lst = Xdata.uns["group_genes"]
    n_gene_group = len(group_genes_lst)
    
    ## GLM model
    for i in range(n_gene_group):
        flag_ = Xdata.var["GeneName"].isin(group_genes_lst[i])
        states_array[flag_, :] = CNV_ratio[i,:]
    return states_array


# Final: CNV Estimation iteration.
## (1) improve the efficiency in each step(completed! 2022-01-20)
## (2) save the iteration results and the estimated optimized CNV ratio.(completed! 2022-03-18)

def check_CNV_ratio(CNV_ratio, CNV_ratio_prev):
    """
    check the nan results in fitted CNV ratio and 
    replace with values in the previous round.
    """
    CNV_ratio_update = CNV_ratio.copy()
    idx_ = np.isnan(CNV_ratio_update).any(axis=1)
    CNV_ratio_update[idx_] = CNV_ratio_prev[idx_]

    return CNV_ratio_update

def combine_BAF_emmprob(Xdata, BAF_Xdata):
    """
    Xdata: RDR
    BAF_Xdata: BAF
    """
    # check BAF_emm_prob
    _gene_flag_ = BAF_Xdata.var["GeneName"].isin(Xdata.var["GeneName"])
    
    Xdata.obs["cell_index"] = Xdata.obs.index
    BAF_Xdata.obs["cell_index"] = BAF_Xdata.obs.index
    _cell_flag_ = BAF_Xdata.obs["cell_index"].isin(Xdata.obs["cell_index"])

    Xdata.uns["BAF_emm_prob_log"] = BAF_Xdata.uns["BAF_emm_prob_log"][:, _gene_flag_,:].copy()
    Xdata.uns["BAF_emm_prob_log"] = Xdata.uns["BAF_emm_prob_log"][_cell_flag_, :, :]
    return Xdata

def CNV_optimazation(Xdata,
                     depth_key = "library_ratio_capped",
                     dispersion_key = "dispersion_capped",
                     init_state_ratio = np.array([0.5, 1.0, 1.5]), 
                     max_iter=20, 
                     min_iter=3, 
                     epsilon_conv=1e-2, 
                     verbose= True,
                     fitCNV_verbose = False,
                     HMM_verbose = False,
                     HMM_brk = "chr_arm",
                     start_prob = None,
                     trans_prob = None,
                     combine_BAF = False,
                     states_weight = np.array([1,2,3]),
                     weights= True,
                     nproc = 1,
                     log_save = False, 
                     **kwargs):
    """
    Function:
    ---------
    iterations of fit_CNV_ratio, the loglikelihood can be used to determine the terminations.

    step1: get the emm_prob_log
    step2: get the posterior_mtx_log
    step3: calculate the log liklihood
    step4: determine the terminations


    Parameters:
    ----------
    Xdata: anndata.
    init_prob:

    **kwargs: parameters for fit_CNV_ratio.

    Return:
    ------
    
    Example:
    -------

    """
    ## time stamp
    start_time_ = datetime.datetime.now()

    CNV_state_num = len(init_state_ratio)

    try:
        gene_group_num = len(Xdata.uns["group_genes"])
    except:
        raise ValueError("[XClone] Xdata do not contain group genes! Pls check and do preprocessing!")

    ## initalization
    ### CNV ratio init for different gene groups
    CNV_ratio_dic = {}
    CNV_ratio_dic[str(0)] = np.tile(init_state_ratio, gene_group_num).reshape(-1, CNV_state_num)
    ### likelihood init
    Logliklihood = np.zeros(max_iter)

    # iteration
    for it in range(max_iter):
        if verbose == True:
            print("[XClone] CNV_optimazation iteration: ", it+1)

        if it == 0: # init round
            specific_states = gene_specific_states(Xdata, CNV_ratio_dic[str(0)])
            Xdata = calculate_Xemm_probTry(Xdata, dispersion_key, 
                                                  states = specific_states, combine_BAF = combine_BAF)
        else:
            ## params setting and fit_CNV_ratio
            fit_cnv_params = {}
            fit_cnv_params['Xdata'] = Xdata
            # fit_cnv_params['init_prob'] = normalize(np.exp(emm_prob_log))
            fit_cnv_params['init_prob_layer'] = "emm_prob_log"

            fit_cnv_params['random_seed'] = 2
            fit_cnv_params['verbose'] = fitCNV_verbose
            fit_cnv_params['depth_key'] = depth_key

            fit_cnv_params.update(**kwargs)
            CNV_ratio, model_results = fit_CNV_ratio(**fit_cnv_params)

            ## check CNV ratio
            CNV_ratio_update  = check_CNV_ratio(CNV_ratio, CNV_ratio_dic[str(it-1)])
            CNV_ratio_dic[str(it)] = CNV_ratio_update

            specific_states = gene_specific_states(Xdata, CNV_ratio_update)
            Xdata = calculate_Xemm_probTry(Xdata, dispersion_key, 
                                                  states = specific_states, combine_BAF = combine_BAF)
        ### posterior_mtx_log
        update_Xdata = XHMM_smoothing(Xdata, brk = HMM_brk, emm_inlayer = "emm_prob_log", 
                                      start_prob = start_prob, trans_prob = trans_prob, 
                                      nproc = nproc, verbose = HMM_verbose)

        if log_save == True:
            CNV_visualization(update_Xdata, states_weight, weights)

        _logLik = cal_log_lik(update_Xdata.layers["emm_prob_log"], update_Xdata.layers["posterior_mtx_log"])
        Logliklihood[it] = _logLik

        if it >= min_iter:
            if Logliklihood[it] < Logliklihood[it - 1]:
                if verbose:
                    print("[XClone] Warning: Lower bound decreases!\n")
                    print("Step %d, loglik decrease from %.2e to %.2e" 
                          %(it, Logliklihood[it-1], Logliklihood[it]))
            elif Logliklihood[it] - Logliklihood[it - 1] <= epsilon_conv:
                break
            elif it == max_iter - 1:
                if verbose:
                    if max_iter > 10:
                        print("[XClone] Warning: CNV ratio optimization did not converge!\n")
                        print("[XClone] Notes: try to increase the max_iter: ", max_iter, "!\n")
    
    Logliklihood = Logliklihood[:it+1]

    if verbose == True:
        print("iteration_end_round: ", it+1)
        print("Logliklihood: ", Logliklihood)
        print("CNV_ratio: ", CNV_ratio_dic)
    
    ## time stamp
    end_time_ = datetime.datetime.now()
    time_used = end_time_ - start_time_
    print("Time used", time_used.seconds, "seconds")
    
    ## Results
    update_Xdata.uns["CNV_ratio"] = CNV_ratio_dic
    update_Xdata.uns["Logliklihood"] = Logliklihood
    return update_Xdata