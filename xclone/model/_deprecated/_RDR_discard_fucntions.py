"""Base functions for XClone RDR processing
deprecated version
"""
# Author: Rongting Huang
# Date: 2022-07-12
# update: 2022-07-12

### discard in _RDR
def V1_rdr_Xdata(Xdata, ref_celltype = None, cell_anno_key = "cell_type", order_cell = True):
    """
    Function:
    deprecated: 2022-04-18
    cell based raw ratio visualization
    version: 0.0.1
    point estimate
    add pseudo count--need prior--[wait for updating]

    pipeline here: filtering ref==0;
    
    Parameters:
    -----------
    Xdata:
    ref_celltype:
    cell_anno_key:
    order_cell: bool.

    Example:
    -------
    
    """
    if ref_celltype is None:
        raise ValueError("ref_celltype should be assigned! Pls check!")
    ref_flag = Xdata.obs[cell_anno_key] == ref_celltype
    ref_Xdata = Xdata[ref_flag,:]
    ref_bulk = ref_Xdata.X.A.sum(axis=0)
    ref_normalization_term = ref_bulk.sum()
    ref_mean = ref_bulk/ref_normalization_term
    print(ref_mean)
    
    gene_flag = ref_mean!=0

    update_Xdata = Xdata[:, gene_flag]
    
    ## update the ref data and obs data
    ref_Xdata = Xdata[ref_flag, gene_flag]
    obs_Xdata = Xdata[~ref_flag, gene_flag]
    ref_mean = ref_mean[gene_flag]

    cell_lib = obs_Xdata.X.A.sum(axis=1, keepdims=True)
    raw_ratio = np.log((obs_Xdata.X.A + 0.1) / (cell_lib*ref_mean)) ## add small value 0.1 for visualization 

    ## anndata
    rr_ad = ad.AnnData(raw_ratio, var=obs_Xdata.var, obs = obs_Xdata.obs)
    
    ## reorder by cell_type for visualization
    if order_cell ==True:
        groups = rr_ad.obs.groupby(cell_anno_key).indices
        rr_ad = ad.concat([rr_ad[inds] for inds in groups.values()], merge="same")
    return rr_ad, update_Xdata

def V2_rdr_Xdata(Xdata, ref_celltype = None, cell_anno_key = "cell_type"):
    """
    Function:
    deprecated: 2022-04-18
    celltype based raw ratio visualization
    version: 0.0.1
    point estimate
    add pseudo count--need prior--[wait for updating]
    
    Parameters:
    ----------

    **kwargs for NB_HMM2(reference, observations, 
           states = np.array([0.5, 1.0]), 
           start_prob = np.array([0.8, 0.2]), 
           trans_prob = np.array([[0.7,0.3],[0.2,0.8]]))
    
    Example:
    -------
    NB_HMM_Xdata(Xdata, ref = "unclassified", **kwargs)
    
    """
    if ref_celltype == None:
        raise ValueError("ref_celltype should be assigned! Pls check!")
    ref_flag = Xdata.obs[cell_anno_key] == ref_celltype
    ref_Xdata = Xdata[ref_flag,:]
    ref_bulk = ref_Xdata.X.A.sum(axis=0)
    ref_normalization_term = ref_bulk.sum()
    ref_mean = ref_bulk/ref_normalization_term
    print(ref_mean)
    
    gene_flag = ref_mean!=0

    update_Xdata = Xdata[:, gene_flag]
    
    ## update the ref data and obs data 1
    ref_Xdata = Xdata[ref_flag, gene_flag]
    obs_Xdata = Xdata[~ref_flag, gene_flag]
    ref_mean = ref_mean[gene_flag]
    ref_bulk = ref_bulk[gene_flag]

    ## step 02-prepare the obs_bulk
    groups = obs_Xdata.obs.groupby(cell_anno_key).indices
    cell_type_order = []
    
    cnt = 0
    for tmp_celltype, inds in groups.items():
        cell_type_order.append(tmp_celltype)
        if cnt==0:
            obs_all_bulk = obs_Xdata[inds,:].X.A.sum(axis=0)
        else:
            obs_all_bulk = np.vstack((obs_all_bulk, obs_Xdata[inds,:].X.A.sum(axis=0)))
        cnt+=1
    
    ## return ref_obs_all_bulk anndata for histogram visualization ## count
    ref_obs_bulk = np.vstack((ref_bulk, obs_all_bulk))
    ref_obs_bulk_anno = cell_type_order.copy()
    ref_obs_bulk_anno.insert(0,"ref_celltype")
    ref_obs_bulk_anno = pd.DataFrame({"cell_type":ref_obs_bulk_anno})
    ref_obs_bulk_ad = ad.AnnData(ref_obs_bulk, var=obs_Xdata.var, obs = ref_obs_bulk_anno)


    ## calculate the raw_logratio

    cell_lib = obs_all_bulk.sum(axis=1, keepdims=True)
    raw_ratio = np.log2((obs_all_bulk+0.1) / (cell_lib*ref_mean)) ## add small value 0.1 for visualization 

    ## anndata
    new_obs =  pd.DataFrame({"cell_type":cell_type_order})
    rr_ad = ad.AnnData(raw_ratio, var=obs_Xdata.var, obs = new_obs)

    return  rr_ad, ref_obs_bulk_ad, update_Xdata



def rr_ad_celltype_processing_x(Xdata, ref_celltype, cell_anno_key):
    """
    version: 0.0.1
    deprecated: 2022-04-18
    """
    ref_flag = Xdata.obs[cell_anno_key] == ref_celltype
    ref_Xdata = Xdata[ref_flag,:]
    # obs_Xdata = Xdata[~ref_flag,:]
    obs_Xdata = Xdata #/ 20220315 for ref_obs_ad2
    
    ref_bulk = ref_Xdata.X.A.sum(axis=0)
    ref_normalization_term = ref_bulk.sum()
    ref_norm = ref_bulk/ref_normalization_term

    # 01-celltype based anndata
    ## prepare the obs_bulk
    groups = obs_Xdata.obs.groupby(cell_anno_key).indices ## dict
    cell_type_order = []
    
    cnt = 0
    for tmp_celltype, inds in groups.items():
        cell_type_order.append(tmp_celltype)
        if cnt==0:
            obs_all_bulk = obs_Xdata[inds,:].X.A.sum(axis=0)
        else:
            obs_all_bulk = np.vstack((obs_all_bulk, obs_Xdata[inds,:].X.A.sum(axis=0)))
        cnt+=1
    
    ## return ref_obs_bulk anndata for histogram visualization
    ref_obs_bulk = np.vstack((ref_bulk, obs_all_bulk))
    ref_obs_bulk_anno = cell_type_order.copy()
    ref_obs_bulk_anno.insert(0,"ref_celltype")
    ref_obs_bulk_anno = pd.DataFrame({"cell_type":ref_obs_bulk_anno})
    ref_obs_bulk_ad = ad.AnnData(ref_obs_bulk, var=obs_Xdata.var.copy(), obs = ref_obs_bulk_anno)

    ## calculate ratio-bulk
    ## return rr_ad_celltype for log raw ratio visualization and choosing normal chrs
    cell_lib = obs_all_bulk.sum(axis=1, keepdims=True) # equivalent normalization term like ref_bulk
    raw_ratio1 = np.log2((obs_all_bulk+0.1) / (cell_lib*ref_norm)) ## add small value 0.1 for visualization 
    celltype_obs =  pd.DataFrame({"cell_type": cell_type_order})
    rr_ad_celltype = ad.AnnData(raw_ratio1, var=obs_Xdata.var.copy(), obs = celltype_obs)

    return ref_obs_bulk_ad, rr_ad_celltype

def rr_ad_cell_processing_x(Xdata, ref_celltype, cell_anno_key, order_cell):
    """
    version: 0.0.1
    deprecated: 2022-04-18
    """
    ref_flag = Xdata.obs[cell_anno_key] == ref_celltype
    ref_Xdata = Xdata[ref_flag,:]
    obs_Xdata = Xdata[~ref_flag,:]
    
    ref_bulk = ref_Xdata.X.A.sum(axis=0)
    ref_normalization_term = ref_bulk.sum()
    ref_norm = ref_bulk/ref_normalization_term
    
    ref_cell_num = ref_flag.sum()
    ref_ave = ref_bulk/ref_cell_num # equal to ref_mean
    
    # ref_mean = ref_Xdata.X.mean(0)

    # Xdata.var['ref_mean'] = ref_Xdata.X.mean(0)
    # Xdata.var['ref_var'] = ref_Xdata.X.var(0)
        
    # 02-cell based anndata
    ## return ref_obs anndata for histogram visualization: ref_obs_ad1, ref_obs_ad2
    ## return rr_ad_cell/rr_ad_cell_ordered for log raw ratio visualization
    
    ### version1 merge all ref cells as ref mean
    if order_cell == True:
        obs_ad_cell = reorder_data_by_cellanno(obs_Xdata, cell_anno_key = cell_anno_key)
    elif order_cell == False:
        obs_ad_cell = obs_Xdata

    ref_obs1 = np.vstack((ref_ave, obs_ad_cell.X.A))
    obs_anno1 = obs_ad_cell.obs.copy()
    obs_anno1.loc["ref_ave"] = "ref_celltype" # insert ref_celltype in the last row
    ref_obs_anno1 = obs_anno1.iloc[np.arange(-1, len(obs_anno1)-1)]

    ref_obs_ad1 = ad.AnnData(ref_obs1, var = obs_Xdata.var.copy(), obs = ref_obs_anno1)


    ### version2 keep all ref cells-original order-add ref_ave in the first rows
    ref_obs2 = np.vstack((ref_ave, Xdata.X.A))
    obs_anno2 = Xdata.obs.copy()
    obs_anno2.loc["ref_ave"] = "ref_celltype" # insert ref_celltype
    ref_obs_anno2 = obs_anno2.iloc[np.arange(-1, len(obs_anno2)-1)]

    ref_obs_ad2 = ad.AnnData(ref_obs2, var=Xdata.var.copy(), obs = ref_obs_anno2)

    ## calculate ratio-cell
    cell_lib = obs_Xdata.X.A.sum(axis=1, keepdims=True)
    raw_ratio2 = np.log((obs_Xdata.X.A + 0.1) / (cell_lib*ref_norm)) ## add small value 0.1 for visualization 

    rr_ad_cell = ad.AnnData(raw_ratio2, var=obs_Xdata.var.copy(), obs = obs_Xdata.obs.copy())
    
    print("output anndata is not sparse matrix.")
    ### reorder by cell_type for visualization
    if order_cell == True:
        rr_ad_cell_ordered = reorder_data_by_cellanno(rr_ad_cell, cell_anno_key = cell_anno_key)
        return ref_obs_ad1, ref_obs_ad2, rr_ad_cell_ordered
    elif order_cell == False:
        return ref_obs_ad1, ref_obs_ad2, rr_ad_cell

## discard in HMM_base
## stratege 1
def NB_HMM1(reference, observations, 
           states = np.array([0.5, 1.0]), 
           start_prob = np.array([0.8, 0.2]), 
           trans_prob = np.array([[0.7,0.3],[0.2,0.8]])):
    """
    test version
    reference = np.array([10, 10, 10,10,10, 10,10])
    observations = np.array([5., 9, 3, 5, 12, 6, 9])
    """
    emm_prob_log = generate_nb_emmission_logprob(reference, observations, states)
    ## check emm_prob
    emm_flag = emm_prob == 0
    emm_flag_ = emm_flag.any(axis=1)
    if emm_flag_.sum() != 0:
        # raise ValueError("[XClone]-pls check emm prob")
        print("ValueError: [XClone]-pls check emm prob")
        return emm_flag_
    ## run fwd_bwd_prob
    states_num = len(states)
    obs_num = len(observations)
    res = fwd_bkw_prob1(obs_num, states_num, start_prob, trans_prob, emm_prob_log)
    return res, emm_prob_log

t = 1e-6
def NB_HMM2(reference, observations, 
           states = np.array([0.5, 1.0, 1.5]), 
           start_prob = np.array([0.15, 0.7, 0.15]), 
           trans_prob = np.array([[1-2*t, t, t],[t, 1-2*t, t],[t, t, 1-2*t]])):
    """
    reference = np.array([10, 10, 10,10,10, 10,10])
    observations = np.array([5., 9, 3, 5, 12, 6, 9])
    """
    emm_prob_log = generate_nb_emmission_logprob(reference, observations, states)
    ## check emm_prob
    emm_flag = emm_prob == 0
    emm_flag_ = emm_flag.any(axis=1)
    if emm_flag_.sum() != 0:
        # raise ValueError("[XClone]-pls check emm prob")
        print("ValueError: [XClone]-pls check emm prob")
        return emm_flag_
    ## run fwd_bwd_prob
    states_num = len(states)
    obs_num = len(observations)
    res = fwd_bkw_prob2(obs_num, states_num, start_prob, trans_prob, emm_prob_log)
    return res, emm_prob

### for xdata format to run the HMM
import anndata as ad
import itertools
import numpy as np
import pandas as pd

## stratege 1
### calculate emm_prob for each celltype

## todo check!
def cal_raw_logratio(ref_all_bulk, obs_all_bulk):
    """
    celltype based bulk
    """
    ref_normalization_term = ref_all_bulk.sum()
    _GEX_ref = ref_all_bulk/ref_normalization_term
    obs_normalization_term = obs_all_bulk.sum(axis=1)
    raw_ratio = obs_all_bulk / (_GEX_ref * obs_normalization_term[:,np.newaxis])
    raw_logratio = np.log(raw_ratio)

    return raw_logratio


### from _RDR_CNVratio
def fit_CNV_ratio_test(Xdata, init_prob, hard_assign = False, n_sample_cells = 10, replace = False,
                  dispersion_set = None, random_seed = None, verbose=True, test_func = False,
                  NB_kwargs={'disp': True, 'skip_hessian': True}):
    """
    Function:
    version 0.0.0 [seems more efficient, can be used]
    deprecated-- no libratio counted in[2022-02-25]
    based on selected group genes and GLM model to get the estimate CNV ratio
    for different state, default state is copy loss, neutral, copy gain.

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
    dispersion_set: 
    random_seed:       
    verbose : bool
            Whether print out log info
    test_func: bool.

    
    Return:
    ------
    Xdata: anndata.


    Example:
    -------

    """

    if random_seed is not None:
        rvgen = np.random.RandomState(random_seed)
    else:
        rvgen = np.random
    
    # time stamp
    start_time_ = datetime.datetime.now()

    if hard_assign == True:
        init_prob = hard_assign_prob(init_prob)
        # if test_func:
        #     print("hard_assign prob",init_prob)

    # check input
    if Xdata.shape[0] - 1 == init_prob.shape[0]:
        if Xdata.shape[1] == init_prob.shape[1]:
            if verbose:
                print("input data is in matched dimension")
        else:
            raise ValueError("[XClone] gene dimension is not matched! Pls check!")
    else:
        raise ValueError("[XClone] cell dimension is not matched! Pls check!")
    
    ## data
    group_genes_lst = Xdata.uns["group_genes"]
    n_gene_group = len(group_genes_lst)
    n_states = init_prob.shape[2]
    n_cells = Xdata.shape[0]-1

    model_results = []
    CNV_ratio = []

    ## GLM model
    for i in range(n_gene_group):
        flag_ = Xdata.var["GeneName"].isin(group_genes_lst[i])

        ref_vec = Xdata.X[0, flag_]
        tmp_mtx = Xdata.X[1:,flag_]
        tmp_prob = init_prob[:,flag_, :]
        if dispersion_set is None:
            dispersion_ = Xdata.var["dispersion"][flag_].median()
            print(dispersion_)
        else:
            dispersion_ = dispersion_set

        for gene_idx in range(tmp_mtx.shape[1]):
            # random cell index
            # rs = np.ceil(rvgen.rand(n_sample_cells) * (n_cells-2)).astype(int)
            if replace == False:
                rs = rvgen.choice(n_cells, n_sample_cells, replace = False)
            else:
                rs = rvgen.randint(0, n_cells, n_sample_cells)

            if gene_idx == 0:
                ref_exp_ = np.repeat(ref_vec[gene_idx], n_sample_cells)
                gene_exp_ = tmp_mtx[rs, gene_idx]
                freq_prob_ = tmp_prob[rs, gene_idx, :]
                if test_func:
                    print("test", freq_prob_.shape)
                    print(freq_prob_)
            else:
                ref_exp_ = np.vstack([ref_exp_, np.repeat(ref_vec[gene_idx], n_sample_cells)])
                gene_exp_ = np.vstack([gene_exp_, tmp_mtx[rs, gene_idx]])
                freq_prob_ = np.vstack([freq_prob_, tmp_prob[rs, gene_idx, :]])

        obs_y = gene_exp_.reshape(-1)
        feature_x = np.ones(len(obs_y))
        exposure_x = ref_exp_.reshape(-1)


        for state_idx in range(n_states):
            freq_w = freq_prob_[:, state_idx]
            print(freq_w)

            try:
            # if True:
                W_NB_glm = sm.GLM(obs_y, 
                                  feature_x,
                                  exposure = exposure_x,
                                  family=sm.families.NegativeBinomial(alpha=dispersion_),
                                  freq_weights = freq_w)
            
                W_NB_res = W_NB_glm.fit(**NB_kwargs)
            except:
                model_results.append("Error in this gene group")
                tmp_ratio = np.NaN
                CNV_ratio.append(tmp_ratio)
            else:
                model_results.append(W_NB_res)
                tmp_ratio = np.exp(W_NB_res.params[0])
                CNV_ratio.append(tmp_ratio)

    CNV_ratio = np.array(CNV_ratio).reshape(-1, n_states)

    # time stamp
    end_time_ = datetime.datetime.now()
    time_used = end_time_ - start_time_
    print("Time used", time_used.seconds, "seconds")

    return CNV_ratio, model_results
    # return Xdata

## from RDR_CNVratio.
def CNV_optimazation2(Xdata, 
                    init_state_ratio = np.array([0.5, 1.0, 1.5]), 
                    max_iter=20, 
                    min_iter=3, 
                    epsilon_conv=1e-2, 
                    verbose=True,  
                    nproc = 1,
                    log_save = False, 
                    **kwargs):
    """
    Function:(todo func)
    iterations of fit_CNV_ratio, the loglikelihood can be used to determine the terminations.

    step1: get the emm_prob_log
    step2: get the posterior_mtx_log
    step3: calculate the log liklihood
    step4: determine the terminations

    # init_CNV_ratio = np.array([[0.5,1.0,1.5],
    #                  [0.5,1.0,1.5],
    #                  [0.5,1.0,1.5],
    #                  [0.5,1.0,1.5],
    #                  [0.5,1.0,1.5]])

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
        gene_group_num = len(Xdata.uns["group_genes"]) # len(ref_obs_ad1.uns["group_genes"])
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
            # specific_states = xclone.model.gene_specific_states(ref_obs_ad2_test, init_CNV_ratio)
            specific_states = gene_specific_states(Xdata, CNV_ratio_dic[str(0)])
            emm_prob_log = calculate_Xemm_prob2(Xdata, 
                                            states = specific_states,
                                            gene_specific = True, overdispersion = Xdata.var["dispersion"],
                                            ref_normalization_term = 1,
                                            obs_normalization_term = Xdata.obs["library_ratio"][1:])
        else:
            ## params setting and fit_CNV_ratio
            fit_cnv_params = {}
            fit_cnv_params['Xdata'] = Xdata
            fit_cnv_params['init_prob'] = normalize(np.exp(emm_prob_log))
            fit_cnv_params['random_seed'] = 2
            fit_cnv_params['verbose'] = verbose
            fit_cnv_params['feature_X'] = Xdata.uns["NMF_pred_obs"]


            fit_cnv_params.update(**kwargs)
            # CNV_ratio,  model_results = xclone.model.fit_CNV_ratio(ref_obs_ad2_test, xclone.model.normalize(np.exp(emm_prob_log)), 
            # n_sample_cells = 10, dispersion_set = None, random_seed = 2, verbose=True)
            # CNV_ratio, model_results = fit_CNV_ratio(**fit_cnv_params)
            CNV_ratio, ref_bio_ratio, model_results = fit_CNV_ratio(**fit_cnv_params)

            ## check CNV ratio
            CNV_ratio_update  = check_CNV_ratio(CNV_ratio, CNV_ratio_dic[str(it-1)])
            CNV_ratio_dic[str(it)] = CNV_ratio_update

            specific_states = gene_specific_states(Xdata, CNV_ratio_update)
            emm_prob_log = calculate_Xemm_prob2(Xdata, 
                                            states = specific_states,
                                            gene_specific = True, overdispersion = Xdata.var["dispersion"],
                                            ref_normalization_term = 1,
                                            obs_normalization_term = Xdata.obs["library_ratio"][1:])

        ### posterior_mtx_log
        update_Xdata = XHMM_smoothing(Xdata, emm_prob_log = emm_prob_log, nproc = nproc)

        if log_save == True:
            CNV_visualization(update_Xdata)

        # check emm_prob_log and res_log shape.
        _logLik = cal_log_lik(update_Xdata.uns["emm_prob_log"], update_Xdata.uns["posterior_mtx_log"])
        # _logLik = cal_log_lik(emm_prob_log, posterior_mtx_log)
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
    
    ## time stamp
    end_time_ = datetime.datetime.now()
    time_used = end_time_ - start_time_
    print("Time used", time_used.seconds, "seconds")
    
    ## Results
    update_Xdata.uns["CNV_ratio"] = CNV_ratio_dic
    update_Xdata.uns["Logliklihood"] = Logliklihood
    return update_Xdata