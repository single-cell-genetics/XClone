"""Base functions for XClone RDR HMM based emm_prob generation
some testing record
deprecated maybe
"""

# Author: Rongting Huang
# Date: 2021/09/20
# update: 2021/12/17
import numpy as np
from scipy.stats import nbinom
from scipy.stats import poisson
from scipy.stats import betabinom

from .base_utils import normalize, loglik_amplify


# source- HMM_base.py
## Emmission_prob
### For RDR
def generate_nb_emmission_logprob_v1(reference, observations, overdispersion = None,
                                  ref_normalization_term = None, obs_normalization_term = None, 
                                  states = np.array([0.5, 1.0])):
    """
    version 0.0.1
    deprecated
    emission prob log format - gene specific distribution-params
    negative binomial distributions[scipy]

    states: default np.array([0.5, 1.0]) 
    # example for copy loss and neutral

    Example
    ------
    reference = np.array([10, 10, 10,10,10, 10,10])
    # reference = np.array([5., 9, 3, 5, 12, 6, 9])
    observations = np.array([5., 9, 3, 5, 12, 6, 9])
    states = np.array([0.5, 1.0])
    emm_prob_log = generate_nb_emmission_logprob(reference, observations, states = states)
    """
    emm_prob_log = np.zeros((len(observations),len(states)))
    
    ## Normally, the normalization term is for cell library, but sometimes we generate the emm_prob for each chromosome
    ## so need specify the normalization term in the cell scale if the input(reference and observations) are chr specific.
    
    ## default None means the input ref and obs are genome scale-can achieve the lib size via sum.
    if ref_normalization_term == None:
        ref_normalization_term = reference.sum()
    if obs_normalization_term == None:
        obs_normalization_term= observations.sum()
    if overdispersion == None:
        overdispersion = 0.1
    
    for idx, obs_ in enumerate(observations):
        _GEX_ref = reference[idx]/ref_normalization_term
        _mu = states * _GEX_ref * obs_normalization_term
        _var = _mu + overdispersion * _mu**2 # example over dispersion # todo
        ## todo if the overdispersion is gene specific, then use overdispersion[idx]
         
        # # reparameterization in tensorflow negative binomial
        # _nb_prob = 1 - _mu / _var
        # _nb_total = _mu * (1 - _nb_prob) / _nb_prob
        
        # reparameterization in scipy negative binomial
        _nb_prob =  _mu / _var
        _nb_total = _mu * _nb_prob / (1 - _nb_prob)
        emm_prob_log[idx, :] = nbinom.logpmf(obs_, _nb_total, _nb_prob)

    print("generate a emm_prob matrix for %d states, matrix shape is %s" %(len(states), str(emm_prob_log.shape)))
    return emm_prob_log


def generate_nb_emmission_logprob(reference, observations, overdispersion = None, gene_specific= False,
                                  ref_normalization_term = None, obs_normalization_term = None, 
                                  states = np.array([0.5, 1.0])):
    """

    update version for numpy array-observations one-dim array
    emission prob log format - gene specific distribution-params
    negative binomial distributions[scipy]

    states: default np.array([0.5, 1.0]) 
    # example for copy loss and neutral

    Example
    ------
    reference = np.array([10, 10, 10,10,10, 10,10])
    # reference = np.array([5., 9, 3, 5, 12, 6, 9])
    observations = np.array([5., 9, 3, 5, 12, 6, 9])
    states = np.array([0.5, 1.0])
    emm_prob_log = generate_nb_emmission_logprob(reference, observations, states = states)

    emm_prob_log_update = generate_nb_emmission_logprob_update(reference, observations, states = states) # default overdispersion 0.1
    emm_prob_log_update = generate_nb_emmission_logprob_update(reference, observations, overdispersion = 0.2, states = states)

    overdispersion = np.array([0.1, 0.2, 0.3, 0.1, 0.2, 0.3, 0.1])
    emm_prob_log_update = generate_nb_emmission_logprob_update(reference, observations, overdispersion, gene_specfic= True, states = states)
    """
    # emm_prob_log = np.zeros((len(observations),len(states)))
    
    ## Normally, the normalization term is for cell library, but sometimes we generate the emm_prob for each chromosome
    ## so need specify the normalization term in the cell scale if the input(reference and observations) are chr specific.
    
    ## default None means the input ref and obs are genome scale-can achieve the lib size via sum.
    if ref_normalization_term == None:
        ref_normalization_term = reference.sum()
    if obs_normalization_term == None:
        obs_normalization_term= observations.sum()
    
    if gene_specific:
        if len(overdispersion) == len(observations):
            pass
        else:
            raise ValueError("[XClone]-pls check overdispersion")
    elif overdispersion == None:
        overdispersion = 0.1
    
    
    _GEX_ref = reference/ref_normalization_term

    _mu = states * _GEX_ref[:,np.newaxis] * obs_normalization_term

    ## if overdispersion is gene specfic, then use overdispersion[:,np.newaxis]
    if gene_specific:
        _var = _mu + overdispersion[:,np.newaxis] * _mu**2
    else:
        _var = _mu + overdispersion * _mu**2
    

    _nb_prob =  _mu / _var
    _nb_total = _mu * _nb_prob / (1 - _nb_prob)

    obs_ = np.repeat(observations, len(states)).reshape(len(observations),len(states))

    emm_prob_log = nbinom.logpmf(obs_, _nb_total, _nb_prob)

    print("generate a emm_prob matrix for %d states, matrix shape is %s" %(len(states), str(emm_prob_log.shape)))
    return emm_prob_log

## 
def generate_nb_emmission_logprob_update(reference, observations, overdispersion = None, gene_specific= False,
                                  ref_normalization_term = None, obs_normalization_term = None, 
                                  states = np.array([0.5, 1.0])):
    """
    Function: 
    negative binomial distributions[scipy]

    update version for numpy array-
    emission prob log format - gene specific distribution-params

    Parameters
    ----------
    reference:
    observations:

    states: default np.array([0.5, 1.0]) # example for copy loss and neutral


    Returns
    -------
    
    Example
    -------
    reference = np.array([10, 10, 10,10,10, 10,10])
    # reference = np.array([5., 9, 3, 5, 12, 6, 9])
    observations = np.array([[5., 9, 3, 5, 12, 6, 9]])
    observations = np.array([[5., 9, 3, 5, 12, 6, 9],[5., 9, 3, 5, 12, 6, 9]])

    states = np.array([0.5, 1.0])
    emm_prob_log = generate_nb_emmission_logprob_update(reference, observations, states = states)

    emm_prob_log_update = generate_nb_emmission_logprob_update(reference, observations, states = states) # default overdispersion 0.1
    emm_prob_log_update = generate_nb_emmission_logprob_update(reference, observations, overdispersion = 0.2, states = states)

    overdispersion = np.array([0.1, 0.2, 0.3, 0.1, 0.2, 0.3, 0.1])
    emm_prob_log_update = generate_nb_emmission_logprob_update(reference, observations, overdispersion, gene_specfic= True, states = states)
    """
    # emm_prob_log = np.zeros((len(observations),len(states)))
    
    ## Normally, the normalization term is for cell library, but sometimes we generate the emm_prob for each chromosome
    ## so need specify the normalization term in the cell scale if the input(reference and observations) are chr specific.

    ## best way is calculate the emm_prob at first, in the cell scale
    
    ## default None means the input ref and obs are genome scale-can achieve the lib size via sum.
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
    ## to do
    for i in range(observations.shape[0]):
        if i == 0:
            _GEX_ref_extend =  _GEX_ref[np.newaxis,:]
        else:
            _GEX_ref_extend =  np.vstack((_GEX_ref_extend, _GEX_ref[np.newaxis,:]))
    
    _mu =  _GEX_ref_extend * obs_normalization_term[:,np.newaxis,np.newaxis]

    ## if overdispersion is gene specfic, then use overdispersion[:,np.newaxis]
    if gene_specific:
        _var = _mu + overdispersion[np.newaxis, :, np.newaxis] * _mu**2
    else:
        _var = _mu + overdispersion * _mu**2
    

    _nb_prob =  _mu / _var
    _nb_total = _mu * _nb_prob / (1 - _nb_prob)


    for i in range(observations.shape[0]):
        obs_tmp = np.repeat(observations[i], len(states)).reshape(observations.shape[1],len(states))
        if i == 0:
            obs_ = obs_tmp[np.newaxis,:]
        else:
            obs_ = np.vstack((obs_, obs_tmp[np.newaxis,:]))

    emm_prob_log = nbinom.logpmf(obs_, _nb_total, _nb_prob)

    print("generate a emm_prob matrix for %d states, matrix shape is %s" %(len(states), str(emm_prob_log.shape)))
    return emm_prob_log


def generate_poisson_emmission_logprob(reference, observations, states = np.array([0.5, 1.0])):
    ## emission prob-gene specific distribution-params
    ## poissondistributions
    emm_prob_log = np.zeros((len(observations),len(states)))
    for idx, obs_ in enumerate(observations):
        _GEX_ref = reference[idx]
        _mu = states * _GEX_ref # example for copy loss and neutral

        emm_prob_log[idx, :] = poisson.logpmf(obs_, _mu)
    print("generate a emm_prob matrix for %d states, matrix shape is %s" %(len(states), str(emm_prob.shape)))
    return emm_prob_log


### for BAF
def generate_bb_emmission_prob(observations_AD, observations_DP, 
                               alpha = np.array([0.5, 20, 60]) , 
                               beta =np.array([0.5, 20, 30]), 
                               states=np.array([0.5, 1.0, 1.5])):
    """
    ## seems not gene/ block specific distributions
    ## to do list[for testing]

    observations  [AD, DP]

    copy loss BAF=0 or 1, 
    copy gain BAF=1/3 or 2/3
    copy neutral BAF=0.5
    # """
    states_num = len(states)
    x = np.repeat(observations_AD, states_num).reshape(-1, states_num)
    n = np.repeat(observations_DP, states_num).reshape(-1, states_num)

    emm_prob = betabinom.pmf(x, n, alpha, beta)
    return emm_prob


    
from ..plot._data import transform_anndata

def NB_HMM_Xdata_typebased(Xdata,
                           ref_celltype = None, 
                           cell_anno_key = "cell_type", 
                           brk = "chr_arm", 
                           states = None,
                           overdispersion = 0.1,
                           gene_specific= False,
                           ref_normalization_term = None, 
                           obs_normalization_term = None, 
                           filter_emm_prob = False, 
                           filter_markergenes = False, 
                           gene_dis_cutoff = 0.8,
                           emm_prob_log = None,
                           **kwargs):
    """
    Funcrtion:
    ---------

    brk can be "chr"

    **kwargs for NB_HMM2(reference, observations, 
           states = np.array([0.5, 1.0]), 
           start_prob = np.array([0.8, 0.2]), 
           trans_prob = np.array([[0.7,0.3],[0.2,0.8]]))
    
    Example:
    -------
    NB_HMM_Xdata(Xdata, ref_celltype = "unclassified", **kwargs)
    
    """
    # default states
    if states == None:
        states = np.array([0.5, 1.0, 1.5])
    
    ref_flag = Xdata.obs[cell_anno_key] == ref_celltype
    ref_Xdata = Xdata[ref_flag,:]
    ref_bulk = ref_Xdata.X.A.sum(axis=0)

    ## step01: filter genes in ref_bulk = 0
    gene_flag = ref_bulk!=0
    ### update the ref data and obs data 1
    ref_Xdata = Xdata[ref_flag, gene_flag]
    obs_Xdata = Xdata[~ref_flag, gene_flag]
    ref_bulk = ref_bulk[gene_flag]
    print(ref_bulk)
    
    # ## filter markergenes use KL divergency/JS divergency [for each celltype]
    # if filter_markergenes:
    #     groups = obs_Xdata.obs.groupby(cell_anno_key).indices
    #     filter_flag_lst = []
    #     for tmp_celltype, inds in groups.items():
    #         tmp_obs = obs_Xdata[inds]
    #         obs_bulk = tmp_obs.X.A.sum(axis=0)
    #         tmp_emm_prob = generate_nb_emmission_prob(ref_bulk, obs_bulk, overdispersion = overdispersion, states=states)
        
    #         ## filter emm_prob 0 value
    #         if filter_emm_prob:
    #             gene_flag1 = (tmp_emm_prob == 0).any(axis=1)
    #         else:
    #         ## not filter 0
    #             gene_flag1 = (tmp_emm_prob != tmp_emm_prob).any(axis=1)
        
    #         ## filter celltype specific marker gene
    #         tmp_gene_distance = JSD_emm_genes(tmp_emm_prob)
    #         gene_flag2 = (tmp_gene_distance > gene_dis_cutoff)
        
    #         filter_flag = np.vstack((gene_flag1,gene_flag2)).any(axis=0)
    #         print("filter %d genes for %s" %(filter_flag.sum(), tmp_celltype))
        
    #         filter_flag_lst.append(filter_flag)
    
    #     gene_filter_flag = np.vstack(filter_flag_lst).any(axis=0)
    
    #     ## update the ref data and obs data 2
    #     ref_Xdata = ref_Xdata[:, ~gene_filter_flag]
    #     obs_Xdata = obs_Xdata[:, ~gene_filter_flag]
    
    
    # groups = obs_Xdata.obs.groupby(cell_anno_key).indices

    # filter_flag_lst = []
    # emm_prob_log_dic = {}
    # for tmp_celltype, inds in groups.items():
    #     print(tmp_celltype)
    #     # calculate emm_prob across whole genome (not the selected chrs)
    #     ref_all_bulk = ref_Xdata.X.A.sum(axis=0)
    #     obs_all_bulk = obs_Xdata[inds,:].X.A.sum(axis=0)
        
    #     ## calculate the emm_prob across the genome (not the chr specific region)
    #     tmp_emm_prob_log = generate_nb_emmission_logprob(ref_all_bulk, obs_all_bulk, overdispersion = overdispersion, states=states)
    #     emm_prob_log_dic[tmp_celltype] = tmp_emm_prob_log

    #     ### filter emm_prob value equal to 0
    #     if filter_emm_prob:
    #         gene_flag_ = (tmp_emm_prob_log == -np.inf).any(axis=1)
    #     else:
    #         ## not filter 0
    #         gene_flag_ = (tmp_emm_prob_log != tmp_emm_prob_log).any(axis=1)
        
    #     filter_flag_lst.append(gene_flag_)
    #     # print(filter_flag_lst.sum())
    
    # gene_filter_flag = np.vstack(filter_flag_lst).any(axis=0)
    # print(gene_filter_flag.sum())    
    
    # ## step04-
    # ## update the ref data and obs data
    # ref_Xdata = ref_Xdata[:, ~gene_filter_flag]
    # obs_Xdata = obs_Xdata[:, ~gene_filter_flag]

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

    if obs_all_bulk.ndim != 2:
        print("ndim of obs bulk:", obs_all_bulk.ndim)
    elif obs_all_bulk.ndim == 1:
        obs_all_bulk = obs_all_bulk[np.newaxis,:]
    
    ## step03-calculate emm_prob across whole genome (not the selected chrs)
    ref_all_bulk = ref_Xdata.X.A.sum(axis=0)

    ## calculate the emm_prob across the genome (not the chr specific region)
    if emm_prob_log is not None:
        print("RDR VB emm_prob_log:", emm_prob_log.shape)
    else:
        emm_prob_log = generate_nb_emmission_logprob_update(ref_all_bulk, obs_all_bulk, 
                                                            overdispersion = overdispersion, 
                                                            gene_specific= gene_specific, 
                                                            ref_normalization_term = ref_normalization_term, 
                                                            obs_normalization_term = obs_normalization_term,
                                                            states=states)
    ## calculate raw ratio
    raw_ratio = cal_raw_logratio(ref_all_bulk, obs_all_bulk)
    new_var = obs_Xdata.var
    new_obs = pd.DataFrame({"cell_type":cell_type_order})
    rr_ad = transform_anndata(raw_ratio, new_var, new_obs)
                           
    ## step04-run HMM ## todo -XC_HMM
    brk_item = obs_Xdata.var[brk].drop_duplicates(keep="first")
    res_dict = {}
    res_log_dict = {}
    for i in range(emm_prob_log.shape[0]):
        tmp_celltype = cell_type_order[i]
        tmp_emm_prob_log = emm_prob_log[i]
        cnt = 0
        for brk_ in brk_item:
            tmp_region_flag = obs_Xdata.var[brk] == brk_
            tmp_res, emm_prob_log_ = XC_HMM(tmp_emm_prob_log[tmp_region_flag,:])
            if cnt == 0:
                res_ = tmp_res[0]
                res_log_ = tmp_res[1]
            else:
                res_ = np.vstack((res_, tmp_res[0]))
                res_log_ = np.vstack((res_log_, tmp_res[1]))
            cnt+=1
        
        ## dic record
        res_dict[tmp_celltype] = res_
        res_log_dict[tmp_celltype] = res_log_

        ## numpy array record for loglik cal
        if i == 0:
            res = res_[np.newaxis,:]
            res_log = res_log_[np.newaxis, :]
        else:
            res = np.vstack((res, res_[np.newaxis,:]))
            res_log = np.vstack((res_log, res_log_[np.newaxis, :]))

        
    ## calculate the likelihood as ref of overdispersion grid search
    log_lik = cal_log_lik(emm_prob_log, res_log)


    ### for output
    v_ad = ad.concat([ref_Xdata, obs_Xdata], merge="same")
    print("used_data:" , v_ad.X.shape)

    
    # brk_item = obs_Xdata.var[brk].drop_duplicates(keep="first")
    # res_dict = {}
    # res_log_dict = {}
    
    # for tmp_celltype, inds in groups.items():
    #     print(tmp_celltype)
    #     tmp_emm_prob_log = emm_prob_log_dic[tmp_celltype][~gene_filter_flag, :]

    #     cnt = 0
    #     for brk_ in brk_item:
    #         tmp_region_flag = obs_Xdata.var[brk] == brk_
    #         tmp_res, emm_prob_log_ = XC_HMM(tmp_emm_prob_log[tmp_region_flag,:])
    #         if cnt == 0:
    #             res_ = tmp_res[0]
    #             res_log = tmp_res[1]

    #         else:
    #             res_ = np.vstack((res_, tmp_res[0]))
    #             res_log = np.vstack((res_log, tmp_res[1]))
    #         cnt+=1
    #     res_dict[tmp_celltype] = res_
    #     res_log_dict[tmp_celltype] = res_log


    ## output anndata for visualization 
    ### convert the prob dict to ad
    cnt = 0
    celltype_lst = [] # in a order
    for celltype, prob_ in res_dict.items():
        celltype_lst.append(celltype)
        if cnt==0:
            res_prob = prob_.T
            res_cnv = np.argmax(prob_.T, axis=0)
        else:
            res_prob = np.vstack([res_prob, prob_.T])
            res_cnv = np.vstack([res_cnv, np.argmax(prob_.T, axis=0)])
        cnt+=1
    
    K = emm_prob_log.shape[2] # states num
    prob_celltype = list(itertools.chain.from_iterable(itertools.repeat(i, K)
                                                    for i in celltype_lst))
    obs_df1 = pd.DataFrame(index = prob_celltype)
    res_prob_ad =  ad.AnnData(res_prob, var=obs_Xdata.var, obs = obs_df1)
    res_prob_ad.obs_names_make_unique()
    res_prob_ad.obs["celltype"] = res_prob_ad.obs.index
    
    obs_df2 = pd.DataFrame(index = celltype_lst)
    res_cnv_ad = ad.AnnData(res_cnv, var=obs_Xdata.var, obs = obs_df2)
    
    ## check nan value
    nan_count = np.isnan(res_prob_ad.X).sum() 
    prob_check_nan = nan_count!=0
    if prob_check_nan:
        print("there are %d nan value in the prob mtx" %(nan_count))
    # return raw_ratio, rr_ad, log_lik, res_prob_ad, res_cnv_ad, res_dict, res_log_dict, emm_prob_log, v_ad
    ## raw_ratio is rr_ad.X
    return rr_ad, log_lik, res_prob_ad, res_cnv_ad, res_dict, res_log_dict, emm_prob_log, v_ad


## stratege 2
### calculate emm_prob for each cell and share the state for each celltype

def NB_HMM_Xdata_cellbased(Xdata, 
                            ref_celltype = "ref_celltype", 
                            cell_anno_key = "cell_type", 
                            brk = "chr_arm", 
                            filter_emm_prob = True, 
                            filter_markergenes = True, 
                            gene_dis_cutoff = 0.8, 
                            states = None,
                            emm_prob_log = None, 
                            **kwargs):
    """

    simple version 20211103
    input the emm_prob_log from RDR VB model
    Function:
    ---------
    
    Example:
    -------
    NB_HMM_Xdata_cellbased(Xdata, ref_celltype = "unclassified", **kwargs)
    
    
    OUTPUT:
    -------
    cellbased for visulization
    """

    ref_flag = Xdata.obs[cell_anno_key] == ref_celltype

    ### update the ref data and obs data
    ref_Xdata = Xdata[ref_flag,:]
    obs_Xdata = Xdata[~ref_flag,:]

    # # print(obs_Xdata.shape)
    # # print(emm_prob_log.shape[0])
    # print(obs_Xdata.obs.index[0])
    # return True

    ## run HMM based on emm_prob_log 
    ## todo -XC_HMM
    brk_item = obs_Xdata.var[brk].drop_duplicates(keep="first")
    res_dict = {}
    res_log_dict = {}
    for i in range(emm_prob_log.shape[0]):
        # tmp_celltype = cell_type_order[i]
        tmp_cell_id = obs_Xdata.obs.index[i]
        tmp_emm_prob_log = emm_prob_log[i]
        cnt = 0
        for brk_ in brk_item:
            tmp_region_flag = obs_Xdata.var[brk] == brk_
            tmp_res, emm_prob_log_ = XC_HMM(tmp_emm_prob_log[tmp_region_flag,:])
            if cnt == 0:
                res_ = tmp_res[0]
                res_log_ = tmp_res[1]
            else:
                res_ = np.vstack((res_, tmp_res[0]))
                res_log_ = np.vstack((res_log_, tmp_res[1]))
            cnt+=1
        
        ## dic record
        res_dict[tmp_cell_id] = res_
        res_log_dict[tmp_cell_id] = res_log_

    #     ## numpy array record for loglik cal
    #     if i == 0:
    #         res = res_[np.newaxis,:]
    #         res_log = res_log_[np.newaxis, :]
    #     else:
    #         res = np.vstack((res, res_[np.newaxis,:]))
    #         res_log = np.vstack((res_log, res_log_[np.newaxis, :]))

        
    # ## calculate the likelihood as ref of overdispersion grid search
    # log_lik = cal_log_lik(emm_prob_log, res_log)
    
    return res_dict, res_log_dict


import scipy as sp
def NB_HMM_Xdata_cellbased_test(Xdata, 
                                ref_celltype = None, 
                                cell_anno_key = "cell_type",
                                gene_specific = True,
                                overdispersion = None,
                                ref_normalization_term = None,
                                obs_normalization_term = None,
                                brk = "chr_arm",
                                states = None,
                                emm_prob_log = None, 
                                **kwargs):
    """
    Version: 20211117
    -------

    Function:
    ---------

    Parameters
    ----------
    Xdata:
    ref_celltype:

    **kwargs for NB_HMM2(reference, observations, 
           states = np.array([0.5, 1.0]), 
           start_prob = np.array([0.8, 0.2]), 
           trans_prob = np.array([[0.7,0.3],[0.2,0.8]]))
    Return
    ------
    cellbased for visulization
    
    Example:
    -------
    NB_HMM_Xdata(Xdata, ref = "unclassified", **kwargs)
    """
    ## default states
    if states == None:
        states = np.array([0.5, 1.0, 1.5])
    
    ## data preparation
    ref_flag = Xdata.obs[cell_anno_key] == ref_celltype
    
    ### update the ref data and obs data

    ref_Xdata = Xdata[ref_flag,:]
    obs_Xdata = Xdata[~ref_flag,:]

    if sp.sparse.issparse(Xdata.X):
        ref_mtx = ref_Xdata.X.A
        obs_mtx = obs_Xdata.X.A
    else:
        ref_mtx = ref_Xdata.X
        obs_mtx = obs_Xdata.X

    ref_bulk = ref_mtx.sum(axis=0)
    ref_cell_num = ref_mtx.shape[0]

    if ref_normalization_term is None:
        ref_normalization_term = 1

    ### prepare the input for generate_nb_emmission_logprob
    ref_ = ref_bulk/ref_cell_num
    obs_ = obs_mtx

    print(ref_)
    print(obs_)

    #### check the shape for the input
    if gene_specific == True:
        if len(overdispersion) != len(ref_):
            raise ValueError("[XClone]overdispersion doesn't match with the gene numbers! Pls check!")

    ## calculate emm_prob_log
    if emm_prob_log is not None:
        print("RDR VB emm_prob_log:", emm_prob_log.shape)
    else:
        emm_prob_log = generate_nb_emmission_logprob_update(ref_, 
                                                            obs_, 
                                                            overdispersion = overdispersion, 
                                                            gene_specific = gene_specific, 
                                                            ref_normalization_term = ref_normalization_term, 
                                                            obs_normalization_term = obs_normalization_term,
                                                            states = states)
    print(emm_prob_log)
    
    ## filter nan emm_prob_log and update the Xdata
    idx_1 = np.where(emm_prob_log == emm_prob_log)
    gene_idx = np.unique(idx_1[1])
    emm_prob_log = emm_prob_log[:,gene_idx,:]

    update_Xdata = Xdata[:,gene_idx]
    obs_Xdata = update_Xdata[~ref_flag,:]
    if len(gene_idx) >= 1:
        print("filter nan emm_prob")

    ## run HMM based on emm_prob_log 
    ## todo -XC_HMM
    brk_item = obs_Xdata.var[brk].drop_duplicates(keep="first")
    print(brk_item)
    res_dict = {}
    res_log_dict = {}
    for i in range(emm_prob_log.shape[0]):
        # tmp_celltype = cell_type_order[i]
        tmp_cell_id = obs_Xdata.obs.index[i]
        tmp_emm_prob_log = emm_prob_log[i]
        cnt = 0
        for brk_ in brk_item:
            tmp_region_flag = obs_Xdata.var[brk] == brk_
            tmp_res, emm_prob_log_ = XC_HMM(tmp_emm_prob_log[tmp_region_flag,:])
            # print(tmp_res)
            if cnt == 0:
                res_ = tmp_res[0]
                res_log_ = tmp_res[1]
            else:
                res_ = np.vstack((res_, tmp_res[0]))
                res_log_ = np.vstack((res_log_, tmp_res[1]))
            cnt+=1
        
        ## dic record
        res_dict[tmp_cell_id] = res_
        res_log_dict[tmp_cell_id] = res_log_

    #     ## numpy array record for loglik cal
    #     if i == 0:
    #         res = res_[np.newaxis,:]
    #         res_log = res_log_[np.newaxis, :]
    #     else:
    #         res = np.vstack((res, res_[np.newaxis,:]))
    #         res_log = np.vstack((res_log, res_log_[np.newaxis, :]))

        
    # ## calculate the likelihood as ref of overdispersion grid search
    # log_lik = cal_log_lik(emm_prob_log, res_log)

    # ## output anndata for visualization 
    # ### convert the prob dict to ad
    # cnt = 0
    # celltype_lst = [] # in a order
    # for celltype, prob_ in res_dict.items():
    #     celltype_lst.append(celltype)
    #     if cnt==0:
    #         res_prob = prob_.T
    #         res_cnv = np.argmax(prob_.T, axis=0)
    #     else:
    #         res_prob = np.vstack([res_prob, prob_.T])
    #         res_cnv = np.vstack([res_cnv, np.argmax(prob_.T, axis=0)])
    #     cnt+=1
    
    # K = emm_prob.shape[1] # states num ## todo
    # prob_celltype = list(itertools.chain.from_iterable(itertools.repeat(i, K)
    #                                                 for i in celltype_lst))
    # obs_df1 = pd.DataFrame(index = prob_celltype)
    # res_prob_ad =  ad.AnnData(res_prob, var=v_ad.var, obs = obs_df1)
    # res_prob_ad.obs_names_make_unique()
    
    # obs_df2 = pd.DataFrame(index = celltype_lst)
    # res_cnv_ad = ad.AnnData(res_cnv, var=v_ad.var, obs = obs_df2)

    # return res_prob_ad, res_cnv_ad, res_dict, v_ads
    
    return emm_prob_log, res_dict, res_log_dict, update_Xdata
 
    
    # brk_item = obs_Xdata.var[brk].drop_duplicates(keep="first")
 
    # groups = obs_Xdata.obs.groupby(cell_anno_key).indices
    # res_dict = {}
    # for tmp_celltype, inds in groups.items():
    #     print(tmp_celltype)
    #     # calculate emm_prob across whole genome (selected chrs)
    #     ref_all_bulk = ref_Xdata.X.A.sum(axis=0)
    #     # obs_all_bulk = obs_Xdata[inds,:].X.A.sum(axis=0)
        
    #     ## calculate the emm_prob across the genome (not the chr specific region)
    #     ### calculate for each cell-- to do---save log value!!! 20210818
    #     cnt=0
    #     for tmp_idx in inds:
    #         if cnt == 0:
    #             tmp_emm_prob = generate_nb_emmission_prob(ref_all_bulk, obs_Xdata[tmp_idx,:].X.A[0], states=states)
    #             # cell_emm_prob = tmp_emm_prob
    #         else:
    #             tmp_emm_prob_ = generate_nb_emmission_prob(ref_all_bulk, obs_Xdata[tmp_idx,:].X.A[0], states=states)
    #             tmp_emm_prob = tmp_emm_prob * tmp_emm_prob_
    #             # cell_emm_prob = np.vstack((cell_emm_prob, tmp_emm_prob_))
    #         cnt+=1
        
    #     cnt = 0
    #     for brk_ in brk_item:
    #         tmp_region_flag = obs_Xdata.var[brk] == brk_
    #         tmp_res, emm_prob = XC_HMM(tmp_emm_prob[tmp_region_flag,:])
    #         if cnt == 0:
    #             res_ = tmp_res[0]
    #         else:
    #             res_ = np.vstack((res_, tmp_res[0]))
    #         cnt+=1
    #     res_dict[tmp_celltype] = res_
    # # return res_dict