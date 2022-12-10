"""Base funcs for data visualization of CNVs in XClone"""
# Author: Rongting Huang
# Date: 2022-05-04
# update: 2022-05-04

## Part IV: Results Visualization
import numpy as np
import pandas as pd
import anndata as ad

def convert_res_to_ann(Xdata, 
                       res_prob_layer, 
                       weights = False, 
                       states_weight = np.array([1,2,3])):
    """
    Function: 
    Convert cell based res to anndata for visualization.
    Check nan value.
    todo: change legend for copy states

    Parameters:
    ----------
    res: np.array. same order with Xdata.X
    Xdata: anndata. Provide the annotation.

    Return:
    ------
    
    Example:
    -------
    """
    res_prob= Xdata.layers[res_prob_layer].copy()
    if weights:
        if states_weight is None:
            raise ValueError("[XClone-convert_res_ann] Pls check states_weight!")
        res_cnv_weights = (res_prob * states_weight).sum(axis=-1)
        res_cnv_weights_ad = Xdata.copy()
        res_cnv_weights_ad.X = res_cnv_weights.copy()

        return res_cnv_weights_ad
    else:
        res_cnv = np.argmax(res_prob, axis=-1)
        res_cnv_ad = Xdata.copy()
        res_cnv_ad.X = res_cnv.copy()

        return res_cnv_ad 

def convert_res_visual(res_dict, ref_obs_ad, states_weight = np.array([1,2,3]), states_K=3, obs_names_make_unique=False):
    """
    Convert cell based res to anndata for visualization.
    Check nan value.

    From dict to anndata.
    """

    obs_Xdata = ref_obs_ad[1:,:]
    ## output anndata for visualization 
    ### convert the prob dict to ad
    cnt = 0
    cell_lst = [] # in a order
    for cell_, prob_ in res_dict.items():
        cell_lst.append(cell_)
        if cnt==0:
            res_prob = prob_.T
            res_cnv = np.argmax(prob_.T, axis=0)
            res_cnv_weights = (prob_ * states_weight).sum(axis=1)
        else:
            res_prob = np.vstack([res_prob, prob_.T])
            res_cnv = np.vstack([res_cnv, np.argmax(prob_.T, axis=0)])
            res_cnv_weights = np.vstack([res_cnv_weights, (prob_ * states_weight).sum(axis=1)]) 
        cnt+=1
    
    # K = emm_prob_log.shape[2] # states num
    K = states_K
    prob_cell = list(itertools.chain.from_iterable(itertools.repeat(i, K)
                                                    for i in cell_lst))
    obs_df1 = pd.DataFrame(index = prob_cell)
    obs_df1["index"] = obs_df1.index
    obs_df1 = pd.merge(obs_df1, obs_Xdata.obs, how = "left", left_on = "index", right_index=True)
    # obs_df1 = pd.merge(obs_df1, obs_Xdata.obs, left_index=True, right_index=True)
    res_prob_ad =  ad.AnnData(res_prob, var=obs_Xdata.var, obs = obs_df1)
    if obs_names_make_unique == True:
        res_prob_ad.obs_names_make_unique()
    
    obs_df2 = pd.DataFrame(index = cell_lst)
    obs_df2 = pd.merge(obs_df2, ref_obs_ad.obs, left_index=True, right_index=True)
    res_cnv_ad = ad.AnnData(res_cnv, var=obs_Xdata.var, obs = obs_df2)

    res_cnv_weights_ad = ad.AnnData(res_cnv_weights, var=obs_Xdata.var, obs = obs_df2)
    res_cnv_weights_ad1 = ad.AnnData(res_cnv_weights-1, var=obs_Xdata.var, obs = obs_df2)
    
    ## check nan value
    nan_count = np.isnan(res_prob_ad.X).sum() 
    prob_check_nan = nan_count!=0
    if prob_check_nan:
        print("there are %d nan value in the prob mtx" %(nan_count))

    return res_prob_ad, res_cnv_ad, res_cnv_weights_ad, res_cnv_weights_ad1


def get_emm_prob_ad(emm_prob_log, cell_type_order, var_df, obs_names_make_unique=True):
    """
    output emm_prob_log for visualization

    cell_type_order: for constructing the anndata's obs
    var_df: for constructing the anndata's var
    --------
    Example:
    emm_prob_normalization_ad, emm_prob_log_ad = xclone.model.get_emm_prob_ad(emm_prob_log, rr_ad.obs["cell_type"], rr_ad.var)
    """
    K = emm_prob_log.shape[2] # states num
    prob_celltype = list(itertools.chain.from_iterable(itertools.repeat(i, K)
                                                    for i in cell_type_order))
    obs_df1 = pd.DataFrame(index = prob_celltype)
    
    ## log scale
    emm_log_prob_mtx = transfer_emmprob_to_ad(emm_prob_log)

    ## normalization in original scale
    emm_log_normalization = normalization_emm_log(emm_prob_log)
    emm_normalization = np.exp(emm_log_normalization)
    
    emm_prob_normalization_mtx = transfer_emmprob_to_ad(emm_normalization)

    ## output ad
    emm_prob_normalization_ad =  ad.AnnData(emm_prob_normalization_mtx, var=var_df, obs = obs_df1)
    if obs_names_make_unique == True:
        emm_prob_normalization_ad.obs_names_make_unique()
    emm_prob_normalization_ad.obs["celltype"] = emm_prob_normalization_ad.obs.index
    
    emm_prob_log_ad =  ad.AnnData(emm_log_prob_mtx, var=var_df, obs = obs_df1)
    if obs_names_make_unique == True:
        emm_prob_log_ad.obs_names_make_unique()
    emm_prob_log_ad.obs["celltype"] = emm_prob_log_ad.obs.index

    # emm_prob_ad.layers['emm_prob_normalization'] = emm_prob_normalization_mtx
    # emm_prob_ad.layers['emm_logprob'] = emm_log_prob_mtx
    return emm_prob_normalization_ad, emm_prob_log_ad


def transfer_emmprob_to_ad(emm_prob_log):
    """
    like (6,8427,3) shape to (18, 8427) shape
    """
    for i in range(emm_prob_log.shape[0]):
        if i == 0:
            emm_log_prob_mtx = emm_prob_log[i].T
        else:
            emm_log_prob_mtx = np.vstack((emm_log_prob_mtx, emm_prob_log[i].T))
    return emm_log_prob_mtx

def transfer_cnvprob_to_ad(cnv_prob, var_, celltype, states = np.array([1,2,3])):
    """
    like (833,3205,3) shape to (3205, 833) shape
    emm_prob_ad, 
    emm_prob_ad1 just for visualization 0 1 2
    
    """
    emm_prob = cnv_prob.transpose((1,0,2))
    res = (emm_prob *states).sum(axis=2)

    obs_df = pd.DataFrame(index = celltype)
    
    emm_prob_ad =  ad.AnnData(res, var=var_, obs = obs_df)
    emm_prob_ad.obs["celltype"] = emm_prob_ad.obs.index

    emm_prob_ad1 = ad.AnnData(res-1, var=emm_prob_ad.var, obs = emm_prob_ad.obs)

    return emm_prob_ad, emm_prob_ad1

def normalization_emm_log(emm_prob_log):
    """
    """
    emm_log_normalization = emm_prob_log - logsumexp(emm_prob_log, axis=2, keepdims=True)
    return emm_log_normalization

## Part V: Results Analysis
def transform_dic_to_array(emm_prob_log_dic, res_log_dict):
    """
    add newaxis
    """

    cnt = 0
    cell_type_order = []
    for cell_type in emm_prob_log_dic.keys():
        cell_type_order.append(cell_type)
        if cnt ==0:
            res_log = res_log_dict[cell_type][np.newaxis,:]
            emm_prob_log = emm_prob_log_dic[cell_type][np.newaxis,:]
        else:
            res_log = np.vstack((res_log, res_log_dict[cell_type][np.newaxis,:]))
            emm_prob_log = np.vstack((emm_prob_log, emm_prob_log_dic[cell_type][np.newaxis,:]))
        cnt+=1
    print(res_log.shape)
    print(emm_prob_log.shape)

    return res_log, emm_prob_log, cell_type_order

## Part VI:

def cal_raw_logratio(ref_all_bulk, obs_all_bulk):
    """
    celltype based bulk
    ## add small value 0.1 for visualization
    """
    ref_normalization_term = ref_all_bulk.sum()
    _GEX_ref = ref_all_bulk/ref_normalization_term
    obs_normalization_term = obs_all_bulk.sum(axis=1)
    raw_ratio = (obs_all_bulk + 0.1) / (_GEX_ref * obs_normalization_term[:,np.newaxis])
    raw_logratio = np.log2(raw_ratio)

    return raw_logratio


