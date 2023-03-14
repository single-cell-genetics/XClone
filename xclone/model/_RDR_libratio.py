"""Base functions for XClone RDR processing
strategies for libratio estimation
"""

# Author: Rongting Huang
# Date: 2022-01-20
# update: 2022-03-15

## RDR
## Part-II libsize ratio estimation part in RDR module [source _RDR]

# strategy notes(pipeline here):
# 1) preprocess the Xdata here for following analysis.
# 2) select 4 normal chrs based on the celltype bulk Xdata.
# 3) calculate libratio using GLM for celltype/cell(use the extended selected chrs).
# 4) get emm_prob from RDR VB model and put back to HMM framework.

import numpy as np
import pandas as pd
from scipy import stats

import datetime
import multiprocessing

from ._RDR_base import list_update

# from .._logging import get_logger

## 1) preprocessing for libratio fitting
def select_normal_CHR1(Xdata, select_chr_num = 4):
    """
    deprecated 
    Function:
    chr based: select 4 normal chrs(based on raw ratio) to 
               get the lib size for rescaling.

    Parameters:
    -----------
    Xdata: celltype_ad from Xdata_preprocess.
    use the layer["raw_ratio"]
    select_chr_num: selected chrs num for following analysis, default: 4.
    
    Return:
    ------
    Xdata: with obsm including t_stat, p_val and select_chr_index

    Example:
    --------
    # rr_ad_celltype = xclone.model.select_normal_CHR(rr_ad_celltype)
    celltype_ad = xclone.model.select_normal_CHR(celltype_ad)
    """
    # sub_logger = get_logger("RDR.libratio.select.normal")
    # sub_logger.info("select_normal_CHR")
    
    chr_lst = Xdata.var["chr"].drop_duplicates(keep="first")
    chr_len = chr_lst.shape[0]

    for i in range(chr_len):
        chr_flag = Xdata.var["chr"] == chr_lst[i]
        tmp_chr = Xdata.layers["raw_ratio"][:, chr_flag]
        pop_mean = np.zeros((1, tmp_chr.shape[0]))
        ## for raw ratio is in log format.
        # t-test
        tmp_t, tmp_p = stats.ttest_1samp(tmp_chr, pop_mean, axis=1, nan_policy='propagate', alternative='two-sided')
        
        if i == 0:
            t_stat = tmp_t
            p_val = tmp_p
        else:
            t_stat = np.vstack((t_stat, tmp_t))
            p_val = np.vstack((p_val, tmp_p))
    
    t_stat = t_stat.T
    p_val = p_val.T
    
    ## select the near normal chrs based on the t-test p_value
    order_index = np.argsort(p_val, axis=1)
    order_index = order_index[:, ::-1] # debugs notes: select higher p_val chr
    select_chr_index = order_index[:,0:select_chr_num]

    Xdata.obsm["chr_index"] = select_chr_index
    Xdata.obsm["t_stat"] = t_stat
    Xdata.obsm["p_val"] = p_val
    print(Xdata.obs)
    print(select_chr_index)
    return Xdata


import numpy as np
from sklearn.mixture import GaussianMixture

def select_normal_CHR(Xdata, select_chr_num = 4):
    """
    problem: maybe not all chr have 2 GMM comp, e.g., chr with copy loss.
    """
    chr_lst = Xdata.var["chr"].drop_duplicates(keep="first")
    chr_len = chr_lst.shape[0]

    for i in range(chr_len):
        chr_flag = Xdata.var["chr"] == chr_lst[i]
        tmp_chr = Xdata.layers["raw_ratio"][:, chr_flag]

        gm = GaussianMixture(n_components=2, random_state=0).fit(tmp_chr.T)
        # gm.means_
        # gm.covariances_
        type_lst = [*range(0,tmp_chr.shape[0])]
        means_ = gm.means_[(np.argsort(gm.means_, axis = 0)[1,:], type_lst)]
        vars_ = gm.covariances_[(np.argsort(gm.means_, axis = 0)[1,:], type_lst, type_lst)]

        import scipy.stats
        #find p-value(z-score)
        tmp_pval = scipy.stats.norm.sf(abs((0-means_)/np.sqrt(vars_)))
        if i == 0:
            p_val = tmp_pval
        else:
            p_val = np.vstack((p_val, tmp_pval))
    p_val = p_val.T
    
    ## select the near normal chrs based on the t-test p_value
    order_index = np.argsort(p_val, axis=1)
    order_index = order_index[:, ::-1]
    select_chr_index = order_index[:,0:select_chr_num]

    Xdata.obsm["chr_index"] = select_chr_index
    Xdata.obsm["p_val"] = p_val
    print(Xdata.obs)
    print(p_val)
    print(select_chr_index)
    return Xdata


def equal_content(a, b):
    """
    FUNC used in check_cell_celltype()

    """
    if a == b:
        return 1
    else:
        return 0
    
def check_cell_celltype(rr_ad_cell_anno, rr_ad_celltype_anno):
    """
    """
    ## build the dataframe 
    
#     df = pd.DataFrame(rr_ad_cell_anno)
#     df["b"] = rr_ad_celltype_anno.values
#     df.columns = ['a','b']
    
    data = {"a": rr_ad_cell_anno.values,
            "b": rr_ad_celltype_anno.values}
    df = pd.DataFrame(data)
    
    ## compare the order of the celltype Xdata and cellbased Xdata
    df["res"] = df.apply(lambda x: equal_content(x['a'], x['b']), axis=1)
    
    if df["res"].sum() == len(df["res"]):
        return True, df
    else:
        return False, df

def extend_chr_index(rr_ad_cell, rr_ad_celltype, select_chr_index):
    """
    First version:(might be deprecated)
    for ordered rr_ad_cell
    Function:
    ---------
    extend the select_chr_index from celltype-based to cell based index
    and then can used the new extended index to get libratio for cell based Xdata.

    Parameters:
    ----------

    Return:
    ------

    """
    rr_ad_cell_anno = rr_ad_cell.obs["cell_type"].drop_duplicates().reset_index(drop=True)
    rr_ad_celltype_anno = rr_ad_celltype.obs["cell_type"]

    res, cell_celltype_anno_df = check_cell_celltype(rr_ad_cell_anno, rr_ad_celltype_anno)
    if res == False:
        raise ValueError("cell based data should be reordered by the celltype based Xdata! Pls check!")
    

    cell_counts = rr_ad_cell.obs.groupby("cell_type").count().values[:,0]


    for idx_ in range(len(cell_counts)):
        if idx_ == 0:
            extend_chr_index = np.tile(select_chr_index[idx_], (cell_counts[idx_],1))
        else:
            extend_chr_index = np.vstack([extend_chr_index, np.tile(select_chr_index[idx_], (cell_counts[idx_],1))])
    return extend_chr_index


def extend_chr_index_to_singlecell(Xdata, celltype_ad, cell_anno_key = "cell_type"):
    """
    Function: extend the select_chr_index from celltype-based to cell based index
    and then can used the new extended index to get libratio for cell based Xdata.

    Here try to find the same celltype in Xdata and then map the chr_index in 
    `celltype_ad` to `extend_chr_index` in `Xdata`.

    Notes: if the celltype in Xdata is not in `celltype_ad`, the value in 
    `extend_chr_index` will be NAN.
    
    Parameters:
    -----------
    Xdata: anndata.
    celltype_ad: anndata.
    cell_anno_key: char.
    
    Return:
    -------
    Xdata: anndata.

    Example:
    -------
    For celltype:
    ref_obs_bulk_ad = xclone.model.extend_chr_index_to_singlecell(ref_obs_bulk_ad, update_rr_ad_celltype, cell_anno_key = "cell_type")
    For cell:
    ref_obs_ad1 = xclone.model.extend_chr_index_to_singlecell(ref_obs_ad1, update_rr_ad_celltype, cell_anno_key = "cell_type")
    """
    cell_anno = Xdata.obs[cell_anno_key]

    rr_ad_celltype_anno = celltype_ad.obs[cell_anno_key]
    select_chr_index = celltype_ad.obsm["chr_index"]

    extend_chr_index = np.empty((Xdata.shape[0], select_chr_index.shape[1]))
    extend_chr_index[:] = np.NaN

    for cell_type_, chr_index_ in zip(rr_ad_celltype_anno, select_chr_index):
        cell_type_flag_ = cell_anno == cell_type_
        extend_chr_index[cell_type_flag_, :] = chr_index_
    
    Xdata.obsm["select_chr_index"] = extend_chr_index
    return Xdata


## 2) libratio estimation

# import statsmodels
import statsmodels.api as sm
import statsmodels as sm

def get_libsize(X_data, select_chr_index, verbose=True, 
    NB_kwargs={'disp': True, 'skip_hessian': True}):
    """
    deprecated. 20230105 Notes
    ---> fit_lib_ratio FUNC
    based on select normal chrs ans GLM model to get the libsize
    here use NB glm
    https://www.statsmodels.org/stable/generated/statsmodels.discrete.discrete_model.NegativeBinomial.html

    X_data: the first row is ref; and the follwing rows are other obs;
    (for both celltype besed bulk count and cell based count)

    Example:
    get_libsize(ref_obs_bulk_ad, select_chr_index) # for celltype

    get_libsize(ref_obs_ad, extend_select_chr_index) # for cell
    """
    ## init the output
    sample_libsize = []
    sample_alpha = []

    sample_chr_total = []
    ref_chr_total = []
    
    model_results = []
    ## use the GLM to get the libsize for each sample
    obs_num = X_data.X.shape[0]-1
    for i in range(obs_num):
        select_chr = select_chr_index[i]
        select_chr_lst = []
        for k in range(len(select_chr)):
            select_chr_lst.append(str(select_chr[k]+1))
        
        chr_flag = X_data.var["chr"].isin(select_chr_lst)
        if verbose:
            print(i, chr_flag.sum())
        obs_y = X_data.X[i+1, chr_flag]
        feature_x = np.ones(len(obs_y))
        exposure_x = X_data.X[0, chr_flag]
        NB_glm = sm.discrete.discrete_model.NegativeBinomial(obs_y, feature_x, 
                                                            loglike_method = 'nb2', offset = None, exposure = exposure_x, missing = 'none', check_rank = True)
        NB_results = NB_glm.fit(**NB_kwargs)
        model_results.append(NB_results)
        libsize_ = np.exp(NB_results.params[0])
        alpha_ = NB_results.params[1]
        sample_libsize.append(libsize_)
        sample_alpha.append(alpha_)
        
        sample_chr_total.append(obs_y.sum())

        # ref total when select different chrs
        ref_chr_total.append(exposure_x.sum())
        sample_chr_total_normalization = [a/b for a, b in zip(sample_chr_total, ref_chr_total)]
        

    return sample_libsize, sample_alpha, sample_chr_total, ref_chr_total, sample_chr_total_normalization, model_results


def cell_filter_by_celltype(Xdata, anno_key = "select_chr_index", verbose = True):
    """
    Function:
    default: filter the cells with no cell annotation. get valiad dataset.
    Here based on the select_chr_index.
    
    Parameters:
    ----------
    Xdata: anndata.

    Return:
    ------
    Xdata: anndata.

    Example:
    --------

    ref_obs_ad1 = cell_filter(ref_obs_ad1)
    """
    flag_ = np.isnan(Xdata.obsm[anno_key][1:,:]).any(axis=1)
    if verbose:
        print("remove cells num:", flag_.sum())
    flag_ = ~flag_
    flag_ = np.insert(flag_, 0, True)
    Xdata = Xdata[flag_, :]
    return Xdata

##  https://www.statsmodels.org/stable/generated/statsmodels.genmod.generalized_linear_model.GLM.html#statsmodels.genmod.generalized_linear_model.GLM
## for celltype/cell maybe use different glm model
    
def fit_lib_ratio(Xdata, select_chr_index = None, verbose=False, avg_key = "ref_avg",
    model_results_path=None, NB_kwargs={'disp': True, 'skip_hessian': True}):
    """
    Function:
    version: 0.0.0
    based on selected normal chrs ans GLM model to get the libsize
    here use NB glm
    https://www.statsmodels.org/stable/generated/statsmodels.discrete.discrete_model.NegativeBinomial.html
    
    Parameter:
    ---------
    Xdata: the first row is ref; and the follwing rows are other obs;
    (for both celltype besed bulk count and cell based count)
    
    Return:
    -------
    Example:
    --------
    >>> fit_lib_ratio(ref_obs_bulk_ad, select_chr_index = ref_obs_bulk_ad.obsm["select_chr_index"][1:,:], verbose=True) # for celltype

    >>> fit_lib_ratio(ref_obs_ad, select_chr_index = ref_obs_ad.obsm["select_chr_index"][1:,:], verbose=True) # for cell
    """
    import statsmodels as sm
    import scipy as sp

    start_time_ = datetime.datetime.now()
    ## default sort index for using
    ## extend chr index when 23, 24--> X, Y
    sorted_chr = ['1', '2','3', '4', '5', '6', '7', '8', '9','10', '11', '12', '13', '14', '15', '16', '17', '18', '19','20', '21', '22', 'X', 'Y']

    if sp.sparse.issparse(Xdata.X):
        Xmtx = Xdata.X.A
    else:
        Xmtx = Xdata.X

    if select_chr_index is None:
        select_chr_index = Xdata.obsm["select_chr_index"]
    
    ## init the output
    sample_libsize = []
    sample_alpha = []

    sample_chr_total = []
    ref_chr_total = []
    
    model_results = []
    removed_cells_lst = []

    ## use the GLM to get the libsize for each sample
    obs_num = Xdata.shape[0]

    for i in range(obs_num):
        cell_barcodes = Xdata.obs.index[i]
        select_chr = select_chr_index[i]
        select_chr_lst = []
        for k in range(len(select_chr)):
            select_chr_lst.append(sorted_chr[int(select_chr[k])])
        if verbose:
            print(select_chr_lst)
        chr_flag = Xdata.var["chr"].isin(select_chr_lst)
        if verbose:
            print(i, chr_flag.sum())

        obs_y = Xmtx[i, chr_flag]
        feature_x = np.ones(len(obs_y))

        exposure_x = Xdata.var[avg_key][chr_flag]
        try:
            NB_glm = sm.discrete.discrete_model.NegativeBinomial(obs_y, feature_x, loglike_method = 'nb2', offset = None, exposure = exposure_x, missing = 'none', check_rank = True)
            NB_results = NB_glm.fit(**NB_kwargs)
        except:
            removed_cells_lst.append(cell_barcodes)
            model_results.append("Error in this cell" +cell_barcodes )

            libsize_ = np.NaN
            alpha_ = np.NaN
            
            sample_libsize.append(libsize_)
            sample_alpha.append(alpha_)
            sample_chr_total.append(obs_y.sum())
            ref_chr_total.append(exposure_x.sum())
            sample_chr_total_normalization = [a/b for a, b in zip(sample_chr_total, ref_chr_total)]

        else:
            model_results.append(NB_results)
            libsize_ = np.exp(NB_results.params[0])
            alpha_ = NB_results.params[1]
            sample_libsize.append(libsize_)
            sample_alpha.append(alpha_)
        
            sample_chr_total.append(obs_y.sum())

            # ref total when select different chrs
            ref_chr_total.append(exposure_x.sum())
            sample_chr_total_normalization = [a/b for a, b in zip(sample_chr_total, ref_chr_total)]
    
    if model_results_path is not None:
        np.save(model_results_path, model_results)

    ## update
    Xdata.obs["library_ratio"] = sample_libsize
    Xdata.obs["library_alpha"] = sample_alpha

    Xdata.obs["sample_chr_total"] = sample_chr_total
    Xdata.obs["ref_chr_total"] = ref_chr_total
    Xdata.obs["sample_chr_total_normalization"] = sample_chr_total_normalization

    # Xdata.uns["fit_lib_ratio_model"] = model_results
    # Xdata.uns["fit_lib_removed_cells"] = removed_cells_lst

    # time stamp
    end_time_= datetime.datetime.now()
    time_used = end_time_ - start_time_
    print("[XClone RDR library ratio fitting] Time used:", time_used.seconds, "seconds")

    return Xdata

def show_progress(RV=None):
    print(RV)
    return RV

def libratio_GLM_base(obs_y, feature_x, exposure_x, cell_barcodes,
    NB_kwargs={'disp': True, 'skip_hessian': True}):
    """
    Function:
    """
    try:
        NB_glm = sm.discrete.discrete_model.NegativeBinomial(obs_y, feature_x, 
                                                            loglike_method = 'nb2', 
                                                            offset = None, 
                                                            exposure = exposure_x, 
                                                            missing = 'none', 
                                                            check_rank = True)
        NB_results = NB_glm.fit(**NB_kwargs)
    except:
        cell_remove = True
        # removed_cells_lst.append(cell_barcodes)
        model_results_= "Error in this cell" + cell_barcodes
        # model_results.append("Error in this cell" + cell_barcodes)

        libsize_ = np.NaN
        alpha_ = np.NaN
        
        # sample_libsize.append(libsize_)
        # sample_alpha.append(alpha_)
        # sample_chr_total.append(obs_y.sum())
        # ref_chr_total.append(exposure_x.sum())
        # sample_chr_total_normalization = [a/b for a, b in zip(sample_chr_total, ref_chr_total)]


    else:
        cell_remove = False
        # model_results.append(NB_results)
        model_results_ = NB_results
        libsize_ = np.exp(NB_results.params[0])
        alpha_ = NB_results.params[1]
        # sample_libsize.append(libsize_)
        # sample_alpha.append(alpha_)
    
    # sample_chr_total = obs_y.sum()
    # # ref total when select different chrs
    # ref_chr_total= exposure_x.sum()
    # sample_chr_total_normalization = [a/b for a, b in zip(sample_chr_total, ref_chr_total)]

    # return model_results_, libsize_, alpha_, cell_remove, sample_chr_total, ref_chr_total, sample_chr_total_normalization
    return model_results_, libsize_, alpha_, cell_remove


## run GLM with single or multiple processes
def fit_lib_ratio_accelerate(Xdata, select_chr_index = None, avg_key = "ref_avg", verbose=True, 
    model_results_path=None, nproc = 1, NB_kwargs={'disp': True}):
    """
    Function:
    version: 0.0.1
    update for GLM multi-processing[need fix the bug in nproc>1]
    based on selected normal chrs ans GLM model to get the libsize
    here use NB glm
    https://www.statsmodels.org/stable/generated/statsmodels.discrete.discrete_model.NegativeBinomial.html
    
    Parameter:
    ---------
    Xdata: the first row is ref; and the follwing rows are other obs;
    (for both celltype besed bulk count and cell based count)
    
    Return:
    -------
    Example:
    --------
    fit_lib_ratio(ref_obs_bulk_ad, select_chr_index = ref_obs_bulk_ad.obsm["select_chr_index"][1:,:], verbose=True) # for celltype

    fit_lib_ratio(ref_obs_ad, select_chr_index = ref_obs_ad.obsm["select_chr_index"][1:,:], verbose=True) # for cell
    """
    # import statsmodels.api as sm
    # import statsmodels as sm
    import scipy as sp

    start_time_ = datetime.datetime.now()
    ## default sort index for using
    ## extend chr index when 23, 24--> X, Y
    sorted_chr = ['1', '2','3', '4', '5', '6', '7', '8', '9','10', '11', '12', '13', '14', '15', '16', '17', '18', '19','20', '21', '22', 'X', 'Y']
    
    if sp.sparse.issparse(Xdata.X):
        Xmtx = Xdata.X.A
    else:
        Xmtx = Xdata.X

    if select_chr_index is None:
        select_chr_index = Xdata.obsm["select_chr_index"][1:,:]
    
    print(select_chr_index)
    ## init the output
    sample_libsize = []
    sample_alpha = []

    sample_chr_total = []
    ref_chr_total = []
    
    model_results = []
    removed_cells_lst = [] 

    ## use the GLM to get the libsize for each sample
    # obs_num = Xdata.shape[0]-1
    obs_num = Xdata.shape[0]

    if nproc == 1:
        for i in range(obs_num):
            # cell_barcodes = Xdata.obs.index[i+1]
            cell_barcodes = Xdata.obs.index[i]
            select_chr = select_chr_index[i]
            select_chr_lst = []
            for k in range(len(select_chr)):
                select_chr_lst.append(sorted_chr[int(select_chr[k])])
            print(select_chr_lst)
            chr_flag = Xdata.var["chr"].isin(select_chr_lst)
            if verbose:
                print(i, chr_flag.sum())
            # obs_y = Xdata.X[i+1, chr_flag]
            obs_y = Xmtx[i, chr_flag]
            feature_x = np.ones(len(obs_y))
            # exposure_x = Xdata.X[0, chr_flag]
            exposure_x = Xdata.var[avg_key][chr_flag]

            model_results_, libsize_, alpha_, cell_remove = libratio_GLM_base(
                obs_y, feature_x, exposure_x, cell_barcodes, NB_kwargs)

            model_results.append(model_results_)
            sample_libsize.append(libsize_)
            sample_alpha.append(alpha_)
            removed_cells_lst.append(cell_remove)
            
            # sample_chr_total.append(obs_y.sum())

    elif nproc > 1:
        print("[XClone] libratio fitting: multiprocessing for each cell")
        print("nproc:", nproc)

        pool = multiprocessing.Pool(processes=nproc)

        fit_result = []

        for i in range(obs_num):
            # cell_barcodes = Xdata.obs.index[i+1]
            cell_barcodes = Xdata.obs.index[i]
            select_chr = select_chr_index[i]
            select_chr_lst = []
            for k in range(len(select_chr)):
                select_chr_lst.append(sorted_chr[int(select_chr[k])])
            if verbose:
                print(select_chr_lst)
            chr_flag = Xdata.var["chr"].isin(select_chr_lst)
            if verbose:
                print(i, chr_flag.sum())
            
            # obs_y = Xdata.X[i+1, chr_flag]
            obs_y = Xmtx[i, chr_flag]
            feature_x = np.ones(len(obs_y))
            # exposure_x = Xdata.X[0, chr_flag]
            exposure_x = Xdata.var[avg_key][chr_flag]

            ## todo- test and add real call back # bug in nproc>1 # no reason # no efficiency improvement
            fit_result.append(pool.apply_async(libratio_GLM_base, 
                (obs_y, feature_x, exposure_x, cell_barcodes, NB_kwargs),
                callback=show_progress
                ))
        pool.close()
        pool.join()

        print(fit_result)

        ## processing result   todo here
        # result = [res.get() for res in result]
        for tmp_res in fit_result:
            model_results.append(tmp_res.get()[0])
            sample_libsize.append(tmp_res.get()[1])
            sample_alpha.append(tmp_res.get()[2])
            removed_cells_lst.append(tmp_res.get()[3])

    if model_results_path is not None:
        np.save(model_results_path, model_results)    

    Xdata.obs["library_ratio"] = sample_libsize
    Xdata.obs["library_alpha"] = sample_alpha

    # Xdata.uns["fit_lib_ratio_model"] = model_results
    # Xdata.uns["fit_lib_removed_cells"] = removed_cells_lst

    # time stamp
    end_time_= datetime.datetime.now()
    time_used = end_time_ - start_time_
    print("Time used", time_used.seconds, "seconds")

    return Xdata
    
def remove_cells(Xdata, mode = "NaN"):
    """
    Fucntion:
    removed cells with nan/inf library_ratio
    actually, when the mode is "GLM", also remove part of nan value.
    for those without results from GLM we assign nan value.    
    Parameters:
    ----------

    Return:
    -------

    Example:
    -------

    """
    # remove cells with nan results
    if mode == "GLM":
        FLAG_ = Xdata.obs.index.isin(Xdata.uns["fit_lib_removed_cells"])
        print("remove no GLM results cells num:", FLAG_.sum())
    if mode == "NaN":
        FLAG_ = np.isnan(Xdata.obs["library_ratio"])
        print("remove nan library_ratio cells num:", FLAG_.sum())
    if mode == "INF":
        FLAG_ = np.isinf(Xdata.obs["library_ratio"])
        print("remove inf library_ratio cells num:", FLAG_.sum())
    
    Xdata = Xdata[~FLAG_, :]
    ## just observations have fitted model
    # update_model_lst = list_update(Xdata.uns["fit_lib_ratio_model"].copy(), ~FLAG_[1:])

    # Xdata.uns["fit_lib_ratio_model"] = update_model_lst

    return Xdata


def total_count_ratio(ref_obs_bulk_ad):
    """
    Function:
    For getting the init libratio.
    whole chromosome here.
    maybe need update for select chrs.#TODO
    """
    ref = ref_obs_bulk_ad.X[0,:].sum()
    obs = ref_obs_bulk_ad.X[1:,:].sum(axis=1)
    
    return obs/ref

## 20210930
## RDR ELBO
## notes global optimization scipy.optimize.basinhopping

### a test version for getting libratio
### deprecated one
from scipy.optimize import minimize
# import numpy as np
                                                                                                                                                                                              
def libratio_elbo_fun(r_counts, ref_counts, alpha_, beta_):
    """
    version: return a func

    depprecated
    test version
    
    Example:
    ----
    libratio_elbo_fun(r_counts, ref_counts, alpha_, beta_)(libratio)
    
    """
    func = lambda x: -np.dot(r_counts, np.log(x)).sum() + (ref_counts * np.dot(alpha_/beta_, x).T).sum()
    return func

def libratio_elbo_func(x, r_counts, ref_counts, alpha_, beta_, verbose=True):
    """
    test version verbose=True for print.
    For multiple cells or for one cell
    here gene i  cell j--according to the derivation

    # r_counts, ref_counts, alpha_, beta_ =  args
    # term1_ = np.dot(r_counts, np.log(phi_)).sum()
    # term2_ = (ref_counts * np.dot(alpha_/beta_, phi_).T).sum()
    # return term1_ + term2_
    ----
    
    Example:
    ----
    
    libratio = np.ones((8,1))
    ref_counts = np.random.randint(10, size = (1,10))
    r_counts = np.random.randint(10, size = (10,8))

    alpha_ =  np.random.rand(10,8)
    beta_ = np.random.rand(10,8)
    ----
    libratio_elbo_func(libratio, r_counts, ref_counts, alpha_, beta_)
    """
    if verbose:
        print(x)
    if r_counts.ndim == 2:
        res = -np.dot(r_counts, np.log(x)).sum() + (ref_counts * np.dot(alpha_/beta_, x).T).sum()
    elif r_counts.ndim == 1:
        res = - (r_counts*np.log(x)).sum() + (ref_counts * (alpha_/beta_) * x).sum()
    else:
        print("[XClone warning]: check the input r_counts dimension!")
    return res



def optmize_libratio(libratio_init, r_counts, ref_counts, alpha_, beta_):
    """
    one format for using the scipy minimize
    """
#     res = minimize(libratio_elbo_fun(r_counts, ref_counts, alpha_, beta_), libratio_init,  method="SLSQP")
    res = minimize(libratio_elbo_fun(r_counts, ref_counts, alpha_, beta_), libratio_init,  method="L-BFGS-B")#BFGS  # tol=1
    return res

def optmize_libratio2(libratio_init, r_counts, ref_counts, alpha_, beta_, **kwargs):
    """
    Function:

    Params:
    ------
    **kwargs: maybe todo-here can be extended to adopt any params in scipy minimize optimization func, like the method:
    "L-BFGS-B", "SLSQP" and any other method.
    ------
    ##  minimize-local optimization function

    Example:
    -------
    """
    if r_counts.ndim == 2:
        cells_num = r_counts.shape[0]
        res = minimize(libratio_elbo_func, libratio_init, args = (r_counts, ref_counts, alpha_, beta_), method="L-BFGS-B", bounds = [(1e-6, None)]*cells_num)#BFGS  # tol=1 #method="SLSQP"
    elif r_counts.ndim == 1:
        cells_num = 1
        res = minimize(libratio_elbo_func, libratio_init, args = (r_counts, ref_counts, alpha_, beta_), method="L-BFGS-B", bounds = [(1e-6, None)]*cells_num)#BFGS  # tol=1
    else:
        print("[XClone warning]: check the input r_counts dimension!")
    return res

## whole genome can not get better ratio, so we try to get ratio from four best(normal) chrs
## because different cell type will get different 4 chrs, so the gene length are different (also the ref
## for different celltype, so we try to process different celltype, not in a structured matrix)

# def init_count_ratio(ref, obs):
#     """
    
#     """
#     ref = ref_obs_bulk_ad.X[0,:].sum()
#     obs = ref_obs_bulk_ad.X[1:,:].sum(axis=1)
    
#     return obs/ref

def get_libsize2(X_data, select_chr_index, verbose=False):
    """
    based on select normal chrs and scipy opt to get the libsize
    here use libratio_elbo_fun2

    X_data: the first row is ref; and the follwing rows are other obs;
    (for both celltype besed bulk count and cell based count)

    Example:
    get_libsize(ref_obs_bulk_ad, select_chr_index)
    """
    ## init the output
    sample_libsize = []
#     sample_alpha = []
    libratio_init = [] # record for comparing
    sample_chr_total = []
    ref_chr_total = []
    
    ## use the GLM to get the libsize for each sample
    obs_num = X_data.X.shape[0]-1
    for i in range(obs_num):
        select_chr = select_chr_index[i]
        select_chr_lst = []
        
        for k in range(len(select_chr)):
            select_chr_lst.append(str(select_chr[k]+1))
        
        chr_flag = X_data.var["chr"].isin(select_chr_lst)
        print(i, chr_flag.sum())
        r_counts = X_data.X[i+1, chr_flag].T
        ref_counts = X_data.X[0, chr_flag].T
        
        libratio_init_ = r_counts.sum()/ref_counts.sum() ## need update maybe for matrix,cell based (celltype block)
        libratio_init.append(libratio_init_)
        ## have transposed the counts
#         if r_counts.ndim == 2:
#             genes_numi = r_counts.shape[0]
#             cells_numj = r_counts.shape[1]
            
#         elif r_counts.ndim == 1:
#             print("r_counts ndim = 1")
#             genes_numi = r_counts.shape[0]
#             cells_numj = 1
#         else:
#             print("[XClone warning]: check the input r_counts dimension!")
        
        alpha_ =  np.ones(r_counts.shape)
        beta_ = np.ones(r_counts.shape)
        
        opt_res = optmize_libratio2(libratio_init_, r_counts, ref_counts, alpha_, beta_, verbose=verbose)
        if verbose:
            print("print opt_res:", opt_res)
        
        libsize_ = opt_res.x
        sample_libsize.append(libsize_[0])

        sample_chr_total.append(r_counts.sum())
        # ref total when select different chrs
        ref_chr_total.append(ref_counts.sum())
    sample_chr_total_normalization = [a/b for a, b in zip(sample_chr_total, ref_chr_total)]

    return sample_libsize, libratio_init, sample_chr_total, ref_chr_total, sample_chr_total_normalization

def check_libratio(Xdata, anno_key = "library_ratio"):
    print("[XClone RDR library ratio]: checking")
    print("max_value:", Xdata.obs[anno_key].max())
    print("min_value:", Xdata.obs[anno_key].min())

    print("qt_0.95_value:", np.quantile(Xdata.obs[anno_key], 0.95))
    print("qt_0.05_value:", np.quantile(Xdata.obs[anno_key], 0.05))
    


def libsize_clip(Xdata, anno_key = "library_ratio", min_threshold = 0, max_threshold = 5):
    """
    YH: one extreme value; better clipped, though not making to much difference
    """
    update_anno_key = anno_key + "_capped"
    
    Xdata.obs[update_anno_key] = np.clip(
    Xdata.obs[anno_key], min_threshold, max_threshold)
    print("[XClone RDR library ratio]: clipping")

    return Xdata

def libsize_select(Xdata, anno_key = "library_ratio", threshold = 10, verbose =True):
    """
    Select the rational estimated libsize/common total counts
    filter or keep it as the threshold value
    
    Parameters:
    ----------

    Return:
    -------

    Example:
    --------
    
    """
    FLAG_ = Xdata.obs[anno_key] > threshold
    if verbose:
        print("removed cells", Xdata.obs[anno_key][FLAG_])
    
    Xdata = Xdata[~FLAG_, :]
    
    ## just observations have fitted model  ~FLAG_[1:]
    # update_model_lst = list_update(Xdata.uns["fit_lib_ratio_model"].copy(), ~FLAG_[1:])
    # Xdata.uns["fit_lib_ratio_model"] = update_model_lst

    return Xdata