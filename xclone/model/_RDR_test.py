"""Base functions for XClone RDR processing
"""

# Author: Rongting Huang
# Date: 2021-12-20
# update: 2021-12-20


import numpy as np
import pandas as pd
import anndata as ad
import scipy as sp
from scipy import stats


import datetime


# import statsmodels
# import statsmodels.api as sm
import statsmodels as sm



def fit_lib_ratio(Xdata, select_chr_index = None, verbose=True, model_results_path=None):
    """
    Function:
    based on select normal chrs ans GLM model to get the libsize
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

    start_time_ = datetime.datetime.now()
    ## default sort index for using
    ## extend chr index when 23, 24--> X, Y
    sorted_chr = ['1', '2','3', '4', '5', '6', '7', '8', '9','10', '11', '12', '13', '14', '15', '16', '17', '18', '19','20', '21', '22', 'X', 'Y']


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
    obs_num = Xdata.shape[0]-1
    for i in range(obs_num):
        cell_barcodes = Xdata.obs.index[i+1]
        select_chr = select_chr_index[i]
        select_chr_lst = []
        for k in range(len(select_chr)):
            select_chr_lst.append(sorted_chr[int(select_chr[k])])
            # select_chr_lst.append(str(int(select_chr[k]+1))) # error
        print(select_chr_lst)
        chr_flag = Xdata.var["chr"].isin(select_chr_lst)
        if verbose:
            print(i, chr_flag.sum())
        obs_y = Xdata.X[i+1, chr_flag]
        feature_x = np.ones(len(obs_y))
        exposure_x = Xdata.X[0, chr_flag]
        try:
            NB_glm = sm.discrete.discrete_model.NegativeBinomial(obs_y, feature_x, loglike_method = 'nb2', offset = None, exposure = exposure_x, missing = 'none', check_rank = True)
            NB_results = NB_glm.fit()
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
    
    # obs_Xdata = Xdata[1:,:].copy()
    # obs_Xdata.obs["library_ratio"] = sample_libsize
    # obs_Xdata.obs["library_alpha"] = sample_alpha
    # obs_Xdata.obs["sample_chr_total"] = sample_chr_total
    # obs_Xdata.obs["ref_chr_total"] = ref_chr_total
    # obs_Xdata.obs["sample_chr_total_normalization"] = sample_chr_total_normalization
    # obs_Xdata.uns["fit_lib_ratio_model"] = model_results
    # obs_Xdata.uns["fit_lib_removed_cells"] = removed_cells_lst

    ref_obs_Xdata = Xdata.copy()
    sample_libsize.insert(0, 1)
    ref_obs_Xdata.obs["library_ratio"] = sample_libsize
    sample_alpha.insert(0, 0)
    ref_obs_Xdata.obs["library_alpha"] = sample_alpha
    
    ref_total = Xdata.X[0, :].sum()
    sample_chr_total.insert(0, ref_total)
    ref_obs_Xdata.obs["sample_chr_total"] = sample_chr_total
    ref_chr_total.insert(0, ref_total)
    ref_obs_Xdata.obs["ref_chr_total"] = ref_chr_total
    sample_chr_total_normalization.insert(0, 1)
    ref_obs_Xdata.obs["sample_chr_total_normalization"] = sample_chr_total_normalization
    ref_obs_Xdata.uns["fit_lib_ratio_model"] = model_results
    ref_obs_Xdata.uns["fit_lib_removed_cells"] = removed_cells_lst

    # time stamp
    end_time_= datetime.datetime.now()
    time_used = end_time_ - start_time_
    print("Time used", time_used.seconds, "seconds")
    
    return ref_obs_Xdata
    # return obs_Xdata
    # return obs_data, sample_libsize, sample_alpha, sample_chr_total, ref_chr_total, sample_chr_total_normalization, model_results
