"""Base functions for XClone RDR processing
strategies for gene-specific dispersion estimation
"""

# Author: Rongting Huang
# Date: 2022-01-20
# update: 2022-03-15

## RDR

## Part-III gene-specific Dispersions estimation part in RDR module [source _RDR]

# strategy notes
# 1. view celltype information to find which celltype have the most cells
# 2. select the celltype with most cell nums
# 3. and estimate gene-specific dispersions based on the selected celltype(do not perform 
# this across wholedataset to avoid counting in the biological difference) 

import numpy as np
import scipy as sp
import pandas as pd
import datetime

from ._RDR_base import list_update

def view_celltype(Xdata, cell_anno_key = "cell_type"):
    """
    """
    celltype_count = Xdata.obs[cell_anno_key].value_counts()
    print(celltype_count)
    return celltype_count


def select_celltype(Xdata, cell_anno_key = "cell_type", select_celltype = "Paneth"):
    """
    """
    Flag_ = Xdata.obs[cell_anno_key] == select_celltype
    
    update_Xdata = Xdata[Flag_, :]
    
    return update_Xdata

def Estimate_Dispersions(Xdata, sample_libsize = None, verbose=False, 
    model_results_path=None, NB_kwargs={'disp': True, 'skip_hessian': True}):
    """
    Function:[deprecated]
    Estimate overdispersion based on mean value for each gene.

    todo
    1. need save var annotations 
    2. multiprocessing for each gene
    
    
    Parameters:
    -----------
    sample_libsize: here can input total counts/learned libsize ratio
    
    Return:
    -------

    Example:
    --------

    Estimate_Dispersions(Xdata, sample_libsize = None, verbose=True, 
    model_results_path = "/storage/yhhuang/users/rthuang/xclone/testing/test_gx109_inferCNV_estDispersion.npy")
    
    test_a = np.load("/storage/yhhuang/users/rthuang/xclone/testing/test_gx109_inferCNV_estDispersion.npy", allow_pickle =True)
    """

    if sp.sparse.issparse(Xdata.X):
        X_mtx = Xdata.X.A
    else:
        X_mtx = Xdata.X
    
    sample_libsize = X_mtx.sum(axis=1)  # cell total counts

    # check total count/libratio
    if (sample_libsize == 0).sum() > 0:
        raise ValueError("[XClone]Zero value in sample_libsize/total counts! Pls check!")
    
    ## init the output
    model_results = []
    gene_alpha = []
    gene_alpha_bse = []
    
    
    # total_counts = X_mtx.sum(axis=1)
    
    cell_num = X_mtx.shape[0]
    gene_num = X_mtx.shape[1]

    for i in range(gene_num):
        obs_y = X_mtx[:, i]
        feature_x = np.ones(cell_num)

        NB_glm = sm.discrete.discrete_model.NegativeBinomial(obs_y, feature_x, 
                                                            loglike_method = 'nb2', 
                                                            exposure = sample_libsize, 
                                                            offset = None, 
                                                            missing = 'none', 
                                                            check_rank = True)
        NB_results = NB_glm.fit(**NB_kwargs)
        model_results.append(NB_results)
        
        alpha_ = NB_results.params[1]
        bse_ = NB_results.bse[1]

        gene_alpha.append(alpha_)
        gene_alpha_bse.append(bse_)
    
    if verbose == True:
        np.save(model_results_path, model_results)

    return gene_alpha, gene_alpha_bse, model_results


def fit_Dispersions(Xdata, libsize_key = "library_ratio", verbose=True, 
    model_results_path=None, NB_kwargs={'disp': True, 'skip_hessian': True}):
    """
    Function:
    Estimate overdispersion based on mean value for each gene.
    
    todo
    2. multiprocessing for each gene?  
    
    Parameters:
    ----------
    libsize_key: sample_libsize, here can input total counts/learned libsize ratio

    Return:
    ------
    Xdata: anndata. with dispersion related items in `var` and model res in `uns`.
    
    Example:
    -------
    obs_ad1 = fit_Dispersions(obs_ad1, verbose=False, model_results_path=None)
    """
    import statsmodels as sm

    start_time_ = datetime.datetime.now()

    if sp.sparse.issparse(Xdata.X):
        X_mtx = Xdata.X.A
    else:
        X_mtx = Xdata.X
    
    if libsize_key is None:
        sample_libsize = X_mtx.sum(axis=1)  # cell total counts
        # total_counts = X_mtx.sum(axis=1) # wrong here-- need counts ratio[todo]
    else:
        sample_libsize = Xdata.obs[libsize_key]

    # check total count/libratio
    if (sample_libsize == 0).sum() > 0:
        raise ValueError("[XClone]Zero value in sample_libsize/total counts! Pls check!")
    
    ## init the output
    model_results = []
    gene_alpha = []
    gene_alpha_bse = []
     
    cell_num = X_mtx.shape[0]
    gene_num = X_mtx.shape[1]

    removed_gene_lst = []

    for i in range(gene_num):
        obs_y = X_mtx[:, i]
        feature_x = np.ones(cell_num)
        if verbose:
            print(cell_num)

        gene_name = Xdata.var["GeneName"][i]

        try:
            NB_glm = sm.discrete.discrete_model.NegativeBinomial(obs_y, feature_x, 
                                                            loglike_method = 'nb2', exposure = sample_libsize, offset = None, missing = 'none', check_rank = True)
            NB_results = NB_glm.fit(**NB_kwargs)

        except:
            
            removed_gene_lst.append(gene_name)
            model_results.append("Error in this gene:" + gene_name)

            alpha_ = np.NaN
            bse_ = np.NaN

            gene_alpha.append(alpha_)
            gene_alpha_bse.append(bse_)
        else:
            model_results.append(NB_results)
        
            alpha_ = NB_results.params[1]
            bse_ = NB_results.bse[1]

            gene_alpha.append(alpha_)
            gene_alpha_bse.append(bse_)
    
    if model_results_path is not None:
        np.save(model_results_path, model_results)
    
    if verbose:
        print(removed_gene_lst)
        print("Removed genes number", len(removed_gene_lst))

    Xdata.var["dispersion"] = gene_alpha
    Xdata.var["gene_dispersion_bse"] = gene_alpha_bse
    # Xdata.uns["fit_Dispersions_model"] = model_results
    Xdata.uns["fit_dispersion_removed_genes"] = removed_gene_lst

    # time stamp
    end_time_= datetime.datetime.now()
    time_used = end_time_ - start_time_
    print("[XClone RDR gene dispersion fitting] Time used:", time_used.seconds, "seconds")
    return Xdata


def map_var_info(Xdata_all, Xdata_specific, specific_celltype = None, var_key_lst = ["dispersion", "gene_dispersion_bse"],
                 uns_key_lst = ["fit_dispersion_removed_genes"]):
    """
    Function:
    map var information from celltype specific anndata back to completed anndata
    
    Parameters:
    -----------
    Xdata_all: anndata.
    Xdata_specific: anndata.

    Return:
    -------
    Xdata_all: updated anndata.
    
    Example:
    --------
    map_var_info(obs_ad1, obs_ad1_Paneth, specific_celltype = "Paneth", var_key_lst = ["dispersion", "gene_dispersion_bse"])
    """
    for anno_key in var_key_lst:
        Xdata_all.var[anno_key] = Xdata_specific.var[anno_key].copy()
    for anno_key in uns_key_lst:
        Xdata_all.uns[anno_key] = Xdata_specific.uns[anno_key].copy()
    
    Xdata_all.uns["dispersion_base_celltype"] = specific_celltype

    return Xdata_all

def remove_genes(Xdata, mode = "GLM"):
    """
    Function:
    removed genes with nan dispersion
    
    Parameters:
    ----------
    mode: NaN
    
    Return:
    -------
    
    Example:
    --------
    """
    # remove genes without GLM dispersion results
    if mode == "GLM":
        FLAG_ = Xdata.var["GeneName"].isin(Xdata.uns["fit_dispersion_removed_genes"])
        print("remove no GLM results genes num:", FLAG_.sum())
    # remove genes with nan results
    if mode == "NaN":
        FLAG_ = np.isnan(Xdata.var["dispersion"])
        print("remove nan dispersion genes num:", FLAG_.sum())
    if mode == "INF":
        FLAG_ = np.isinf(Xdata.var["dispersion"])
        print("remove inf dispersion genes num:", FLAG_.sum())
    
    Xdata = Xdata[:, ~FLAG_]
    
    # update_model_lst = list_update(Xdata.uns["fit_Dispersions_model"].copy(), ~FLAG_)
    # Xdata.uns["fit_Dispersions_model"] = update_model_lst
    return Xdata

def check_dispersion(Xdata, anno_key = "dispersion"):
    print("[XClone RDR gene-specific dispersion]: checking")
    print("max_value:", Xdata.var[anno_key].max())
    print("min_value:", Xdata.var[anno_key].min())

    print("qt_0.95_value:", np.quantile(Xdata.var[anno_key], 0.95))
    print("qt_0.05_value:", np.quantile(Xdata.var[anno_key], 0.05))


def dispersion_clip(Xdata, anno_key = "dispersion", 
                    min_threshold = 0, max_threshold = 20,
                    qt_low = 0.05, qt_up =0.95, verbose = True):
    """
    """
    update_anno_key = anno_key + "_capped"

    select_array = Xdata.var[anno_key]
    q_maxvalue = np.quantile(select_array, qt_up)
    q_minvalue = np.quantile(select_array, qt_low)

    if max_threshold is None:
        max_ = q_maxvalue
    else:
        max_ = min(q_maxvalue, max_threshold)
    
    if min_threshold is None:
        min_ = q_minvalue
    else:
        min_ = max(q_minvalue, min_threshold)
    
    Xdata.var[update_anno_key] = np.clip(
    Xdata.var[anno_key], min_, max_)
    print("[XClone RDR dispersion]: clipping")

    return Xdata


def dispersion_select(Xdata, anno_key = "dispersion", min_value = 0.01, max_value = None, 
                      qt_low = 0.05, qt_up =0.95, mode = "narrow", verbose = True):
    """
    Function:
    Select the rational estimated gene-specific dispersion
    
    Parameters:
    ----------
    Xdata: anndata.
    anno_key: char.
    min_value: float. predefined Lower limit.
    max_value: float. predefined up limits.
    qt_low: float. quantile percent.
    qt_up: float. quantile percent.
    mode: default mode is `narrow`. narrow down the boundary of dispersion 
          used in the following analysis.
          another mode is `filter`, which filter out the genes with dispersion 
          not in the specified boundary.
    verbose: bool.
    
    Return:
    -------
    Xdata: anndata. with updated var["dispersion"].

    Example:
    -------
    """
    update_anno_key = anno_key + "_capped"
    Xdata.var[update_anno_key] = Xdata.var[anno_key].copy()

    select_array = Xdata.var[anno_key]
    q_maxvalue = np.quantile(select_array, qt_up)
    q_minvalue = np.quantile(select_array, qt_low)

    if max_value is None:
        max_ = q_maxvalue
    else:
        max_ = min(q_maxvalue, max_value)
    
    if min_value is None:
        min_ = q_minvalue
    else:
        min_ = max(q_minvalue, min_value)

    FLAG_1 = select_array > max_
    FLAG_2 = select_array < min_

    if mode == "narrow":
        Xdata.var[update_anno_key][FLAG_1] = max_
        Xdata.var[update_anno_key][FLAG_2] = min_
        if verbose:
            print("max_value:", max_)
            print("min_value:", min_)
            print("change %d terms to be max" %FLAG_1.sum())
            print("change %d terms to be min" %FLAG_2.sum())

    elif mode == "filter":
        df_ = pd.DataFrame({"max":FLAG_1,
                            "min":FLAG_2})
        df_["final"] = df_[["max", "min"]].any(1)
        FLAG_ = df_["final"]
        Xdata = Xdata[:,~FLAG_]
        # update_model_lst = list_update(Xdata.uns["fit_Dispersions_model"].copy(), ~FLAG_)
        # Xdata.uns["fit_Dispersions_model"] = update_model_lst

    return Xdata