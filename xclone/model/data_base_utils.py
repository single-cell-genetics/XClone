"""Utility functions for XClone data format."""

# Author: Rongting Huang
# Date: 2021-12-13
# update: 2021-03-10

import anndata as ad
import numpy as np

import matplotlib.pylab as plt

# preprocessing
def Xdata_mapping(Xdata, ref_Xdata):
    """
    Function:
    Actually, after many steps of preprocessing, cells and features will be filtered.
    And new dataset will be stored.
    Here can mapping the filterings on the original dataset for the follwing analysis.

    todo: However, maybe add new layer at the first is better.? not sure. but mapping is 
    another strategy. Based on the obs index and var cols.
    Parameters:
    ----------

    Return:
    -------

    Example:
    -------
    gx109_res = Xdata_mapping(gx109_res, ref_obs_ad1)

    """
    is_gene = Xdata.var["GeneName"].isin(ref_Xdata.var["GeneName"])
    is_cell =  Xdata.obs.index.isin(ref_Xdata.obs.index)

    Xdata = Xdata[is_cell, is_gene]

    return Xdata

# Part I: Exploring the data
def mtx_describe(Xdata, quantile_lst = None, low_lst = None, up_lst = None):
    """
    Function:
    basic description for Xdata.
    """
    # data
    if type(Xdata) == ad._core.anndata.AnnData:
        mtx_ = Xdata.copy().X
    if type(Xdata) == np.ndarray:
        mtx_ = Xdata.copy()

    array_ = mtx_.reshape(-1)
    

    print("mean value: ", array_.mean())
    print("max value: ", array_.max())
    print("min value: ", array_.min())
    # quantile
    if quantile_lst is None:
        quantile_lst = [0.25, 0.5, 0.75, 0.85, 0.90, 0.99]
    
    for qt_ in quantile_lst:
        print("quantile ", qt_, "value: ", np.quantile(array_, qt_))
    

    ## visualization
    figsize = (8, 8)
    sharex = False
    sharey = False
    
    if low_lst is None:
        low_lst = [0.01, 0.01, 0.01]
    
    if up_lst is None:
        up_lst = [2, 2.5, 3]
    
    plots_num = len(low_lst) + 1
    plots_row = 2
    plots_col = int(plots_num/2)


    fig, axs = plt.subplots(plots_row, plots_col, figsize = figsize, tight_layout = True, sharex=sharex, sharey=sharey)

    axs[0, 0].hist(array_)
    i_row_idx = []
    j_col_idx = []
    for i_row_ in range(plots_row):
        for j_col_ in range(plots_col):
            i_row_idx.append(i_row_)
            j_col_idx.append(j_col_)

    for i_row, j_col, low_, up_ in zip(i_row_idx[1:], j_col_idx[1:], low_lst, up_lst):
        print("boundary:", low_, up_)
        anno_ = "boundary:" + str(low_) + "-" + str(up_)
        flag_ = (array_ > low_) & (array_ <= up_)
        # for i_row in range(plots_row):
        #     for j_col in range(1, plots_col):
        axs[i_row, j_col].hist(array_[flag_])
        axs[i_row, j_col].set_title(anno_)
    
    plt.show()
    return array_


# Part II: preprocessing for generating valid data
def mtx_bound(Xdata, upper_lim = None, lower_lim = None, 
              up_quantile = None, low_quantile = None):
    """
    Function:
    Detect the outliers in the Xdata matrix and set boundary.

    Parameters:
    ----------
    upper_lim: set upper limit for the mtx.
    lower_lim: set lower limit for the mtx.

    up_quantile:
    low_quantile:


    Return:
    -------
    Xdata_update: XClone anndata.

    Example:
    -------
    Xdata_update = mtx_bound(Xdata, upper_lim = None, lower_lim = None, 
              up_quantile = None, low_quantile = None)

    Xdata_update = mtx_bound(Xdata, upper_lim = None, lower_lim = None, 
              up_quantile = None, low_quantile = None)

    """

    mtx_ = Xdata.copy().X ## notes:use copy and maintain the original Xdata.
    if up_quantile is None:
        pass
    else:
        mtx_array = mtx_.reshape(-1)
        upper_lim = np.quantile(mtx_array, up_quantile)
        lower_lim = np.quantile(mtx_array, low_quantile)

    flag1_ = mtx_ > upper_lim
    flag2_ = mtx_ < lower_lim
    mtx_[flag1_] = upper_lim
    mtx_[flag2_] = lower_lim

    Xdata_update = Xdata.copy()
    Xdata_update.X = mtx_

    return Xdata_update

def mtx_bound1(mtx_, upper_lim = None, lower_lim = None, up_quantile = None, low_quantile = None):
    """
    """

#     mtx_ = Xdata.X
    if up_quantile is None:
        pass
    else:
        mtx_array = mtx_.reshape(-1)
        upper_lim = np.quantile(mtx_array, up_quantile)
        lower_lim = np.quantile(mtx_array, low_quantile)

    flag1_ = mtx_ > upper_lim
    flag2_ = mtx_ < lower_lim
    mtx_[flag1_] = upper_lim
    mtx_[flag2_] = lower_lim

#     Xdata.X = mtx_

    return mtx_

