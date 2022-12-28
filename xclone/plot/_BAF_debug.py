"""Base functions for XClone BAF HMM_BB debug
"""

# Author: Rongting Huang
# Date: 2022/05/25
# update: 2022/06/09

import datetime
import numpy as np
import pandas as pd
import scipy as sp

## debug analysis part

## debug plotting

def bulk_AF_chr(Xdata, chr_select = "10", 
                cell_anno_key = "Clone_ID", cluster_ID = 1,
                ad_key = "AD",
                dp_key = "DP",
                fig_title = "bulk AD/DP"):
                
    """
    debug BCH869 chr10p AF.
    """
    import matplotlib.pylab as plt
    BAF_adata = Xdata[:, Xdata.var["chr"] == chr_select]
    cell_flag = Xdata.obs[cell_anno_key] == cluster_ID
    BAF_adata = BAF_adata[cell_flag, :].copy()
    
    BAF_adata_p = BAF_adata[:, BAF_adata.var["arm"] == "p"].copy()
    # BAF_adata_q = BAF_adata[:, BAF_adata.var["arm"] == "q"].copy()
    
    plt.plot(BAF_adata.layers[dp_key].sum(axis=0).reshape(-1, 1),alpha = 0.6)
    plt.plot(BAF_adata_p.layers[dp_key].sum(axis=0).reshape(-1, 1), alpha = 0.6)

    plt.plot(BAF_adata.layers[ad_key].sum(axis=0).reshape(-1, 1), alpha = 0.8)
    plt.plot(BAF_adata_p.layers[ad_key].sum(axis=0).reshape(-1, 1), alpha = 0.8)
    plt.title(fig_title)


def state_specific_df(state_prob, Xdata, CNV_states_type = "copy_loss"):
    """
    """
    # init the matrix to dataframe
    data_mtx = state_prob.copy()
    row_index = Xdata.obs.index.copy() 
    col_ = Xdata.var.index.copy()
    data_df = pd.DataFrame(data = data_mtx, index = row_index, columns = col_)
    
    data_df.index.name = "cell_barcodes"
    data_df.columns.name = "bins_loc"
    # stack the data to get the prob value series
    data_df = data_df.stack()
    # reset index [Multiindex will be columns]
    data_df = data_df.reset_index()
    
    data_df.columns = ("cell_barcodes", "bins_loc", "prob")
    
    data_df["CNV_states"] = CNV_states_type
    
    return data_df
    


def transfer_prob_to_df(emm_prob_log, Xdata, CNV_state_order = ["copy loss(0.01)", "copy loss(0.99)", "copy neutral", "copy gain(1/3)", "copy gain(2/3)"]):
    """
    
    """
    prob_use = emm_prob_log.copy()
    update_Xdata = Xdata.copy()
    
    ## check cnv states
    if len(CNV_state_order) == prob_use.shape[-1]:
        pass
    else:
        raise ValueError("[XClone]-pls check the given params: CNV_state_order.")
    
    df_lst = []
    for k in range(prob_use.shape[-1]):
        state_prob = prob_use[:,:,k]
        tmp_df = state_specific_df(state_prob, update_Xdata, CNV_state_order[k])
        df_lst.append(tmp_df)
        
    concated_df = pd.concat(df_lst)
    return concated_df

##debug-exploring research    
## try to do other pivot tries and do the visualization

def visualize_emm_prob(Xdata, 
                       uns_key = "posterior_mtx",
                       chr_lst = ["1q"], 
                       CNV_state_order = ["copy loss(0.01)", "copy loss(0.99)", "copy neutral", "copy gain(1/3)", "copy gain(2/3)"]):
    """
    """
    print("[XClone]visualization emm_prob of chr: ", chr_lst )
    bin_flag = Xdata.var["chr_arm"].isin(chr_lst)
    
    prob_used = Xdata.uns[uns_key][:, bin_flag,:].copy()
    Xdata_used = Xdata[:, bin_flag].copy()
    
    transferrred_df = transfer_prob_to_df(prob_used, Xdata_used, CNV_state_order = CNV_state_order)
    
    
    data_used = transferrred_df.copy()
    obs_df = Xdata_used.obs.copy()
    data_used1 = pd.merge(data_used, obs_df["cell_type"], left_on = "cell_barcodes", right_index = True)
    
    import seaborn as sns
    import matplotlib.pylab as plt
    # %matplotlib inline
    
    ## for all cells(no clear information)
    fig, ax = plt.subplots(figsize=(40, 10))
    ax = sns.boxplot(x="bins_loc", y="prob", hue="CNV_states", data=data_used1, palette="Set3", fliersize=1,
                     linewidth=1)
    ax.set_title('Box plot of beta binomial prob for different CNV states')
    
    ## for each celltype
    sns.catplot(x="bins_loc", y="prob", hue="CNV_states", data=data_used1, palette="Set3", fliersize=1,
                     linewidth=1, col="cell_type", kind="box",
                height=4, aspect=5, orient = "v", col_wrap =1)
    plt.show()
    return None
## some color setting test
# import palettable
# sns.catplot(x="bins_loc", y="prob", hue="CNV_states", data=data_used1, palette=palettable.lightbartlein.diverging.BlueDarkRed12_5.mpl_colors, fliersize=1,
#                      linewidth=1, col="cell_type", kind="box",
#                 height=6, aspect=0.8, orient = "v", col_wrap =4)
# # ax.set_title('Box plot of beta binomial prob for different CNV states')
# # plot.fig.suptitle("Value of Tips Given to Waiters, by Days of the Week and Sex",
# #                   fontsize=24, fontdict={"weight": "bold"})
# sns.catplot(x="bins_loc", y="prob", hue="CNV_states", data=data_used1, palette=palettable.scientific.diverging.Vik_5.mpl_colors, fliersize=1,
#                      linewidth=1, col="cell_type", kind="box",
#                 height=6, aspect=0.8, orient = "v", col_wrap =4)

# sns.catplot(x="bins_loc", y="prob", hue="CNV_states", data=data_used1, palette="husl", fliersize=1,
#                      linewidth=1, col="cell_type", kind="box",
#                 height=6, aspect=0.8, orient = "v", col_wrap =4)


def get_count_df(Xdata, Xlayer, cell_anno_key = "cell_type", 
                 region_key = "chr_arm", chr_lst = None):
    """
    """
    bin_flag = Xdata.var[region_key].isin(chr_lst)
    Xdata_used = Xdata[:, bin_flag].copy()

    # init the matrix to dataframe
    if sp.sparse.issparse(Xdata_used.layers[Xlayer]):
        data_mtx = Xdata_used.layers[Xlayer].A.copy()
    else:
        data_mtx = Xdata_used.layers[Xlayer].copy()
    
    row_index = Xdata_used.obs.index.copy() 
    col_ = Xdata_used.var.index.copy()
    data_df = pd.DataFrame(data = data_mtx, index = row_index, columns = col_)

    data_df.index.name = "cell_barcodes"
    data_df.columns.name = "bins_loc"

    # stack the data to get the prob value series
    data_df = data_df.stack()
    # reset index [Multiindex will be columns]
    data_df = data_df.reset_index()

    data_df.columns = ("cell_barcodes", "bins_loc", "count")

    df_used = pd.merge(data_df, Xdata.obs[cell_anno_key], left_on = "cell_barcodes", right_index = True)

    return df_used

def visualize_DP(Xdata, Xlayer = "dp_bin", chr_lst = ["18p", "18q"]):
    """
    """
    data_used = get_count_df(Xdata, Xlayer, chr_lst)

    import seaborn as sns
    import matplotlib.pylab as plt
    ## for all cells(no clear information)
    fig, ax = plt.subplots(figsize=(40, 10))
    ax = sns.boxplot(x="bins_loc", y="count",  data=data_used, palette="Set3", fliersize=1,
                     linewidth=1)
    ax.set_title('Box plot of DP count-all cells and for different cell type')

    ## for each celltype
    sns.catplot(x="bins_loc", y="count", data=data_used, palette="Set3", fliersize=1,
                     linewidth=1, col="cell_type", kind="box",whis=2, 
                height=4, aspect=5, orient = "v", col_wrap =1)
    plt.show()   

    return None

def BAF_to_df(Xdata, Xlayer = "fill_BAF_phased"):
    """
    """
    data_mtx = Xdata.layers[Xlayer].copy()
    row_index = Xdata.obs.index.copy()
    col_ = Xdata.var.index.copy()
    data_df = pd.DataFrame(data = data_mtx, index = row_index, columns = col_)
    
    data_df.index.name = "cell_barcodes"
    data_df.columns.name = "bins_loc"
    
    data_df = data_df.stack()
    # reset index [Multiindex will be columns]
    data_df = data_df.reset_index()
    
    data_df.columns = ("cell_barcodes", "bins_loc", "BAF")
    
    return data_df

def visualize_BAF_distribution(Xdata, Xlayer = "fill_BAF_phased", select_anno_key = "chr_arm", chr_lst = ["1q"]):
    """
    """
    bin_flag = Xdata.var[select_anno_key].isin(chr_lst)
    Xdata_used = Xdata[:, bin_flag].copy()
    
    transferrred_df = BAF_to_df(Xdata_used, Xlayer)
    ## add celltype info
    data_used = transferrred_df.copy()
    obs_df = Xdata_used.obs.copy()
    data_used1 = pd.merge(data_used, obs_df["cell_type"], left_on = "cell_barcodes", right_index = True)
    
    import seaborn as sns
    import matplotlib.pylab as plt
    # %matplotlib inline
    
    fig_waitTriage, ax = plt.subplots(figsize=(40, 10))
    ax = sns.boxplot(x="bins_loc", y="BAF", hue="cell_type", data=data_used1, palette="Set3", fliersize=1,
                     linewidth=1)
    #place legend outside top right corner of plot
    plt.legend(bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0)
    ax.set_title('Box plot of BAF')
    plt.show()
    
    return data_used1