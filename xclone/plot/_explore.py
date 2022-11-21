"""Base functions for XClone debug visualization.
"""
# Author: Rongting Huang
# Date: 2022/08/29
# update: 2022/08/29

# smart-seq BCH dataset

import numpy as np
import random
import pandas as pd
from ._BAF_debug import transfer_prob_to_df

def visualize_prob(Xdata, 
                   Xlayer = "posterior_mtx",
                   cell_anno_key = "cell_type",
                   region_key = "chr_arm",
                   chr_lst = ["1q"], 
                   random_lst = None,
                   CNV_state_order = ["copy loss", "copy neutral", "copy gain"]):
    """
    """
    print("[XClone] visualization probability of chr: ", chr_lst)
    
    region_flag = Xdata.var[region_key].isin(chr_lst)
    prob_used = Xdata.layers[Xlayer][:, region_flag,:].copy()
    Xdata_used = Xdata[:, region_flag].copy()

    if random_lst is not None:
        prob_used = prob_used[:, random_lst, :]
        Xdata_used = Xdata_used[:, random_lst]

    transferrred_df = transfer_prob_to_df(prob_used, Xdata_used, 
                                          CNV_state_order = CNV_state_order)
    
    data_used = transferrred_df.copy()
    obs_df = Xdata_used.obs.copy()
    data_used1 = pd.merge(data_used, obs_df[cell_anno_key], 
                          left_on = "cell_barcodes", right_index = True)
    
    ## visualization
    import seaborn as sns
    import matplotlib.pylab as plt

    ## for all cells(no clear information)
    fig, ax = plt.subplots(figsize=(40, 10))
    ax = sns.boxplot(x="bins_loc", y="prob", hue="CNV_states", data=data_used1, 
                    palette="Set3", fliersize=1,
                    linewidth=1)
    ax.set_title('Box plot of beta binomial prob for different CNV states')
    
    ## for each celltype
    sns.catplot(x="bins_loc", y="prob", hue="CNV_states", data=data_used1, 
                palette="Set3", fliersize=1,
                linewidth=1, col="cell_type", kind="box",
                height=4, aspect=5, orient = "v", col_wrap =1)
    plt.show()

    return Xdata_used

from ._BAF_debug import get_count_df

def visualize_count(Xdata,
                    Xlayer = "raw_expr",
                    cell_anno_key = "cell_type",
                    region_key = "chr_arm",
                    chr_lst = ["1q"]
                    ):
    """
    """
    data_used = get_count_df(Xdata, Xlayer, cell_anno_key, region_key, chr_lst)
    
    # visualization
    import seaborn as sns
    import matplotlib.pylab as plt
    ## for all cells(no clear information)
    fig, ax = plt.subplots(figsize=(40, 10))
    ax = sns.boxplot(x="bins_loc", y="count", data=data_used,
                    palette="Set3", fliersize=1,linewidth=1)
    ax.set_title('Box plot of DP count-all cells and for different cell type')

    ## for each celltype
    sns.catplot(x="bins_loc", y="count", data=data_used, 
                palette="Set3", fliersize=1,
                linewidth=1, col=cell_anno_key, kind="box",whis=2, 
                height=4, aspect=5, orient = "v", col_wrap =1)
    plt.show()

def compare_prob(Xdata, 
                 pre_Xlayer = "emm_prob_noHMM",
                 after_Xlayer = "posterior_mtx",
                 cell_anno_key = "cell_type",
                 region_key = "chr_arm",
                 chr_lst = ["1q"], 
                 random_select = True,
                 CNV_state_order = ["copy loss", "copy neutral", "copy gain"]):
    """
    """
    region_flag = Xdata.var[region_key].isin(chr_lst)
    if random_select == True:
        length_ = region_flag.sum()
        range_lst = list(range(length_))
        random_lst = random.sample(range_lst, 8)
        print(random_lst)
        
        visualize_prob(Xdata, pre_Xlayer, cell_anno_key = cell_anno_key,
                 region_key = region_key,
                 chr_lst = chr_lst, 
                 random_lst = random_lst,
                 CNV_state_order = CNV_state_order)
        
        Xdata_used = visualize_prob(Xdata, after_Xlayer, cell_anno_key = cell_anno_key,
                 region_key = region_key,
                 chr_lst = chr_lst, 
                 random_lst = random_lst,
                 CNV_state_order = CNV_state_order)
        
        visualize_count(Xdata_used, "raw_expr", 
                        cell_anno_key = cell_anno_key,
                        region_key = region_key,
                        chr_lst = chr_lst)


