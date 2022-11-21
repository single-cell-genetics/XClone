"""Base functions for XClone RDR preprocessing
mainly for genes selection before RDR model.
"""

# Author: Rongting Huang
# Date: 2021/12/16
# update: 2022/01/20

import numpy as np
import pandas as pd
import scanpy as sc


## Part I: filter marker genes for celltype
def get_markers(Xdata, 
                target_sum=1e4, 
                rank_group = "cell_type",
                test_method='wilcoxon',
                top_n =15, 
                marker_genes_key = "rank_marker_genes"
                ):
    """
    scanpy pipeline.

    Example:

    """
    Xdata_use = Xdata.copy()
    Xdata_use.var.index = Xdata_use.var["GeneName"]

    sc.pp.normalize_total(Xdata_use, target_sum=target_sum)
    sc.pp.log1p(Xdata_use)

    sc.tl.rank_genes_groups(Xdata_use, rank_group, method=test_method)
    top_marker_genes = pd.DataFrame(Xdata_use.uns['rank_genes_groups']['names']).head(top_n).values
    marker_genes_rank = np.unique(top_marker_genes)
    
    Xdata.uns[marker_genes_key] = marker_genes_rank
    return marker_genes_rank

def filter_markers(Xdata, anno_key="GeneName", marker_genes_key = None, top_n = 15, marker_lst = None):
    """
    based on the marker list.

    marker_genes_key: marker genes saved in Xdata.uns.  'rank_genes_groups'
    marker_lst: marker genes lst provided by users.

    """
    update_Xdata = Xdata.copy()
    if marker_genes_key is not None:
        marker_genes = Xdata.uns[marker_genes_key]['names'].head(top_n).values
        print("[XClone] use marker genes in uns:", marker_lst)
        is_marker = Xdata.var[anno_key].isin(marker_genes)
    elif marker_lst is not None:
        print("[XClone] use marker genes provided by users: \n", marker_lst)
        is_marker = Xdata.var[anno_key].isin(marker_lst)

    filter_genes_num = is_marker.sum()
    used_genes_num = (~is_marker).sum()
    print("filter_genes_num:", filter_genes_num)
    print("used_genes_num:", used_genes_num)

    adata_NOmarkers = update_Xdata[:, ~is_marker]
    return adata_NOmarkers


## Part II:
# FUNC--JSD_emm_genes  in HMM_base

from scipy.spatial import distance
import matplotlib.pyplot as plt
import seaborn as sns

def JSD_emm_genes(emm_prob, show_plot=True):
    """
    jensenshannon
    """
    start_ = distance.jensenshannon(emm_prob[0,:], emm_prob[1,:])
    end_ = distance.jensenshannon(emm_prob[-2,:], emm_prob[-1,:])
    
    # one step back
    emm_prob_bak = emm_prob[2:,:]
    # one step forward
    emm_prob_fwd = emm_prob[:-2,:]
    d1 = distance.jensenshannon(emm_prob[1:-1,:], emm_prob_bak, axis=1) # scipy >=1.7.0
    d2 = distance.jensenshannon(emm_prob[1:-1,:], emm_prob_fwd, axis=1)
    d_ = (d1+d2)/2
    
    gene_distance = np.insert(d_, 0, start_)
    gene_distance = np.append(gene_distance, end_)
    if show_plot:
        sns.displot(gene_distance, kde=True)
        plt.show()
    return gene_distance