"""Denoise functions for XClone probability.
"""

# Author: Rongting Huang

import numpy as np
import anndata as ad
from sklearn.mixture import GaussianMixture

from ..analysis._clustering import XClustering

def denoise_gene_scale(Xdata, Xlayer,
                       neutral_index = 2, 
                       cnv_index = [0,1,3,4], 
                       GMM_detection = True,
                       gmm_comp = 2,
                       cell_prop_cutoff = 0.05,
                       out_layer = "denoised_gene_prob"):
    """
    default method: GMM detect the cell proportion identified as CNV states for each gene.
    alternative method: set cell_prop_cutoff.
    """
    Xmtx = Xdata.layers[Xlayer]
    argmax_ = np.argmax(Xmtx, axis = -1)

    neutral_sum = (argmax_ == neutral_index).sum(axis = 0)
    total_sum = argmax_.shape[0]
    proportion_ = (total_sum - neutral_sum)/total_sum
    if GMM_detection:
        gmm = GaussianMixture(n_components = gmm_comp).fit(np.array(proportion_).reshape(-1,1))
        lower_index = np.argmin(gmm.means_[:,0])
        print(gmm.means_)

        label_ = gmm.predict(proportion_.reshape(-1,1))
        mask_flag = label_ == lower_index 
    else:
        mask_flag = proportion_ <= cell_prop_cutoff
    
    denoise_mtx  = Xmtx.copy()
    denoise_mtx[:,mask_flag, neutral_index] = 1
    for i in cnv_index:
        denoise_mtx[:,mask_flag, i] = 0
     
    Xdata.layers[out_layer] = denoise_mtx
    return Xdata

def denoise_rdr_by_baf(Xdata, 
                       RDR_layer, 
                       BAF_layer, 
                       out_RDR_layer,
                       cell_prop_cutoff = 0.999
                       ):
    """
    default method: GMM detect the cell proportion identified as allele bias states, and use the states for rdr denoise.
    """
    Xmtx = Xdata.layers[BAF_layer]
    if Xmtx.shape[-1] == 3:
        b_neutral_index = 1
    elif Xmtx.shape[-1] == 5:
        b_neutral_index = 2


    argmax_ = np.argmax(Xmtx, axis = -1)

    neutral_sum = (argmax_ == b_neutral_index).sum(axis = 0)
    total_sum = argmax_.shape[0]
    proportion_ = neutral_sum/total_sum
    mask_flag = proportion_ >= cell_prop_cutoff
    
    denoise_mtx  = Xdata.layers[RDR_layer].copy()
    if denoise_mtx.shape[-1] == 3:
        r_neutral_index = 1
        r_cnv_index = [0,2]
    r_argmax_ = np.argmax(denoise_mtx, axis = -1)
    r_neutral_sum = (r_argmax_ == r_neutral_index).sum(axis = 0)
    r_total_sum = r_argmax_.shape[0]
    r_proportion_ = r_neutral_sum/r_total_sum
    r_mask_flag = r_proportion_ <=0.2

    mask_flag = mask_flag*r_mask_flag

    denoise_mtx[:,mask_flag, r_neutral_index] = 1
    for i in r_cnv_index:
        denoise_mtx[:,mask_flag, i] = 0
    
    Xdata.layers[out_RDR_layer] = denoise_mtx
    return Xdata

def denoise_by_cluster(Xdata, 
                       Xlayer,
                       neutral_index = 2, 
                       cnv_index = [0,1,3],
                       GMM_detection = True,
                       gmm_comp = 2,
                       cell_prop_cutoff = 0.3,
                       out_layer = "XC_denoise",
                       Clustering_X = True,
                       Clustering_layer = None,
                       N_Cluster = 2,
                       cell_prop_threshold = 0.1,
                       cluster_anno_key = None):
    """
    Function:
    ---------
    develop version.

    Params:
    -------
    Xdata: anndata.
    Xlayer: in_layer to perform `denoise_by_cluster` FUNC.
    neutral_index: `denoise_by_cluster` params.
    cnv_index: `denoise_by_cluster` params.
    GMM_detection: `denoise_by_cluster` params.
    gmm_comp: `denoise_by_cluster` params.
    cell_prop_cutoff: `denoise_by_cluster` params.
    out_layer: out_layer to save `denoise_by_cluster` results.

    Clustering_X: bool. default True: perform Xclustering to get cluster ID.
    Clustering_layer: in_layer to perform `XClustering`. 
                      default use `Xlayer` from combine module, also can try different 
                      module's layer, e.g., prob from BAF.
    N_Cluster: Specify cluster numbers for `XClustering`.
    cell_prop_threshold: cells proportion less than cutoff, skip denoise strategy.

    cluster_anno_key: if `Clustering_X` is False, need specify an existed cluster_anno_key. 

    Example:
    --------


    """
    ## Step1: Perform XClustering before denoise_by_cluster.
    ## (if Clustering_X is default True.)
    if Clustering_X:
        if Clustering_layer is None:
            Clustering_layer = Xlayer
        
        Z, embedding = XClustering(Xdata, Clustering_layer, n_clusters = N_Cluster)

        Xdata.obs["pre_denoise_XClone_ID"] = Z
        cluster_anno_key = "pre_denoise_XClone_ID"
    else:
        if cluster_anno_key is None:
            raise ValueError("[XClone] warning: please specify a cluster anno key.")

    ## Step2: Perform denoise_by_cluster.
    Cluster_ID = np.unique(Xdata.obs[cluster_anno_key])
    denoise_Xdata = []
    for XCluster_ID in Cluster_ID:
        flag_ = Xdata.obs[cluster_anno_key] == XCluster_ID
        tmp_Xdata = Xdata[flag_, :].copy()

        if tmp_Xdata.shape[0]/Xdata.shape[0] < cell_prop_threshold:
            tmp_Xdata.layers[out_layer] = tmp_Xdata.layers[Xlayer].copy()
        else:
            tmp_Xdata = denoise_gene_scale(tmp_Xdata, Xlayer,
                                           neutral_index = neutral_index, 
                                           cnv_index = cnv_index, 
                                           GMM_detection = GMM_detection,
                                           gmm_comp = gmm_comp,
                                           cell_prop_cutoff = cell_prop_cutoff,
                                           out_layer = out_layer)
        denoise_Xdata.append(tmp_Xdata)

    ordered_Xdata = ad.concat(denoise_Xdata, merge="same")
    return ordered_Xdata

def merge_denoise_prob(Xdata,
              in_layer,
              out_layer):
    """
    default: merge_index = [0,1] for BCH869 dataset. 
    todo:extend to other merge situation
    
    Example:
    --------
    >>> combine_Xdata = merge_denoise_prob(combine_Xdata, 
                           in_layer = "denoised_gene_prob",
                           out_layer = "denoised_gene_prob_merge")
    """
    Xdata.layers[out_layer] = Xdata.layers[in_layer].copy()
    Xdata.layers[out_layer][:,:,0] = Xdata.layers[out_layer][:,:,0] + Xdata.layers[out_layer][:,:,1]
    Xdata.layers[out_layer] = np.delete(Xdata.layers[out_layer], [1], axis = -1)
    return Xdata