"""Denoise functions for XClone probability.
"""

# Author: Rongting Huang

def denoise_gene_scale(Xdata, Xlayer,
                       neutral_index = 2, 
                       cnv_index = [0,1,3], 
                       GMM_detection = True,
                       gmm_comp = 2,
                       cell_prop_cutoof = 0.05,
                       out_layer = "denoised_gene_prob"):
    """
    default method: GMM detect the cell proportion identified as CNV states for each gene.
    alternative method: set cell_prop_cutoof.[not test yet, todo]
    """
    import numpy as np
    from sklearn.mixture import GaussianMixture
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
        mask_flag = proportion_ <= cell_prop_cutoof
    
    denoise_mtx  = Xmtx.copy()
    denoise_mtx[:,mask_flag, neutral_index] = 1
    for i in cnv_index:
        denoise_mtx[:,mask_flag, i] = 0
     
    Xdata.layers[out_layer] = denoise_mtx
    return Xdata

import numpy as np
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