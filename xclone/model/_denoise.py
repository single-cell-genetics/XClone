"""Denoise functions for XClone probability.
"""

# Author: Rongting Huang

def denoise(Xdata, Xlayer, 
            brk = "chr_arm", 
            neutral_index = 2, 
            cnv_index = [0,1,3], 
            gmm_comp = 3,
            out_layer = "denoised_prob"):
    """
    
    """
    import numpy as np
    from sklearn.mixture import GaussianMixture
    Xmtx = Xdata.layers[Xlayer]
    argmax_ = np.argmax(Xmtx, axis = -1)
    
    brk_item = Xdata.var[brk].drop_duplicates(keep="first")
    
    proportion_ = []
    for brk_ in brk_item:
        tmp_region_flag = Xdata.var[brk] == brk_
        tmp_mtx = argmax_[:,tmp_region_flag]
        neutral_sum =  (tmp_mtx == neutral_index).sum()
        total_sum = tmp_mtx.shape[0] * tmp_mtx.shape[1]
        prop_ = (total_sum - neutral_sum)/total_sum
        proportion_.append(prop_)
#     return brk_item, proportion_
    if gene_denoise:

        neutral_sum = (argmax_ == neutral_index).sum(axis = 0)
        total_sum = argmax_.shape[0]
        prop_ = (total_sum - neutral_sum)/total_sum


    gmm = GaussianMixture(n_components =gmm_comp).fit(np.array(proportion_).reshape(-1,1))
    comp_means = gmm.means_.reshape(1,-1)
    sort_means = np.sort(comp_means)[0,:]
    threshold = sort_means[0] + 0.05
    
    flag_ = proportion_ <= threshold
    
    mask_brk_item = brk_item[flag_]
    print(mask_brk_item)
    mask_flag = Xdata.var[brk].isin(mask_brk_item)
    
    denoise_mtx  = Xmtx.copy()
    denoise_mtx[:,mask_flag, neutral_index] = 1
    for i in cnv_index:
        denoise_mtx[:,mask_flag, i] = 0
     
    Xdata.layers[out_layer] = denoise_mtx
    return Xdata

def denoise_gene_scale(Xdata, Xlayer,
            neutral_index = 2, 
            cnv_index = [0,1,3], 
            gmm_comp = 2,
            out_layer = "denoised_gene_prob"):
    """
    """
    import numpy as np
    from sklearn.mixture import GaussianMixture
    Xmtx = Xdata.layers[Xlayer]
    argmax_ = np.argmax(Xmtx, axis = -1)

    neutral_sum = (argmax_ == neutral_index).sum(axis = 0)
    total_sum = argmax_.shape[0]
    proportion_ = (total_sum - neutral_sum)/total_sum

    gmm = GaussianMixture(n_components = gmm_comp).fit(np.array(proportion_).reshape(-1,1))
    # comp_means = gmm.means_.reshape(1,-1)
    # sort_means = np.sort(comp_means)[0,:]
    # print(sort_means)
    lower_index = np.argmin(gmm.means_[:,0])
    print(gmm.means_)
    # if gmm.means_[0][0] < gmm.means_[1][0]:
    #     lower_index = 0
    # else:
    #     lower_index = 1
    
    label_ = gmm.predict(proportion_.reshape(-1,1))
    mask_flag = label_ == lower_index
    
    denoise_mtx  = Xmtx.copy()
    denoise_mtx[:,mask_flag, neutral_index] = 1
    for i in cnv_index:
        denoise_mtx[:,mask_flag, i] = 0
     
    Xdata.layers[out_layer] = denoise_mtx
    return Xdata
