"""Base functions for XClone Gaussian model.
"""

# Author: Rongting Huang
# Date: 2021-05-24
# update: 2022-05-24

import numpy as np
import pandas as pd
import scipy as sp

import anndata as ad
import itertools
import datetime

## init the framework

import numpy as np
from sklearn.mixture import GaussianMixture
from timeit import default_timer as timer

def get_initial_means(X, init_params, r):
    # Run a GaussianMixture with max_iter=0 to output the initalization means
    gmm = GaussianMixture(
        n_components=5, init_params=init_params, tol=1e-9, max_iter=0, random_state=r
    ).fit(X)
    return gmm.means_

def get_gmm_means(X, init_params, r, max_iter=50):
    gmm = GaussianMixture(
        n_components=4, init_params=init_params, tol=1e-9, max_iter=max_iter, random_state=r
    ).fit(X)
    return gmm.means_


def get_CNV_states(Xdata, Xlayer = "BAF_phased", 
                   gene_specific = False,
                   n_components = 5, 
                   means_init = None,
                   max_iter = 50, 
                   **kwargs):
    """
    """
    start = timer()

    r = np.random.RandomState(seed=1234)
    if gene_specific:
        X = Xdata.layers[Xlayer].copy()
    else:
        X = Xdata.layers[Xlayer].copy().reshape(-1,1)

    params = {}
    params["n_components"] = n_components
    params["max_iter"] = max_iter
    params["random_state"] = r
    params["tol"] = 1e-9
    if means_init is not None:
        params["means_init"] = means_init
    params.update(**kwargs)
    gmm = GaussianMixture(**params).fit(X)

    end = timer()
    elapsed_sec = end - start
    print("[XClone get_CNV_states] time_used: " + "{:.2f}".format(elapsed_sec) + "seconds")
    return gmm.means_

def guide_BAF_theo_states(CNV_states, threshold = 0.1, theo_neutral = 0.5):
    """
    """
    sort_array = np.sort(CNV_states[:,0])
    
    states_num = len(sort_array)
    
    if states_num == 3:
        guide_cnv_ratio = sort_array[[0,2]]
        guide_cnv_ratio = correct_guide_BAF_cnv(guide_cnv_ratio, threshold = threshold, theo_neutral = theo_neutral)
    elif states_num == 5:
        guide_cnv_ratio = sort_array[[0,1,3,4]]
        guide_cnv_ratio = correct_guide_BAF_cnv(guide_cnv_ratio,  threshold = threshold, theo_neutral = theo_neutral)
    
    return guide_cnv_ratio

def correct_guide_BAF_cnv(guide_cnv_ratio, threshold = 0.1, theo_neutral = 0.5):
    """
    if guide cnv ratio values are similar,(then maybe no cnv states exists here) 
    return default theoratical value

    todo: may change threshold.
    """
    states_num = len(guide_cnv_ratio)
    if states_num == 2:
        if abs(theo_neutral - guide_cnv_ratio[0]) < threshold:
            guide_cnv_ratio[0] = 0.3
            print("correct BAF CNV guiding copy loss-B ratio")
        if abs(guide_cnv_ratio[1] - theo_neutral) < threshold:
            guide_cnv_ratio[1] = 0.7
            print("correct BAF CNV guiding copy loss-A ratio")
    if states_num == 4:
        ## todo need correct again maybe.
        if abs(theo_neutral - guide_cnv_ratio[0]) < 2*threshold:
            guide_cnv_ratio[0] = 0.15
            print("correct BAF CNV guiding copy loss-B ratio")
        if abs(guide_cnv_ratio[3] - theo_neutral) < 2*threshold:
            guide_cnv_ratio[3] = 0.85
            print("correct BAF CNV guiding copy loss-A ratio")
        
        if abs(theo_neutral - guide_cnv_ratio[1]) < threshold:
            guide_cnv_ratio[1] = 1/3
            print("correct BAF CNV guiding copy gain-A ratio")
        if abs(guide_cnv_ratio[2] - theo_neutral) < threshold:
            guide_cnv_ratio[2] = 2/3
            print("correct BAF CNV guiding copy gain-B ratio")

    return guide_cnv_ratio