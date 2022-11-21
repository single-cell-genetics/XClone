"""Base functions for XClone plotting test.
"""
# Author: Rongting Huang
# Date: 2021/07/23
# update: 2021/07/23


import random

def random_assign_celltype(Xdata, n_cluster, col_idxname = "random_celltype", celltype_prefix="celltype"):
    """
    generate cell/cluster type for Xheatmap testing.

    Examples:
    -----
    ## generate 5  randome cluster type for test_adata
    test_annotation = random_assign_celltype(test_adata, 5,"random_cluster","cluster")
    """
    n_obs = Xdata.n_obs
    Xdata.obs[col_idxname] = 'init'
    slice = random.sample(range(n_obs), n_cluster-1)
    slice.sort()
    slice.insert(0,0)
    slice.append(n_obs)
    for i in range(len(slice)-1):
        Xdata.obs[col_idxname][slice[i]:slice[i+1]] = celltype_prefix + str(i)
    return Xdata.obs[col_idxname]