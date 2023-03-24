"""Base functions for XClone RDR processing
strategies for reference biological information removing
"""

# Author: Rongting Huang
# Date: 2021-03-10
# update: 2021-03-18

## Part-V reference biological ratio prediction based on NMF

import numpy as np
import scanpy as sc


from .base_utils import normalize

def extra_preprocess(adata, ref_celltype, cluster_key='cell_type',
                     avg_key = "ref_avg", depth_key='counts_ratio', 
                     low_dim=False, run_KNN=False, KNN_neighbors = 10,
                     copy=False):
    """
    """
    import scipy as sp
    adata = adata.copy() if copy else adata
    
    if sp.sparse.issparse(adata.X):
        Xmtx = adata.X.A
    else:
        Xmtx = adata.X
    
    if ref_celltype not in list(adata.obs[cluster_key]):
        print("Error: %s not exist as ref cell type" %(ref_celltype))
        return None
        
    _is_ref = adata.obs[cluster_key] == ref_celltype
    adata.var['ref_avg'] = Xmtx[_is_ref, :].mean(axis=0)
    adata.obs['counts_ratio'] = Xmtx.sum(axis=1) / adata.var['ref_avg'].sum()

    ## check before generate expected layer
    if depth_key == "library_ratio_capped":
        adata.obs[depth_key] = np.where(np.isinf(adata.obs["library_alpha"]), adata.obs["counts_ratio"], adata.obs[depth_key])
        ## todo capped again; add in later version[not test yet] 20230106
        # adata.obs[depth_key] = np.where(adata.obs[depth_key] < 0.001*adata.obs["counts_ratio"], adata.obs["counts_ratio"], adata.obs[depth_key])

    ## generate normalised
    X_baseline = (adata.obs[depth_key].values.reshape(-1, 1) *
                  adata.var[avg_key].values.reshape(1, -1))
    adata.layers['ref_normalized'] = (Xmtx / X_baseline)
    
    ## generate lower dimension on residues
    if low_dim:
        _off_set = 0.3 # offset of normalized value; can be further tuned
        X_residual = np.log(adata.layers['ref_normalized'] + _off_set)

        ## issue: negative values
        _model = PCA(n_components=20)
        _model.fit(X_residual[_is_ref, :])

        W = _model.transform(X_residual)
        H = _model.components_

        adata.obsm['Residual_W'] = W
        adata.varm['Residual_H'] = H.T
        X_predict = np.clip(np.dot(W, H), -3, 3) - _off_set
        adata.layers['expected'] = np.exp(X_predict) * X_baseline
    else:
        adata.layers['expected'] = X_baseline

    ## generate KNN graph for smoothing across cell neighbourhood
    if run_KNN:
        import scanpy as sc
        raw_X = adata.X.copy()
        adata.X = np.log(adata.layers['ref_normalized'] + 0.3)
        sc.pp.pca(adata)
        sc.pp.neighbors(adata, n_neighbors = KNN_neighbors, n_pcs=40)
        ## Notes: connectivities and distances can be slightly different rvery run
        ## even the random_state = 0 (default).
        adata.X = raw_X
        ## connectivities normalization
        adata.obsp['connectivities'] = normalize(adata.obsp['connectivities'])
    
    return adata if copy else None


## (1) modelling
from sklearn.decomposition import NMF, PCA

def NMF_confounder(Xdata, 
                   anno_key = "cell_type",
                   normalization_method = "sc",
                   log_transform = True,
                   ref_celltype = "unclassified",
                   NMF_components = 20,
                   **kwargs):
    """
    Function:
    Reference cells used for biological ratio prediction.

    Parameters:
    -----------
    Xdata: anndata.
    anno_key: str.
    normalization_method:"sc" or "xc", scanpy/xclone method, same results.
    log_transform: bool.
    ref_celltype: str.
    NMF_components: int. `n_components` in sklearn.decomposition.NMF
    
    **kwargs: params in sklearn.decomposition.NMF
    https://scikit-learn.org/stable/modules/generated/
    sklearn.decomposition.NMF.html#sklearn.decomposition.NMF
    max_iter: int
    random_state: int

    Example:
    -------
    gx109_pred = xclone.model.NMF_confounder(gx109_adata, 
    normalization_method = "xc", max_iter = 300, random_state = 2)
    """
    ## filtering data without annotaion
    valid_cells = Xdata.obs[anno_key] == Xdata.obs[anno_key]
    Xdata = Xdata[valid_cells, :]
    
    ## normalization and log transformation
    if normalization_method == "xc":
        Xdata_norm = Xdata.copy()
        if log_transform:
            Xdata_norm.X = np.log(Xdata_norm.X/Xdata_norm.X.sum(1) * 10000 + 1)
        else:
            Xdata_norm.X = Xdata_norm.X/Xdata_norm.X.sum(1) * 10000

    if normalization_method == "sc":
        sc.pp.normalize_total(Xdata, target_sum=1e4)
        if log_transform:
            sc.pp.log1p(Xdata)
        Xdata_norm = Xdata.copy()
    
    ## NMF modelling
    _is_ref = Xdata_norm.obs[anno_key] == ref_celltype

    _X_norm = np.array(Xdata_norm.X) + 0.0
    
    ## NMF decomposition params setting
    NMF_params = {}
    NMF_params["n_components"] = NMF_components
    NMF_params.update(**kwargs)

    _model = NMF(**NMF_params)
    _model.fit(_X_norm[_is_ref, :])

    W = _model.transform(_X_norm)
    H = _model.components_
    _X_predicted = np.dot(W, H)

    Xdata_pred = Xdata_norm.copy()
    Xdata_pred.X = _X_predicted

    return Xdata_pred
