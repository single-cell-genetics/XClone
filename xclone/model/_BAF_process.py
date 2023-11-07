## todo get_KNN_connectivities_from_expr and use different name.
from .base_utils import normalize
import scipy as sp
import scanpy as sc

def extra_preprocess_BAF(adata, Xlayer = "fill_BAF_phased",
                         KNN_neighbors = 10,
                         run_KNN = False, 
                         copy=False):
    """s
    """
    adata = adata.copy() if copy else adata

    ## generate KNN graph for smoothing across cell neighbourhood
    if run_KNN:
        raw_X = adata.X.copy()
        adata.X = adata.layers[Xlayer].copy()
        sc.pp.pca(adata)
        sc.pp.neighbors(adata, n_neighbors= KNN_neighbors, n_pcs=40)
        adata.X = raw_X
        adata.obsp['connectivities'] = normalize(adata.obsp['connectivities'])

    return adata if copy else None


def get_KNN_connectivities_from_expr(Xdata, expr_adata):
    """
    """
    ## check data cell order the same
    cell_barcodes_BAF =  Xdata.obs.index.values
    cell_barcodes_RDR =  expr_adata.obs.index.values

    flag_ = (cell_barcodes_BAF == cell_barcodes_RDR).sum() == Xdata.obs.shape[0]
    if flag_:
        Xdata.obsp["connectivities_expr"] = expr_adata.obsp["connectivities"].copy()
    else:
        raise ValueError("[XClone] Xdata cell order are not matched! Pls check and do preprocessing!")
    return Xdata

def BAF_prob_correct(Xdata, Xlayer1, Xlayer2):
    """
    todo: maybe add function to use 3 states layer 
    to correct 5 states layer BAF prob noise.
    """
    pass