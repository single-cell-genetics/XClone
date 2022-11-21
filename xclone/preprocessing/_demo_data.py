"""Base functions for XClone demo data loading.
"""

# Author: Rongting Huang
# Date: 2021/07/15
# update: 2022/07/13

import pkg_resources
import pandas as pd
import scipy as sp

from ._data import xclonedata
# import anndata as an
import warnings
import xclone


## demo datasets
def load_TNBC1_RDR():
    """Return demo TNBC1_RDR adata (cells by genes).
    
    Example
    --------
    >>> import xclone
    >>> RDR_adata = xclone.pp.load_TNBC1_RDR()
    """
    # stream = pkg_resources.resource_stream(__name__, '../data/demo_datasets/TNBC1/TNBC1_RDR_adata_demo.h5ad')
    # RDR_adata = an.read_h5ad(stream)
    ## the saved one with no layers for saving storage.
    # RDR_adata.layers["raw_expr"] = RDR_adata.X.copy()
    warnings.warn(
        "`xclone.pp.load_TNBC1_RDR` is deprecated since XClone "
        "v0.3.0 and will be removed in a future version. Please use "
        "`xclone.data.tnbc1_rdr` instead.",
        DeprecationWarning,
        stacklevel=2,
    )

    return xclone.data.tnbc1_rdr()
    # return RDR_adata

def load_TNBC1_BAF():
    """Return demo TNBC1_BAF adata (cells by genes).

    Example
    -------
    >>> import xclone
    >>> BAF_adata = xclone.pp.load_TNBC1_BAF()
    """
    # stream = pkg_resources.resource_stream(__name__, '../data/demo_datasets/TNBC1/TNBC1_BAF_adata_demo.h5ad')
    # return an.read_h5ad(stream)

    warnings.warn(
        "`xclone.pp.load_TNBC1_BAF` is deprecated since XClone "
        "v0.3.0 and will be removed in a future version. Please use "
        "`xclone.data.tnbc1_baf` instead.",
        DeprecationWarning,
        stacklevel=2,
    )

    return xclone.data.tnbc1_baf()

def load_TNBC1_demo():
    """
    deprecated-2022-11
    Function: load the demo datasets for testing
    ------
    example:
    > import xclone
    > BAF_adata, RDR_adata = xclone.pp.load_TNBC1_demo()
    > RDR_adata
    AnnData object with n_obs × n_vars = 1097 × 33472
    var: 'GeneName', 'GeneID', 'chr', 'start', 'stop', 'arm', 'chr_arm', 'band'
    uns: 'data_mode', 'data_notes'
    > BAF_adata
    AnnData object with n_obs × n_vars = 1097 × 61775
    var: 'chr', 'start', 'stop', 'arm'
    uns: 'data_mode', 'data_notes'
    layers: 'AD', 'DP'
    # ###--------save the data in h5ad
    # > BAF_adata.write("./TNBC1_BAF_adata.h5ad")
    # > RDR_adata.write("./TNBC1_RDR_adata.h5ad")
    # ... storing 'chr' as categorical
    # ... storing 'arm' as categorical
    # ... storing 'chr_arm' as categorical
    # ... storing 'band' as categorical
    """
    mtx_barcodes_file = pkg_resources.resource_stream(__name__, '../data/demo_datasets/TNBC1/cellSNP.samples.tsv')
    mtx_barcodes = pd.read_table(mtx_barcodes_file, header = None, index_col=0)
    AD_file = pkg_resources.resource_stream(__name__, '../data/demo_datasets/TNBC1/kat-csp-post.50kb.block.AD.mtx')
    DP_file = pkg_resources.resource_stream(__name__, '../data/demo_datasets/TNBC1/kat-csp-post.50kb.block.DP.mtx')
    RDR_file = pkg_resources.resource_stream(__name__, '../data/demo_datasets/TNBC1/kat-rdr-feature-matrix.mtx')
    AD = sp.io.mmread(AD_file).tocsr()
    DP = sp.io.mmread(DP_file).tocsr()
    RDR = sp.io.mmread(RDR_file).tocsr()
    BAF_adata = xclonedata([AD,DP],'BAF', mtx_barcodes, "hg38_blocks", "[demo] TNBC1 scRNA-seq data in copyKAT")
    RDR_adata = xclonedata(RDR, 'RDR', mtx_barcodes, "hg38_genes", "[demo] TNBC1 scRNA-seq data in copyKAT")
    return BAF_adata,RDR_adata


