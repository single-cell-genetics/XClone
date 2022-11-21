"""Builtin datasets loadding.
"""
from pathlib import Path
from typing import Optional, Union

import numpy as np
import pandas as pd
# import anndata as ad

from scanpy import read
# from scanpy import read, read_10x_h5, read_10x_mtx, write, read_visium

url_datadir = "https://github.com/Rongtingting/xclone-data/raw/main/"

def tnbc1_rdr(file_path: Union[str, Path] = "data/TNBC1/TNBC1_RDR_adata_demo.h5ad"):
    """Triple-negative breast cancer.
    From `Gao R et al. (2021) <https://doi.org/10.1038%2Fs41587-020-00795-2>`__.
    scRNA-seq Source: GSE148673
    Platform: 10x scRNA
    Number of cells: 1097
    Genome version: refdata-cellranger-GRCh38-3.0.0
    Preprocessd: xcltk v0.1.15 RDR pipeline.

    Arguments
    ---------
    file_path
        Path where to save dataset and read it from.
    
    Returns
    -------
    Returns `adata` object
    
    Example
    -------
    >>> import xclone
    >>> RDR_adata = xclone.data.tnbc1_rdr()
    """

    url = f"{url_datadir}TNBC1_scRNA/TNBC1_RDR_adata_demo.h5ad"
    adata = read(file_path, backup_url=url)
    adata.layers["raw_expr"] = adata.X.copy()

    return adata

def tnbc1_baf(file_path: Union[str, Path] = "data/TNBC1/TNBC1_BAF_adata.h5ad"):
    """Triple-negative breast cancer.
    From `Gao R et al. (2021) <https://doi.org/10.1038%2Fs41587-020-00795-2>`__.
    scRNA-seq Source: GSE148673
    Platform: 10x scRNA
    Number of cells: 1097
    Genome version: refdata-cellranger-GRCh38-3.0.0
    Preprocessd: xcltk v0.1.15 BAF pipeline.

    Arguments
    ---------
    file_path
        Path where to save dataset and read it from.
    
    Returns
    -------
    Returns `adata` object
    
    Example
    -------
    >>> import xclone
    >>> BAF_adata = xclone.data.tnbc1_baf()
    """
    url = f"{url_datadir}TNBC1_scRNA/TNBC1_BAF_adata.h5ad"
    adata = read(file_path, backup_url=url)

    return adata

def bch869_rdr(file_path: Union[str, Path] = "data/BCH869/BCH869_RDR_adata_demo.h5ad"):
    """A glioma sample BCH869 with histone H3 lysine27-to-methionine mutations (H3K27M-glioma), 
    where 489 malignant cells and 3 non-tumour cells were probed by smart-seq2.
    From `Filbin MG et al. (2018) <https://doi.org/10.1126%2Fscience.aao4750>`__.
    scRNA-seq Source: GSE102130
    Platform: SMART-seq2
    Number of cells: 960 in total; 489 in detected clones
    annotation: Single Cell Portal Study: single-cell analysis in pediatric 
    midline gliomas with histone H3K27M mutation
    Genome version: refdata-cellranger-hg19-3.0.0
    Preprocessd: xcltk v0.1.15 RDR pipeline.

    Arguments
    ---------
    file_path
        Path where to save dataset and read it from.
    
    Returns
    -------
    Returns `adata` object

    Example
    -------
    >>> import xclone
    >>> RDR_adata = xclone.data.bch869_rdr()
    """

    url = f"{url_datadir}BCH869_scRNA/BCH869_RDR_adata_demo.h5ad"
    adata = read(file_path, backup_url=url)
    adata.layers["raw_expr"] = adata.X.copy()

    return adata

def bch869_baf(file_path: Union[str, Path] = "data/BCH869/BCH869_BAF_adata.h5ad"):
    """A glioma sample BCH869 with histone H3 lysine27-to-methionine mutations (H3K27M-glioma), 
    where 489 malignant cells and 3 non-tumour cells were probed by smart-seq2.
    From `Filbin MG et al. (2018) <https://doi.org/10.1126%2Fscience.aao4750>`__.
    scRNA-seq Source: GSE102130
    Platform: SMART-seq2
    Number of cells: 960 in total; 489 in detected clones
    annotation: Single Cell Portal Study: single-cell analysis in pediatric 
    midline gliomas with histone H3K27M mutation
    Genome version: refdata-cellranger-hg19-3.0.0
    Preprocessd: xcltk v0.1.15 BAF pipeline.

    Arguments
    ---------
    file_path
        Path where to save dataset and read it from.
    
    Returns
    -------
    Returns `adata` object

    Example
    -------
    >>> import xclones
    >>> BAF_adata = xclone.data.bch869_baf()
    """
    url = f"{url_datadir}BCH869_scRNA/BCH869_BAF_adata.h5ad"
    adata = read(file_path, backup_url=url)

    return adata