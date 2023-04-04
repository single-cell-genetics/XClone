"""Builtin datasets loadding.
"""
from pathlib import Path
from typing import Optional, Union

import numpy as np
import pandas as pd
# import anndata as ad
from scanpy import read
# from scanpy import read, read_10x_h5, read_10x_mtx, write, read_visium

## global url_datadir
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

    gene mode

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

def bch869_baf_50kb(file_path: Union[str, Path] = "data/BCH869/BCH869_BAF_adata_50kb.h5ad"):
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

    50kb block mode

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
    >>> BAF_adata = xclone.data.bch869_baf_50kb()
    """
    url = f"{url_datadir}BCH869_scRNA/BCH869_BAF_adata_50kb.h5ad"
    adata = read(file_path, backup_url=url)

    return adata

from ...analysis.extract import dir_make
import requests
import zipfile
import os

def download_unzip_file(url, file_path, extract_path):
    """
    Example:
    --------
    url = "http://example.com/example.zip"
    file_path = "/path/to/save/example.zip"
    extract_path = "/path/to/extract"
    download_unzip_file(url, file_path, extract_path)
    """
    # Download the ZIP file from the URL
    response = requests.get(url)

    # Save the ZIP file to the specified file path
    with open(file_path, "wb") as f:
        f.write(response.content)

    # Extract the contents of the ZIP file to the specified extract path
    with zipfile.ZipFile(file_path, "r") as zip_ref:
        zip_ref.extractall(extract_path)

def gx109_rdr(file_path: Union[str, Path] = "/data/GX109-T1c/GX109-T1c_RDR_adata.zip",
              extract_path: Union[str, Path] = "/data/GX109-T1c/",
              raw_expr_layer: bool = True):
    """Gastric cancer tissue.
    scRNA-seq Source: GEOXXXX (will provide shortly)
    Platform: 5' 10x Genomics (scRNA-seq)
    Number of cells: 5400
    Genome version: refdata-cellranger-GRCh38-3.0.0
    cellranger: (v2.2.0)
    Preprocessd: xcltk v0.1.15 RDR pipeline.
    Notes: Released by XClone paper, there is a curated annotation 
    that contains 5400 cells from 4 different cell types.

    Arguments
    ---------
    file_path
        Path where to save dataset and read it from.
    extract_path
        Path to extract zip files to h5ad file.
    raw_expr_layer
        bool, add raw_expr layer or not.
    
    Returns
    -------
    Returns `adata` object
    
    Example
    -------
    >>> import xclone
    >>> RDR_adata = xclone.data.gx109_rdr()
    """
    h5ad_file = f"./data/GX109-T1c/GX109-T1c_RDR_adata.h5ad"
    
    if os.path.exists(h5ad_file):
        print("load the GX109 rdr data from downloaded file.")
    else:
        print("load the GX109 rdr data from xclone-data github.")
        url = f"{url_datadir}GX109-T1c_scRNA/GX109-T1c_RDR_adata.zip"
        cwd = os.getcwd()
        extract_path = f'{cwd}{extract_path}'
        if os.path.exists(extract_path):
            pass
        else:
            dir_make(extract_path)
        file_path = f'{cwd}{file_path}'
        download_unzip_file(url, file_path, extract_path)

    adata = read(h5ad_file)
    if raw_expr_layer:
        adata.layers["raw_expr"] = adata.X.copy()
    return adata