"""Base functions for XClone base anno data loading.
"""
# Author: Rongting Huang
# Date: 2021/07/15
# update: 2021/07/21

## Part I: load annotation data
## load package data

import pkg_resources
import pandas as pd

def load_anno(genome_mode):
    """Return a dataframe about the anno data
    for hg19/hg38.
    gene_based or block based
    # usage:
    # from xclone.utils import load_anno
    """
    # This is a stream-like object. If you want the actual info, call
    # stream.read()
    if genome_mode == "hg38_genes":
        stream = pkg_resources.resource_stream(__name__, '../data/anno_data/annotate_genes_hg38_update.txt')
    if genome_mode == "hg38_blocks":
        stream = pkg_resources.resource_stream(__name__, '../data/anno_data/annotate_blocks_hg38_update.txt')
    if genome_mode == "hg19_genes":
        stream = pkg_resources.resource_stream(__name__, '../data/anno_data/annotate_genes_hg19_update.txt')
    if genome_mode == "hg19_blocks":
        stream = pkg_resources.resource_stream(__name__, '../data/anno_data/annotate_blocks_hg19_update.txt')
    return pd.read_table(stream)#encoding='latin-1'


def load_hg38_genes():
    """Return a genes list of hg38.
    ## example usage of load_anno
    usage: from xclone.utils import load_hg38_genes
    load_hg38_genes()
    0         MIR1302-2HG
    1             FAM138A
    2               OR4F5
    3          AL627309.1
    4          AL627309.3
             ...     
    33467          TTTY4C
    33468         TTTY17C
    33469    LINC00266-4P
    33470            CDY1
    33471           TTTY3
    Name: GeneName, Length: 33472, dtype: object
    """
    return load_anno("hg38_genes")['GeneName']

def load_hg19_genes():
    """Return a genes list of hg19.
    ## example usage of load_anno
    usage: from xclone.utils import load_hg19_genes
    load_hg19_genes()
    """
    return load_anno("hg19_genes")['GeneName']

def load_cc_genes():
    """Return a list of cell cycle genes.
    """
    stream = pkg_resources.resource_stream(__name__, '../data/anno_data/cellcycle_genes.txt')
    return pd.read_table(stream)


def load_hk_genes():
    """Return a list of house keeping genes.
    """
    stream = pkg_resources.resource_stream(__name__, '../data/anno_data/housekeeping_genes.txt')
    return pd.read_table(stream)