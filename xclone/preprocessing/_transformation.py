"""Base functions for XClone data transformation preprocessing.
init for smart-seq scRNA data.
"""
# Author: Rongting Huang
# Date: 2022/09/28
# update: 2022/09/28


from scipy import sparse
import numpy as np

## scRNA-seq data; smart-seq
def Xtransformation(Xdata, transform = True, Xlayers = ["raw_expr"]):
    """
    Expression transformation for smart-seq scRNA count.
    Xlayers for transformation.
    """
    if transform:
        Xdata.X = sparse.csr_matrix(np.round(np.log(Xdata.X.A + 1)))
        for layer_ in Xlayers:
            Xdata.layers[layer_] = Xdata.X
    return Xdata

## scATAC-seq data 