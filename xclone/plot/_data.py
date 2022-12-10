"""Base functions for XClone RDR plotting
data preprocessing for **visualization** and *analysis*.
"""

# Author: Rongting Huang
# Date: 2021/12/17
# update: 2022/12/10

## obs[‘viz_index’] = index ## Todo1: improve efficiency.

import anndata as ad

def reorder_data_by_cellanno(Xdata, cell_anno_key):
    """
    Function:
    order the anndata by cell_annotation(e.g., celltype) for visualization.
    https://anndata.readthedocs.io/en/latest/concatenation.html
    
    Parameters:
    ----------
    Xdata: anndata.
    cell_anno_key: char.
    
    Return:
    ------
    ordered_Xdata: anndata
    """
    groups = Xdata.obs.groupby(cell_anno_key).indices # dict, index can be saved in Xdata.uns
    ordered_Xdata = ad.concat([Xdata[inds] for inds in groups.values()], merge="same")
    return ordered_Xdata


def transform_anndata(raw_ratio, new_var, new_obs):
    """
    transform the log raw ratio to anndata for visualization.
    """
    rr_ad = ad.AnnData(raw_ratio, var=new_var, obs = new_obs)
    return rr_ad
