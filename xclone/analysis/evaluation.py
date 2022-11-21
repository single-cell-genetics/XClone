"""Base functions for XClone performance evaluation.
"""
import numpy as np
import pandas as pd

def extract_Xdata(Xdata, use_cell_file, use_gene_file):
    """
    Notes:
    use_cell_file, use_gene_file from xianjie
    
    Example:
    --------
    >>> data_dir = "..."
    >>> use_cell_file = dat_dir + "GX109.isec.cell.anno.tsv"
    >>> use_gene_file = dat_dir + "GX109.isec.gene.anno.tsv"
    >>> RDR_adata = extract_Xdata(RDR_adata,use_cell_file, use_gene_file)
    """
    
    use_cell = pd.read_csv(use_cell_file, sep = "\t")
    use_gene = pd.read_csv(use_gene_file, sep = "\t")
    
    use_cell_lst = use_cell["cell"].tolist()
    use_gene_lst = use_gene["Gene"].tolist()
    
    cell_flag = Xdata.obs.index.isin(use_cell_lst)
    gene_flag = Xdata.var["GeneName"].isin(use_gene_lst)
    
    update_Xdata = Xdata.copy()
    update_Xdata = update_Xdata[cell_flag,:].copy()
    update_Xdata = update_Xdata[:, gene_flag].copy()
    
    return update_Xdata

def Ground_truth_mtx_generate(Xdata, genes, cells_type = ['Paneth', 'TA']):
    """
    Function:
    --------
    Generate binary matrix: region with CNV state is 1, otherwise 0.
    
    Example:
    Generate ground_truth_matrix for copy loss evaluation
    --------
    >>> cells_type = ["Paneth", "TA", "Chief", "M"]
    >>> copyloss_genes = pd.read_csv(loss_truth_file, sep = "\t")
    >>> Groud_truth_mtx = Groud_truth_mtx_generate(RDR_adata, copyloss_genes, cells_type)
    """
    
    Ground_truth_mtx = np.zeros(Xdata.shape)
    
    cell_flag = Xdata.obs["cell_type"].isin(cells_type)
    gene_flag = Xdata.var["GeneName"].isin(genes["Gene"])
    
    gene_index = np.flatnonzero(gene_flag.values)
    cell_index = np.flatnonzero(cell_flag.values)
    
    Ground_truth_mtx[np.repeat(cell_index,len(gene_index)),  np.tile(gene_index, len(cell_index))] = 1
    
    return Ground_truth_mtx


def get_confusion(ids1, ids2):
    """Get confusion matrix
    
    Parameters
    ----------
    ids1: numpy.array or list
        id list in the first annotation
    ids2: numpy.array or list
        id list in the second annotation
        
    Return
    ------
    (confuse_mat, ids1_uniq, ids2_uniq)
    confuse_mat[i, j]: 
        number of samples have ids1 == ids1_uniq[i]
        and ids2 == id2_uniq[j]
    """
    if type(ids1) == list: ids1 = np.array(ids1)
    if type(ids2) == list: ids2 = np.array(ids2)
    
    ids1_uniq = np.unique(ids1)
    ids2_uniq = np.unique(ids2)
    
    confuse_mat = np.zeros((len(ids1_uniq), len(ids2_uniq)), dtype=int)
    for i, _id1 in enumerate(ids1_uniq):
        for j, _id2 in enumerate(ids2_uniq):
            confuse_mat[i, j] = np.sum((ids1 == _id1) * (ids2 == _id2))
            
    return confuse_mat, ids1_uniq, ids2_uniq

def get_confuse_mat_df(confuse_mat,index_names, columns_names):
    import pandas as pd
    df = pd.DataFrame(confuse_mat, index = index_names, columns = columns_names)
    return df