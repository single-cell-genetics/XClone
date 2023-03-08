"""Base functions for XClone performance evaluation.
"""
import numpy as np
import pandas as pd
import anndata as ad

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


def extract_Xdata1(Xdata, GT_df):
    """
    Notes:
    update for BCH869 dataset.
    
    Example:
    --------
    >>> data_dir = "..."
    >>> use_cell_file = dat_dir + "GX109.isec.cell.anno.tsv"
    >>> use_gene_file = dat_dir + "GX109.isec.gene.anno.tsv"
    >>> RDR_adata = extract_Xdata(RDR_adata,use_cell_file, use_gene_file)
    """


    use_cell_lst = GT_df.index.tolist()
    use_gene_lst = GT_df.columns.tolist()
    
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
    develop base on GX109 example.
    
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

def load_ground_truth(file_path,  out_df = True, BCHdataset = True):
    """
    support for load benchmarking ground truth matrix from .rds file.
    develop base on BCH869 example.
    """
    import pyreadr
    result = pyreadr.read_r(file_path) 
    GT_df = result[None]
    if BCHdataset:
        GT_df.index = GT_df.index.str.replace("BCH", "BT_", regex = False).str.replace(".", "-", regex = False)
    Ground_truth_mtx = GT_df.values
    print("Loading Ground_truth_mtx, shape (cells, genes): ", Ground_truth_mtx.shape)
    if out_df:
        return GT_df
    else:
        return Ground_truth_mtx

def resort_Ground_truth_mtx(GT_df, Xdata, verbose = True):
    """
    resort Ground_truth_mtx by Xdata order.
    develop base on BCH869 example.
    """
    cell_barcodes = pd.DataFrame(Xdata.obs.index, columns = ["cell_barcodes"])
    gene_names = pd.DataFrame(Xdata.var["GeneName"])

    GT_df = gene_names.merge(GT_df.T, left_on = "GeneName", right_index=True, how = "left")
    GT_df.set_index(keys = "GeneName", drop = True, inplace = True)
    if verbose:
        print("Ground_truth_mtx reorder genes")
    GT_df = cell_barcodes.merge(GT_df.T, left_on = "cell_barcodes", right_index=True, how = "left")
    GT_df.set_index(keys = "cell_barcodes", drop = True, inplace = True)
    if verbose:
        print("Ground_truth_mtx reorder cells")
    Ground_truth_mtx = GT_df.values
    if verbose:
        print("Ground_truth_mtx shape (cells, genes): ", Ground_truth_mtx.shape)
    return Ground_truth_mtx


def base_evaluate_map(adata, ground_truth_file):
    """
    """
    GT_df = load_ground_truth(ground_truth_file,  out_df = True, BCHdataset = True)
    Xdata = extract_Xdata1(adata, GT_df)
    Ground_truth_mtx = resort_Ground_truth_mtx(GT_df, Xdata)
    
    return Ground_truth_mtx, Xdata

def roc_auc_evaluate(Ground_truth_mtx, Xdata, Xlayer = "XC_denoise", state_idx = 0, ROC_show = False):
    """
    ref: # https://github.com/huangyh09/foundation-data-science/blob/main/w10-classification/P3_ROC_NaiveBayes.ipynb
    
    state_idx = 0  copy_loss
    state_idx = 1 loh
    state_idx = 3 copy_gain
    """
   
    import numpy as np
    from sklearn import metrics
    fpr, tpr, thresholds = metrics.roc_curve(Ground_truth_mtx.flatten(), Xdata.layers[Xlayer][:,:,state_idx].flatten())
    roc_auc = metrics.roc_auc_score(Ground_truth_mtx.flatten(), Xdata.layers[Xlayer][:,:,state_idx].flatten())
    print("AUC value: ", roc_auc)
    
    if ROC_show:
        _query_threshold1 = 0.2
        idx1 = np.argmin(np.abs(thresholds - _query_threshold1))
        _fpr1, _tpr1, _threshold1 = fpr[idx1], tpr[idx1], thresholds[idx1]

        _query_threshold2 = 0.5
        idx2 = np.argmin(np.abs(thresholds - _query_threshold2))
        _fpr2, _tpr2, _threshold2 = fpr[idx2], tpr[idx2], thresholds[idx2]

        _query_threshold3 = 0.7
        idx3 = np.argmin(np.abs(thresholds - _query_threshold3))
        _fpr3, _tpr3, _threshold3 = fpr[idx3], tpr[idx3], thresholds[idx3]


        import matplotlib.pylab as plt
        fig = plt.figure(dpi=150)
        plt.plot(fpr, tpr, label='AUC=%.4f' %(roc_auc))

        plt.plot(_fpr1, _tpr1, 'o', label='threshold: %.2f' %(_threshold1))
        plt.plot(_fpr2, _tpr2, 'o', label='threshold: %.2f' %(_threshold2))
        plt.plot(_fpr3, _tpr3, 'o', label='threshold: %.2f' %(_threshold3))

        plt.plot([0, 1], [0, 1], '--', color='grey')
        plt.plot([0, 0, 1], [0, 1, 1], '--', color='purple')

        plt.grid(alpha=0.4)

        plt.xlabel("False positive rate (1 - specificity)")
        plt.ylabel("True positive rate (sensitivity)")
        plt.title("ROC curve")
        plt.legend()
        plt.show()

    return roc_auc

def fast_evaluate(GT_dat_dir, Combine_adata, Xlayer = "XC_denoise", dataset_name = "BCH869", update = False):
    AUC_record = {}
    ## copy gain
    state_name = "copy_gain"
    ground_truth_file = f"{GT_dat_dir}{state_name}/truth/{dataset_name}.isec.copy.gain.binary.matrix.rds"
    if update:
        ground_truth_file = f"{GT_dat_dir}{state_name}/result/s4_truth/{dataset_name}.{state_name}.gene_scale.truth.cell_x_gene.binary.mtx.rds"

    Ground_truth_mtx, Xdata = base_evaluate_map(Combine_adata, ground_truth_file)
    roc_auc = roc_auc_evaluate(Ground_truth_mtx, Xdata, Xlayer, state_idx = 3, ROC_show = True)

    AUC_record[state_name] = roc_auc
    
    ## copy loss
    state_name = "copy_loss"
    ground_truth_file = f"{GT_dat_dir}{state_name}/truth/{dataset_name}.isec.copy.loss.binary.matrix.rds"
    if update:
        ground_truth_file = f"{GT_dat_dir}{state_name}/result/s4_truth/{dataset_name}.{state_name}.gene_scale.truth.cell_x_gene.binary.mtx.rds"

    Ground_truth_mtx, Xdata = base_evaluate_map(Combine_adata, ground_truth_file)
    roc_auc = roc_auc_evaluate(Ground_truth_mtx, Xdata, Xlayer, state_idx = 0, ROC_show = True)

    AUC_record[state_name] = roc_auc
    
    ## loh
    state_name = "loh"
    ground_truth_file = f"{GT_dat_dir}{state_name}/truth/{dataset_name}.isec.loh.binary.matrix.rds"
    if update:
        ground_truth_file = f"{GT_dat_dir}{state_name}/result/s4_truth/{dataset_name}.{state_name}.gene_scale.truth.cell_x_gene.binary.mtx.rds"

    Ground_truth_mtx, Xdata = base_evaluate_map(Combine_adata, ground_truth_file)
    roc_auc = roc_auc_evaluate(Ground_truth_mtx, Xdata, Xlayer, state_idx = 1, ROC_show = True)

    AUC_record[state_name] = roc_auc

    print("Xlayer used: ", Xlayer)
    print("Results: ", AUC_record)

    return AUC_record
    

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

def change_res_resolution(Xdata, Xlayer, outlayer = "prob_merge", brk = "chr_arm"):
    """
    change result (probability) resolution, e.g., from gene to chr_arms or chrs.
    For methods benchmarking.
    """
    merge_results = []
    
    brk_item = Xdata.var[brk].drop_duplicates(keep="first")
    for brk_ in brk_item:
        flag_= Xdata.var[brk] == brk_
        tmp_res = np.median(Xdata.layers[Xlayer][:,flag_,:], axis = 1)
        
        ## normalized merged probability
        # from sklearn.preprocessing import normalize
        # tmp_res = normalize(tmp_res, norm="l1")
        tmp_res = tmp_res/tmp_res.sum(axis = -1, keepdims = True)
        
        merge_results.append(np.expand_dims(tmp_res, axis = 1))
    res_ = np.hstack(merge_results)
    
    # ## check res_
    # print((abs(res_.sum(axis=-1) - 1) > 1e-15).sum())
    
    ## generate merge_Xdata var
    merge_var1 = Xdata.var.copy()
    merge_var1.drop_duplicates(brk, keep="first", inplace=True)
    # keep_cols = ["chr", "start", "stop", "arm", "chr_arm", "band"]
    keep_cols = ["chr", "start", "stop", "arm", "chr_arm"]
    merge_var1 = merge_var1[keep_cols]
    merge_var1.rename(columns={'stop': 'gene1_stop'}, inplace=True)

    merge_var2 = Xdata.var.copy()
    merge_var2.drop_duplicates(brk, keep="last", inplace=True)
    # keep_cols = ["chr", "start", "stop", "arm", "chr_arm", "band"]
    keep_cols = ["chr", "start", "stop", "arm", "chr_arm"]
    merge_var2 = merge_var2[keep_cols]
    merge_var2.rename(columns={'start': 'lastgene_start'}, inplace=True)
    
    merge_var = merge_var1.merge(merge_var2, left_on=brk, right_on=brk, suffixes=('_start', '_stop'))
    merge_var.drop(columns=['chr_stop', 'arm_stop'], inplace  = True)
    merge_var.rename(columns={'chr_start': 'chr', "arm_start": "arm"}, inplace=True)
    ## generate new merged Xdata with the merge res in new layer
    x_array = np.zeros((res_.shape[0], res_.shape[1]))
    merge_Xdata = ad.AnnData(x_array, var = merge_var, obs = Xdata.obs.copy())
    merge_Xdata.layers[outlayer] = res_
    
    return merge_Xdata