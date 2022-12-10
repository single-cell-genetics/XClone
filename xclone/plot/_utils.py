"""Utils functions for XClone data plotting data
include some relative analysis
"""
# Author: Rongting Huang
# Date: 2021/07/24
# update: 2021/07/24

import numpy as np

def base_heatmap(
    Xdata: AnnData,
    **kwargs):
    """
    Function: base Heatmap for XClone confusion matrix visualization.
    maybe not for anndata
    """
    pass


## For extracting BAF cell cluster 
def get_cluster_barcodes_lst(cluster_id, cluster_index_lst, mtx_barcodes_lst):
    """
    Function:
    
    Example:
    uni_cluster_barcodes = get_cluster_barcodes_lst(1, cluster_index_lst, mtx_barcodes_lst)
    """
    if len(cluster_index_lst) != len(mtx_barcodes_lst):
        print("error! No matching in cluster_index and mtx_barcodes")
    flag_ = cluster_index_lst ==  cluster_id
    barcodes_ = mtx_barcodes_lst[flag_]
    return list(barcodes_[0])


def identified_node_lst(mtx_barcodes_lst, barcodes_node_lst, node_ID="node_10050"):
    """
    Function: generate list for confusion matrix building
    mtx_barcodes_lst: ref list with the same order as the matrix used
    barcodes_node_lst: barcodes for specified node, like 10050
    
    Example:
    baf_dat_dir = '/storage/yhhuang/users/rthuang/processed_data/xcltk/xianjie-cpos/mkn45_021021/mkn45-fulldepth/mkn45-fulldepth-csp-post/phase-snp-even/'
    mtx_barcodes_file = baf_dat_dir + "cellSNP.samples.tsv"
    mtx_barcodes_lst = pd.read_table(mtx_barcodes_file,header = None)

    mkn45_10x_dloupe = "/storage/yhhuang/research/mito/mkn45/fulldepth/mode2/mkn45_5k_dloupe-group_10438-region_1_1-5120000_Y_56320001-59373566-heatmap_copy number.csv"
    barcodes_10050_lst = get_node_barcodeslst(mkn45_10x_dloupe, 10050)
    node_lst = identified_node_lst(mtx_barcodes_lst, barcodes_10050_lst, node_ID="node_10050")
    """
    bool_node_lst = mtx_barcodes_lst.isin(barcodes_node_lst)
    print("node_cells_num",bool_node_lst.sum())
    node_lst_ = bool_node_lst.applymap(lambda x: node_ID if x == True else "non_" + node_ID)[0].values
    return node_lst_

## For confusion heatmap
### from vireoSNP
### https://github.com/single-cell-genetics/vireo/blob/master/vireoSNP/utils/base_utils.py#L3
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

def get_confuse_mat_df(confuse_mat,index_names=None,clolums_names=None):
    """
    Function: construct df for confuse_mat for heatmap plot
    ========================
    Example:
    confuse_mat_df = get_confuse_mat_df(confuse_mat,index_names=["cluster0", "cluster1"],clolums_names=["node_10050", "non node_10050"])
    """
    confuse_mat_df = pd.DataFrame(confuse_mat)
    if index_names:
        confuse_mat_df.index = index_names
    if clolums_names:
        confuse_mat_df.columns = clolums_names
    return confuse_mat_df
    

## comparison analysis-confusion matrix and visualization
# import pandas as pd
# import numpy as np
def compare_clusters(model_proba, mtx_barcodes_file, node_barcodes_lst):
    """
    Function:
    
    Example:
    #### data
    filename = "mkn45-500k-BAF-chr2/_model_baf_sub2_6.model"
    _model_baf_sub2_6 = joblib.load(filename)
    
    baf_dat_dir = '/storage/yhhuang/users/rthuang/processed_data/xcltk/xianjie-cpos/mkn45_020821/mkn45-500k/mkn45-500k-csp-post/phase-snp-even/'
    mtx_barcodes_file = baf_dat_dir + "cellSNP.samples.tsv"
    
    mkn45_10x_dloupe = "/storage/yhhuang/research/mito/mkn45/fulldepth/mode2/mkn45_5k_dloupe-group_10438-region_1_1-5120000_Y_56320001-59373566-heatmap_copy number.csv"
    barcodes_8805_lst = get_node_barcodeslst(mkn45_10x_dloupe, 8805)
    
    ## usage
    compare_clusters(_model_baf_sub2_6.proba_, mtx_barcodes_file, barcodes_8805_lst)
    
    """
    cluster_index_lst = np.argmax(model_proba, axis=1)
    uni_cluster = np.unique(cluster_index_lst)
    mtx_barcodes_lst = pd.read_table(mtx_barcodes_file,header = None)
#     print(mtx_barcodes_lst)
    for uni_id in uni_cluster:
        print("cluster",uni_id, ":", np.sum(cluster_index_lst == uni_id))
        uni_cluster_barcodes = get_cluster_barcodes_lst(uni_id, cluster_index_lst, mtx_barcodes_lst)
#         print(uni_cluster_barcodes)
        inter_set = (set(uni_cluster_barcodes)&set(node_barcodes_lst))
#         print(inter_set)
        print("inter_set", len(inter_set))




# def confuse_heatmap(confuse_mat_df, cmap = "Blues", version1 = True, version2 =True, save_pdf = None):
#     """
#     Function:
#     confusion matrix heatmap
    
#     confuse_mat_df: count- confusion matrix in pd.DataFrame
    
#     Example：
#     confuse_heatmap(confuse_mat_df, version1 = False, save_pdf = True)
#     """
#     import matplotlib.pyplot as plt
#     import seaborn as sns
#     ## setting
#     font = {'family' : 'DejaVu Sans',
#             'size'   : 15}
#     plt.rc('font', **font)
#     if version1:
#         sns.heatmap(confuse_mat_df, annot=True, fmt="d", cmap=cmap)
#     if version2:
#         ## data normalization
#         plot_df = (confuse_mat_df.T/confuse_mat_df.sum(axis=1).values).T
#         ## plot percentage
#         plt.figure(figsize=[8,6])
#         ax = sns.heatmap(plot_df, cmap=cmap) # annot=True,
#         ## annot original count
#         height, width = np.array(confuse_mat_df).shape
#         text_colors = ['black', 'white']
#         for x in range(width):
#             for y in range(height):
#                 ax.annotate(str(np.array(confuse_mat_df)[y][x]), xy=(x+0.5, y+0.5), 
#                             horizontalalignment='center',
#                             verticalalignment='center', color=text_colors[int(np.array(plot_df)[y][x] > 0.5)], fontsize = 15)
#         plt.xlabel('CellRanger nodes')
#         plt.ylabel('CNV clones')
#         # plt.yticks(rotation=45)
#         # plt.xticks(rotation=45)
#         # plt.xticks(range(3), confusion_matrix.columns) #, rotation=315)
#         # plt.yticks(range(len(norm_conf)), set(clone_id))
#         plt.tight_layout()
#         plt.title('Confusion matrix for comparing CNV clones with node_10050', fontsize = 18)
#         if save_pdf:
#             plt.savefig('conf_mat_CNV_2clones-node10050_update.pdf', dpi=200, bbox_inches='tight')

## confusion matrix heatmap
def confuse_heatmap(confuse_mat_df, cmap = "Blues", version1 = True, version2 =True, plot_xlabel = None,
                    plot_ylabel = None, plt_title = None, save_file_name = None):
    """
    Function:
    ---------
    confusion matrix heatmap
    confuse_mat_df: count- confusion matrix in pd.DataFrame
    
    Example:
    --------           
    >>> confuse_mat, ids1_uniq, ids2_uniq = get_confusion(CopyKAT_lst, expression_lst)
    >>> confuse_mat_df = get_confuse_mat_df(confuse_mat,
                                        index_names=["cloneA", "cloneB","Normal"],
                                        clolums_names=["cloneA", "cloneB","Normal"])
    >>> confuse_heatmap(confuse_mat_df,  
                    plot_xlabel = "CopyKAT",
                    plot_ylabel = "Ground Truth", 
                    plt_title = "Concordance in subclone identification", 
                    save_file_name = 'CopyKAT_vs_Groudtruth1.pdf')
    """
    import matplotlib.pyplot as plt
    import seaborn as sns
    ## setting
    font = {'family' : 'DejaVu Sans',
            'size'   : 15}
    plt.rc('font', **font)
    if version1:
        plt.figure(figsize=[8,6])
        sns.heatmap(confuse_mat_df, annot=True, fmt="d", cmap=cmap)
    if version2:
        ## data normalization
        plot_df = (confuse_mat_df.T/confuse_mat_df.sum(axis=1).values).T
        ## plot percentage
        plt.figure(figsize=[8,6])
        ax = sns.heatmap(plot_df, cmap=cmap) # annot=True,
        ## annot original count
        height, width = np.array(confuse_mat_df).shape
        text_colors = ['black', 'white']
        for x in range(width):
            for y in range(height):
                ax.annotate(str(np.array(confuse_mat_df)[y][x]), xy=(x+0.5, y+0.5), 
                            horizontalalignment='center',
                            verticalalignment='center', color=text_colors[int(np.array(plot_df)[y][x] > 0.5)], fontsize = 15)
        if plot_xlabel is None:
            pass
        else:
            plt.xlabel(plot_xlabel)
        if plot_ylabel is None:
            pass
        else:
            plt.ylabel(plot_ylabel)
        # plt.yticks(rotation=45)
        # plt.xticks(rotation=45)
        # plt.xticks(range(3), confusion_matrix.columns) #, rotation=315)
        # plt.yticks(range(len(norm_conf)), set(clone_id))
        plt.tight_layout()
        if plt.title is None:
            pass
        else:
            plt.title(plt_title, fontsize = 18) 
            #plt.title('Concordance in subclone identification', fontsize = 18)
        if save_file_name is None:
            pass
        else:
            plt.savefig(save_file_name, dpi=300, bbox_inches='tight')
              
