"""Utils functions for XClone data preprocessing/analysis
"""
# Author: Rongting Huang
# Date: 2021/07/22
# update: 2021/07/22

import pandas as pd
import matplotlib.pyplot as plt

def get_node_barcodeslst(dloupe_node_file, node_id):
    """
    Function: get the barcodes list for specific node (cell cluster)
    
    dloupe_node_file: path of the file, defalut for dloupe output csv file,
    which contains node_id, barcodes, num_cells, num_noisy and blabla...
    node_id: int, the ID of the node need to be extrcted.
    
    example:
    mkn45_10x_dloupe = "/storage/yhhuang/research/mito/mkn45/fulldepth/mode2/
    mkn45_5k_dloupe-group_10438-region_1_1-5120000_Y_56320001-59373566-heatmap_copy number.csv"
    barcodes_10050_lst = get_node_barcodeslst(mkn45_10x_dloupe, 10050)
    """
    data_10x_cluster = pd.read_table(dloupe_node_file,sep=",")
    node_data = data_10x_cluster[data_10x_cluster["node_id"] == node_id]
    node_barcodes = node_data["barcodes"].str.split(";")
    node_barcodes_lst = list(node_barcodes)[0]
    print("node:", node_id,"num_cells:",len(node_barcodes_lst))
    return node_barcodes_lst

def get_non_node_barcodeslst(dloupe_node_file, node_id):
    """
    Function: get the barcodes list for excluded node (cell cluster)
    specific for CellRanger dloupe output csv file.
    
    dloupe_node_file: path of the file, defalut for dloupe output csv file,
    which contains node_id, barcodes, num_cells, num_noisy and blabla...
    node_id: int, the ID of the excluded node need to be extrcted.
    
    example:
    mkn45_10x_dloupe = "/storage/yhhuang/research/mito/mkn45/fulldepth/mode2/
    mkn45_5k_dloupe-group_10438-region_1_1-5120000_Y_56320001-59373566-heatmap_copy number.csv"
    barcodes_non_node_lst = get_non_node_barcodeslst(mkn45_10x_dloupe, 10050)
    """
    data_10x_cluster = pd.read_table(dloupe_node_file,sep=",")
    Flag_ = data_10x_cluster["node_id"] != node_id
    node_data = data_10x_cluster[Flag_]
    node_barcodes = node_data["barcodes"].str.split(";")
    node_barcodes_lst = []
    for idx_ in node_data.index:
        for i in range(len(node_barcodes[idx_])):
            node_barcodes_lst.append(node_barcodes[idx_][i])     
    print("non_node:", node_id,"num_cells:",len(node_barcodes_lst))
    return node_barcodes_lst

def processing_dlouple_file(dloupe_node_file, save_status = False, save_file = "xxxx_scDNA_cellranger_node.txt"):
    """
    Function: get the barcodes list from the cellranger scDNA dloupe heatmap csv file
    for each node (cell cluster)
    
    dloupe_node_file: path of the file, defalut for dloupe output csv file,
    which contains node_id, barcodes, num_cells, num_noisy and blabla...
    save_status: default False for saving the result for future processing
    
    return:
    barcodes and corresponding node_id, order is the same with the heatmap in loupe scDNA Browser.
    
    example:
    mkn45_10x_dloupe = "/storage/yhhuang/research/mito/mkn45/fulldepth/mode2/
    mkn45_5k_dloupe-group_10438-region_1_1-5120000_Y_56320001-59373566-heatmap_copy number.csv"
    cellranger_resutls = processing_dlouple_file(mkn45_10x_dloupe, 10050)
    """
    data_10x = pd.read_table(dloupe_node_file, sep=",")
    data_10x["node_barcodes_lst"] = data_10x["barcodes"].str.split(";")
    node_barcodes_lst = []
    node_lst = []
    for idx_ in data_10x.index:
        node_barcodes_lst_tmp = data_10x["node_barcodes_lst"][idx_]
        node_tmp = data_10x["node_id"][idx_]
        for i in range(len(node_barcodes_lst_tmp)):
            node_barcodes_lst.append(node_barcodes_lst_tmp[i])
            node_lst.append(node_tmp)
    cellranger_data = {'barcodes':node_barcodes_lst,
                     'node_id': node_lst}
    cellranger_results = pd.DataFrame(cellranger_data)
    if save_status:
        cellranger_results.to_csv(save_file, sep="\t", index=False)
    return cellranger_results

def get_barcodeslst(cellranger_results_file, node_id_lst, include=True, return_index=False):
    """
    Function:
    cellranger_results_file: processed and saved cellbarcodes and corresponding node_id
    
    usage:
    node_barcodes_lst = get_barcodeslst(cellranger_results_file,[2710,2701],include=False)
    """
    cellranger_results = pd.read_csv(cellranger_results_file, sep="\t")
    if include:
        flag_ = cellranger_results["node_id"].isin(node_id_lst)
    else:
        flag_ = ~(cellranger_results["node_id"].isin(node_id_lst))
    barcodes_lst = cellranger_results_test["barcodes"][flag_].tolist()
    print("node_id_lst: ", node_id_lst)
    print("include_node_status: ", include)
    print("get barcodeslst cell numbers: ", len(barcodes_lst))
    barcodes_lst_index = flag_.index
    if return_index:
        return barcodes_lst_index, barcodes_lst
    return barcodes_lst

def annotate_node(cellranger_results, node_cut_id):
    """
    Function: annotate the cellranger clusters with the node_cut_id
    cellranger_results: the dataframe with the cellbarcodes and correspoding node_id;
    node_cut_id: the node_id lst for cutting the clusters, [ ) format
    return:
    new dataframe with annotated cell clusters in the new column "cluster"
    
    usage:
    cellranger_cut_results = annotate_node(cellranger_results, [2536])
    cellranger_cut_results = annotate_node(cellranger_results, [2536,2543])
    
    """
    flag_ = cellranger_results["node_id"].drop_duplicates().isin(node_cut_id)
    idx_ = cellranger_results["node_id"].drop_duplicates().isin(node_cut_id).index
    cut_idx_ = idx_[flag_]
    cluster_num = len(cut_idx_)
    cellranger_results["cluster"] = cluster_num
    print(cluster_num)
    for i in reversed(range(cluster_num)):
        tmp_cut = cut_idx_[i]
        cellranger_results["cluster"][cellranger_results.index < tmp_cut ] = i
    return cellranger_results


def get_region_position(regions,chr):
    """
    Function: get the start position for chromosome in the genome coordinate
    regions:chr\t block_start \t block_end for each line
    chr: int, specific chr num
    Example:
    baf_dat_dir = '/storage/yhhuang/users/rthuang/processed_data/xcltk/xianjie-cpos/mkn45_020821/mkn45-500k/mkn45-500k-csp-post/phase-snp-even/'
    regions = np.genfromtxt(baf_dat_dir + "blocks.50kb.tsv", delimiter="\t")
    get_region_position(regions,4)
    """
    flag_regions = list(regions[:,0] == chr)
    first_index = flag_regions.index(True)
    ## TODO THERE SHOULD BE ANOTHER FUNCTION TO GET THE FIRST INDEX DIRECTLY #rt
    print("Start position for chr",chr,": ",first_index)
    return first_index



