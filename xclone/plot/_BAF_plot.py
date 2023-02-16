"""Base functions for XClone BAF plotting
"""

# Author: Rongting Huang
# Date: 2021/03/15
# update: 2022/04/21

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from ._base_manhattan import manhattan
from ._base_xanndata import Xheatmap
from ._data import reorder_data_by_cellanno


## Part I: bulk BAF visualization
def prepare_manhattan_df(Xdata, AD_key = "AD", DP_key = "DP", filter_nan = True, fillna_value = 0.5):
    """
    Function: prepare pd.DataFrame from anndata for manhattan visualization.
    Default: for BAF pesudo bulk allele frequency visualization
    """
    # import scipy as sp
    from scipy.sparse import issparse
    
    ## regions information
    df_plot = Xdata.var[["chr", "stop"]]
    df_plot = df_plot.rename(columns={"chr":"chrom", "stop": "pos"})

    ### merge AD and DP for cells to get the bulk level
    if issparse(Xdata.layers[AD_key]):
        AD = Xdata.layers[AD_key].A
    else:
        AD = Xdata.layers[AD_key]
    
    if issparse(Xdata.layers[DP_key]):
        DP = Xdata.layers[DP_key].A
    else:
        DP = Xdata.layers[DP_key]
    
    AD_sum = AD.sum(axis=0)
    DP_sum = DP.sum(axis=0) 

    df_plot["AD"] = AD_sum
    df_plot["DP"] = DP_sum

    if filter_nan:
        FLAG_ = ~(df_plot['DP'] == 0)
        df_plot = df_plot[FLAG_]
        df_plot["BAF"] = df_plot["AD"]/df_plot["DP"]
    else:
        df_plot["BAF"] = df_plot["AD"]/df_plot["DP"]
        df_plot["BAF"] = df_plot["BAF"].fillna(fillna_value)

    ### get the BAF_value
    # df_plot["BAF"] = df_plot["AD"]/df_plot["DP"]
    ## return the dataframe for plot
    return df_plot[["chrom","pos", "BAF"]]


def visualize_bulk_BAF(Xdata,
                       AD_key = "AD",
                       DP_key = "DP",
                       filter_nan = True,
                       chr_lst = None, 
                       plot_title = "visualize_bulk_BAF test", 
                       save_fig = None, 
                       **kwargs):
    """
    Function:
    visualize_bulk_BAF in manhattan plot.
    params:
    -------
    df_plot = prepare_manhattan_df(Xdata), pd.DataFrame with 3 cols named ["chrom","pos", "BAF"]
    plot_title: title for the fig
    chr_lst = ['2','4'], default None -> plot the whole genome in data df_plot
    save_fig: path+filename for saving the fig, default None means just show the plot
    **kwargs: params for manhattan plot, like figsize=(18,6)... 

    Example:
    -------
    xclone.pl.visualize_bulk_BAF(df_plot, plot_title = "test")
    xclone.pl.visualize_bulk_BAF(df_plot, plot_title = "test", figsize=(16,5))
    xclone.pl.visualize_bulk_BAF(df_plot, plot_title = "test", save_fig="./test.pdf",figsize=(16,5))

    xclone.pl.visualize_bulk_BAF(df_plot, chr_lst = ["2","4", "6"], plot_title = "test")
    """
    import matplotlib.pyplot as plt
    ## prepare data
    df_plot = prepare_manhattan_df(Xdata, AD_key, DP_key, filter_nan)

    if chr_lst==None:
        # plt = lp.get_pyplot()
        manhattan(df_plot,colora="#F79489", colorb="#2E8BC0",**kwargs)
        plt.axhline(0.5, color='red')
    else:
        df_plot = df_plot[df_plot["chrom"].isin(chr_lst)]
        manhattan(df_plot, colora="#75E6DA", colorb="#EF7C8E", pts_kws = {"markersize": 8, "marker":'.'}, **kwargs)
        plt.axhline(0.5, color='#167D7F')
    # plt.xlabel(xlabel)
    # plt.ylabel(ylabel)
    plt.title(plot_title, fontsize=18)
    ## savefig
    if save_fig == None:
        plt.show()
    else:
        plt.savefig(save_fig)

def compare_visualize_bulk_BAF(Xdata, 
                               chr_lst = None,
                               AD_key = "AD",
                               DP_key = "DP",
                               filter_nan = True,   
                               cell_anno_key = "cell_type",
                               save_fig = False,
                               **kwargs):
    """
    For each celltype

    ## todo 
    # merge multi-plots and save
    """

    ## check Xdata annotation
    anno_is_na_ = Xdata.obs[cell_anno_key].isna()
    update_Xdata = Xdata.copy()[~anno_is_na_, :]

    ## prepare diffrent df for each celltype and plot multiplots
    _clusters = update_Xdata.obs[cell_anno_key].astype("str")
    for uni_cluster_ in _clusters.unique():
        is_cluster = _clusters == uni_cluster_
        df_plot = prepare_manhattan_df(update_Xdata[is_cluster, :], AD_key, DP_key, filter_nan)
        if chr_lst==None:
            manhattan(df_plot,colora="#F79489", colorb="#2E8BC0", **kwargs)
            plt.axhline(0.5, color='red')
            
        else:
            df_plot = df_plot[df_plot["chrom"].isin(chr_lst)]
            manhattan(df_plot, colora="#75E6DA", colorb="#EF7C8E", pts_kws = {"markersize": 8, "marker":'.'}, **kwargs)
            plt.axhline(0.5, color='#167D7F')
        plt.title(uni_cluster_, fontsize=18)
        plt.show()

## Part II: cell BAF visualization
def calculate_cell_BAF(Xdata, AD_key = "AD", DP_key = "DP", BAF_key = "BAF"):
    """
    Xdata: anndata.
    AD_key: Layer name for AD used.
    DP_key: Layer name for DP used.
    BAF_key: Layer name for the calculated BAF.

    """
    import scipy as sp
    # check data format
    if sp.sparse.issparse(Xdata.layers[AD_key]):
        AD = Xdata.layers[AD_key].A
    else:
        AD = Xdata.layers[AD_key]
    if sp.sparse.issparse(Xdata.layers[DP_key]):
        DP = Xdata.layers[DP_key].A
    else:
        DP = Xdata.layers[DP_key]
    
    ## calculate_cell_BAF
    Xdata.layers[BAF_key] = AD/DP

    update_Xdata = Xdata.copy()
    return update_Xdata

def visualize_cell_BAF(Xdata, 
                       Xlayer = "BAF",
                       cell_anno_key = "cell_type",
                       chr_anno_key = "chr",
                       chr_lst = None, 
                       shrink_BAF = False,
                       default_show = False,
                        **kwargs):
    """
    BAF visualization.
    
    """
    if chr_lst is None:
        pass
    else:
        FLAG_ = Xdata.var[chr_anno_key].isin(chr_lst)
        Xdata = Xdata.copy()[:,FLAG_]
    
    ## reorder for celltype
    Xdata_re = reorder_data_by_cellanno(Xdata, cell_anno_key=cell_anno_key)
    if default_show:
        Xheatmap(Xdata_re, Xlayer, center = 0.5, fillna_value = 0.5, 
                cell_anno_key = cell_anno_key, cmap="vlag",
                colorbar_ticks = None, colorbar_label = None,
                colorbar_name = "BAF values",**kwargs)
    else:
        Xheatmap(Xdata_re, Xlayer, center = 0.5, fillna_value = 0.5, 
                cell_anno_key = cell_anno_key, cmap="vlag",
                colorbar_ticks = [0, 0.5, 1], colorbar_label = [0, 0.5, 1],
                colorbar_name = "BAF values",**kwargs)

    ## change value around 0.5
    if shrink_BAF:
        shrink_Xdata_re = shrink_BAF_center(Xdata_re, Xlayer, **kwargs)
        Xheatmap(shrink_Xdata_re, Xlayer, center = 0.5, fillna_value = 0.5, 
                 cell_anno_key = cell_anno_key, cmap="vlag",
                 colorbar_ticks = None, colorbar_label = None,
                 colorbar_name = "BAF values",**kwargs)

    return None

def shrink_BAF_center(Xdata, Xlayer = "BAF", 
                      center_low = 0.45, center_high = 0.55, center_value = 0.5):
    """
    """
    update_Xdata = Xdata.copy()
    flag_matrix = (update_Xdata.layers[Xlayer] >= center_low) * (update_Xdata.layers[Xlayer] <= center_high)
    update_Xdata.layers[Xlayer][flag_matrix] = center_value
    return update_Xdata


## bulk plot setting

## -------------------------------------------------
### analysis func for MQuad
# -used for cellranger output[node cluster]
## visualization version1
def BAF_plot(data,xlabel,ylabel,title,display_state="show"):
    """
    Function: plot BAF-heatmap
    display_state: default:show  alternative: "save",saved as jpg
    Example:
    BAF_plot(data=AD.T/DP.T,xlabel="blocks in genome",ylabel="cells",title="node 8805",display_state="show")
    """
    fig = plt.figure(figsize=(6, 5))
    im = plt.imshow(data, interpolation='none', 
                aspect='auto', cmap='coolwarm') #seismic, bwr
    plt.colorbar()
    plt.tight_layout()
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)
    if display_state == "show":
        plt.show()
    if display_state == "save":
        plt.savefig(title + ".jpg")

from scipy.sparse import issparse

def plot_BAF_scatter(AD,DP,xlabel,ylabel,title,display_state="show"):
    """
    Function: plot BAF for subcluster-scatter plot
    
    Example1:
    plot_BAF_scatter(AD_sub4_node,DP_sub4_node,xlabel="chr4",ylabel="BAF",title="BAF_node8805_chr4")
    Example2:
    plot_BAF_scatter(AD_node,DP_node,xlabel="genome",ylabel="BAF",title="BAF_node8805_genome")
    ### TODO-get the accurate position and then display it(maybe try to use the pd.Series's index)
    
    """
    ## preprocessing-merge all cells and filter nan
    print("AD shape:", AD.shape)
    print("DP shape:", DP.shape)
    if type(AD) == np.matrix or issparse(AD):
        AD = AD.A
    if type(DP) == np.matrix or issparse(DP):
        DP = DP.A
    print("AD type:", type(AD))
    print("DP type:", type(DP))
    ### merge all cells to get BAF
    AD_sum = AD.sum(axis=1)
    DP_sum = DP.sum(axis=1)
    ### Filter nan value
    AD_FLAG = ~(AD_sum == 0)
    DP_FLAG = ~(DP_sum == 0)
    FLAG_ = AD_FLAG & DP_FLAG
    AD_sum_filter = AD_sum[FLAG_]
    DP_sum_filter = DP_sum[FLAG_]
    print("filter AD",len(AD_sum_filter))
    print("filter DP",len(DP_sum_filter))
    
    x_len = len(AD_sum_filter)
    y_data = AD_sum_filter/DP_sum_filter
    plt.figure(figsize=(24, 6))
    plt.scatter(range(x_len), y_data, alpha=0.5)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)


def compare_plot_BAF(AD,DP, model_proba, mtx_barcodes_file,regions_file, chr_lst,xlabel1="genome",xlabel2="chr4", plot_title="BAF_cluster_chr4(6-clusters)"):
    """
    Example:
    compare_plot_BAF(AD,DP, model_proba, mtx_barcodes_file,regions_file, [4])
    
    """
    cluster_index_lst = np.argmax(model_proba, axis=1)
    uni_cluster = np.unique(cluster_index_lst)
    mtx_barcodes_lst = pd.read_table(mtx_barcodes_file,header = None)
    regions = np.genfromtxt(regions_file, delimiter="\t")
    for uni_id in uni_cluster:
        uni_cluster_barcodes = get_cluster_barcodes_lst(uni_id, cluster_index_lst, mtx_barcodes_lst)
        AD_cluster_tmp = sub_CellCluster(AD, mtx_barcodes_file, uni_cluster_barcodes)
        DP_cluster_tmp = sub_CellCluster(DP, mtx_barcodes_file, uni_cluster_barcodes)
        AD_sub_cluster = sub_chr(AD_cluster_tmp, regions,"AD", chr_lst)
        DP_sub_cluster = sub_chr(DP_cluster_tmp, regions,"DP", chr_lst)
        plot_BAF_scatter(AD_cluster_tmp, DP_cluster_tmp, xlabel=xlabel1, ylabel="BAF", title= plot_title)
        plot_BAF_scatter(AD_sub_cluster, DP_sub_cluster, xlabel=xlabel2, ylabel="BAF", title= plot_title)



## visualization version2
# def get_BAF_plot_data(regions_file, AD, DP,filter_nan=True):
#     """
#     deprecated 2021-07-31
#     Function: get the fitted format for plot BAF scatter
#     regions_file:contain the chr and also the blocks info for each chr
#     AD: AD information or all cells/subcluster cells, mtx
#     DP: DP information or all cells/subcluster cells, mtx
    
#     Example:
    
#     """
#     regions_df = pd.read_table(regions_file, header = None)
#     regions_df.columns = ["chrom","chr_block_start","chr_block_stop"]
# #     regions_df["chrom"] = regions_df["chrom"].astype('category')
#     regions_df["genome_pos"] = regions_df.index
#     regions_df["pos_in_chr"] = regions_df["chr_block_stop"]/50000
#     regions_df["pos_in_chr"] = regions_df["pos_in_chr"].astype(np.int64)
#     ## add AD and DP
#     ### merge AD and DP for cells to get the bulk level
#     AD_sum = AD.sum(axis=1)
#     DP_sum = DP.sum(axis=1)
#     regions_df["AD"] = pd.Series(AD_sum.A[:,0])
#     regions_df["DP"] = pd.Series(DP_sum.A[:,0])
#     ## filter tha nan
#     if filter_nan:
# #         AD_flag = ~(regions_df['AD'] == 0)
#         DP_flag = ~(regions_df['DP'] == 0)
# #         FLAG_ = AD_flag & DP_flag
#         FLAG_ = DP_flag
#         regions_df = regions_df[FLAG_]
#     ### get the BAF_value
#     regions_df.loc[:,("BAF")] = regions_df.loc[:,("AD")] / regions_df.loc[:,("DP")]
#     ## return the dataframe for plot
#     return regions_df

# from ._base_manhattan import get_pyplot

def BAF_scatter_plot2(df_plot,plot_title, chr_lst=['2','4'], sub_chr_plot=False, filenotes = "filename"):
    # import limix_plot as lp
    import matplotlib.pyplot as plt
    ## todo- can change code in limix_plot@Rongtingting
#     df_plot = BAF_df[["chrom","genome_pos", "BAF"]]
#     df_plot = df_plot.rename(columns={"genome_pos": "pos"})
    # df_plot = BAF_df[["chrom","pos_in_chr", "BAF"]]
    # df_plot = df_plot.rename(columns={"pos_in_chr": "pos"})
    if sub_chr_plot:
        df_plot = df_plot[df_plot["chrom"].isin(chr_lst)]
        # plt = lp.get_pyplot()
        # plt = get_pyplot()
        manhattan(df_plot,colora="#75E6DA", colorb="#EF7C8E", pts_kws = {"markersize": 8, "marker":'.'})  
        plt.axhline(0.5, color='#167D7F')
        plt.title(plot_title, fontsize=18)
        plt.savefig(filenotes + "sub_chr.pdf")
#         plt.xlabel(xlabel)
#         plt.ylabel(ylabel)
    else:
        # plt = lp.get_pyplot()
        # plt = get_pyplot()
        manhattan(df_plot,colora="#F79489", colorb="#2E8BC0")
        plt.axhline(0.5, color='red')
        plt.title(plot_title, fontsize=18)
        plt.savefig(filenotes + "genome.pdf")


def compare_BAF_plot2(AD,DP, model_proba, mtx_barcodes_file,regions_file, chr_lst, plot_title="BAF_cluster_chr2+4(2-clusters)", cluster_note=["Clone1 (xx cells)","Clone2 (xx cells)"]):
    """
    Function:
    integration of BAF scatter plot for different clusters for specific regions and genome
    Example:
    compare_BAF_plot2(AD,DP, model_proba, mtx_barcodes_file,regions_file, ['2','4'],)
    
    """
    cluster_index_lst = np.argmax(model_proba, axis=1)
    uni_cluster = np.unique(cluster_index_lst)
    mtx_barcodes_lst = pd.read_table(mtx_barcodes_file,header = None)
    regions = np.genfromtxt(regions_file, delimiter="\t")
    cluster_cnt = 0
    for uni_id in uni_cluster:
        uni_cluster_barcodes = get_cluster_barcodes_lst(uni_id, cluster_index_lst, mtx_barcodes_lst)
        AD_cluster_tmp = sub_CellCluster(AD, mtx_barcodes_file, uni_cluster_barcodes)
        DP_cluster_tmp = sub_CellCluster(DP, mtx_barcodes_file, uni_cluster_barcodes)
        ## prepare data in dataframe for plot
        BAF_df = get_BAF_plot_data(regions_file, AD_cluster_tmp, DP_cluster_tmp)
        BAF_scatter_plot2(BAF_df, plot_title= plot_title + cluster_note[cluster_cnt], filenotes = str(uni_id))
        BAF_scatter_plot2(BAF_df, plot_title= plot_title + cluster_note[cluster_cnt], chr_lst=chr_lst, sub_chr_plot=True, filenotes = str(uni_id))
        cluster_cnt +=1

        
def node_BAF_plot2(AD,DP,mtx_barcodes_file,node, regions_file, chr_lst, plot_title1="BAF_cluster_chr2+4(node-xxx)", plot_title2="BAF_cluster_chr2+4(non_node-xxx)"):
    mkn45_10x_dloupe = "/storage/yhhuang/research/mito/mkn45/fulldepth/mode2/mkn45_5k_dloupe-group_10438-region_1_1-5120000_Y_56320001-59373566-heatmap_copy number.csv"
#     barcodes_10050_lst = get_node_barcodeslst(mkn45_10x_dloupe, 10050)
    barcodes_node_lst = get_node_barcodeslst(mkn45_10x_dloupe, node)
    barcodes_non_node_lst = get_non_node_barcodeslst(mkn45_10x_dloupe, node)

    AD_node = sub_CellCluster(AD, mtx_barcodes_file,barcodes_node_lst)
    DP_node = sub_CellCluster(DP, mtx_barcodes_file,barcodes_node_lst)
    BAF_df = get_BAF_plot_data(regions_file, AD_node, DP_node)
    BAF_scatter_plot2(BAF_df, plot_title= plot_title1, filenotes = "node10050")
    BAF_scatter_plot2(BAF_df, plot_title= plot_title1, chr_lst=chr_lst, sub_chr_plot=True, filenotes = "node10050")
    ## 
    AD_non_node = sub_CellCluster(AD, mtx_barcodes_file,barcodes_non_node_lst)
    DP_non_node = sub_CellCluster(DP, mtx_barcodes_file,barcodes_non_node_lst)
    BAF_df_ = get_BAF_plot_data(regions_file, AD_non_node, DP_non_node)
    BAF_scatter_plot2(BAF_df_, plot_title= plot_title2, filenotes = "non_node")
    BAF_scatter_plot2(BAF_df_, plot_title= plot_title2, chr_lst=chr_lst, sub_chr_plot=True, filenotes = "non_node")
    