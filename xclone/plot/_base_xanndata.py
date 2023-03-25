"""[XClone] 
Base heatmap plotting functions for XClone CNV calling.
"""
# Author: Rongting Huang
# rthuang@connect.hku.hk
# Date: 2021/07/23
# update: 2021/07/23

import os
import warnings

from anndata import AnnData

from itertools import cycle


import numpy as np
import pandas as pd

from matplotlib import pyplot as plt
# from matplotlib.patches import Patch
import matplotlib.patches as mpatches
import matplotlib.image as mpimg
import seaborn as sns

from scipy.sparse import issparse

def prepare_Xheatmap_df(Xdata, Xlayer):
    """
    Function: prepare pd.DataFrame from anndata for heatmap visualization.
    """
    if Xlayer == None:
        if issparse(Xdata.X):
            X_df = pd.DataFrame(Xdata.X.A, columns = Xdata.var["chr"], index= Xdata.obs.index)
        else:
            X_df = pd.DataFrame(Xdata.X, columns = Xdata.var["chr"], index= Xdata.obs.index)
    else:
        if issparse(Xdata.layers[Xlayer]):
            X_df = pd.DataFrame(Xdata.layers[Xlayer].A, columns = Xdata.var["chr"], index= Xdata.obs.index)
        else:
            X_df = pd.DataFrame(Xdata.layers[Xlayer], columns = Xdata.var["chr"], index= Xdata.obs.index)
    X_df.columns.names = [""] 
    # not display the df.columns name in the plot, eg. chr
    return X_df


def get_chr_breakpoints(chr_df):
    """
    """
    chr_start_idxs = chr_df["chr"].drop_duplicates(keep="first").index
    chr_end_idxs = chr_df["chr"].drop_duplicates(keep="last").index
    break_tuples = []  
    for start_idx, end_idx in zip(chr_start_idxs, chr_end_idxs):
        break_tuples.append((start_idx, end_idx))
    return break_tuples

def add_chrlabel(g, chr_df, set_grid = False, labelcolor="k"):
    """
    add chr information in sns.clustermap.

    update: use major and minor tick to add chr labels in the middle of the chr
    and add ticks at the breakpoints.

    """
    break_tuples = get_chr_breakpoints(chr_df)
    # break_tuples format: [(0,122),(122,280)...]
    ax = g.ax_heatmap
    ticks_label = []
    for break_points in break_tuples:
        ## breakpoint-major ticks
        major_ticks = np.append(ax.get_xticks(), int(break_points[0]))
        ax.set_xticks(major_ticks)
        ax.tick_params(which = "major", direction='out', length=3, width=1, colors='k',
               labelcolor = labelcolor)
        ## chrlabel-minor ticks
        minor_ticks = np.append(ax.get_xticks(minor=True), int(float(break_points[1] + break_points[0] + 1) / 2.0))
        ax.set_xticks(minor_ticks, minor = True)
        ax.tick_params(which = "minor", direction='out', length=2, width=0.5, colors='grey',
               labelcolor = labelcolor, labeltop = False)
        # ticks_label.append(chr_df["chr"][break_points[0]]) # 1,2,3...
        ticks_label.append("chr" + chr_df["chr"][break_points[0]]) # chr1, chr2...
    ax.set_xticklabels(ticks_label, rotation=45, ha='center', minor = True)
    # ax.set_xticklabels(ax.xaxis.get_majorticklabels(), rotation=45, ha='center')
    # ax.set_yticklabels(ax.get_yticklabels(), rotation=0)
    if set_grid:
        ax.grid(which="major", color="black", linestyle='-', linewidth=1)


def set_Xheatmap_args(legend_pos=2):
    """
    set default params for Xheatmap
    """
    paras = {}
    ## axes setting
    paras['row_cluster'] = False
    paras['col_cluster'] = False
    paras['yticklabels'] = False
    paras['xticklabels'] = False
    ## data setting
    paras['cmap'] = 'vlag'
    paras['center'] = False
    paras['rasterized'] = True
    ## figure setting
    paras['figsize'] = (16,8) # publication
    ### webplot
    # paras['figsize'] = (12,8)

    ## legend
    legend_paras = {}
    legend_paras["pos"] = legend_pos
    return paras, legend_paras

def Xheatmap(
    Xdata: AnnData,
    Xlayer = None,
    fillna_value = 0,
    cell_anno_key = "Cluster",
    clusters_display_name = "Clone",
    title = 'XClone',
    X_xlabel = "Chromosomes",
    X_ylabel = "Cells",
    legend_pos = 2,
    add_chr: bool = True,
    change_colorbar: bool = True,
    colorbar_ticks = [0, 1, 2], 
    colorbar_label = ["copy loss", "copy neutral", "copy gain"], 
    colorbar_name = "CNV states",
    set_grid: bool = False, 
    save_file: bool = False,
    set_dpi = 150,
    out_file = 'XClone_Xheatmap_test.pdf',
    **kwargs):
    """
    Function: base Heatmap for XClone RDR and BAF analysis.
    params:
    ------
    Xdata: AnnData in XClone;
    cell_anno_key: col(colnames) in Xdata.obs used for annotate clusters,like group_by[in common];
    palette: color palette in seaborn, default is 'vlag';
    **kwargs:
    can assign any params in sns.clustermap; can also update the default paras 
    in set_Xheatmap_args()

    """
    ## check Xdata annotation
    anno_is_na_ = Xdata.obs[cell_anno_key].isna()
    update_Xdata = Xdata.copy()[~anno_is_na_, :]

    ## pd.DataFrame format for sns.clustermap
    X_df = prepare_Xheatmap_df(update_Xdata, Xlayer)
        
    ## fillna
    X_df = X_df.fillna(fillna_value)
    
    ## _clusters for obs(rows) and _chrs for cols
    _clusters = update_Xdata.obs[cell_anno_key].astype("str")
    _clusters.name = clusters_display_name ## not display the name for cell clusters eg. random cluster
    _chrs = X_df.columns
    
    ## color setting-Lookup Table(lut) for color mapping
    cell_palette = cycle(sns.color_palette("colorblind", len(set(_clusters))))
    lut1 = dict(zip(_clusters.unique(), cell_palette))
    row_colors = _clusters.map(lut1)
    
    chr_palette = cycle(['#525252', '#cccccc'])
    lut2 = dict(zip(_chrs.unique(), chr_palette))
    col_colors = _chrs.map(lut2)
    
    ##params setting
    paras = {}
    paras['data'] = X_df
    paras['row_colors'] = row_colors
    paras['col_colors'] = col_colors
    
    ## colorbar setting(1)
    if change_colorbar:
        if colorbar_ticks is not None:
            paras['cbar_kws'] = {"ticks": colorbar_ticks}
            #  "orientation": "horizontal"

    ### update the default paras
    args, legend_paras = set_Xheatmap_args(legend_pos=legend_pos)
    paras.update(**args)
    ### update the user assign paras 
    paras.update(**kwargs)

    # paras["cmap"].set_bad(color='gray')

    ## plot the heatmap with sns.clustermap
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        g = sns.clustermap(**paras)

        if add_chr:
            ## prepare for chr_label-in add_chrlabel
            chr_df = update_Xdata.var["chr"].reset_index()
            del chr_df["index"]
            add_chrlabel(g, chr_df, set_grid = set_grid)
        
        ## colorbar setting(2)
        if change_colorbar:
            # print(g.ax_cbar.get_yticks(minor = True))
            # print(g.ax_cbar.get_yticks())
            # print(g.ax_cbar.get_position())

            # g.ax_cbar.set_yticks(colorbar_ticks)
            ## should set at `cbar_kws`
            if colorbar_label is not None:
                g.ax_cbar.set_yticklabels(colorbar_label)

        g.ax_cbar.set_xlabel(colorbar_name)

            # g.cax.set_yticks([0.125, 0.375, 0.625,0.875])
            # g.cax.set_yticklabels(colorbar_label)
        g.fig.suptitle(title)
        g.ax_heatmap.set_xlabel(X_xlabel)
        g.ax_heatmap.set_ylabel(X_ylabel)
    
    ## add clusters legend
    cluster_handles = [mpatches.Patch(facecolor=lut1[name]) for name in lut1]
    if legend_paras["pos"] == 1:
        plt.legend(cluster_handles, lut1, title = cell_anno_key,
                bbox_to_anchor=(1, 1), bbox_transform=plt.gcf().transFigure, loc='upper right')
    if legend_paras["pos"] == 2:
        plt.legend(cluster_handles, lut1, title = cell_anno_key,
                bbox_to_anchor=(0, 0), bbox_transform=plt.gcf().transFigure, loc='lower left')
    if save_file == False:
        plt.show()
    else:
        plt.savefig(out_file, bbox_inches='tight', dpi=set_dpi)
        plt.show()
        plt.close()



def set_XXheatmap_args(legend_pos=2, legend_mode=1):
    """
    set default params for XXheatmap
    """
    paras = {}
    ## axes setting
    paras['row_cluster'] = False
    paras['col_cluster'] = False
    paras['yticklabels'] = False
    paras['xticklabels'] = False
    ## data setting
    paras['cmap'] = 'vlag'
    paras['center'] = False
    paras['rasterized'] = True
    ## figure setting
    paras['figsize'] = (16,8) # publication
    ### webplot
    # paras['figsize'] = (12,8)


    ## legend
    legend_paras = {}
    legend_paras["pos"] = legend_pos
    legend_paras["mode"] = legend_mode
    return paras, legend_paras

def multi_lut_set(Xdata, cell_anno_key, clusters_display_name, sns_palette_lst):
    """
    Function: for multi cluster mappings
    paras
    ------
    anno_key, 
    anno_display_name, 
    sns_palette,

    return
    ------
    row_colors, 
    row_luts, dict for one legend
    lut_dic_lst, dict lst for multiple legend
    """
    lut_dic_lst = []  ## multi-luts
    for anno_key, anno_display_name, sns_palette in zip(cell_anno_key,clusters_display_name, sns_palette_lst):
        idx = cell_anno_key.index(anno_key)
        _clusters  = Xdata.obs[anno_key]
        _clusters.name = anno_display_name ## not display the name for cell clusters eg. random cluster
        ## color setting-Lookup Table(lut) for color mapping
        cell_palette = cycle(sns.color_palette(sns_palette, len(set(_clusters))))
        tmp_lut = dict(zip(_clusters.unique(), cell_palette))
        tmp_lut_copy = tmp_lut.copy()
        lut_dic_lst.append(tmp_lut_copy) ## important
        # https://www.kite.com/python/answers/how-to-append-a-dictionary-to-a-list-in-python
        tmp_row_colors = _clusters.map(tmp_lut)
        # for add separate cell annotation legends
        if idx == 0:
            row_colors = tmp_row_colors
            row_luts = tmp_lut
        else:
            row_colors = pd.concat([row_colors,tmp_row_colors],axis=1)
            row_luts.update(**tmp_lut) # for add united cell annotation legends

    return row_colors, row_luts, lut_dic_lst

def XXheatmap(
    Xdata: AnnData,
    Xlayer = None,
    cell_anno_key = ["Cluster", "Celltype"],
    clusters_display_name = ["Clone", "Celltype"],
    sns_palette_lst = ["colorblind","muted"],
    cell_legend_title = "Cell anno",
    title = 'XClone',
    X_xlabel = "Chromosomes",
    X_ylabel = "Cells",
    legend_pos = 2,
    legend_mode = 2,
    add_chr: bool = True,
    change_colorbar: bool = True,
    colorbar_ticks = [0, 1, 2], 
    colorbar_label = ["copy loss", "copy neutral", "copy gain"], 
    colorbar_name = "CNV states",
    save_file: bool = False,
    set_dpi = 150,
    out_file = 'XClone_XXheatmap_test.pdf',
    **kwargs):
    """
    Function: base Heatmap for XClone RDR and BAF analysis.
    params:
    ------
    Xdata: AnnData in XClone;
    cell_anno_key: col(colnames) in Xdata.obs used for annotate clusters,like group_by[in common];
    update for list to contain mulituple annotations, construct color mapping dataframe for multi-cols.
    can be more than two cols;
    clusters_display_name:
    sns_palette_lst: color palette in seaborn, default is 'vlag';
    **kwargs:
    can assign any params in sns.clustermap; can also update the default paras 
    in set_Xheatmap_args()

    """
    ## pd.DataFrame format for sns.clustermap
    X_df = prepare_Xheatmap_df(Xdata, Xlayer)
    
    ## _clusters for obs(rows) and _chrs for cols
    ### obs
    row_colors, row_luts, lut_dic_lst = multi_lut_set(Xdata, cell_anno_key, clusters_display_name, sns_palette_lst)
    
    ### vars
    _chrs = X_df.columns
    #### color setting-Lookup Table(lut) for color mapping
    chr_palette = cycle(['#525252', '#cccccc'])
    lut_chr = dict(zip(_chrs.unique(), chr_palette))
    col_colors = _chrs.map(lut_chr)
    
    ## params setting
    paras = {}
    paras['data'] = X_df
    paras['row_colors'] = row_colors
    paras['col_colors'] = col_colors

    ## colorbar setting(1)
    if change_colorbar:
        if colorbar_ticks is not None:
            paras['cbar_kws'] = {"ticks": colorbar_ticks}
            #  "orientation": "horizontal"

    ### update the default paras
    args, legend_paras = set_XXheatmap_args(legend_pos, legend_mode)
    paras.update(**args)
    ### update the user assign paras
    paras.update(**kwargs)

    ## plot the heatmap with sns.clustermap
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        g = sns.clustermap(**paras)
        if add_chr:
            ## prepare for chr_label-in add_chrlabel
            chr_df = Xdata.var["chr"].reset_index()
            del chr_df["index"]
            add_chrlabel(g, chr_df)
        
        ## colorbar setting(2)
        if change_colorbar:
            if colorbar_label is not None:
                g.ax_cbar.set_yticklabels(colorbar_label)
        
        g.ax_cbar.set_xlabel(colorbar_name)
        
        g.fig.suptitle(title)
        g.ax_heatmap.set_xlabel(X_xlabel)
        g.ax_heatmap.set_ylabel(X_ylabel)
        # ax = g.fig.gca()  ## seem around legend
        # ax.set_xlabel('testx')
        # ax.set_ylabel('testy')
    
    ## add legend
    ### add clusters legend version 1
    if legend_paras["mode"] == 1:
        cluster_handles = [mpatches.Patch(facecolor=row_luts[name]) for name in row_luts]
        if legend_paras["pos"] == 1:
            plt.legend(cluster_handles, row_luts, title = cell_legend_title,
                    bbox_to_anchor=(1, 1), bbox_transform=plt.gcf().transFigure, loc='upper right')
        if legend_paras["pos"] == 2:
            plt.legend(cluster_handles, row_luts, title = cell_legend_title,
                    bbox_to_anchor=(0, 0), bbox_transform=plt.gcf().transFigure, loc='lower left')

    ### add clusters legend version 2
    if legend_paras["mode"] == 2:
        for tmp_lut in lut_dic_lst:
            idx = lut_dic_lst.index(tmp_lut)
            # Create a legend for the first line.
            # Create another legend for the second line...
            if idx < len(lut_dic_lst)-1: 
                cluster_handles = [mpatches.Patch(facecolor=tmp_lut[name]) for name in tmp_lut]
                tmp_legend = plt.legend(handles=cluster_handles, labels = tmp_lut, 
                    title = clusters_display_name[idx],
                    bbox_to_anchor=(1, 1), bbox_transform=plt.gcf().transFigure, loc='upper right')
                # Add the legend manually to the current Axes.
                plt.gca().add_artist(tmp_legend)
            else:
                cluster_handles2 = [mpatches.Patch(facecolor=tmp_lut[name]) for name in tmp_lut]
                plt.legend(cluster_handles2, tmp_lut, title = clusters_display_name[idx],
                    bbox_to_anchor=(1, 1), bbox_transform=plt.gcf().transFigure, loc='upper left')
    ## show or save the file
    if save_file == False:
        plt.show()
    else:
        plt.savefig(out_file, bbox_inches='tight', dpi=set_dpi)
        plt.show()
        plt.close()


def Xheatmap_addref(
    Xdata: AnnData,
    ref_Xdata: AnnData,
    Xlayer = None,
    ref_Xlayer = None,
    cell_anno_key = "Cluster",
    clusters_display_name = "",
    plot_title = 'XClone with ref',
    ref_legend_pos = 1,
    add_ref_chr: bool = False,
    figsize=(14,12),
    remove_tmp: bool = False,
    save_file: bool = False,
    merge_out_file = 'XClone_Xheatmap_with_reference_test.pdf',
    dpi = 300,
    **kwargs):
    """
    add subplot Xheatmap for ref dataset
    Here, we just subset the ref_Xdata and Xdata 
    from the original whole Xdata, assumes the two have same chrs info

    **kwargs: for Xheatmap
    """
    ## 1.params setting
    ref_paras = {}
    ref_paras["Xdata"] = ref_Xdata
    ref_paras["Xlayer"] = ref_Xlayer
    ref_paras["cell_anno_key"] = cell_anno_key
    ref_paras["add_chr"] = add_ref_chr
    ref_paras["legend_pos"] = ref_legend_pos
    ref_paras["clusters_display_name"] = clusters_display_name
    ref_paras["X_xlabel"] = ""
    ref_paras["X_ylabel"] = "Refernece cells"
    ref_paras["save_file"] = True # default save
    ref_paras["out_file"] = 'tmp_XClone_ref.png'
    ref_paras["figsize"] = (12,6)
    ref_paras.update(**kwargs)

    obs_paras = {}
    obs_paras["Xdata"] = Xdata
    obs_paras["Xlayer"] = Xlayer
    obs_paras["cell_anno_key"] = cell_anno_key
    obs_paras["title"] = ""
    obs_paras["X_ylabel"] = "Cells"
    obs_paras["save_file"] = True
    obs_paras["out_file"] = 'tmp_XClone_obs.png'
    obs_paras["figsize"] = (12,6)
    obs_paras.update(**kwargs)

    ## 2. plot and save for merge
    Xheatmap(**ref_paras)
    Xheatmap(**obs_paras)

    ## 3.create your suplots from temporal images
    fig, axarr = plt.subplots(2, 1, figsize=figsize)

    axarr[0].imshow(mpimg.imread(ref_paras["out_file"]))
    axarr[1].imshow(mpimg.imread(obs_paras["out_file"]))
    
    ### setting for the merge plot
    ### turn off x and y axis
    [ax.set_axis_off() for ax in axarr.ravel()]
    fig.suptitle(plot_title) # fontsize=16
    plt.tight_layout()

    ## save and remove temporal images
    if save_file == False:
        plt.show()
    else:
        plt.savefig(merge_out_file, bbox_inches='tight', dpi=dpi) # dpi=600
        plt.close()
        if remove_tmp:
            os.remove(ref_paras["out_file"])
            os.remove(obs_paras["out_file"])

    ## develop notes
    ## previous version-try to plot two sns.clustermap meanwhile
    ## but failed--deprecated one
    # ## pd.DataFrame format for sns.clustermap
    # X_df = prepare_Xheatmap_df(Xdata, Xlayer)
    # refX_df = prepare_Xheatmap_df(ref_Xdata, ref_Xlayer)
    # .......
    # ## plot the heatmap with sns.clustermap

    # fig, axes = plt.subplots(2,1,figsize=(12,16))
    # paras.update(ax=axes[1])
    # ref_paras.update(ax=axes[0])

    # with warnings.catch_warnings():
    #     warnings.simplefilter("ignore")
    #     ## ref
    #     sns.clustermap(**ref_paras)
    #     # sns.clustermap(ax=axes[1], **paras)
    #     sns.clustermap(**paras)
    #     if add_chr:
    #         ## prepare for chr_label-in add_chrlabel
    #         chr_df = Xdata.var["chr"].reset_index()
    #         del chr_df["index"]
    #         add_chrlabel(axes[1], chr_df)
    # fig.suptitle(title)
    #####--problems-- sns.clustermap have complex axes, can not use ax here#######
    

## to do list----- defeine params func for both obs and ref setting
def XXheatmap_addref(
    Xdata: AnnData,
    ref_Xdata: AnnData,
    Xlayer = None,
    ref_Xlayer = None,
    cell_anno_key = ["Cluster","Celltype"],
    clusters_display_name = ["Clone","Celltype"],
    sns_palette_lst = ["colorblind","muted"],
    cell_legend_title = "Cell anno",
    legend_pos = 2,
    legend_mode = 1,
    ref_cell_anno_key = ["Cluster","Celltype"],
    ref_clusters_display_name = ["Clone","Celltype"],
    ref_sns_palette_lst = ["colorblind","muted"],
    ref_cell_legend_title = "Cell anno",
    ref_legend_pos = 2,
    ref_legend_mode = 1, 
    plot_title = 'XClone XXheatmap with ref',
    add_ref_chr: bool = False,
    figsize=(14,12),
    remove_tmp: bool = False,
    save_file: bool = False,
    merge_out_file = 'XClone_XXheatmap_with_reference_test.pdf',
    dpi=300,
    **kwargs):
    """
    add subplot XXheatmap for ref dataset
    Here, we just subset the ref_Xdata and Xdata 
    from the original whole Xdata, assumes the two have same chrs info

    **kwargs: for XXheatmap
    """
    ## 1.params setting
    ref_paras = {}
    ref_paras["Xdata"] = ref_Xdata
    ref_paras["Xlayer"] = ref_Xlayer
    ref_paras["cell_anno_key"] = ref_cell_anno_key
    ref_paras["clusters_display_name"] = ref_clusters_display_name
    ref_paras["sns_palette_lst"] = ref_sns_palette_lst
    ref_paras["cell_legend_title"] = ref_cell_legend_title
    ref_paras["add_chr"] = add_ref_chr
    ref_paras["legend_pos"] = ref_legend_pos
    ref_paras["legend_mode"] = ref_legend_mode
    ref_paras["X_xlabel"] = ""
    ref_paras["X_ylabel"] = "Refernece cells"
    ref_paras["save_file"] = True # default save
    ref_paras["out_file"] = 'tmp_XClone_XXheatmap_ref.png'
    ref_paras["figsize"] = (12,6)
    ref_paras.update(**kwargs)

    obs_paras = {}
    obs_paras["Xdata"] = Xdata
    obs_paras["Xlayer"] = Xlayer
    obs_paras["cell_anno_key"] = cell_anno_key
    obs_paras["clusters_display_name"] = clusters_display_name
    obs_paras["sns_palette_lst"] = sns_palette_lst
    obs_paras["cell_legend_title"] = cell_legend_title
    obs_paras["legend_pos"] = legend_pos
    obs_paras["legend_mode"] = legend_mode
    obs_paras["title"] = ""
    obs_paras["X_ylabel"] = "Cells"
    obs_paras["save_file"] = True
    obs_paras["out_file"] = 'tmp_XClone_XXheatmap_obs.png'
    obs_paras["figsize"] = (12,6)
    obs_paras.update(**kwargs)

    ## 2. plot and save for merge
    XXheatmap(**ref_paras)
    XXheatmap(**obs_paras)

    ## 3.create your suplots from temporal images
    fig, axarr = plt.subplots(2, 1, figsize=figsize)

    axarr[0].imshow(mpimg.imread(ref_paras["out_file"]))
    axarr[1].imshow(mpimg.imread(obs_paras["out_file"]))
    
    ### setting for the merge plot
    ### turn off x and y axis
    [ax.set_axis_off() for ax in axarr.ravel()]
    fig.suptitle(plot_title) # fontsize=16
    plt.tight_layout()

    ## save and remove temporal images
    if save_file == False:
        plt.show()
    else:
        plt.savefig(merge_out_file, bbox_inches='tight', dpi=dpi) # dpi=600
        plt.close()
        if remove_tmp:
            os.remove(ref_paras["out_file"])
            os.remove(obs_paras["out_file"])


## Notes: to do list
## R heatmap-- multi layers annotation
## https://github.com/jokergoo/ComplexHeatmap
## or the new ggheatmap