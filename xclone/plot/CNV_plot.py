"""Base functions for XClone CNV plotting
"""
# Author: Rongting Huang
# Date: 2022/06
# Update: 2022/09

import numpy as np
from ._visualize  import convert_res_to_ann
from ._data import reorder_data_by_cellanno
from ._base_xanndata import Xheatmap, XXheatmap

## utils
def remove_Xdata_layers(Xdata, copy = True):
    """
    copy: bool. default, True: return copy of Xdata but remove layers.
    save memory. False: remove the layers of the Xdata.
    """
    def remove_layers(Xdata):
        dict_keys = []
        for key_ in Xdata.layers.keys():
            dict_keys.append(key_)
        for i in range(len(dict_keys)):
            Xdata.layers.pop(dict_keys[i])
        return Xdata
    
    if copy:
        Xdata_copy = Xdata.copy()
        Xdata_copy = remove_layers(Xdata_copy)
        return Xdata_copy
    else:
        Xdata = remove_layers(Xdata)
        return Xdata

def color_mapping(name = "haplotyped_cmap"):
    """
    color mapping for combination results visualization.
    """
    import palettable
    white_c = palettable.colorbrewer.diverging.PuOr_3.mpl_colors[1]
    loh1_c = palettable.colorbrewer.diverging.PuOr_3.mpl_colors[0]
    loh2_c = palettable.colorbrewer.diverging.PuOr_3.mpl_colors[2]
    # loss1_c = palettable.colorbrewer.sequential.YlGnBu_9.mpl_colors[-2]
    # loss2_c = palettable.colorbrewer.sequential.YlGnBu_9.mpl_colors[-4]
    # loss1_c = palettable.colorbrewer.qualitative.Pastel1_3.mpl_colors[0]
    # loss2_c = palettable.colorbrewer.qualitative.Pastel1_3.mpl_colors[1]
    loss1_c = palettable.colorbrewer.qualitative.Paired_5.mpl_colors[-1]
    loss2_c = palettable.colorbrewer.qualitative.Paired_8.mpl_colors[1]
    loss_c = palettable.colorbrewer.diverging.RdBu_8.mpl_colors[-1]
    gain_c = palettable.colorbrewer.diverging.RdBu_8.mpl_colors[0]

    loh_c = palettable.colorbrewer.diverging.BrBG_4.mpl_colors[-2]

    if name == "loh_haplotyped_cmap":
        color_map = [loh1_c, white_c,  loh2_c]
    if name == "loh_cmap":
        color_map = [loh_c, white_c]

    if name == "haplotyped_cmap":
        color_map = [white_c, loss_c, loh1_c, loh2_c, gain_c]
    if name == "no_haplotyped_cmap":
        color_map = [white_c, loss_c, loh_c, gain_c]
    
    if name == "haplotyped_cmap1":
        color_map = [white_c, loss1_c, loss2_c, loh1_c, loh2_c, gain_c]
    if name == "haplotyped_cmap2":
        color_map = [white_c, loss1_c, loss2_c, loh_c, gain_c]
    if name == "combine_cmap":
        color_map = [loss_c, loh_c, white_c,gain_c]
    if name == "combine_cmap2":
        color_map = [loss1_c, loss2_c, loh_c, white_c, gain_c]
    if name == "combine_cmap3":
        color_map = [loss1_c, loss2_c, loh1_c, loh2_c, white_c, gain_c]
    if name == "combine_cmap4":
        color_map = [loss_c, loh1_c, loh2_c, white_c, gain_c]
    if name == "combine_cmap5":
        color_map = [loss1_c, loss2_c, white_c, gain_c]
    return color_map

def CNV_visualization_combine(Xdata,
                              Xlayer,
                              states_num = 4,
                              cell_anno_key = "cell_type", 
                              color_map_name = None,
                              colorbar_ticks = None,
                              colorbar_label = None,
                              **kwargs):
    """
    4 states: copy loss, loh, copy neutral, copy gain
    """
    import palettable
    import seaborn as sns

    res_cnv_ad = convert_res_to_ann(Xdata, Xlayer)
    ## reorder the cells based on annotation
    res_cnv_ad_re = reorder_data_by_cellanno(res_cnv_ad, cell_anno_key=cell_anno_key)

    ## plot the Xheatmap
    center_value = float((states_num-1)/2)
    if color_map_name is None:
        color_map = color_mapping(name = "combine_cmap")
    else:
        color_map = color_mapping(name = color_map_name)

    if colorbar_ticks is None:
        colorbar_ticks = [0,1,2,3]
    if colorbar_label is None:
        colorbar_label = ["copy loss", "loh", "copy neutral", "copy gain"]

    Xheatmap(res_cnv_ad_re, cell_anno_key = cell_anno_key, center=center_value, 
    colorbar_ticks = colorbar_ticks,
    colorbar_label = colorbar_label,
    cmap = color_map, **kwargs)


def complex_CNV_visualization_combine(Xdata, 
                                      Xlayer,
                                      states_num = 4,
                                      cell_anno_key = ["cluster", "cell_type"], 
                                      clusters_display_name = ["Clone", "Celltype"], 
                                      color_map_name = None,
                                      colorbar_ticks = None,
                                      colorbar_label = None,
                                      **kwargs):
    """
    """
    import palettable
    import seaborn as sns

    res_cnv_ad = convert_res_to_ann(Xdata, Xlayer)
    ## reorder the cells based on annotation
    res_cnv_ad_re = reorder_data_by_cellanno(res_cnv_ad, cell_anno_key=cell_anno_key)

    ## plot the Xheatmap
    center_value = float((states_num-1)/2)

    if color_map_name is None:
        color_map = color_mapping(name = "combine_cmap")
    else:
        color_map = color_mapping(name = color_map_name)
    
    if colorbar_ticks is None:
        colorbar_ticks = [0,1,2,3]
    if colorbar_label is None:
        colorbar_label = ["copy loss", "loh", "copy neutral", "copy gain"]

    XXheatmap(res_cnv_ad_re, 
    cell_anno_key = cell_anno_key, 
    clusters_display_name = clusters_display_name,
    center = center_value,
    colorbar_ticks = colorbar_ticks,
    colorbar_label = colorbar_label, 
    cmap = color_map,
    **kwargs)


## Visualization part in CNVratio estimation performance
## RDR CNV
def CNV_visualization(update_Xdata, Xlayer = "posterior_mtx", weights = True, 
                      states_weight = np.array([1,2,3]), states_num = 3,
                      cell_anno_key = "cell_type", **kwargs):
    """
    updated: 2022-08-05
    For RDR CNV visualization.
    default using .layers["posterior_mtx"] for visualization
    RDR 3 states: copy gain, copy neutral, copy loss.
    """
    # https://jiffyclub.github.io/palettable/
    import palettable
    import seaborn as sns
    if weights:
        ## transfer data to anndata for visualization
        res_cnv_weights_ad = convert_res_to_ann(update_Xdata, Xlayer, 
                                                weights = weights, states_weight = states_weight)
        center_value = 2.0
        ## reorder the cells based on annotation
        res_cnv_weights_ad1_re = reorder_data_by_cellanno(res_cnv_weights_ad, 
                                                        cell_anno_key=cell_anno_key)
        ## plot the Xheatmap
        Xheatmap(res_cnv_weights_ad1_re, cell_anno_key = cell_anno_key, center = center_value,
        colorbar_ticks = [1.1,2,3],
        colorbar_label = ["copy loss", "copy neutral", "copy gain"], **kwargs)
    else:
        # category visualization
        ## transfer data to anndata for visualization
        res_cnv_ad = convert_res_to_ann(update_Xdata, Xlayer)
        center_value = float(states_num//2)
        ## reorder the cells based on annotation
        res_cnv_ad_re = reorder_data_by_cellanno(res_cnv_ad, cell_anno_key=cell_anno_key)
        ## plot the Xheatmap
        white_c = palettable.colorbrewer.diverging.PuOr_3.mpl_colors[1]
        loss_c = palettable.colorbrewer.diverging.RdBu_8.mpl_colors[-1]
        gain_c = palettable.colorbrewer.diverging.RdBu_8.mpl_colors[0]
        color_map = [loss_c, white_c, gain_c]
        Xheatmap(res_cnv_ad_re, cell_anno_key = cell_anno_key, center=center_value, 
        colorbar_ticks = [0,1,2],
        colorbar_label = ["copy loss",  "copy neutral", "copy gain"],
        cmap = color_map, **kwargs)

def CNV_visualization_complex(update_Xdata, Xlayer = "posterior_mtx", weights = True, 
                      states_weight = np.array([1,2,3]), states_num = 3,
                      cell_anno_key = ["cluster", "cell_type"],
                      clusters_display_name = ["Clone", "Celltype"], **kwargs):
    """
    updated: 2022-08-23
    For RDR CNV visualization.
    default using .layers["posterior_mtx"] for visualization
    RDR 3 states: copy gain, copy neutral, copy loss.
    """
    # https://jiffyclub.github.io/palettable/
    import palettable
    import seaborn as sns
    if weights:
        ## transfer data to anndata for visualization
        res_cnv_weights_ad = convert_res_to_ann(update_Xdata, Xlayer, 
                                                weights = weights, states_weight = states_weight)
        center_value = 2.0
        ## reorder the cells based on annotation
        res_cnv_weights_ad1_re = reorder_data_by_cellanno(res_cnv_weights_ad, 
                                                        cell_anno_key=cell_anno_key)
        ## plot the Xheatmap
        XXheatmap(res_cnv_weights_ad1_re, cell_anno_key = cell_anno_key, clusters_display_name = clusters_display_name,
        center = center_value,
        colorbar_ticks = [1,2,3],
        colorbar_label = ["copy loss", "copy neutral", "copy gain"],**kwargs)
    else:
        # category visualization
        ## transfer data to anndata for visualization
        res_cnv_ad = convert_res_to_ann(update_Xdata, Xlayer)
        center_value = float(states_num//2)
        ## reorder the cells based on annotation
        res_cnv_ad_re = reorder_data_by_cellanno(res_cnv_ad, cell_anno_key=cell_anno_key)
        ## plot the Xheatmap
        white_c = palettable.colorbrewer.diverging.PuOr_3.mpl_colors[1]
        loss_c = palettable.colorbrewer.diverging.RdBu_8.mpl_colors[-1]
        gain_c = palettable.colorbrewer.diverging.RdBu_8.mpl_colors[0]
        color_map = [loss_c, white_c, gain_c]
        XXheatmap(res_cnv_ad_re, cell_anno_key = cell_anno_key, clusters_display_name = clusters_display_name,
        center=center_value, 
        colorbar_ticks = [0,1,2],
        colorbar_label = ["copy loss",  "copy neutral", "copy gain"],
        cmap = color_map, **kwargs)


## BAF CNV 
def CNV_visualization2(update_Xdata, Xlayer = "posterior_mtx", weights = True, 
                       states_weight = np.array([1,2,3]), states_num = 3,
                       cell_anno_key = "cell_type", **kwargs):
    """
    updated: 2022-08-05
    default using .layers["posterior_mtx"] for visualization
    BAF 3 states-  2 copy loss(haplotype) and one copy neutral
    """
    from seaborn import palettes
    # https://jiffyclub.github.io/palettable/
    import palettable
    import seaborn as sns
    ## transfer data to anndata for visualization
    if weights:
        ## transfer data to anndata for visualization
        res_cnv_weights_ad = convert_res_to_ann(update_Xdata, Xlayer, weights = weights, states_weight = states_weight)
        center_value = 2.0
        ## reorder the cells based on annotation
        res_cnv_weights_ad1_re = reorder_data_by_cellanno(res_cnv_weights_ad, cell_anno_key=cell_anno_key)
        ## plot the Xheatmap
        Xheatmap(res_cnv_weights_ad1_re, cell_anno_key = cell_anno_key, center = center_value,
        colorbar_ticks = [1, 2, 3], colorbar_label = ["copy loss-A",  "copy neutral",  "copy loss-B"],
        cmap = palettable.colorbrewer.diverging.PuOr_3.mpl_colormap, **kwargs)
    else:
        # category visualization
        ## transfer data to anndata for visualization
        res_cnv_ad = convert_res_to_ann(update_Xdata, Xlayer)
        center_value = float(states_num//2)
        ## reorder the cells based on annotation
        res_cnv_ad_re = reorder_data_by_cellanno(res_cnv_ad, cell_anno_key=cell_anno_key)
        ## plot the Xheatmap
        Xheatmap(res_cnv_ad_re, cell_anno_key = cell_anno_key, center = center_value, 
        colorbar_ticks = [0.3, 1, 1.7], colorbar_label = ["copy loss-A",  "copy neutral",  "copy loss-B"],
        cmap = palettable.colorbrewer.diverging.PuOr_3.mpl_colors, **kwargs)
        # colorbar_ticks = [0, 1, 2]

def CNV_visualization2_complex(update_Xdata, Xlayer = "posterior_mtx", weights = True, 
                       states_weight = np.array([1,2,3]), states_num = 3,
                       cell_anno_key = ["cluster", "cell_type"],
                       clusters_display_name = ["Clone", "Celltype"], **kwargs
                       ):
    """
    updated: 2022-08-23
    default using .layers["posterior_mtx"] for visualization
    BAF 3 states-  2 copy loss(haplotype) and one copy neutral
    """
    from seaborn import palettes
    # https://jiffyclub.github.io/palettable/
    import palettable
    import seaborn as sns
    ## transfer data to anndata for visualization
    if weights:
        ## transfer data to anndata for visualization
        res_cnv_weights_ad = convert_res_to_ann(update_Xdata, Xlayer, weights = weights, states_weight = states_weight)
        center_value = 2.0
        ## reorder the cells based on annotation
        res_cnv_weights_ad1_re = reorder_data_by_cellanno(res_cnv_weights_ad, cell_anno_key=cell_anno_key)
        ## plot the Xheatmap
        XXheatmap(res_cnv_weights_ad1_re, cell_anno_key = cell_anno_key, clusters_display_name = clusters_display_name,
        center = center_value,
        colorbar_ticks = [1, 2, 3], colorbar_label = ["copy loss",  "copy neutral",  "copy loss"],
        cmap = palettable.colorbrewer.diverging.PuOr_3.mpl_colormap, **kwargs)
    else:
        # category visualization
        ## transfer data to anndata for visualization
        res_cnv_ad = convert_res_to_ann(update_Xdata, Xlayer)
        center_value = float(states_num//2)
        ## reorder the cells based on annotation
        res_cnv_ad_re = reorder_data_by_cellanno(res_cnv_ad, cell_anno_key=cell_anno_key)
        ## plot the Xheatmap
        XXheatmap(res_cnv_ad_re, cell_anno_key = cell_anno_key, clusters_display_name = clusters_display_name,
        center = center_value, 
        colorbar_ticks = [0, 1, 2], colorbar_label = ["copy loss",  "copy neutral",  "copy loss"],
        cmap = palettable.colorbrewer.diverging.PuOr_3.mpl_colors, **kwargs)




def CNV_LOH_visualization(Xdata, 
                          haplotype = False,
                          Xlayer = None,
                          cell_anno_key = "cell_type",
                          clusters_display_name = "Celltype", **kwargs):
    """
    """
    res_cnv_ad = remove_Xdata_layers(Xdata)
    if Xlayer is not None:
        pass
    else:
        if haplotype == True:
            Xlayer = "lohprob_corrected"
        else:
            Xlayer = "lohprob_corrected_merge"

    res_cnv_ad.layers[Xlayer] = Xdata.layers[Xlayer].copy()
    
    res_cnv_ad1 = convert_res_to_ann(res_cnv_ad, Xlayer)
    res_cnv_ad_re = reorder_data_by_cellanno(res_cnv_ad1, cell_anno_key=cell_anno_key)

    if haplotype == True:
        ## support 3 states
        ## copy neutral, copy loss, loh1, loh2, copy gain
        ## loh1, copy neutral, loh2
        color_map = color_mapping("loh_haplotyped_cmap")
        Xheatmap(res_cnv_ad_re, cell_anno_key = cell_anno_key, 
        clusters_display_name = clusters_display_name,
        center = 1,
        # colorbar_ticks = [0, 1, 2], 
        colorbar_ticks = [0.25, 1, 1.75], 
        colorbar_label = [ "LOH1", "copy neutral", "LOH2"],
        cmap = color_map, **kwargs)

    else:
        ## support 2 states
        ##  loh, copy neutral
        color_map = color_mapping("loh_cmap")
        Xheatmap(res_cnv_ad_re, cell_anno_key = cell_anno_key, 
        clusters_display_name = clusters_display_name,
        center = 0.5,
        # colorbar_ticks = [0, 1], 
        colorbar_ticks = [0.25, 0.75], 
        colorbar_label = ["LOH", "copy neutral"],
        cmap = color_map, **kwargs)


def CNV_combine_visualization(Xdata, 
                              Xlayer = "combined_states", 
                              mode = 0,
                              cell_anno_key = "cell_type",
                              **kwargs):
    """
    XClone RDR and BAF combination visualization.
    """
    res_cnv_ad = remove_Xdata_layers(Xdata)
    res_cnv_ad.X = Xdata.layers[Xlayer].copy()
    res_cnv_ad_re = reorder_data_by_cellanno(res_cnv_ad, cell_anno_key=cell_anno_key)
    
    if mode == 0:
        ## merge_copyloss = False, merge_loh = False
        ## support 6 states
        color_map = color_mapping("haplotyped_cmap1")
        Xheatmap(res_cnv_ad_re, cell_anno_key = cell_anno_key, 
            center = 2.5,
            colorbar_ticks = [0, 1, 2, 3, 4, 5], 
            colorbar_label = ["copy neutral", "copy lossA", "copy lossB", "LOH-A", "LOH-B", "copy gain"],
            cmap = color_map, **kwargs)
    
    if mode == 1:
        ## merge_copyloss = False, merge_loh = True
        ## support 5 states
        ## copy neutral, copy lossA, copy lossB, LOH, copy gain
        color_map = color_mapping("haplotyped_cmap2")
        Xheatmap(res_cnv_ad_re, cell_anno_key = cell_anno_key, 
            center = 2,
            colorbar_ticks = [0, 1, 2, 3, 4], 
            colorbar_label = ["copy neutral", "copy lossA", "copy lossB", "LOH", "copy gain"],
            cmap = color_map, **kwargs)

    if mode == 2:
        ## merge_copyloss = True, merge_loh = False
        ## support 5 states
        ## copy neutral, copy loss, loh1, loh2, copy gain
        color_map = color_mapping("haplotyped_cmap")
        Xheatmap(res_cnv_ad_re, cell_anno_key = cell_anno_key, 
            center = 2,
            colorbar_ticks = [0, 1, 2, 3, 4], 
            colorbar_label = ["copy neutral", "copy loss", "LOH-A", "LOH-B", "copy gain"],
            cmap = color_map, **kwargs)
            
    if mode == 3:
        ## merge_copyloss = True, merge_loh = True
        ## support 4 states
        ## copy neutral, copy loss, loh, copy gain
        color_map = color_mapping("no_haplotyped_cmap")
        Xheatmap(res_cnv_ad_re, cell_anno_key = cell_anno_key, 
            center = 1.5, 
            colorbar_ticks = [0, 1, 2, 3], 
            colorbar_label = ["copy neutral", "copy loss", "LOH", "copy gain"],
            cmap = color_map, **kwargs)


def CNV_combine_visualization_complex(Xdata, Xlayer = "combined_states", haplotype = False,
                              cell_anno_key = ["cluster", "cell_type"],
                              clusters_display_name = ["Clone", "Celltype"], **kwargs):
    """
    XClone RDR and BAF combination visualization.
    need update according to the FUNC `CNV_combine_visualization`.
    """
    res_cnv_ad = remove_Xdata_layers(Xdata)
    
    res_cnv_ad.X = Xdata.layers[Xlayer].copy()
    res_cnv_ad_re = reorder_data_by_cellanno(res_cnv_ad, cell_anno_key = cell_anno_key)

    if haplotype == True:
        ## support 5 states
        ## copy neutral, copy loss, loh1, loh2, copy gain
        color_map = color_mapping("haplotyped_cmap")
        XXheatmap(res_cnv_ad_re, cell_anno_key = cell_anno_key, clusters_display_name = clusters_display_name,
        center = 2,
        colorbar_ticks = [0, 1, 2, 3, 4], 
        colorbar_label = ["copy neutral", "copy loss", "LOH1", "LOH2", "copy gain"],
        cmap = color_map, **kwargs)

    else:
        ## support 4 states
        ## copy neutral, copy loss, loh, copy gain
        color_map = color_mapping("no_haplotyped_cmap")
        XXheatmap(res_cnv_ad_re, cell_anno_key = cell_anno_key, clusters_display_name = clusters_display_name,
        center = 1.5, 
        colorbar_ticks = [0, 1, 2, 3], 
        colorbar_label = ["copy neutral", "copy loss", "LOH", "copy gain"],
        cmap = color_map, **kwargs)