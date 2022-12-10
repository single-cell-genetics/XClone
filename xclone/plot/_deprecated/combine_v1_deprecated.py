"""deprecated functions for XClone combine version1.
record.

Rongting Huang
"""

from .._visualize  import convert_res_to_ann
from .._data import reorder_data_by_cellanno
from .._base_xanndata import Xheatmap, XXheatmap

from ..CNV_plot import remove_Xdata_layers, color_mapping

#=====================================
### combined version1: deprecated part
## from CNV_plot.py
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