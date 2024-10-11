"""Base functions for XClone CNV plotting.
Visualization part in CNVratio estimation performance.
"""
# Author: Rongting Huang
# Date: 2022/06
# Update: 2022/12

import numpy as np

from ._visualize  import convert_res_to_ann
from ._data import reorder_data_by_cellanno
from ._base_xanndata import Xheatmap, XXheatmap

import gc

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
        
        del Xdata
        gc.collect()

        return Xdata_copy
        
    else:
        Xdata = remove_layers(Xdata)
        return Xdata

def color_mapping(cm_name = "haplotyped_cmap"):
    """
    color mapping for combination results visualization.
    For more information on color palettes, see Palettable: https://jiffyclub.github.io/palettable/
    """
    import palettable
    # color definition for XClone detected CNV states.
    # -see combine strategy.[Supp figs]
    ## neutral state
    white_c = palettable.colorbrewer.diverging.PuOr_3.mpl_colors[1]
    ## loh state
    loha_c = palettable.colorbrewer.diverging.PuOr_3.mpl_colors[0] ## orange yellow
    lohb_c = palettable.colorbrewer.diverging.PuOr_3.mpl_colors[2] ## purple
    loh_c = palettable.colorbrewer.diverging.BrBG_4.mpl_colors[-2] ## cyan-blue
    ## loss state
    lossa_c = palettable.colorbrewer.qualitative.Paired_5.mpl_colors[-1] ## pink
    lossb_c = palettable.colorbrewer.qualitative.Paired_8.mpl_colors[1] ## blue
    loss_c = palettable.colorbrewer.diverging.RdBu_8.mpl_colors[-1] ## dark blue
    ## gain state
    gain_c = palettable.colorbrewer.diverging.RdBu_8.mpl_colors[0] ## red
    gaina_c = palettable.colorbrewer.diverging.PiYG_8.mpl_colors[0] ## rose
    gainb_c = palettable.colorbrewer.diverging.Spectral_8.mpl_colors[1] ## orange
    ## RDR cmap
    if cm_name == "RDR_cmap":
        color_map = [loss_c, white_c, gain_c]
    if cm_name == "RDR_cmap1":
        color_map = [loss_c, white_c, gaina_c, gain_c]
    ## BAF cmap
    if cm_name == "BAF_cmap":
        color_map = [lohb_c, white_c, loha_c]
    if cm_name == "BAF_cmap1":
        color_map = [lohb_c, gaina_c, white_c, gainb_c,loha_c]
    if cm_name == "BAF_cmap_continuous":
        color_map = palettable.colorbrewer.diverging.PuOr_3.mpl_colormap.reversed()
    ## Combine cmap
    if cm_name == "combine_cmap":
        color_map = [loss_c, loh_c, white_c, gain_c]
    ### merge loh
    if cm_name == "combine_cmap2":
        color_map = [lossa_c, lossb_c, loh_c, white_c, gain_c]
    if cm_name == "combine_cmap3":
        color_map = [lossa_c, lossb_c, loha_c, lohb_c, white_c, gain_c]
    ### merge loss
    if cm_name == "combine_cmap4":
        color_map = [loss_c, loha_c, lohb_c, white_c, gain_c]
    if cm_name == "combine_cmap5":
        color_map = [lossa_c, lossb_c, white_c, gain_c]
    
    ## will be deprecated
    if cm_name == "loh_haplotyped_cmap":
        color_map = [lohb_c, white_c,  loha_c]
    if cm_name == "loh_cmap":
        color_map = [loh_c, white_c]
    ## will be deprecated
    if cm_name == "haplotyped_cmap":
        color_map = [white_c, loss_c, loha_c, lohb_c, gain_c]
    if cm_name == "no_haplotyped_cmap":
        color_map = [white_c, loss_c, loh_c, gain_c]
    ## will be deprecated
    if cm_name == "haplotyped_cmap1":
        color_map = [white_c, lossa_c, lossb_c, loha_c, lohb_c, gain_c]
    if cm_name == "haplotyped_cmap2":
        color_map = [white_c, lossa_c, lossb_c, loh_c, gain_c]
    return color_map



## Visualization part in CNVratio estimation performance
## RDR CNV
def CNV_visualization(Xdata, Xlayer = "posterior_mtx", weights = True, 
                      states_weight = np.array([1,2,3]), states_num = 3,
                      colorbar_name = "RDR states", 
                      cell_anno_key = "cell_type", **kwargs):
    """
    Visualize Copy Number Alterations (CNA) for the provided data.

    This function visualizes inferred CNA states based on the Read Depth Ratio (RDR). 
    By default, it uses the `.layers["posterior_mtx"]` for visualization. The 
    RDR states typically represent copy loss, copy neutral, and copy gain.

    Parameters
    ----------

        Xdata : anndata.AnnData
            The annotated data matrix containing the CNA probability inferred from RDR.
        Xlayer : str, optional
            The layer in `Xdata` to be used for visualization. Default is "posterior_mtx".
        weights : bool, optional
            If True, use weighted visualization. 
            If False, use category visualization. Default is True.
        states_weight : numpy.ndarray, optional
            The weights for the CNA states. Default is `np.array([1, 2, 3])`.
        states_num : int, optional
            The number of CNA states. Default is 3 in RDR module.
        colorbar_name : str, optional
            The name to be displayed on the colorbar. Default is "RDR states".
        cell_anno_key : str, optional
            The key for cell annotations used to reorder cells for visualization. 
            Default is "cell_type".
        **kwargs : dict
            Additional keyword arguments passed to the `Xheatmap` function.

    Returns
    -------

        None

    Example
    -------

        .. code-block:: python

            import xclone
            import anndata
            import numpy as np

            # Visualize CNA with default settings
            xclone.pl.CNV_visualization(Xdata)

            # Visualize CNA with custom settings
            xclone.pl.CNV_visualization(Xdata, Xlayer="custom_layer", 
                                        colorbar_name="Custom RDR states", 
                                        cell_anno_key="custom_cell_type")


    """
    if weights:
        ## transfer data to anndata for visualization
        res_cnv_weights_ad = convert_res_to_ann(Xdata, Xlayer, 
                                                weights = weights, states_weight = states_weight)
        center_value = 2.0
        ## reorder the cells based on annotation
        res_cnv_weights_ad1_re = reorder_data_by_cellanno(res_cnv_weights_ad, 
                                                          cell_anno_key = cell_anno_key)
        ## plot the Xheatmap
        Xheatmap(res_cnv_weights_ad1_re, cell_anno_key = cell_anno_key, center = center_value,
        colorbar_ticks = [1.1,2,3],
        colorbar_name = colorbar_name, 
        colorbar_label = ["copy loss", "copy neutral", "copy gain"], **kwargs)

        del res_cnv_weights_ad
        del res_cnv_weights_ad1_re
        gc.collect()

    else:
        # category visualization
        ## transfer data to anndata for visualization
        res_cnv_ad = convert_res_to_ann(Xdata, Xlayer)
        center_value = float(states_num//2)
        ## reorder the cells based on annotation
        res_cnv_ad_re = reorder_data_by_cellanno(res_cnv_ad, cell_anno_key=cell_anno_key)
        ## plot the Xheatmap
        color_map = color_mapping(cm_name = "RDR_cmap")
        Xheatmap(res_cnv_ad_re, cell_anno_key = cell_anno_key, center=center_value, 
        colorbar_ticks = [0,1,2],
        colorbar_name = colorbar_name, 
        colorbar_label = ["copy loss",  "copy neutral", "copy gain"],
        cmap = color_map, **kwargs)

        del res_cnv_ad
        del res_cnv_ad_re
        gc.collect()

    return None

def Complex_CNV_visualization(Xdata, Xlayer = "posterior_mtx", weights = True, 
                      states_weight = np.array([1,2,3]), states_num = 3,
                      cell_anno_key = ["cluster", "cell_type"],
                      colorbar_name = "RDR states", 
                      clusters_display_name = ["Clone", "Celltype"], **kwargs):
    """
    Visualize complex Copy Number Alterations (CNA) for the provided data with multiple annotations.

    This function visualizes CNA based on the Read Depth Ratio (RDR) states, using multiple cell annotations.
    By default, it uses the `.layers["posterior_mtx"]` for visualization. The RDR states typically represent 
    copy loss, copy neutral, and copy gain.

    Parameters
    ----------

        Xdata : anndata.AnnData
            The annotated data matrix containing the CNA probability inferred from RDR.
        Xlayer : str, optional
            The layer in `Xdata` to be used for visualization. Default is "posterior_mtx".
        weights : bool, optional
            If True, use weighted visualization. 
            If False, use category visualization. Default is True.
        states_weight : numpy.ndarray, optional
            The weights for the CNA states. Default is `np.array([1, 2, 3])`.
        states_num : int, optional
            The number of CNA states. Default is 3.
        cell_anno_key : list of str, optional
            The keys for cell annotations used to reorder cells for visualization. 
            Default is ["cluster", "cell_type"].
        colorbar_name : str, optional
            The name to be displayed on the colorbar. Default is "RDR states".
        clusters_display_name : list of str, optional
            The display names for the clusters annotation in the visualization. 
            Default is ["Clone", "Celltype"].
        **kwargs : dict
            Additional keyword arguments passed to the `XXheatmap` function.

    Returns
    -------

        None

    Example
    -------

        .. code-block:: python

            import xclone
            # Visualize complex CNV with default settings
            xclone.pl.Complex_CNV_visualization(Xdata)

            # Visualize complex CNV with custom settings
            xclone.pl.Complex_CNV_visualization(Xdata, Xlayer="custom_layer", weights=False, 
                                                cell_anno_key=["custom_anno1", "custom_anno2"],
                                                colorbar_name="Custom RDR states", 
                                                clusters_display_name=["Custom Clone", "Custom Celltype"])

                                                
    """
    if weights:
        ## transfer data to anndata for visualization
        res_cnv_weights_ad = convert_res_to_ann(Xdata, Xlayer, 
                                                weights = weights, states_weight = states_weight)
        center_value = 2.0
        ## reorder the cells based on annotation
        res_cnv_weights_ad1_re = reorder_data_by_cellanno(res_cnv_weights_ad, 
                                                        cell_anno_key=cell_anno_key)
        ## plot the Xheatmap
        XXheatmap(res_cnv_weights_ad1_re, cell_anno_key = cell_anno_key, clusters_display_name = clusters_display_name,
        center = center_value,
        colorbar_ticks = [1,2,3],
        colorbar_name = colorbar_name, 
        colorbar_label = ["copy loss", "copy neutral", "copy gain"],**kwargs)

        del res_cnv_weights_ad
        del res_cnv_weights_ad1_re
        gc.collect()

    else:
        # category visualization
        ## transfer data to anndata for visualization
        res_cnv_ad = convert_res_to_ann(Xdata, Xlayer)
        center_value = float(states_num//2)
        ## reorder the cells based on annotation
        res_cnv_ad_re = reorder_data_by_cellanno(res_cnv_ad, cell_anno_key=cell_anno_key)
        ## plot the Xheatmap
        color_map = color_mapping(cm_name = "RDR_cmap")
        XXheatmap(res_cnv_ad_re, cell_anno_key = cell_anno_key, clusters_display_name = clusters_display_name,
        center=center_value, 
        colorbar_ticks = [0,1,2],
        colorbar_name = colorbar_name, 
        colorbar_label = ["copy loss",  "copy neutral", "copy gain"],
        cmap = color_map, **kwargs)

        del res_cnv_ad
        del res_cnv_ad_re
        gc.collect()
    return None

## BAF CNV 
def BAF_CNV_visualization(Xdata, Xlayer = "posterior_mtx", weights = False, 
                          states_weight = np.array([1,2,3]),
                          colorbar_name = "BAF states", 
                          cell_anno_key = "cell_type", **kwargs):
    """
    Visualize Allele bias state inferred from B-allele frequency (BAF) for the provided data.

    By default, it uses the `.layers["posterior_mtx"]` for visualization. The BAF states typically represent 
    "allele-A bias", "allele balance", and "allele-B bias" (copy loss (haplotype specific) and copy neutral). 
    With higher resolution, it can represent "allele-A bias (++)", "allele-A bias (+)", "allele balance",  
    "allele-B bias (+)" and "allele-B bias (++)".
    Allele bias 5 states are only supported for category visualization.
    State number to be displayed are detected automatically from the provided `Xlayer`.


    Parameters
    ----------
    
        Xdata : anndata.AnnData
            The annotated data matrix containing allele bia probability inferred from BAF.
        Xlayer : str, optional
            The layer in `Xdata` to be used for visualization. Default is "posterior_mtx".
        weights : bool, optional
            If True, use weighted visualization. 
            If False, use category visualization. Default is False.
        states_weight : numpy.ndarray, optional
            The weights for the allele bia states. Default is `np.array([1, 2, 3])`.
        colorbar_name : str, optional
            The name to be displayed on the colorbar. Default is "BAF states".
        cell_anno_key : str, optional
            The key for the cell annotation used to reorder cells for visualization. 
            Default is "cell_type".
        **kwargs : dict
            Additional keyword arguments passed to the `Xheatmap` function.

    Returns
    -------

        None

    Example
    -------

        .. code-block:: python

            import xclone
            # Visualize BAF CNV with default settings
            xclone.pl.BAF_CNV_visualization(Xdata)

            # Visualize BAF CNV with custom settings
            xclone.pl.BAF_CNV_visualization(Xdata, Xlayer="custom_layer", weights=True, 
                                            cell_anno_key="custom_anno",
                                            colorbar_name="Custom BAF states")


    """
    ## transfer data to anndata for visualization
    states_num = Xdata.layers[Xlayer].shape[-1]
    if weights:
        ## transfer data to anndata for visualization
        res_cnv_weights_ad = convert_res_to_ann(Xdata, Xlayer, weights = weights, states_weight = states_weight)
        center_value = 2.0
        ## reorder the cells based on annotation
        res_cnv_weights_ad1_re = reorder_data_by_cellanno(res_cnv_weights_ad, cell_anno_key=cell_anno_key)
        ## plot the Xheatmap
        color_map = color_mapping(cm_name = "BAF_cmap_continuous")
        Xheatmap(res_cnv_weights_ad1_re, cell_anno_key = cell_anno_key, center = center_value,
        colorbar_ticks = [1, 2, 3], 
        colorbar_name = colorbar_name, 
        colorbar_label = ["allele-A bias",  "allele balance",  "allele-B bias"],
        cmap = color_map, **kwargs)

        del res_cnv_weights_ad
        del res_cnv_weights_ad1_re
        gc.collect()

    else:
        # category visualization
        ## transfer data to anndata for visualization
        res_cnv_ad = convert_res_to_ann(Xdata, Xlayer)
        center_value = float(states_num//2)
        ## reorder the cells based on annotation
        res_cnv_ad_re = reorder_data_by_cellanno(res_cnv_ad, cell_anno_key=cell_anno_key)
        ## plot the Xheatmap
        if states_num == 3:
            color_map = color_mapping(cm_name = "BAF_cmap")
            Xheatmap(res_cnv_ad_re, cell_anno_key = cell_anno_key, center = center_value, 
            colorbar_ticks = [0.3, 1, 1.7], 
            colorbar_name = colorbar_name, 
            colorbar_label = ["allele-A bias",  "allele balance",  "allele-B bias"],
            cmap = color_map, **kwargs)
            # colorbar_ticks = [0, 1, 2]
        elif states_num == 5:
            color_map = color_mapping(cm_name = "BAF_cmap1")
            Xheatmap(res_cnv_ad_re, cell_anno_key = cell_anno_key, center = center_value, 
            colorbar_ticks = [0.3, 1, 1.9, 2.8, 3.5], 
            colorbar_name = colorbar_name, 
            colorbar_label = ["allele-A bias (++)", "allele-A bias (+)", "allele balance",  "allele-B bias (+)", "allele-B bias (++)"],
            cmap = color_map, **kwargs)
        
        del res_cnv_ad
        del res_cnv_ad_re
        gc.collect()

    return None


def Complex_BAF_CNV_visualization(Xdata, Xlayer = "posterior_mtx", weights = False, 
                       states_weight = np.array([1,2,3]),
                       colorbar_name = "BAF states", 
                       cell_anno_key = ["cluster", "cell_type"],
                       clusters_display_name = ["Clone", "Celltype"], 
                       **kwargs):
    """
    Visualize Complex Allele bias states inferred from B-allele frequency (BAF) for the provided data with multiple annotations.

    By default, it uses the `.layers["posterior_mtx"]` for visualization. The BAF states typically represent
    "allele-A bias", "allele balance", and "allele-B bias" for 3 states, and can represent "allele-A bias (++)",
    "allele-A bias (+)", "allele balance", "allele-B bias (+)", and "allele-B bias (++)" for 5 states.
    State number to be displayed are detected automatically from the provided `Xlayer`.

    Parameters
    ----------

        Xdata : anndata.AnnData
            The annotated data matrix containing allele bias probability inferred from BAF.
        Xlayer : str, optional
            The layer in `Xdata` to be used for visualization. Default is "posterior_mtx".
        weights : bool, optional
            If True, use weighted visualization. 
            If False, use category visualization. Default is False.
        states_weight : numpy.ndarray, optional
            The weights for the allele bias states. Default is `np.array([1, 2, 3])`.
        colorbar_name : str, optional
            The name to be displayed on the colorbar. Default is "BAF states".
        cell_anno_key : list of str, optional
            The keys for the cell annotations used to reorder cells for visualization. Default is ["cluster", "cell_type"].
        clusters_display_name : list of str, optional
            The display names for the clusters. Default is ["Clone", "Celltype"].
        **kwargs : dict
            Additional keyword arguments passed to the `XXheatmap` function.

    Returns
    -------

        None

    Example
    -------

        .. code-block:: python

            import xclone
            # Visualize complex BAF CNV with default settings
            xclone.pl.Complex_BAF_CNV_visualization(Xdata)

            # Visualize complex BAF CNV with custom settings
            xclone.pl.Complex_BAF_CNV_visualization(Xdata, Xlayer="custom_layer", weights=True, 
                                                    cell_anno_key=["custom_cluster", "custom_cell_type"],
                                                    clusters_display_name=["Custom Clone", "Custom Celltype"],
                                                    colorbar_name="Custom BAF states")

 
    """
    ## transfer data to anndata for visualization
    states_num = Xdata.layers[Xlayer].shape[-1]
    if weights:
        ## transfer data to anndata for visualization
        res_cnv_weights_ad = convert_res_to_ann(Xdata, Xlayer, weights = weights, states_weight = states_weight)
        center_value = 2.0
        ## reorder the cells based on annotation
        res_cnv_weights_ad1_re = reorder_data_by_cellanno(res_cnv_weights_ad, cell_anno_key=cell_anno_key)
        ## plot the Xheatmap
        color_map = color_mapping(cm_name = "BAF_cmap_continuous")
        XXheatmap(res_cnv_weights_ad1_re, cell_anno_key = cell_anno_key, clusters_display_name = clusters_display_name,
        center = center_value,
        colorbar_ticks = [1, 2, 3], 
        colorbar_name = colorbar_name, 
        colorbar_label = ["allele-A bias",  "allele balance",  "allele-B bias"],
        cmap = color_map, **kwargs)

        del res_cnv_weights_ad
        del res_cnv_weights_ad1_re
        gc.collect()

    else:
        # category visualization
        ## transfer data to anndata for visualization
        res_cnv_ad = convert_res_to_ann(Xdata, Xlayer)
        center_value = float(states_num//2)
        ## reorder the cells based on annotation
        res_cnv_ad_re = reorder_data_by_cellanno(res_cnv_ad, cell_anno_key=cell_anno_key)
        ## plot the Xheatmap
        if states_num == 3:
            color_map = color_mapping(cm_name = "BAF_cmap")
            XXheatmap(res_cnv_ad_re, cell_anno_key = cell_anno_key, clusters_display_name = clusters_display_name,
            center = center_value, 
            colorbar_ticks = [0, 1, 2], 
            colorbar_name = colorbar_name, 
            colorbar_label = ["allele-A bias",  "allele balance",  "allele-B bias"],
            cmap = color_map, **kwargs)
        elif states_num == 5:
            color_map = color_mapping(cm_name = "BAF_cmap1")
            XXheatmap(res_cnv_ad_re, cell_anno_key = cell_anno_key, clusters_display_name = clusters_display_name,
            center = center_value, 
            colorbar_ticks = [0.3, 1, 1.9, 2.8, 3.5], 
            colorbar_name = colorbar_name, 
            colorbar_label = ["allele-A bias (++)", "allele-A bias (+)", "allele balance",  "allele-B bias (+)", "allele-B bias (++)"],
            cmap = color_map, **kwargs)
        
        del res_cnv_ad
        del res_cnv_ad_re
        gc.collect()

    return None
        

## Visualization part in CNVratio estimation performance
## combined CNV
def Combine_CNV_visualization(Xdata,
                              Xlayer,
                              states_num = 4,
                              cell_anno_key = "cell_type", 
                              color_map_name = None,
                              colorbar_ticks = None,
                              colorbar_label = None,
                              **kwargs):
    """
    Visualize precise CNA states inferred from combination of RDR and BAF information.

    By default, it visualizes 4 states: copy loss, loss of heterozygosity (loh), copy neutral, and copy gain.

    Parameters
    ----------

        Xdata : anndata.AnnData
            The annotated data matrix containing CNA states.
        Xlayer : str
            The layer in `Xdata` to be used for visualization.
        states_num : int, optional
            The number of CNA states to visualize. Default is 4.
        cell_anno_key : str, optional
            The key for the cell annotation used to reorder cells for visualization. Default is "cell_type".
        color_map_name : str, optional
            The name of the color map to use for visualization. Default is None, which uses "combine_cmap".
        colorbar_ticks : list of int, optional
            The positions of the ticks on the colorbar. Default is None, which uses [0, 1, 2, 3].
        colorbar_label : list of str, optional
            The labels for the ticks on the colorbar. Default is None, which uses ["copy loss", "loh", "copy neutral", "copy gain"].
        **kwargs : dict
            Additional keyword arguments passed to the `Xheatmap` function.

    Returns
    -------

        None

    Example
    -------

        .. code-block:: python

            import xclone

            # Visualize combined CNA with custom settings
            colorbar_ticks = [0.25,1,2,2.75]
            colorbar_label = ["copy loss","loh", "copy neutral", "copy gain"]
            xclone.pl.Combine_CNV_visualization(combine_Xdata, Xlayer = "plot_prob_merge1", 
                                        cell_anno_key = plot_cell_anno_key, 
                                        color_map_name = "combine_cmap", 
                                        states_num = 4, 
                                        colorbar_ticks = colorbar_ticks,
                                        colorbar_label = colorbar_label,
                                        title = fig_title,
                                        save_file = True, 
                                        out_file = combine_res_select_fig,
                                        **kwargs)

    """
    res_cnv_ad = convert_res_to_ann(Xdata, Xlayer)
    ## reorder the cells based on annotation
    res_cnv_ad_re = reorder_data_by_cellanno(res_cnv_ad, cell_anno_key = cell_anno_key)

    ## plot the Xheatmap
    center_value = float((states_num-1)/2)
    if color_map_name is None:
        color_map = color_mapping(cm_name = "combine_cmap")
    else:
        color_map = color_mapping(cm_name = color_map_name)

    if colorbar_ticks is None:
        colorbar_ticks = [0,1,2,3]
    if colorbar_label is None:
        colorbar_label = ["copy loss", "loh", "copy neutral", "copy gain"]

    Xheatmap(res_cnv_ad_re, cell_anno_key = cell_anno_key, center=center_value, 
    colorbar_ticks = colorbar_ticks,
    colorbar_label = colorbar_label,
    cmap = color_map, **kwargs)

    del res_cnv_ad
    del res_cnv_ad_re
    gc.collect()

    return None


def Complex_Combine_CNV_visualization(Xdata, 
                                      Xlayer,
                                      states_num = 4,
                                      cell_anno_key = ["cluster", "cell_type"], 
                                      clusters_display_name = ["Clone", "Celltype"], 
                                      color_map_name = None,
                                      colorbar_ticks = None,
                                      colorbar_label = None,
                                      **kwargs):
    """
    Adapted from `Combine_CNV_visualization` to visualize multiple annotations.

    Parameters
    ----------

        Xdata : anndata.AnnData
            The annotated data matrix containing CNA states.
        Xlayer : str
            The layer in `Xdata` to be used for visualization.
        states_num : int, optional
            The number of CNA states to visualize. Default is 4.
        cell_anno_key : list of str, optional
            The keys for the cell annotations used to reorder cells for visualization. Default is ["cluster", "cell_type"].
        clusters_display_name : list of str, optional
            The display names for the clusters. Default is ["Clone", "Celltype"].
        color_map_name : str, optional
            The name of the color map to use for visualization. Default is None, which uses "combine_cmap".
        colorbar_ticks : list of int, optional
            The positions of the ticks on the colorbar. Default is None, which uses [0, 1, 2, 3].
        colorbar_label : list of str, optional
            The labels for the ticks on the colorbar. Default is None, which uses ["copy loss", "loh", "copy neutral", "copy gain"].
        **kwargs : dict
            Additional keyword arguments passed to the `XXheatmap` function.

    Returns
    -------

        None

    Example
    -------

        .. code-block:: python

            import xclone

            # Visualize combined CNA with custom settings
            colorbar_ticks = [0.25,1,2,2.75]
            colorbar_label = ["copy loss","loh", "copy neutral", "copy gain"]
            xclone.pl.Complex_Combine_CNV_visualization(combine_Xdata, Xlayer = "plot_prob_merge1", 
                                        cell_anno_key = plot_cell_anno_key, 
                                        color_map_name = "combine_cmap", 
                                        states_num = 4, 
                                        cell_anno_key = ["cluster", "cell_type"], 
                                        clusters_display_name = ["Clone", "Celltype"], 
                                        colorbar_ticks = colorbar_ticks,
                                        colorbar_label = colorbar_label,
                                        title = fig_title,
                                        save_file = True, 
                                        out_file = combine_res_select_fig,
                                        **kwargs)




    """
    res_cnv_ad = convert_res_to_ann(Xdata, Xlayer)
    ## reorder the cells based on annotation
    res_cnv_ad_re = reorder_data_by_cellanno(res_cnv_ad, cell_anno_key=cell_anno_key)

    ## plot the Xheatmap
    center_value = float((states_num-1)/2)

    if color_map_name is None:
        color_map = color_mapping(cm_name = "combine_cmap")
    else:
        color_map = color_mapping(cm_name = color_map_name)
    
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

    del res_cnv_ad
    del res_cnv_ad_re
    gc.collect()
    
    return None