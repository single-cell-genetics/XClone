

def CNV_combine_visualization_xxx(Xdata, Xlayer, haplotype = False):
    """
    DEPRECATED 20220805
    Function:
    ---------

    Parameters:
    ------------
    Xdata: anndata.
    Xlayer: posterior probability for visualization.
    haplotype: bool, False: copy loss, neutral, copy gain
    haplotype: bool, True: copy loss1, copy loss2, neutral, copy gain


    todo:
    1) for visualization, need transfer to non log version. 
    2) select color mapping.
    """
    import palettable
    if haplotype:
        ## 4 color
        # category visualization
        res_cnv_ad = convert_res_to_ann(Xdata, Xlayer)
        center_value = float(2)
        res_cnv_ad_re = reorder_data_by_cellanno(res_cnv_ad, 
                                             cell_anno_key="cell_type")
        white_color = palettable.colorbrewer.diverging.PuOr_3.mpl_colors[1]
        loss1_color = palettable.colorbrewer.diverging.PuOr_3.mpl_colors[0]
        loss2_color = palettable.colorbrewer.diverging.PuOr_3.mpl_colors[2]
        gain_color = palettable.colorbrewer.diverging.RdBu_4.mpl_colors[0]
        color_lst = [loss1_color, loss2_color, white_color, gain_color, gain_color]
        # Xheatmap(res_cnv_ad_re, cell_anno_key = 'cell_type', 
            # center=center_value, cmap = palettable.colorbrewer.diverging.PuOr_5.mpl_colors)
        Xheatmap(res_cnv_ad_re, cell_anno_key = 'cell_type', center=center_value, 
        colorbar_ticks = [0, 1, 2, 3], colorbar_label = ["copy loss", "copy loss", "copy neutral", "copy gain"],
        cmap = color_lst)
    else:
        ## 3 color
        res_cnv_ad = convert_res_to_ann(Xdata, Xlayer)
        center_value = float(1)
        res_cnv_ad_re = reorder_data_by_cellanno(res_cnv_ad, 
                                             cell_anno_key="cell_type")
                                             
        Xheatmap(res_cnv_ad_re, fillna_value = center_value, cell_anno_key = 'cell_type', center=center_value,
        colorbar_ticks = [0, 1, 2], colorbar_label = ["copy loss", "", "copy neutral", "", "copy gain"],
        cmap = palettable.colorbrewer.diverging.RdBu_3.mpl_colors)