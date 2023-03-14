"""Base pipeline for XClone combination module"""

# Author: Rongting Huang
# Date: 2022-11-17
# update: 2022-11-17

import os
import xclone
import numpy as np
from .._logging import get_logger
from datetime import datetime, timezone


def run_combine(RDR_Xdata,
                BAF_merge_Xdata,
                verbose = True,
                run_verbose = True,
                config_file = None):
    """
    """
    ## settings
    from .._config import XCloneConfig
    
    if config_file == None:
        print (
            f'Model configuration file not specified.\n'
            f'Default settings in XClone-combine will be used.'
        )
        config = XCloneConfig(module = "Combine")

    else:
        config = config_file
    ## base settings
    warninig_ignore = config.warninig_ignore
    if warninig_ignore:
        import warnings
        warnings.filterwarnings('ignore')
    ## general settings
    out_dir = config.outdir
    dataset_name = config.dataset_name
    
    cell_anno_key = config.cell_anno_key
    exclude_XY = config.exclude_XY
    
    ## RDR BAF settings
    BAF_denoise = config.BAF_denoise
    ## combine settings
    copyloss_correct = config.copyloss_correct
    copyloss_correct_mode = config.copyloss_correct_mode
    copygain_correct= config.copygain_correct
    copygain_correct_mode = config.copygain_correct_mode
    
    ## plot settings
    xclone_plot = config.xclone_plot
    plot_cell_anno_key = config.plot_cell_anno_key
    merge_loss = config.merge_loss
    merge_loh = config.merge_loh
    
    ## Result output prepare
    if out_dir is None:
        cwd = os.getcwd()
        out_dir = cwd + "/XCLONE_OUT/"
    
    out_data_dir = str(out_dir) + "/data/"
    xclone.al.dir_make(out_data_dir)

    RDR_combine_corrected_file = out_data_dir + "combine_adata_corrected.h5ad"
    combine_final_file = out_data_dir + "combined_final.h5ad"
    
    ##------------------------
    main_logger = get_logger("Main combine module")
    main_logger.info("XClone combine module Started")
    start_time = datetime.now(timezone.utc)

    if run_verbose:
        print("[XClone Combination module running]************************")
    
    if exclude_XY:
        RDR_Xdata = xclone.pp.exclude_XY_adata(RDR_Xdata)
        BAF_merge_Xdata = xclone.pp.exclude_XY_adata(BAF_merge_Xdata)

        print("[XClone warning] Combine module excelude chr XY analysis.")

    # map BAF to RDR
    if BAF_denoise:
        if "denoised_posterior_mtx" in BAF_merge_Xdata.layers:
            BAF_use_Xlayer = "denoised_posterior_mtx"
        else:
            print("[XClone] warning: please proveide denoised BAF posterior_mtx layer")
            BAF_use_Xlayer = "posterior_mtx"
    else:
        BAF_use_Xlayer = "posterior_mtx"
    
    combine_Xdata = xclone.model.bin_to_gene_mapping(BAF_merge_Xdata,
                        RDR_Xdata,
                        Xlayer = BAF_use_Xlayer,
                        extend_layer = "BAF_extend_post_prob",
                        return_prob = False)
    
    combine_Xdata = xclone.model.CNV_prob_combination(combine_Xdata,
                         RDR_layer = "posterior_mtx",
                         BAF_layer = "BAF_extend_post_prob",
                         copyloss_correct = copyloss_correct,
                         copyloss_correct_mode = copyloss_correct_mode,
                         copygain_correct = copygain_correct,
                         copygain_correct_mode = copygain_correct_mode)
    try:
        combine_Xdata.write(RDR_combine_corrected_file)
    except Exception as e:
        print("[XClone Warning]", e)
    else:
        print("[XClone hint] combine_corrected_file saved in %s." %(out_data_dir))
    

    combine_Xdata = xclone.model.CNV_prob_merge_for_plot(combine_Xdata, Xlayer = "corrected_prob")
    try:
        combine_Xdata.write(combine_final_file)
    except Exception as e:
        print("[XClone Warning]", e)
    else:
        print("[XClone hint] combine_final_file saved in %s." %(out_data_dir))
    

    end_time = datetime.now(timezone.utc)
    time_passed = end_time - start_time
    main_logger.info("XClone combine module finished (%d seconds)" % (time_passed.total_seconds()))

    if xclone_plot:
        set_figtitle = config.set_figtitle
        if run_verbose:
            print("[XClone plotting]")
        if plot_cell_anno_key is None:
            plot_cell_anno_key = cell_anno_key
        run_combine_plot(combine_Xdata, dataset_name, 
                         plot_cell_anno_key, 
                         merge_loss, merge_loh, 
                         set_figtitle,
                         out_dir)
    
    return combine_Xdata


def run_combine_plot(combine_Xdata,
            dataset_name,
            plot_cell_anno_key,
            merge_loss = True,
            merge_loh = True,
            set_figtitle = True,
            out_dir = None,
            **kwargs):
    """
    """

    ## Result output prepare
    if out_dir is None:
        cwd = os.getcwd()
        out_dir = cwd + "/XCLONE_OUT/"
    
    out_plot_dir = str(out_dir) + "/plot/"
    xclone.al.dir_make(out_plot_dir)

    fig_title = ""
    combine_res_base_fig = out_plot_dir + dataset_name + "_combine_base.png"
    combine_res_select_fig = out_plot_dir + dataset_name + "_combine_select.png"

    sub_logger = get_logger("Combine plot module")
    sub_logger.info("Combine plot module started")
    start_time = datetime.now(timezone.utc)
    if set_figtitle:
        fig_title = dataset_name + " XClone Combine module"


    ## BASE PLOT
    xclone.pl.Combine_CNV_visualization(combine_Xdata, Xlayer = "prob1_merge", 
        cell_anno_key = plot_cell_anno_key,  
        title = fig_title,
        save_file = True, out_file = combine_res_base_fig,
        **kwargs)
    
    ## SELECT PLOT
    if merge_loh:
        if merge_loss:
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
        else:
            colorbar_ticks = [0,1,2,3,4]
            colorbar_label = ["copy lossA", "copy lossB", "LOH", "copy neutral", "copy gain"]
            xclone.pl.Combine_CNV_visualization(combine_Xdata, Xlayer = "plot_prob_merge2", 
                                        cell_anno_key = plot_cell_anno_key, 
                                        color_map_name = "combine_cmap2", 
                                        states_num = 5,
                                        colorbar_ticks = colorbar_ticks,
                                        colorbar_label = colorbar_label,
                                        title = fig_title,
                                        save_file = True, 
                                        out_file = combine_res_select_fig,
                                        **kwargs)
    elif merge_loss:
        colorbar_ticks = [0,1,2,3,4]
        colorbar_label = ["copy loss","LOH-A", "LOH-B",  "copy neutral", "copy gain"]
        xclone.pl.Combine_CNV_visualization(combine_Xdata, Xlayer = "plot_prob_merge4", 
                                        cell_anno_key = plot_cell_anno_key, 
                                        color_map_name = "combine_cmap4", 
                                        states_num = 5,
                                        colorbar_ticks = colorbar_ticks,
                                        colorbar_label = colorbar_label,
                                        title = fig_title,
                                        save_file = True, 
                                        out_file = combine_res_select_fig,
                                        **kwargs)
        
    else:
        colorbar_ticks = [0,1,2,3,4,5]
        colorbar_label = ["copy lossA", "copy lossB","LOH-A", "LOH-B", "copy neutral", "copy gain"]
        xclone.pl.Combine_CNV_visualization(combine_Xdata, Xlayer = "plot_prob_merge3", 
                                        cell_anno_key = plot_cell_anno_key, 
                                        color_map_name = "combine_cmap3", 
                                        states_num = 6,
                                        colorbar_ticks = colorbar_ticks,
                                        colorbar_label = colorbar_label,
                                        title = fig_title,
                                        save_file = True, 
                                        out_file = combine_res_select_fig,
                                        **kwargs)
    
    end_time = datetime.now(timezone.utc)
    time_passed = end_time - start_time
    sub_logger.info("Combine plot module finished (%d seconds)" % (time_passed.total_seconds()))
    return None

def plot_Combine_module():
    """
    all plots for Combine module.
    """
    pass