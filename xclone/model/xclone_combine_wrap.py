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

    out_dir = config.outdir
    dataset_name = config.dataset_name
    cell_anno_key = config.cell_anno_key

    xclone_plot = config.xclone_plot
    plot_cell_anno_key = config.plot_cell_anno_key
    merge_loss = config.merge_loss
    merge_loh = config.merge_loh
    
    ## Result output prepare
    if out_dir is None:
        cwd = os.getcwd()
        out_dir = cwd + "/XCLONE_OUT/"
    
    out_data_dir = str(out_dir) + "data/"
    xclone.al.dir_make(out_data_dir)


    ## use_file
    # output after CNV calling
    # BAF_final_file = out_dir + "BAF_merge_Xdata_KNN_HMM_post.h5ad"
    # ## use RDR connectivities
    # RDR_final_file = out_dir + "RDR_adata_KNN_HMM_post.h5ad"

    RDR_combine_corrected_file = out_data_dir + "combine_adata_corrected.h5ad"
    combine_final_file = out_data_dir + "combined_final.h5ad"
    
    ##------------------------
    main_logger = get_logger("Main combine module")
    main_logger.info("XClone combine module Started")
    start_time = datetime.now(timezone.utc)

    if run_verbose:
        print("[XClone Combination module running]************************")

    # map BAF to RDR
    combine_Xdata = xclone.model.bin_to_gene_mapping(BAF_merge_Xdata,
                        RDR_Xdata,
                        Xlayer = "posterior_mtx",
                        extend_layer = "BAF_extend_post_prob",
                        return_prob = False)
    
    combine_Xdata = xclone.model.CNV_prob_combination(combine_Xdata,
                         RDR_layer = "posterior_mtx",
                         BAF_layer = "BAF_extend_post_prob")

    combine_Xdata.write(RDR_combine_corrected_file)

    combine_Xdata = xclone.model.CNV_prob_merge_for_plot(combine_Xdata, Xlayer = "loss_corrected_prob1")
    combine_Xdata.write(combine_final_file)

    end_time = datetime.now(timezone.utc)
    time_passed = end_time - start_time
    main_logger.info("XClone combine module finished (%d seconds)" % (time_passed.total_seconds()))

    if xclone_plot:
        if run_verbose:
            print("[XClone plotting]")
        if plot_cell_anno_key is None:
            plot_cell_anno_key = cell_anno_key
        run_combine_plot(combine_Xdata, dataset_name, plot_cell_anno_key, merge_loss, merge_loh, out_dir)
    
    return combine_Xdata


def run_combine_plot(combine_Xdata,
            dataset_name,
            plot_cell_anno_key,
            merge_loss = True,
            merge_loh = True,
            out_dir = None):
    """
    """

    ## Result output prepare
    if out_dir is None:
        cwd = os.getcwd()
        out_dir = cwd + "/XCLONE_OUT/"
    
    out_plot_dir = str(out_dir) + "plot/"
    xclone.al.dir_make(out_plot_dir)


    combine_res_base_fig = out_plot_dir + dataset_name + "combine_base.png"
    combine_res_select_fig = out_plot_dir + dataset_name + "combine_select.png"

    sub_logger = get_logger("Combine plot module")
    sub_logger.info("Combine plot module started")
    start_time = datetime.now(timezone.utc)

    ## BASE PLOT
    xclone.pl.CNV_visualization_combine(combine_Xdata, Xlayer = "prob1_merge", 
        cell_anno_key = plot_cell_anno_key,  save_file = True, out_file = combine_res_base_fig)
    
    ## SELECT PLOT
    if merge_loh:
        if merge_loss:
            colorbar_ticks = [0.25,1,2,2.75]
            colorbar_label = ["copy loss","loh", "copy neutral", "copy gain"]
            xclone.pl.CNV_visualization_combine(combine_Xdata, Xlayer = "plot_prob_merge1", 
                                        cell_anno_key = plot_cell_anno_key, 
                                        color_map_name = "combine_cmap", 
                                        states_num = 4, 
                                        colorbar_ticks = colorbar_ticks,
                                        colorbar_label = colorbar_label,
                                        save_file = True, 
                                        out_file = combine_res_select_fig)
        else:
            colorbar_ticks = [0,1,2,3,4]
            colorbar_label = ["copy lossA", "copy lossB", "LOH", "copy neutral", "copy gain"]
            xclone.pl.CNV_visualization_combine(combine_Xdata, Xlayer = "plot_prob_merge2", 
                                        cell_anno_key = plot_cell_anno_key, 
                                        color_map_name = "combine_cmap2", 
                                        states_num = 5,
                                        colorbar_ticks = colorbar_ticks,
                                        colorbar_label = colorbar_label,
                                        save_file = True, 
                                        out_file = combine_res_select_fig)
    elif merge_loss:
        colorbar_ticks = [0,1,2,3,4]
        colorbar_label = ["copy loss","LOH-A", "LOH-B",  "copy neutral", "copy gain"]
        xclone.pl.CNV_visualization_combine(combine_Xdata, Xlayer = "plot_prob_merge4", 
                                        cell_anno_key = plot_cell_anno_key, 
                                        color_map_name = "combine_cmap4", 
                                        states_num = 5,
                                        colorbar_ticks = colorbar_ticks,
                                        colorbar_label = colorbar_label,
                                        save_file = True, 
                                        out_file = combine_res_select_fig)
        
    else:
        colorbar_ticks = [0,1,2,3,4,5]
        colorbar_label = ["copy lossA", "copy lossB","LOH-A", "LOH-B", "copy neutral", "copy gain"]
        xclone.pl.CNV_visualization_combine(combine_Xdata, Xlayer = "plot_prob_merge3", 
                                        cell_anno_key = plot_cell_anno_key, 
                                        color_map_name = "combine_cmap3", 
                                        states_num = 6,
                                        colorbar_ticks = colorbar_ticks,
                                        colorbar_label = colorbar_label,
                                        save_file = True, 
                                        out_file = combine_res_select_fig)
    
    end_time = datetime.now(timezone.utc)
    time_passed = end_time - start_time
    sub_logger.info("Combine plot module finished (%d seconds)" % (time_passed.total_seconds()))
    return None