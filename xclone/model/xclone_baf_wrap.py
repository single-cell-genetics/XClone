"""Base pipeline for XClone BAF module"""

# Author: Rongting Huang
# Date: 2022-10-21
# update: 2023-01-07


import os
import xclone
import anndata as an
import numpy as np
from .._logging import get_logger
from datetime import datetime, timezone

def preview_BAF(BAF_adata,
                cell_anno_key,
                ref_celltype,
                ):
    """
    Function:
    some precheck for params setting.
    """
    # xclone.model.view_celltype(BAF_adata, cell_anno_key = cell_anno_key)

    ## smooth plot?
    return None


def run_BAF(BAF_adata,
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
            f'Default settings in XClone-BAF will be used.'
        )
        config = XCloneConfig(module = "BAF")

    else:
        config = config_file
    
    # base settings
    warninig_ignore = config.warninig_ignore
    if warninig_ignore:
        import warnings
        warnings.filterwarnings('ignore')
    
    # general settings
    dataset_name = config.dataset_name
    out_dir = config.outdir

    cell_anno_key = config.cell_anno_key
    ref_celltype = config.ref_celltype
    exclude_XY = config.exclude_XY

    ## RDR related
    RDR_file = config.RDR_file
    remove_marker_genes = config.remove_marker_genes
    
    # BAF settings
    theo_neutral_BAF = config.theo_neutral_BAF
    ref_BAF_clip = config.ref_BAF_clip
    BAF_states_num = config.CNV_N_components
    CNV_N_components = config.CNV_N_components
    gene_specific_concentration = config.gene_specific_concentration
    concentration = config.concentration
    guide_theo_CNV_states = config.guide_theo_CNV_states

    extreme_count_cap = config.extreme_count_cap
    BAF_add = config.BAF_add
    BAF_denoise = config.BAF_denoise
    GMM_detect = config.BAF_denoise_GMM_detection
    BAF_denoise_GMM_comp = config.BAF_denoise_GMM_comp
    BAF_denoise_cellprop_cutoff = config.BAF_denoise_cellprop_cutoff
    
    ## phasing
    feature_mode = config.feature_mode
    phasing_region_key = config.phasing_region_key
    phasing_len = config.phasing_len
    bin_nproc = config.bin_nproc

    ## smoothing
    WMA_window_size = config.WMA_window_size
    WMA_smooth_key = config.WMA_smooth_key

    # HMM settings
    start_prob = config.start_prob
    trans_t = config.trans_t
    trans_prob = config.trans_prob
    HMM_brk = config.HMM_brk

    # plot settings
    xclone_plot = config.xclone_plot
    plot_cell_anno_key = config.plot_cell_anno_key
    

    ## Result output prepare
    if out_dir is None:
        cwd = os.getcwd()
        out_dir = cwd + "/XCLONE_OUT/"
    out_data_dir = str(out_dir) + "/data/"
    xclone.al.dir_make(out_data_dir)
    
    ### output before CNV calling
    BAF_base_file = out_data_dir + "BAF_base_Xdata.h5ad"
    BAF_merge_base_file = out_data_dir + "BAF_merge_base_Xdata.h5ad"
    ### output after CNV calling
    BAF_final_file = out_data_dir + "BAF_merge_Xdata_KNN_HMM_post.h5ad"
    
    ### RDR file for BAF module (load preprocessed dataset)
    ### use RDR connectivities, marker genes, etc.
    if RDR_file is None:
        RDR_final_file = out_data_dir + "RDR_adata_KNN_HMM_post.h5ad"
    else:
        RDR_final_file = RDR_file

    ##------------------------
    main_logger = get_logger("Main BAF module")
    main_logger.info("XClone BAF module Started")
    start_time = datetime.now(timezone.utc)

    if run_verbose:
        print("[XClone BAF module running]************************")
    
    if exclude_XY:
        BAF_adata = xclone.pp.exclude_XY_adata(BAF_adata)
        print("[XClone warning] BAF module excelude chr XY analysis.")
    BAF_adata = xclone.pp.check_BAF(BAF_adata, cell_anno_key = cell_anno_key, verbose = verbose)

    RDR_adata = an.read_h5ad(RDR_final_file)
    ## check validated RDR and BAF
    xclone.pp.check_RDR_BAF_cellorder(RDR_adata, BAF_adata)
    ## Remove marker genes
    if remove_marker_genes:
        marker_genes = RDR_adata.uns["rank_marker_genes"]
        BAF_adata = xclone.pp.Xdata_region_selection(BAF_adata,
                                                    select = False,
                                                    chr_anno_key = "GeneName",
                                                    chr_lst = marker_genes,
                                                    update_uns = False,
                                                    uns_anno_key = None)
    ## BAF Phasing
    if HMM_brk in ["chr", "chr_arm"]:
        BAF_var_add = None
    else:
        BAF_var_add = HMM_brk
    BAF_adata, merge_Xdata =  xclone.model.BAF_Local_phasing(BAF_adata, 
                                                             region_key = phasing_region_key, 
                                                             phasing_len = phasing_len, 
                                                             bin_nproc = bin_nproc,
                                                             feature_mode = feature_mode,
                                                             var_add = BAF_var_add)
    BAF_adata, merge_Xdata = xclone.model.BAF_Global_phasing(BAF_adata, merge_Xdata)
    
    ### check coverage for bins
    merge_Xdata.var[(merge_Xdata.layers["dp_bin"].A.sum(axis=0) == 0)]
    if extreme_count_cap:
        merge_Xdata = xclone.model.extrme_count_capping(merge_Xdata)

    merge_Xdata = xclone.pl.calculate_cell_BAF(merge_Xdata, 
                                               AD_key = "ad_bin", DP_key = "dp_bin", BAF_key = "BAF")
    merge_Xdata = xclone.pl.calculate_cell_BAF(merge_Xdata, 
                                               AD_key = "ad_bin_phased", DP_key = "dp_bin", BAF_key = "BAF_phased")
    xclone.model.BAF_fillna(merge_Xdata, Xlayer = "BAF_phased", out_layer = "fill_BAF_phased")

    ## smoothing
    merge_Xdata = xclone.model.get_KNN_connectivities_from_expr(merge_Xdata, RDR_adata)

    merge_Xdata = xclone.model.KNN_smooth(merge_Xdata, 
                                          run_KNN = False, 
                                          KNN_Xlayer = None, 
                                          KNN_connect_use = "connectivities_expr",
                                          layer = "fill_BAF_phased", 
                                          out_layer='BAF_phased_KNN')
    merge_Xdata = xclone.model.WMA_smooth(merge_Xdata, 
                                          layer="BAF_phased_KNN", 
                                          out_layer='BAF_phased_KNN_WMA', 
                                          window_size = WMA_window_size,
                                          chrom_key = WMA_smooth_key, 
                                          verbose=False)
    
    merge_Xdata = xclone.model.WMA_smooth(merge_Xdata, 
                                          layer="fill_BAF_phased", 
                                          out_layer='BAF_phased_WMA', 
                                          window_size = WMA_window_size, 
                                          chrom_key = WMA_smooth_key,
                                          verbose=False)
    # merge_Xdata = xclone.model.KNN_smooth(merge_Xdata, 
    #                                       KNN_connect_use = "connectivities_expr", 
    #                                       layer="BAF_phased_WMA", 
    #                                       out_layer='BAF_phased_WMA_KNN')
    ## after phasing & smoothing
    try:
        BAF_adata.write(BAF_base_file)
        merge_Xdata.write(BAF_merge_base_file)
    except Exception as e:
        print("[XClone Warning]", e)
    else:
        print("[XClone hint] BAF_base_file and merged_file saved in %s." %(out_data_dir))
    
    ## HMM smoothing for CNV states calling
    if guide_theo_CNV_states is not None:
        guide_theo_states = guide_theo_CNV_states
        print("User specify guide_theo_states: ", guide_theo_states)
    else:
        CNV_states = xclone.model.get_CNV_states(merge_Xdata, Xlayer = "BAF_phased_WMA",
                                                n_components = CNV_N_components,
                                                means_init = None,
                                                max_iter = 200)
        guide_theo_states = xclone.model.guide_BAF_theo_states(CNV_states)

    ## if you have limited ref cells, you can set theo_baf as 0.5   *important
    if theo_neutral_BAF is not None:
        merge_Xdata.var["theo_neutral_BAF"] = theo_neutral_BAF
        used_specific_states = xclone.model.gene_specific_BAF(merge_Xdata, 
                            theo_states= guide_theo_states, specific_BAF = "theo_neutral_BAF")
    else:
        merge_Xdata = xclone.model.get_BAF_ref(merge_Xdata, 
                                           Xlayer = "fill_BAF_phased", 
                                           out_anno = "ref_BAF_phased",
                                           anno_key = cell_anno_key, 
                                           ref_cell = ref_celltype,
                                           clipping = ref_BAF_clip)

        used_specific_states = xclone.model.gene_specific_BAF(merge_Xdata, 
                               theo_states= guide_theo_states, specific_BAF = "ref_BAF_phased")

    if gene_specific_concentration:
        concentration_lower = config.concentration_lower
        concentration_upper = config.concentration_upper
        merge_Xdata = xclone.model.concentration_mapping(merge_Xdata, concentration_lower, concentration_upper)

    merge_Xdata = xclone.model.calculate_Xemm_prob_bb(merge_Xdata, 
                                                      AD_key = "ad_bin_phased", DP_key = "dp_bin",
                                                      outlayer = "bin_phased_BAF_specific_center_emm_prob_log", 
                                                      states = used_specific_states,
                                                      states_num = BAF_states_num,
                                                      gene_specific = gene_specific_concentration, 
                                                      concentration = concentration)

    merge_Xdata = xclone.model.BAF_smoothing(merge_Xdata,
                                             inlayer = "bin_phased_BAF_specific_center_emm_prob_log",
                                             outlayer = "bin_phased_BAF_specific_center_emm_prob_log_KNN",
                                             KNN_connectivities_key = "connectivities_expr",
                                             KNN_smooth = True)
    t = trans_t
    if trans_prob is None:
        if CNV_N_components == 3:
            trans_prob = np.array([[1-2*t, t, t],[t, 1-2*t, t],[t, t, 1-2*t]])
        if CNV_N_components == 5:
            trans_prob = np.array([[1-4*t, t, t, t,t],[t, 1-4*t, t, t,t],[t, t, 1-4*t, t,t], [t, t, t, 1-4*t, t], [t, t, t, t, 1-4*t]])
    
    merge_Xdata = xclone.model.XHMM_smoothing(merge_Xdata,
                                              brk = HMM_brk, 
                                              start_prob = start_prob,  
                                              trans_prob = trans_prob, 
                                              emm_inlayer = "bin_phased_BAF_specific_center_emm_prob_log_KNN", 
                                              nproc = 80, 
                                              verbose = False)
    
    # merge_Xdata = xclone.model.BAF_smoothing(merge_Xdata,
    #                                          inlayer = "posterior_mtx_log",
    #                                          outlayer = "posterior_mtx_log_KNN",
    #                                          KNN_connectivities_key = "connectivities_expr",
    #                                          KNN_smooth = True)
    
    # merge_Xdata.layers["posterior_mtx"] = np.exp(merge_Xdata.layers["posterior_mtx_log"])
    
    
    ## try to add 3 states layer for comparasion or denoise or correction. [testing version]
    ### only support BAF 5 states to add 3 states visualization for comparasion now.[testing version] 
    if CNV_N_components == 5:
        if BAF_add is None:
            BAF_add = True 
            # default for BAF module 5 states mode to have 3 states for comparasion.
    else:
        BAF_add = False
    
    if BAF_add:
        merge_Xdata_copy = merge_Xdata.copy()
        ## HMM smoothing for CNV states calling
        if guide_theo_CNV_states is not None:
            guide_theo_states = guide_theo_CNV_states[[0,3]]
            print("User specify guide_theo_states: ", guide_theo_states)
        else:
            CNV_states = xclone.model.get_CNV_states(merge_Xdata_copy, Xlayer = "BAF_phased_WMA",
                                                n_components = 3,
                                                means_init = None,
                                                max_iter = 200)
            guide_theo_states = xclone.model.guide_BAF_theo_states(CNV_states)

        ## if you have limited ref cells, you can set theo_baf as 0.5   *important
        if theo_neutral_BAF is not None:
            merge_Xdata_copy.var["theo_neutral_BAF"] = theo_neutral_BAF
            used_specific_states = xclone.model.gene_specific_BAF(merge_Xdata_copy, 
                            theo_states= guide_theo_states, specific_BAF = "theo_neutral_BAF")
        else:
            merge_Xdata_copy = xclone.model.get_BAF_ref(merge_Xdata_copy, 
                                           Xlayer = "fill_BAF_phased", 
                                           out_anno = "ref_BAF_phased",
                                           anno_key = cell_anno_key, 
                                           ref_cell = ref_celltype,
                                           clipping = ref_BAF_clip)

            used_specific_states = xclone.model.gene_specific_BAF(merge_Xdata_copy, 
                               theo_states= guide_theo_states, specific_BAF = "ref_BAF_phased")
        
        merge_Xdata_copy = xclone.model.calculate_Xemm_prob_bb(merge_Xdata_copy, 
                                                      AD_key = "ad_bin_phased", DP_key = "dp_bin",
                                                      outlayer = "correct_emm_prob_log", 
                                                      states = used_specific_states,
                                                      gene_specific = gene_specific_concentration, 
                                                      concentration = concentration)
        merge_Xdata_copy = xclone.model.BAF_smoothing(merge_Xdata_copy,
                                             inlayer = "correct_emm_prob_log",
                                             outlayer = "correct_emm_prob_log_KNN",
                                             KNN_connectivities_key = "connectivities_expr",
                                             KNN_smooth = True)
        start_prob = np.array([0.3, 0.4, 0.3])                                 
        trans_prob = np.array([[1-2*t, t, t],[t, 1-2*t, t],[t, t, 1-2*t]])
        merge_Xdata_copy = xclone.model.XHMM_smoothing(merge_Xdata_copy,
                                              brk = HMM_brk, 
                                              start_prob = start_prob,  
                                              trans_prob = trans_prob, 
                                              emm_inlayer = "correct_emm_prob_log_KNN", 
                                              nproc = 80, 
                                              verbose = False)
        neutral_index = 1
        cnv_index = [0,2]
        merge_Xdata_copy = xclone.model.denoise_gene_scale(merge_Xdata_copy, Xlayer = "posterior_mtx",
                       neutral_index = neutral_index, 
                       cnv_index = cnv_index, 
                       GMM_detection = True,
                       gmm_comp = 2,
                       cell_prop_cutoff = 0.05,
                       out_layer = "denoised_posterior_mtx")

        merge_Xdata.layers["add_posterior_mtx"] = merge_Xdata_copy.layers["posterior_mtx"]
        merge_Xdata.layers["denoised_add_posterior_mtx"] = merge_Xdata_copy.layers["denoised_posterior_mtx"]
    
    if BAF_denoise:
        if CNV_N_components == 3:
            neutral_index = 1
            cnv_index = [0,2]
        elif CNV_N_components == 5:
            neutral_index = 2
            cnv_index = [0,1,3,4]
        merge_Xdata = xclone.model.denoise_gene_scale(merge_Xdata, Xlayer = "posterior_mtx",
                       neutral_index = neutral_index, 
                       cnv_index = cnv_index, 
                       GMM_detection = GMM_detect,
                       gmm_comp = BAF_denoise_GMM_comp,
                       cell_prop_cutoff = BAF_denoise_cellprop_cutoff,
                       out_layer = "denoised_posterior_mtx")
    try:
        merge_Xdata.write(BAF_final_file)
    except Exception as e:
        print("[XClone Warning]", e)
    else:
        print("[XClone hint] BAF_final_file saved in %s." %(out_data_dir))

    end_time = datetime.now(timezone.utc)
    time_passed = end_time - start_time
    main_logger.info("XClone BAF module finished (%d seconds)" % (time_passed.total_seconds()))

    if xclone_plot:
        set_figtitle = config.set_figtitle
        if plot_cell_anno_key is None:
            plot_cell_anno_key = cell_anno_key
        run_BAF_plot(merge_Xdata, dataset_name, 
                     plot_cell_anno_key, 
                     set_figtitle,
                     out_dir)

    return merge_Xdata

def run_BAF_plot(merge_Xdata,
                 dataset_name,
                 plot_cell_anno_key,
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
    baf_smooth_fig = out_plot_dir + dataset_name + "_BAF_smooth.png"
    baf_final_fig1 = out_plot_dir + dataset_name + "_BAF_CNV.png"
    baf_final_fig2 = out_plot_dir + dataset_name + "_BAF_CNV_denoise.png"
    baf_final_fig3 = out_plot_dir + dataset_name + "_BAF_CNV_3states.png"
    baf_final_fig4 = out_plot_dir + dataset_name + "_BAF_CNV_3states_denoise.png"

    sub_logger = get_logger("BAF plot module")
    sub_logger.info("BAF plot module started")
    start_time = datetime.now(timezone.utc)

    if set_figtitle:
        fig_title = dataset_name + " BAF smooth"
    xclone.pl.BAF_smooth_visualization(merge_Xdata, Xlayer = "BAF_phased_KNN_WMA", 
                                       cell_anno_key = plot_cell_anno_key, 
                                       change_colorbar = False, 
                                       colorbar_name = "BAF",
                                       title = fig_title,
                                       save_file = True, 
                                       out_file = baf_smooth_fig,
                                       **kwargs)
    if set_figtitle:
        fig_title = dataset_name + " BAF module"
    xclone.pl.BAF_CNV_visualization(merge_Xdata, weights = False, 
                                    cell_anno_key = plot_cell_anno_key, 
                                    title = fig_title,
                                    save_file = True, 
                                    out_file = baf_final_fig1,
                                    **kwargs)
                                    
    
    if "denoised_posterior_mtx" in merge_Xdata.layers:
        xclone.pl.BAF_CNV_visualization(merge_Xdata, Xlayer = "denoised_posterior_mtx",
                                        weights = False, 
                                        cell_anno_key = plot_cell_anno_key, 
                                        title = fig_title + " (denoise)",
                                        save_file = True, 
                                        out_file = baf_final_fig2,
                                        **kwargs)
    
    if "add_posterior_mtx" in merge_Xdata.layers:
        xclone.pl.BAF_CNV_visualization(merge_Xdata, Xlayer = "add_posterior_mtx",
                                        weights = False, 
                                        cell_anno_key = plot_cell_anno_key, 
                                        title = fig_title,
                                        save_file = True, 
                                        out_file = baf_final_fig3,
                                        **kwargs)
    if "denoised_add_posterior_mtx" in merge_Xdata.layers:
        xclone.pl.BAF_CNV_visualization(merge_Xdata, Xlayer = "denoised_add_posterior_mtx",
                                        weights = False, 
                                        cell_anno_key = plot_cell_anno_key, 
                                        title = fig_title + " (denoise)",
                                        save_file = True, 
                                        out_file = baf_final_fig4,
                                        **kwargs)

    
    end_time = datetime.now(timezone.utc)
    time_passed = end_time - start_time
    sub_logger.info("BAF plot module finished (%d seconds)" % (time_passed.total_seconds()))
    return None

def plot_BAF_module():
    """
    all plots for BAF module.
    """
    pass

def plot_processed_BAF(merge_Xdata, 
                       dataset_name,
                       plot_cell_anno_key,
                       set_figtitle = True,
                       set_dpi = 300,
                       out_dir = None,
                       **kwargs):
    """
    BAF phasing and smoothing plots.
    """
    if out_dir is None:
        cwd = os.getcwd()
        out_dir = cwd + "/XCLONE_OUT/"
    
    out_plot_dir = str(out_dir) + "/plot/"
    xclone.al.dir_make(out_plot_dir)

    fig_title = ""
    
    baf_fig1 = out_plot_dir + dataset_name + "_BAF.png"
    baf_fig2 = out_plot_dir + dataset_name + "_BAF_phased.png"
    baf_fig3 = out_plot_dir + dataset_name + "_BAF_phased_KNN.png"
    baf_fig4 = out_plot_dir + dataset_name + "_BAF_phased_KNN_WMA.png"
    
    if set_figtitle:
        fig_title = dataset_name + " BAF"
    xclone.pl.visualize_cell_BAF(merge_Xdata, Xlayer = "BAF", cell_anno_key = plot_cell_anno_key, set_dpi =set_dpi,
                                title = fig_title, save_file = True, out_file = baf_fig1, **kwargs)
    if set_figtitle:
        fig_title = dataset_name + " Phased BAF"
    xclone.pl.visualize_cell_BAF(merge_Xdata, Xlayer = "BAF_phased", cell_anno_key = plot_cell_anno_key, set_dpi =set_dpi,
                                title = fig_title, save_file = True, out_file = baf_fig2, **kwargs)
    if set_figtitle:
        fig_title = dataset_name + " Phased BAF after KNN smoothing"
    xclone.pl.visualize_cell_BAF(merge_Xdata, Xlayer = "BAF_phased_KNN", cell_anno_key = plot_cell_anno_key, set_dpi =set_dpi,
                                title = fig_title, save_file = True, out_file = baf_fig3, **kwargs)
    if set_figtitle:
        fig_title = dataset_name + " Phased BAF after KNN smoothing and WMA smoothing"
#     xclone.pl.visualize_cell_BAF(merge_Xdata, Xlayer = "BAF_phased_KNN_WMA", cell_anno_key = plot_cell_anno_key, set_dpi =set_dpi,
#                                 title = fig_title, save_file = True, out_file = baf_fig4)
    xclone.pl.BAF_smooth_visualization(merge_Xdata, Xlayer = "BAF_phased_KNN_WMA", cell_anno_key=plot_cell_anno_key, change_colorbar = True,  
                                       colorbar_ticks = [0.1, 0.5, 0.9], colorbar_label = ["0",  "0.5",  "1"], colorbar_name = "BAF values",set_dpi =set_dpi,
                                       title = fig_title, save_file = True, out_file = baf_fig4, **kwargs)