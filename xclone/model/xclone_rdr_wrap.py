"""Base pipeline for XClone RDR module"""

# Author: Rongting Huang
# Date: 2022-03-19
# update: 2022-11-17


import os
import xclone
import numpy as np
from .._logging import get_logger
from datetime import datetime, timezone

def preview_RDR(RDR_adata,
                cell_anno_key,):
    """
    Function:
    some precheck for params setting.
    """
    xclone.model.view_celltype(RDR_adata, cell_anno_key = cell_anno_key)

    ## smooth plot?
    return None


def run_RDR(RDR_adata,
            verbose = True,
            run_verbose = True,
            config_file = None
            ):
    """
    """
    ## settings
    from .._config import XCloneConfig
    
    if config_file == None:
        print (
            f'Model configuration file not specified.\n'
            f'Default settings in XClone-RDR will be used.'
        )
        config = XCloneConfig(module = "RDR")

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

    # HMM settings
    start_prob = config.start_prob
    trans_t = config.trans_t
    trans_prob = config.trans_prob
    HMM_brk = config.HMM_brk
    ## optimize 
    max_iter = config.max_iter
    min_iter = config.min_iter

    # RDR settings
    filter_ref_ave = config.filter_ref_ave
    smart_transform = config.smart_transform
    marker_group_anno_key = config.marker_group_anno_key
    top_n_marker = config.top_n_marker
    remove_marker = config.remove_marker
    fit_GLM_libratio = config.fit_GLM_libratio
    dispersion_celltype = config.dispersion_celltype
    gene_exp_group = config.gene_exp_group
    gene_exp_ref_log = config.gene_exp_ref_log
    remove_guide_XY = config.remove_guide_XY
    guide_cnv_ratio = config.guide_cnv_ratio
    ## smoothing settings
    KNN_neighbors = config.KNN_neighbors
    WMA_window_size = config.WMA_window_size
    WMA_smooth_key = config.WMA_smooth_key
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
    RDR_base_file = out_data_dir + "RDR_base_adata.h5ad"
    RDR_bulk_file = out_data_dir + "RDR_bulk_adata.h5ad"
    ### output after CNV calling
    RDR_final_file = out_data_dir + "RDR_adata_KNN_HMM_post.h5ad"
    
    ##------------------------
    main_logger = get_logger("Main RDR module")
    main_logger.info(dataset_name)
    main_logger.info("XClone RDR module Started")
    start_time = datetime.now(timezone.utc)

    if run_verbose:
        print("[XClone RDR module running]************************")
    
    if exclude_XY:
        RDR_adata = xclone.pp.exclude_XY_adata(RDR_adata)
        print("[XClone warning] RDR module excelude chr XY analysis.")
    ## RDR data preprocessing
    RDR_adata = xclone.pp.check_RDR(RDR_adata, 
    cell_anno_key = cell_anno_key, 
    cell_detection_rate = 0.05, verbose = verbose)
    
    ## Transformation for smart-seq data
    RDR_adata = xclone.pp.Xtransformation(RDR_adata, transform = smart_transform, Xlayers = ["raw_expr"])
    RDR_adata = xclone.pp.Xdata_RDR_preprocess(RDR_adata, filter_ref_ave = filter_ref_ave, 
    cell_anno_key = cell_anno_key,  ref_celltype = ref_celltype, mode = "FILTER")
    
    ## marker genes
    if marker_group_anno_key is None:
        marker_group_anno_key = cell_anno_key
    marker_genes = xclone.pp.get_markers(RDR_adata,
                                         top_n = top_n_marker,
                                         target_sum=1e4,
                                         rank_group = marker_group_anno_key,
                                         marker_genes_key = "rank_marker_genes",
                                         test_method='wilcoxon')
    ## remove marker genes
    if remove_marker:
        RDR_adata = xclone.pp.filter_markers(RDR_adata, marker_lst = marker_genes)

    ## get base bulk and single cell adata
    RDR_adata, RDR_adata_bulk = xclone.pp.Xdata_RDR_preprocess(RDR_adata, filter_ref_ave = None, 
    cell_anno_key = cell_anno_key, ref_celltype = ref_celltype)

    if fit_GLM_libratio:
        ## fit libratio for each cell based on select normal chrs*
        select_normal_chr_num = config.select_normal_chr_num
        RDR_adata_bulk = xclone.model.select_normal_CHR(RDR_adata_bulk, select_chr_num = select_normal_chr_num)
        ### extend the chr index to whole dataset for libratio fitting
        RDR_adata_bulk= xclone.model.extend_chr_index_to_singlecell(RDR_adata_bulk, 
                                                                    RDR_adata_bulk, 
                                                                    cell_anno_key = cell_anno_key)
        RDR_adata = xclone.model.extend_chr_index_to_singlecell(RDR_adata, 
                                                                RDR_adata_bulk, 
                                                                cell_anno_key = cell_anno_key)

        RDR_adata_bulk = xclone.model.fit_lib_ratio(RDR_adata_bulk, verbose = True, 
                                                    NB_kwargs={'disp': False, 'skip_hessian': True})
        RDR_adata = xclone.model.fit_lib_ratio(RDR_adata, verbose = False, 
                                                NB_kwargs={'disp': False, 'skip_hessian': True})
    
        xclone.model.check_libratio(RDR_adata, anno_key = "library_ratio")
        RDR_adata = xclone.model.libsize_clip(RDR_adata, max_threshold = 5)
        xclone.model.check_libratio(RDR_adata, anno_key = "library_ratio_capped")

        libsize_key = "library_ratio_capped"
        depth_key = "library_ratio_capped"
    
    else:
        libsize_key = "counts_ratio"
        depth_key = "counts_ratio"
    ## fit gene-specific dispersion based on reference celltype
    if dispersion_celltype is None:
        dispersion_celltype = ref_celltype
    RDR_adata_REF = xclone.model.select_celltype(RDR_adata, 
                                                 cell_anno_key = cell_anno_key, 
                                                 select_celltype = dispersion_celltype)
    RDR_adata_REF = xclone.model.fit_Dispersions(RDR_adata_REF, 
                                                 libsize_key = libsize_key, 
                                                 NB_kwargs={'disp': False, 'skip_hessian': True},
                                                 verbose = False, model_results_path = None)
    xclone.model.check_dispersion(RDR_adata_REF, anno_key = "dispersion")
    RDR_adata = xclone.model.map_var_info(RDR_adata, RDR_adata_REF, specific_celltype = dispersion_celltype)
    RDR_adata = xclone.model.remove_genes(RDR_adata)
    RDR_adata = xclone.model.remove_genes(RDR_adata, mode="INF")
    RDR_adata = xclone.model.dispersion_clip(RDR_adata, qt_low = 0.03, qt_up =0.93,  
                                             min_threshold = None, max_threshold = None)
    xclone.model.check_dispersion(RDR_adata, anno_key = "dispersion_capped")
    
    ## output before CNV calling, save data with fitted libratio and dispersion.
    try:
        RDR_adata.write(RDR_base_file)
        RDR_adata_bulk.write(RDR_bulk_file)
    except Exception as e:
        print("[XClone Warning]", e)
    else:
        print("[XClone hint] RDR_base_file and bulk_file saved in %s." %(out_data_dir))

    RDR_adata = xclone.model.extra_preprocess(RDR_adata, cluster_key = cell_anno_key,
                                              ref_celltype = ref_celltype,
                                              depth_key=depth_key,
                                              low_dim=False, run_KNN=True, 
                                              KNN_neighbors = KNN_neighbors,
                                              copy=True)

    RDR_adata = xclone.model.RDR_smoothing_base(RDR_adata,
                                                clip = True,
                                                outlayer = "RDR_smooth",
                                                cell_anno_key = cell_anno_key,
                                                ref_celltype = ref_celltype,
                                                WMA_window_size = WMA_window_size,
                                                chrom_key = WMA_smooth_key,
                                                KNN_sm = True,
                                                KNN_connect_use = "connectivities")
    
    if guide_cnv_ratio is None:
        chr_anno_key = config.guide_chr_anno_key
        guide_chr_lst = xclone.model.guide_CNV_chrs(RDR_adata, 
                                                    Xlayer = "RDR_smooth", 
                                                    anno_key = chr_anno_key,
                                                    remove_XY = remove_guide_XY)
        guide_qt_lst = config.guide_qt_lst
        guide_cnv_ratio = xclone.model.guide_CNV_states(RDR_adata, 
                                                        Xlayer = "RDR_smooth", 
                                                        chr_lst = guide_chr_lst, 
                                                        anno_key = chr_anno_key, 
                                                        qt_lst = guide_qt_lst, 
                                                        show_boxplot = False)
        
    RDR_adata = xclone.model.gene_exp_group(RDR_adata, 
                                            n_group = gene_exp_group,
                                            ref_log =  gene_exp_ref_log,
                                            verbose = verbose)
    
    t = trans_t
    if trans_prob is None:
        trans_prob = np.array([[1-2*t, t, t],[t, 1-2*t, t],[t, t, 1-2*t]])
    RDR_adata = xclone.model.CNV_optimazation(RDR_adata, depth_key = depth_key, init_state_ratio = guide_cnv_ratio,
                    max_iter = max_iter,
                    min_iter = min_iter,
                    HMM_brk = HMM_brk,
                    start_prob = start_prob,
                    trans_prob = trans_prob,
                    verbose = True)
    
    ## output after CNV calling, save data with CNV posterior.
    try:
        RDR_adata.write(RDR_final_file)
    except Exception as e:
        print("[XClone Warning]", e)
    else:
        print("[XClone hint] RDR_final_file saved in %s." %(out_data_dir))

    end_time = datetime.now(timezone.utc)
    time_passed = end_time - start_time
    main_logger.info("XClone RDR module finished (%d seconds)" % (time_passed.total_seconds()))
    
    if xclone_plot:
        rdr_plot_vmin = config.rdr_plot_vmin
        rdr_plot_vmax = config.rdr_plot_vmax
        set_figtitle = config.set_figtitle
        if run_verbose:
            print("[XClone plotting]")
        if plot_cell_anno_key is None:
            plot_cell_anno_key = cell_anno_key
        
        run_RDR_plot(RDR_adata, dataset_name, 
                     plot_cell_anno_key, 
                     set_figtitle,
                     rdr_plot_vmin, rdr_plot_vmax, 
                     out_dir)
    return RDR_adata                                              
 
def run_RDR_plot(RDR_adata,
            dataset_name,
            plot_cell_anno_key,
            set_figtitle = True,
            rdr_plot_vmin = -0.7,
            rdr_plot_vmax = 0.7,
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

    # default:XClone
    fig_title = ""
    rdr_smooth_fig = out_plot_dir + dataset_name + "_RDR_smooth.png"
    rdr_final_fig = out_plot_dir + dataset_name + "_RDR_CNV.png"

    sub_logger = get_logger("RDR plot module")
    sub_logger.info("RDR plot module started")
    start_time = datetime.now(timezone.utc)
    
    if set_figtitle:
        fig_title = dataset_name + " RDR_smooth (log scale ratio)"

    xclone.pl.RDR_smooth_visualization(RDR_adata, 
                                       Xlayer = "RDR_smooth", 
                                       cell_anno_key = plot_cell_anno_key,
                                       change_colorbar = False,
                                       vmin = rdr_plot_vmin, vmax = rdr_plot_vmax,
                                       title = fig_title,
                                       save_file = True, 
                                       out_file = rdr_smooth_fig,
                                       **kwargs)
    if set_figtitle:
        fig_title = dataset_name + " RDR module"
    
    xclone.pl.CNV_visualization(RDR_adata, 
                                states_weight = np.array([1,2,3]), 
                                weights = True, 
                                cell_anno_key = plot_cell_anno_key, 
                                title = fig_title,
                                save_file = True, 
                                out_file = rdr_final_fig,
                                **kwargs)
    
    end_time = datetime.now(timezone.utc)
    time_passed = end_time - start_time
    sub_logger.info("RDR plot module finished (%d seconds)" % (time_passed.total_seconds()))
    return None

def plot_RDR_module():
    """
    all plots for RDR module.
    """
    pass

def plot_processed_RDR(RDR_Xdata, 
                       dataset_name,
                       plot_cell_anno_key,
                       set_figtitle = True,
                       set_dpi = 300,
                       out_dir = None,
                       **kwargs):
    """
    RDR raw ratio and smoothing plots.
    """
    if out_dir is None:
        cwd = os.getcwd()
        out_dir = cwd + "/XCLONE_OUT/"
    
    out_plot_dir = str(out_dir) + "/plot/"
    xclone.al.dir_make(out_plot_dir)

    fig_title = ""
    
    rdr_fig1 = out_plot_dir + dataset_name + "_RDR_raw_ratio.png"
    rdr_fig2 = out_plot_dir + dataset_name + "_RDR_WMA.png"
    rdr_fig3 = out_plot_dir + dataset_name + "_RDR_WMA_KNN.png"
    
    if set_figtitle:
        fig_title = dataset_name + " RDR raw ratio"
    xclone.pl.raw_ratio_visualization(RDR_Xdata, Xlayer = "raw_ratio", cell_anno_key = plot_cell_anno_key,
                                      set_dpi =set_dpi, change_colorbar = False,vmin = -0.7, vmax = 0.7,
                                      title = fig_title, save_file = True, out_file = rdr_fig1, **kwargs)
    if set_figtitle:
        fig_title = dataset_name + " RDR raw ratio after WMA smoothing"
    xclone.pl.RDR_smooth_visualization(RDR_Xdata, Xlayer = 'WMA_smoothed', cell_anno_key = plot_cell_anno_key, set_dpi =set_dpi,
                                       change_colorbar = False, vmin = -0.7, vmax = 0.7,
                                       title = fig_title, save_file = True, out_file = rdr_fig2, **kwargs)
    if set_figtitle:
        fig_title = dataset_name + " RDR raw ratio after WMA smoothing and KNN smoothing"
    xclone.pl.RDR_smooth_visualization(RDR_Xdata, Xlayer = 'RDR_smooth', cell_anno_key = plot_cell_anno_key, set_dpi =set_dpi,
                                       change_colorbar = False, vmin = -0.7, vmax = 0.7,
                                       title = fig_title, save_file = True, out_file = rdr_fig3, **kwargs)