"""Base pipeline for XClone RDR module"""

# Author: Jiamu James Qiao, adapted from xclone_rdr_wrap.py
# Date: 2025-06-10
# update: 


import os
import xclone
import numpy as np
from .._logging import get_logger
from datetime import datetime, timezone

from .smoothing import make_WMA_connectivity

import gc

def run_RDR(RDR_adata, verbose = True, config_file = None):
    """
    Run the RDR (Read Depth Ratio) analysis on the provided annotated data.

    This function is an updated version of the function `run_RDR` in the XClone package.
    It uses a Gaussian mixture model to process log-transformed read depth ratios
    and performs CNV (Copy Number Variation) analysis based on the provided configuration.

    This function performs the RDR analysis using the provided annotated data (`RDR_adata`).
    It allows for verbose output and the use of a custom configuration object. If no configuration
    object is specified, default settings from XClone's RDR module will be used.

    Parameters
    ----------

        RDR_adata : anndata.AnnData
            The annotated data matrix on which the RDR analysis will be performed.
        verbose : bool, optional
            If True, prints detailed information about the process. Default is True.
        config_file : xclone.XCloneConfig or None, optional
            The XClone configuration object. If None, the default settings in XClone-RDR will be used.
            Default is None.
            The configuration can be created as follows:
            
            .. code-block:: python

                config_file = xclone.XCloneConfig(dataset_name="your_dataset_name", module="RDR_gaussian")

    Returns
    -------

        RDR_adata : anndata.AnnData
            The finalized output with multiple layers of information in the `anndata.AnnData` from RDR module.


    Example
    -------

        .. code-block:: python

            import xclone
            
            # Run RDR analysis with default settings
            RDR_Xdata = xclone.model.run_RDR_gaussian(RDR_adata, verbose=True)
            
            # Run RDR analysis with a custom configuration object
            xconfig = xclone.XCloneConfig(dataset_name = dataset_name, module = "RDR_gaussian")
            xconfig.outdir = out_dir
            #... other specified parameters in `xconfig`
            xconfig.display()
            RDR_Xdata = xclone.model.run_RDR_gaussian(RDR_adata, verbose=True, config_file=xconfig)
    
    """
    ## settings
    from .._config import XCloneConfig
    
    if config_file == None:
        print (
            f'Model configuration file not specified.\n'
            f'Default settings in XClone-RDR will be used.'
        )
        config = XCloneConfig(module = "RDR_gaussian")

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
    min_gene_keep_num = config.min_gene_keep_num
    multi_refcelltype = config.multi_refcelltype
    smart_transform = config.smart_transform
    marker_group_anno_key = config.marker_group_anno_key
    get_marker_genes = config.get_marker_genes
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
    KNN_npcs = config.KNN_npcs
    WMA_window_size = config.WMA_window_size
    WMA_smooth_key = config.WMA_smooth_key
    # plot settings
    xclone_plot = config.xclone_plot
    plot_cell_anno_key = config.plot_cell_anno_key
    # develop mode settings
    develop_mode = config.develop_mode

    ## adaptive baseline settings
    ab_k_neighbors = config.ab_k_neighbors
    ab_pseudo_count = config.ab_pseudo_count

    ## denoise settings
    denoise_sd_amplifier = config.denoise_sd_amplifier

    ## GMM settings
    c_k = config.c_k

    ## low rank settings
    low_rank = config.low_rank
    low_rank_n_components = config.low_rank_n_components


    ## Result output prepare
    if out_dir is None:
        cwd = os.getcwd()
        out_dir = cwd + "/XCLONE_OUT/"
    out_data_dir = str(out_dir) + "/data/"
    xclone.al.dir_make(out_data_dir)
    
    ### output after CNV calling
    RDR_final_file = out_data_dir + "RDR_adata_KNN_HMM_post.h5ad"
    
    ##------------------------
    main_logger = get_logger("Main RDR module")
    main_logger.info(dataset_name)
    main_logger.info("XClone RDR module Started")
    start_time = datetime.now(timezone.utc)

    if verbose:
        print("[XClone RDR module running]************************")
    
    if exclude_XY:
        RDR_adata = xclone.pp.exclude_XY_adata(RDR_adata)
        print("[XClone warning] RDR module excelude chr XY analysis.")
    ## RDR data preprocessing
    RDR_adata = xclone.pp.check_RDR(RDR_adata, 
    cell_anno_key = cell_anno_key, 
    cell_detection_rate = 0.05, verbose = verbose)

    ## check ref_celltype
    # if ref_celltype in RDR_adata.obs[cell_anno_key].values:
        # pass
    # else:
        # raise ValueError(f"[XClone error] Item '{ref_celltype}' not found in the RDR_adata's annotation.")
    # modified for multiple ref_celltype
    if isinstance(ref_celltype, list):
        missing_items = [item for item in ref_celltype if item not in RDR_adata.obs[cell_anno_key].values]
        if missing_items:
            raise ValueError(f"[XClone error] Items {missing_items} not found in the RDR_adata's annotation.")
    else:
        if ref_celltype not in RDR_adata.obs[cell_anno_key].values:
            raise ValueError(f"[XClone error] Item '{ref_celltype}' not found in the RDR_adata's annotation.")
    
    ## Transformation for smart-seq data
    RDR_adata = xclone.pp.Xtransformation(RDR_adata, transform = smart_transform, Xlayers = ["raw_expr"])
    
    if config.set_spatial:
        spatial_transform = config.spatial_transform
        RDR_adata = xclone.pp.Spatial_Xtransformation(RDR_adata, transform = spatial_transform, Xlayers = ["raw_expr"])
    
    ## Consider data coverage before genes filtering (update)
    RDR_adata = xclone.pp.Xdata_RDR_preprocess(RDR_adata, 
                                               filter_ref_ave = filter_ref_ave, 
                                               min_gene_keep_num = min_gene_keep_num,
                                               cell_anno_key = cell_anno_key,  
                                               ref_celltype = ref_celltype, 
                                               mode = "FILTER")
    
    ## marker genes
    if marker_group_anno_key is None:
        marker_group_anno_key = cell_anno_key
    if get_marker_genes:
        marker_genes = xclone.pp.get_markers(RDR_adata,
                                         top_n = top_n_marker,
                                         target_sum=1e4,
                                         rank_group = marker_group_anno_key,
                                         marker_genes_key = "rank_marker_genes",
                                         test_method='wilcoxon')
    ## remove marker genes
    if remove_marker:
        RDR_adata = xclone.pp.filter_markers(RDR_adata, marker_lst = marker_genes)



    ## generate adaptive baseline
    ## out layer: 'log_ratio_ab'
    RDR_adata = xclone.model.preprocess_adaptive_baseline(RDR_adata, 
                                                          cell_anno_key=cell_anno_key,
                                                          ref_celltype=ref_celltype,
                                                          k=ab_k_neighbors, 
                                                          pseudo_count=ab_pseudo_count, 
                                                          verbose=False, 
                                                          plot=False)

    ## WMA smoothing
    RDR_adata = xclone.model.WMA_smooth(RDR_adata, 
                                          layer='log_ratio_ab', 
                                          out_layer='WMA_smoothed_log_ratio_ab', 
                                          window_size = WMA_window_size,
                                          chrom_key = WMA_smooth_key, 
                                          verbose=False)
    # default start layer is 'log_ratio_ab'
    # default out layer is 'WMA_smoothed_log_ratio_ab'
    print("[XClone RDR module] WMA smoothing applied, smoothed layer saved to 'WMA_smoothed_log_ratio_ab'.")
    # WMA_smooth in smoothing.py

    ## denoise
    # default start layer is 'WMA_smoothed_log_ratio_ab'
    RDR_adata, layer_name = xclone.model.denoise(RDR_adata, 
                                                 method='dynamic',
                                                 cell_anno_key=cell_anno_key,
                                                 ref_celltype=ref_celltype,
                                                 sd_amplifier=denoise_sd_amplifier)    
    denoised_layer = layer_name

    ## GMM
    RDR_adata, layer_name = xclone.model.compute_gaussian_probabilities(RDR_adata, 
                                                                        layer=layer_name,
                                                                        cell_anno_key=cell_anno_key,
                                                                        ref_celltype=ref_celltype,
                                                                        c_k=c_k)
    if verbose:
        print("[XClone RDR module] GMM probabilities calculated.")
        print("[XClone RDR module] Layer '%s' contains the GMM probabilities." % layer_name)

    # HMM
    z_post_all = xclone.model.calculate_z_post_parallel(RDR_adata, start_prob, trans_prob)
    RDR_adata.layers['post_HMM_prob'] = z_post_all

    # 2nd iteration
    z_post_all = xclone.model.calculate_z_post_parallel(RDR_adata, start_prob, trans_prob, layer='post_HMM_prob')
    RDR_adata.layers['post_HMM_prob_2'] = z_post_all

    if verbose:
        print("[XClone RDR module] HMM posterior probabilities calculated.")
        print("[XClone RDR module] Layer 'post_HMM_prob' contains the HMM posterior probabilities.")
        print("[XClone RDR module] Layer 'post_HMM_prob_2' contains the 2nd iteration of HMM posterior probabilities.")

    ## denoise
    RDR_adata, layer_name = xclone.model.smooth_anndata_layer(RDR_adata, 
                                                  layer='post_HMM_prob_2',
                                                  method='gauss')
    if verbose:
        print("[XClone RDR module] HMM posterior probabilities smoothed.")
        print("[XClone RDR module] Layer '%s' contains the smoothed HMM posterior probabilities." % layer_name)

    ## recommended low rank
    if low_rank:
        RDR_adata, layer_name = xclone.model.low_rank_approximation(RDR_adata, 
                                                        layer=layer_name,
                                                        n_components=low_rank_n_components)
        ## HMM
        z_post_all = xclone.model.calculate_z_post_parallel(RDR_adata, start_prob, trans_prob, layer=layer_name)
        layer_name = 'post_HMM_prob_2_gauss_lowrank'
        RDR_adata.layers[layer_name] = z_post_all
        print("[XClone hint] Low rank approximation applied to the HMM posterior probabilities.")
        print("[XClone hint] Layer '%s' contains the low rank approximation of HMM posterior probabilities." % layer_name)

    RDR_adata.layers['posterior_mtx'] = RDR_adata.layers[layer_name]
    
    ## output after CNV calling, save data with CNV posterior.
    try:
        if develop_mode:
            pass
        else:
            # layers_to_keep = ['raw_expr', 'log_ratio_ab', 'post_HMM_prob_2_gauss', 'post_HMM_prob_2_gauss_lowrank']
            # RDR_adata = xclone.pp.keep_layers(RDR_adata, layers_to_keep)
            pass
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

        if verbose:
            print("[XClone plotting]")
        if plot_cell_anno_key is None:
            plot_cell_anno_key = cell_anno_key
        
        plot_remove_immune = config.plot_remove_immune
        plot_immune_celltype = config.plot_immune_celltype
        if plot_remove_immune:
            pass

        xclone.model.run_RDR_gaussian_plot(RDR_adata, dataset_name, 
                     plot_cell_anno_key, 
                     log_ratio_layer = denoised_layer,
                     prob_layer = layer_name,
                     set_figtitle= set_figtitle,
                     rdr_plot_vmin=rdr_plot_vmin, 
                     rdr_plot_vmax= rdr_plot_vmax, 
                     out_dir= out_dir)
    return RDR_adata   