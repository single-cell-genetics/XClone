"""Base pipeline for XClone BAF module"""

# Author: Rongting Huang
# Date: 2022-10-21
# update: 2022-11-17


import os
import xclone
import anndata as an
import numpy as np
from .._logging import get_logger
from datetime import datetime, timezone

def load_BAF(BAF_file, anno_file):
    """
    """
    pass
    # return BAF_adata


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
    
    # BAF settings
    RDR_file = config.RDR_file
    theo_neutral_BAF = config.theo_neutral_BAF
    WMA_window = config.WMA_window

    # HMM settings
    start_prob = config.start_prob
    trans_t = config.trans_t

    # plot settings
    xclone_plot = config.xclone_plot
    plot_cell_anno_key = config.plot_cell_anno_key
    

    ## Result output prepare
    if out_dir is None:
        cwd = os.getcwd()
        out_dir = cwd + "/XCLONE_OUT/"
    out_data_dir = str(out_dir) + "data/"
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
    BAF_adata = xclone.pp.check_BAF(BAF_adata, cell_anno_key = cell_anno_key, verbose = verbose)

    RDR_adata = an.read_h5ad(RDR_final_file)
    ## check validated RDR and BAF
    xclone.pp.check_RDR_BAF_cellorder(RDR_adata, BAF_adata)
    ## Remove marker genes
    marker_genes = RDR_adata.uns["rank_marker_genes"]
    BAF_adata = xclone.pp.Xdata_region_selection(BAF_adata,
                                                 select = False,
                                                 chr_anno_key = "GeneName",
                                                 chr_lst = marker_genes,
                                                 update_uns = False,
                                                 uns_anno_key = None)
    ## BAF Phasing
    BAF_adata, merge_Xdata =  xclone.model.BAF_Local_phasing(BAF_adata, 
                                                             region_key = "chr", 
                                                             phasing_len = 100, 
                                                             bin_nproc=20)
    BAF_adata, merge_Xdata = xclone.model.BAF_Global_phasing(BAF_adata, merge_Xdata)

    BAF_adata.write(BAF_base_file)
    merge_Xdata.write(BAF_merge_base_file)
    ### check coverage for bins
    merge_Xdata.var[(merge_Xdata.layers["dp_bin"].A.sum(axis=0) == 0)]

    merge_Xdata = xclone.pl.calculate_cell_BAF(merge_Xdata, 
                                               AD_key = "ad_bin1", DP_key = "dp_bin", BAF_key = "BAF")
    merge_Xdata = xclone.pl.calculate_cell_BAF(merge_Xdata, 
                                               AD_key = "ad_bin1_phased", DP_key = "dp_bin", BAF_key = "BAF_phased")
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
                                          window_size = WMA_window, 
                                          verbose=False)
    
    merge_Xdata = xclone.model.WMA_smooth(merge_Xdata, 
                                          layer="fill_BAF_phased", 
                                          out_layer='BAF_phased_WMA', 
                                          window_size = WMA_window, 
                                          verbose=False)
    # merge_Xdata = xclone.model.KNN_smooth(merge_Xdata, 
    #                                       KNN_connect_use = "connectivities_expr", 
    #                                       layer="BAF_phased_WMA", 
    #                                       out_layer='BAF_phased_WMA_KNN')
    
    BAF_adata.write(BAF_base_file)
    merge_Xdata.write(BAF_merge_base_file)
    ## HMM smoothing for CNV states calling
    CNV_states = xclone.model.get_CNV_states(merge_Xdata, Xlayer = "BAF_phased_WMA",
                   n_components = 3,
                   means_init = None,
                   max_iter = 200)
    guide_theo_states = xclone.model.guide_BAF_theo_states(CNV_states)

    merge_Xdata = xclone.model.get_BAF_ref(merge_Xdata, 
                                           Xlayer = "fill_BAF_phased", 
                                           out_anno = "ref_BAF_phased",
                                           anno_key = cell_anno_key, 
                                           ref_cell = ref_celltype)
    merge_Xdata = xclone.model.get_BAF_ref(merge_Xdata, 
                                           Xlayer = "fill_BAF_phased", 
                                           out_anno = "ref_BAF_phased_clipped",
                                           anno_key = cell_anno_key, 
                                           ref_cell = ref_celltype, 
                                           clipping = True)

    ## if you have limited ref cells, you can set theo_baf as 0.5   *important
    if theo_neutral_BAF is not None:
        merge_Xdata.var["theo_neutral_BAF"] = theo_neutral_BAF
        used_specific_states = xclone.model.gene_specific_BAF(merge_Xdata, 
                            theo_states= guide_theo_states, specific_BAF = "theo_neutral_BAF")
    else:
        used_specific_states = xclone.model.gene_specific_BAF(merge_Xdata, 
                            theo_states= guide_theo_states, specific_BAF = "ref_BAF_phased")

    merge_Xdata = xclone.model.calculate_Xemm_prob_bb(merge_Xdata, 
                                                      AD_key = "ad_bin1_phased", DP_key = "dp_bin", 
                                                      concentration = 100,
                                                      outlayer = "bin_phased_BAF_specific_center_emm_prob_log", 
                                                      states = used_specific_states)

    merge_Xdata = xclone.model.BAF_smoothing(merge_Xdata,
                                             inlayer = "bin_phased_BAF_specific_center_emm_prob_log",
                                             outlayer = "bin_phased_BAF_specific_center_emm_prob_log_KNN",
                                             KNN_connectivities_key = "connectivities_expr",
                                             KNN_smooth = True)

    t = trans_t
    trans_prob = np.array([[1-2*t, t, t],[t, 1-2*t, t],[t, t, 1-2*t]])
    
    merge_Xdata = xclone.model.XHMM_smoothing(merge_Xdata, 
                                              start_prob = start_prob,  
                                              trans_prob = trans_prob, 
                                              emm_inlayer = "bin_phased_BAF_specific_center_emm_prob_log_KNN", 
                                              nproc = 80, 
                                              verbose = False)
    merge_Xdata.write(BAF_final_file)

    end_time = datetime.now(timezone.utc)
    time_passed = end_time - start_time
    main_logger.info("XClone BAF module finished (%d seconds)" % (time_passed.total_seconds()))

    if xclone_plot:
        if plot_cell_anno_key is None:
            plot_cell_anno_key = cell_anno_key
        run_BAF_plot(merge_Xdata, dataset_name, plot_cell_anno_key, out_dir)

    return merge_Xdata

def run_BAF_plot(merge_Xdata,
                 dataset_name,
                 plot_cell_anno_key,
                 out_dir = None):
    """
    """
    ## Result output prepare
    if out_dir is None:
        cwd = os.getcwd()
        out_dir = cwd + "/XCLONE_OUT/"
    
    out_plot_dir = str(out_dir) + "plot/"
    xclone.al.dir_make(out_plot_dir)

    # fig_title = ""
    baf_smooth_fig = out_plot_dir + dataset_name + "_BAF_smooth.png"
    baf_final_fig = out_plot_dir + dataset_name + "_BAF_CNV.png"

    sub_logger = get_logger("BAF plot module")
    sub_logger.info("BAF plot module started")
    start_time = datetime.now(timezone.utc)
    
    xclone.pl.BAF_smooth_visualization(merge_Xdata, Xlayer = "BAF_phased_KNN_WMA", 
                                       cell_anno_key = plot_cell_anno_key, 
                                       change_colorbar = False, 
                                       colorbar_name = "BAF",
                                       save_file = True, 
                                       out_file = baf_smooth_fig)
    xclone.pl.CNV_visualization2(merge_Xdata, weights = False, 
                                 cell_anno_key = plot_cell_anno_key, 
                                 save_file = True, 
                                 out_file = baf_final_fig)
    
    end_time = datetime.now(timezone.utc)
    time_passed = end_time - start_time
    sub_logger.info("BAF plot module finished (%d seconds)" % (time_passed.total_seconds()))
    return None