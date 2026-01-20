"""Base pipeline for XClone BAF module."""

# Author: Rongting Huang
# Date: 2022-10-21
# Update: 2026-01-13 by Jiamu James Qiao

import gc
import os
from datetime import datetime, timezone
from typing import Optional, Tuple

import anndata as an
import numpy as np
import xclone
from .._logging import get_logger
from ._pipeline_utils import (  # type: ignore
    configure_warnings,
    load_config,
    log_duration,
    resolve_output_dirs,
    write_adata_safe,
)


def _validate_ref_celltype(adata, cell_anno_key: str, ref_celltype) -> None:
    """Validate that ref_celltype exists in the annotation."""
    if isinstance(ref_celltype, list):
        missing_items = [
            item for item in ref_celltype if item not in adata.obs[cell_anno_key].values
        ]
        if missing_items:
            raise ValueError(
                f"[XClone error] Items {missing_items} not found in the adata's annotation."
            )
    else:
        if ref_celltype not in adata.obs[cell_anno_key].values:
            raise ValueError(
                f"[XClone error] Item '{ref_celltype}' not found in the adata's annotation."
            )


def _setup_rdr_integration(
    config, out_data_dir: str, merge_Xdata=None, RDR_adata=None
) -> Tuple[Optional[str], bool, bool, str, Optional[any]]:
    """Setup RDR integration parameters and handle KNN connectivities.
    
    Determines RDR file path and, if merge_Xdata and RDR_adata are provided, checks for 
    connectivities in RDR_adata.obsp and updates merge_Xdata with expression-based KNN 
    connectivities if available.
    
    Parameters
    ----------
    config : Config object
        Configuration object with RDR settings
    out_data_dir : str
        Output data directory path
    merge_Xdata : Optional[AnnData]
        BAF merge data (optional, used for connectivity setup)
    RDR_adata : Optional[AnnData]
        RDR annotated data (optional, used for connectivity setup)
    
    Returns
    -------
    Tuple containing:
        - RDR file path (or None)
        - remove_marker_genes flag
        - get_BAF_KNN_connectivities flag
        - KNN_connect_use_key (updated if connectivities found)
        - merge_Xdata (optionally updated with expression-based connectivities)
    """
    rdr_file = None
    remove_marker_genes = False
    knn_key = "connectivities"
    
    if config.update_info_from_rdr or config.RDR_file is not None:
        # Determine RDR file path
        if config.RDR_file is None:
            rdr_file = os.path.join(out_data_dir, "RDR_adata_KNN_HMM_post.h5ad")
        else:
            rdr_file = config.RDR_file
        
        remove_marker_genes = True
        
        # Check for connectivities in RDR_adata if both data objects are provided
        if merge_Xdata is not None and RDR_adata is not None:
            if "connectivities" in RDR_adata.obsp:
                knn_key = "connectivities_expr"
                merge_Xdata = xclone.model.get_KNN_connectivities_from_expr(merge_Xdata, RDR_adata)
    
    return rdr_file, remove_marker_genes, True, knn_key, merge_Xdata


def _perform_baf_phasing(
    BAF_adata,
    phasing_region_key: str,
    phasing_len: int,
    bin_nproc: int,
    feature_mode: str,
    HMM_brk: str,
    extreme_count_cap: bool,
) -> Tuple:
    """Perform BAF local and global phasing."""
    BAF_var_add = None if HMM_brk in ["chr", "chr_arm"] else HMM_brk
    BAF_adata, merge_Xdata = xclone.model.BAF_Local_phasing(
        BAF_adata,
        region_key=phasing_region_key,
        phasing_len=phasing_len,
        bin_nproc=bin_nproc,
        feature_mode=feature_mode,
        var_add=BAF_var_add,
    )

    BAF_adata, merge_Xdata = xclone.model.BAF_Global_phasing(BAF_adata, merge_Xdata)
    BAF_adata, merge_Xdata = xclone.model.BAF_Global_phasing_rev(BAF_adata, merge_Xdata)

    if extreme_count_cap:
        merge_Xdata = xclone.model.extrme_count_capping(merge_Xdata)

    merge_Xdata = xclone.pl.calculate_cell_BAF(
        merge_Xdata, AD_key="ad_bin", DP_key="dp_bin", BAF_key="BAF"
    )
    merge_Xdata = xclone.pl.calculate_cell_BAF(
        merge_Xdata,
        AD_key="ad_bin_phased_rev",
        DP_key="dp_bin",
        BAF_key="BAF_phased",
    )
    xclone.model.BAF_fillna(
        merge_Xdata, Xlayer="BAF_phased", out_layer="fill_BAF_phased"
    )
    return BAF_adata, merge_Xdata


def _perform_baf_smoothing(
    merge_Xdata,
    KNN_Xlayer: str,
    KNN_neighbors: int,
    KNN_npcs: int,
    get_BAF_KNN_connectivities: bool,
    KNN_connect_use_key: str,
    WMA_window_size: int,
    WMA_smooth_key: str,
):
    """Perform KNN and WMA smoothing on BAF data."""
    merge_Xdata = xclone.model.extra_preprocess_BAF(
        merge_Xdata,
        Xlayer=KNN_Xlayer,
        KNN_neighbors=KNN_neighbors,
        KNN_npcs=KNN_npcs,
        run_KNN=get_BAF_KNN_connectivities,
        copy=True,
    )

    merge_Xdata = xclone.model.KNN_smooth(
        merge_Xdata,
        run_KNN=False,
        KNN_Xlayer=None,
        KNN_connect_use=KNN_connect_use_key,
        layer="fill_BAF_phased",
        out_layer="BAF_phased_KNN",
    )
    merge_Xdata = xclone.model.WMA_smooth(
        merge_Xdata,
        layer="BAF_phased_KNN",
        out_layer="BAF_phased_KNN_WMA",
        window_size=WMA_window_size,
        chrom_key=WMA_smooth_key,
        verbose=False,
    )

    merge_Xdata = xclone.model.WMA_smooth(
        merge_Xdata,
        layer="fill_BAF_phased",
        out_layer="BAF_phased_WMA",
        window_size=WMA_window_size,
        chrom_key=WMA_smooth_key,
        verbose=False,
    )
    return merge_Xdata


def _get_guide_theo_states(
    merge_Xdata, guide_theo_CNV_states, CNV_N_components: int, Xlayer: str = "BAF_phased_WMA"
):
    """Get guide theoretical CNV states."""
    if guide_theo_CNV_states is not None:
        return guide_theo_CNV_states
    CNV_states = xclone.model.get_CNV_states(
        merge_Xdata,
        Xlayer=Xlayer,
        n_components=CNV_N_components,
        means_init=None,
        max_iter=200,
    )
    return xclone.model.guide_BAF_theo_states(CNV_states)


def _calculate_ref_prop(merge_Xdata, cell_anno_key: str, ref_celltype) -> float:
    """Calculate reference cell proportion."""
    if isinstance(ref_celltype, list):
        ref_cell_num = merge_Xdata.obs[cell_anno_key].isin(ref_celltype).sum()
    else:
        ref_cell_num = (merge_Xdata.obs[cell_anno_key] == ref_celltype).sum()
    total_cell_num = merge_Xdata.obs.shape[0]
    return ref_cell_num / total_cell_num


def _setup_baf_ref_states(
    merge_Xdata,
    guide_theo_states,
    theo_neutral_BAF,
    cell_anno_key: str,
    ref_celltype,
    ref_BAF_clip: bool,
) -> np.ndarray:
    """Setup BAF reference states and return used specific states."""
    if theo_neutral_BAF is not None:
        merge_Xdata.var["theo_neutral_BAF"] = theo_neutral_BAF
        return xclone.model.gene_specific_BAF(
            merge_Xdata, theo_states=guide_theo_states, specific_BAF="theo_neutral_BAF"
        )

    ref_prop = _calculate_ref_prop(merge_Xdata, cell_anno_key, ref_celltype)
    if ref_prop <= 0.01:
        merge_Xdata = xclone.model.get_BAF_ref_limited(
            merge_Xdata,
            Xlayer="BAF_phased_KNN_WMA",
            out_anno="ref_BAF_phased",
            anno_key=cell_anno_key,
            ref_cell=ref_celltype,
            clipping=ref_BAF_clip,
        )
    else:
        merge_Xdata = xclone.model.get_BAF_ref(
            merge_Xdata,
            Xlayer="fill_BAF_phased",
            out_anno="ref_BAF_phased",
            anno_key=cell_anno_key,
            ref_cell=ref_celltype,
            clipping=ref_BAF_clip,
        )

    return xclone.model.gene_specific_BAF(
        merge_Xdata, theo_states=guide_theo_states, specific_BAF="ref_BAF_phased"
    )


def _compute_baf_emm_and_hmm(
    merge_Xdata,
    used_specific_states: np.ndarray,
    BAF_states_num: int,
    gene_specific_concentration: bool,
    concentration: float,
    KNN_connect_use_key: str,
    trans_t: float,
    trans_prob,
    CNV_N_components: int,
    HMM_brk: str,
    start_prob,
    HMM_nproc: int,
):
    """Compute BAF emission probabilities and apply HMM smoothing."""
    merge_Xdata = xclone.model.calculate_Xemm_prob_bb(
        merge_Xdata,
        AD_key="ad_bin_phased_rev",
        DP_key="dp_bin",
        outlayer="bin_phased_BAF_specific_center_emm_prob_log",
        states=used_specific_states,
        states_num=BAF_states_num,
        gene_specific=gene_specific_concentration,
        concentration=concentration,
    )

    merge_Xdata = xclone.model.BAF_smoothing(
        merge_Xdata,
        inlayer="bin_phased_BAF_specific_center_emm_prob_log",
        outlayer="bin_phased_BAF_specific_center_emm_prob_log_KNN",
        KNN_connectivities_key=KNN_connect_use_key,
        KNN_smooth=True,
    )

    if trans_prob is None:
        if CNV_N_components == 3:
            trans_prob = np.array([
                [1 - 2 * trans_t, trans_t, trans_t],
                [trans_t, 1 - 2 * trans_t, trans_t],
                [trans_t, trans_t, 1 - 2 * trans_t],
            ])
        elif CNV_N_components == 5:
            trans_prob = np.array([
                [1 - 4 * trans_t, trans_t, trans_t, trans_t, trans_t],
                [trans_t, 1 - 4 * trans_t, trans_t, trans_t, trans_t],
                [trans_t, trans_t, 1 - 4 * trans_t, trans_t, trans_t],
                [trans_t, trans_t, trans_t, 1 - 4 * trans_t, trans_t],
                [trans_t, trans_t, trans_t, trans_t, 1 - 4 * trans_t],
            ])

    merge_Xdata = xclone.model.XHMM_smoothing(
        merge_Xdata,
        brk=HMM_brk,
        start_prob=start_prob,
        trans_prob=trans_prob,
        emm_inlayer="bin_phased_BAF_specific_center_emm_prob_log_KNN",
        nproc=HMM_nproc,
        verbose=False,
    )
    return merge_Xdata


def _get_denoise_indices(CNV_N_components: int) -> Tuple[int, list]:
    """Get neutral and CNV indices for denoising based on CNV components."""
    if CNV_N_components == 3:
        return 1, [0, 2]
    elif CNV_N_components == 5:
        return 2, [0, 1, 3, 4]
    raise ValueError(f"Unsupported CNV_N_components: {CNV_N_components}")


def _compute_baf_add_states(
    merge_Xdata,
    guide_theo_CNV_states,
    theo_neutral_BAF,
    cell_anno_key: str,
    ref_celltype,
    ref_BAF_clip: bool,
    gene_specific_concentration: bool,
    concentration: float,
    KNN_connect_use_key: str,
    trans_t: float,
    HMM_brk: str,
    HMM_nproc: int,
) -> None:
    """Compute additional 3-state BAF analysis for comparison (5-state mode only)."""
    merge_Xdata_copy = merge_Xdata.copy()

    guide_theo_states = _get_guide_theo_states(merge_Xdata_copy, guide_theo_CNV_states, 3)
    if guide_theo_CNV_states is not None:
        guide_theo_states = guide_theo_CNV_states[[0, 3]]
        print("User specify guide_theo_states: ", guide_theo_states)

    used_specific_states = _setup_baf_ref_states(
        merge_Xdata_copy,
        guide_theo_states,
        theo_neutral_BAF,
                cell_anno_key,
                ref_celltype,
        ref_BAF_clip,
    )

    merge_Xdata_copy = xclone.model.calculate_Xemm_prob_bb(
        merge_Xdata_copy,
        AD_key="ad_bin_phased_rev",
        DP_key="dp_bin",
        outlayer="correct_emm_prob_log",
        states=used_specific_states,
        gene_specific=gene_specific_concentration,
        concentration=concentration,
    )
    merge_Xdata_copy = xclone.model.BAF_smoothing(
        merge_Xdata_copy,
        inlayer="correct_emm_prob_log",
        outlayer="correct_emm_prob_log_KNN",
        KNN_connectivities_key=KNN_connect_use_key,
        KNN_smooth=True,
    )

    start_prob_3state = np.array([0.3, 0.4, 0.3])
    trans_prob_3state = np.array([
        [1 - 2 * trans_t, trans_t, trans_t],
        [trans_t, 1 - 2 * trans_t, trans_t],
        [trans_t, trans_t, 1 - 2 * trans_t],
    ])
    merge_Xdata_copy = xclone.model.XHMM_smoothing(
        merge_Xdata_copy,
        brk=HMM_brk,
        start_prob=start_prob_3state,
        trans_prob=trans_prob_3state,
        emm_inlayer="correct_emm_prob_log_KNN",
        nproc=HMM_nproc,
        verbose=False,
    )

    neutral_index, cnv_index = _get_denoise_indices(3)
    merge_Xdata_copy = xclone.model.denoise_gene_scale(
        merge_Xdata_copy,
        Xlayer="posterior_mtx",
        neutral_index=neutral_index,
        cnv_index=cnv_index,
        GMM_detection=True,
        gmm_comp=2,
        cell_prop_cutoff=0.05,
        out_layer="denoised_posterior_mtx",
    )

    merge_Xdata.layers["add_posterior_mtx"] = merge_Xdata_copy.layers["posterior_mtx"]
    merge_Xdata.layers["denoised_add_posterior_mtx"] = merge_Xdata_copy.layers[
        "denoised_posterior_mtx"
    ]


def preview_BAF(BAF_adata, cell_anno_key, ref_celltype):
    """
    Function:
    some precheck for params setting.
    """
    # xclone.model.view_celltype(BAF_adata, cell_anno_key = cell_anno_key)

    ## smooth plot?
    return None


def run_BAF(BAF_adata, verbose = True, run_verbose = True, config_file = None):
    """
    Run the BAF (B Allele Frequency) analysis on the provided annotated data.

    This function performs the BAF analysis using the provided annotated data (`BAF_adata`).
    It allows for verbose output and the use of a custom configuration object. If no configuration
    object is specified, default settings from XClone's BAF module will be used.

    Parameters
    ----------

        BAF_adata : anndata.AnnData
            The annotated data matrix on which the BAF analysis will be performed.
        verbose : bool, optional
            If True, prints detailed information about the process. Default is True.
        run_verbose : bool, optional
            If True, provides verbose output during the run. Default is True.
        config_file : xclone.XCloneConfig or None, optional
            The XClone configuration object. If None, the default settings in XClone-BAF will be used.
            Default is None.
            The configuration can be created as follows:
            
            .. code-block:: python

                config_file = xclone.XCloneConfig(dataset_name="your_dataset_name", module="BAF")

    Returns
    -------

        merge_Xdata : anndata.AnnData
            The finalized output with multiple layers of information in the `anndata.AnnData` from BAF module.
            It is at bin level.


    Example
    -------

        .. code-block:: python

            import xclone
            
            # Run BAF analysis with default settings
            BAF_merge_Xdata = xclone.model.run_BAF(BAF_adata, verbose=True, run_verbose=True)
            
            # Run BAF analysis with a custom configuration object
            xconfig = xclone.XCloneConfig(dataset_name = dataset_name, module = "BAF")
            xconfig.outdir = out_dir
            #... other specified parameters in `xconfig`
            xconfig.display()
            BAF_merge_Xdata = xclone.model.run_BAF(BAF_adata, verbose=True, run_verbose=True, config_file=xconfig)
    
    """
    config = load_config("BAF", config_file)
    configure_warnings(config.warninig_ignore)
    
    # general settings
    dataset_name = config.dataset_name
    out_dir, out_data_dir, _ = resolve_output_dirs(config.outdir)

    cell_anno_key = config.cell_anno_key
    ref_celltype = config.ref_celltype
    exclude_XY = config.exclude_XY

    ## RDR related
    update_info_from_rdr = config.update_info_from_rdr
    RDR_file = config.RDR_file
    remove_marker_genes = config.remove_marker_genes
    # ## KNN related
    KNN_connect_use_key = config.KNN_connect_use_key
    get_BAF_KNN_connectivities = config.get_BAF_KNN_connectivities
    KNN_Xlayer = config.KNN_Xlayer
    KNN_neighbors = config.KNN_neighbors
    KNN_npcs = config.KNN_npcs
    
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

    HMM_nproc = config.HMM_nproc

    # HMM settings
    start_prob = config.start_prob
    trans_t = config.trans_t
    trans_prob = config.trans_prob
    HMM_brk = config.HMM_brk
    

    # plot settings
    xclone_plot = config.xclone_plot
    plot_cell_anno_key = config.plot_cell_anno_key

    # develop mode settings
    develop_mode = config.develop_mode
    
    
    ### output before CNV calling
    BAF_base_file = os.path.join(out_data_dir, "BAF_base_Xdata.h5ad")
    BAF_merge_base_file = os.path.join(out_data_dir, "BAF_merge_base_Xdata.h5ad")
    ### output after CNV calling
    BAF_final_file = os.path.join(out_data_dir, "BAF_merge_Xdata_KNN_HMM_post.h5ad")
        

    ##------------------------
    main_logger = get_logger("Main BAF module")
    main_logger.info("XClone BAF module Started")
    start_time = datetime.now(timezone.utc)

    if run_verbose:
        print("[XClone BAF module running]************************")

    _validate_ref_celltype(BAF_adata, cell_anno_key, ref_celltype)

    if exclude_XY:
        BAF_adata = xclone.pp.exclude_XY_adata(BAF_adata)
        print("[XClone warning] BAF module exclude chr XY analysis.")
    BAF_adata = xclone.pp.check_BAF(BAF_adata, cell_anno_key=cell_anno_key, verbose=verbose)

    # Load RDR data if RDR integration is enabled
    RDR_adata = None
    remove_marker_genes = False
    if config.update_info_from_rdr:
        # Get RDR file path from _setup_rdr_integration
        rdr_file, remove_marker_genes, _, _, _ = _setup_rdr_integration(config, out_data_dir)
        if rdr_file is not None:
            RDR_adata = an.read_h5ad(rdr_file)
            BAF_adata = BAF_adata[BAF_adata.obs.index.isin(RDR_adata.obs.index), :]
            xclone.pp.check_RDR_BAF_cellorder(RDR_adata, BAF_adata)

    if remove_marker_genes and RDR_adata is not None:
        marker_genes = RDR_adata.uns["rank_marker_genes"]
        BAF_adata = xclone.pp.Xdata_region_selection(
            BAF_adata,
            select=False,
            chr_anno_key="GeneName",
            chr_lst=marker_genes,
            update_uns=False,
            uns_anno_key=None,
        )

    BAF_adata, merge_Xdata = _perform_baf_phasing(
        BAF_adata,
        phasing_region_key=phasing_region_key,
        phasing_len=phasing_len,
        bin_nproc=bin_nproc,
        feature_mode=feature_mode,
        HMM_brk=HMM_brk,
        extreme_count_cap=extreme_count_cap,
    )

    # Setup RDR integration and handle KNN connectivities
    rdr_file, remove_marker_genes, get_BAF_KNN_connectivities, KNN_connect_use_key, merge_Xdata = _setup_rdr_integration(
        config, out_data_dir, merge_Xdata, RDR_adata
    )
    ## clean up memory
    if RDR_adata is not None:
        del RDR_adata
        gc.collect()

    merge_Xdata = _perform_baf_smoothing(
        merge_Xdata,
        KNN_Xlayer,
        KNN_neighbors,
        KNN_npcs,
        get_BAF_KNN_connectivities,
        KNN_connect_use_key,
        WMA_window_size,
        WMA_smooth_key,
    )

    if develop_mode:
        write_adata_safe(BAF_adata, BAF_base_file, "BAF_base_file")
        write_adata_safe(merge_Xdata, BAF_merge_base_file, "BAF_merge_base_file")

    del BAF_adata
    gc.collect()

    guide_theo_states = _get_guide_theo_states(
        merge_Xdata, guide_theo_CNV_states, CNV_N_components
    )
    if guide_theo_CNV_states is not None:
        print("User specify guide_theo_states: ", guide_theo_states)

    used_specific_states = _setup_baf_ref_states(
        merge_Xdata,
        guide_theo_states,
        theo_neutral_BAF,
        cell_anno_key,
        ref_celltype,
        ref_BAF_clip,
    )

    if gene_specific_concentration:
        concentration_lower = config.concentration_lower
        concentration_upper = config.concentration_upper
        merge_Xdata = xclone.model.concentration_mapping(
            merge_Xdata, concentration_lower, concentration_upper
        )

    merge_Xdata = _compute_baf_emm_and_hmm(
        merge_Xdata,
        used_specific_states,
        BAF_states_num,
        gene_specific_concentration,
        concentration,
        KNN_connect_use_key,
        trans_t,
        trans_prob,
        CNV_N_components,
        HMM_brk,
        start_prob,
        HMM_nproc,
    )

    if CNV_N_components == 5:
        if BAF_add is None:
            BAF_add = develop_mode
    else:
        BAF_add = False

    if BAF_add:
        _compute_baf_add_states(
            merge_Xdata,
            guide_theo_CNV_states,
            theo_neutral_BAF,
            cell_anno_key,
            ref_celltype,
            ref_BAF_clip,
            gene_specific_concentration,
            concentration,
            KNN_connect_use_key,
            trans_t,
            HMM_brk,
            HMM_nproc,
        )

    if BAF_denoise:
        neutral_index, cnv_index = _get_denoise_indices(CNV_N_components)
        merge_Xdata = xclone.model.denoise_gene_scale(
            merge_Xdata,
            Xlayer="posterior_mtx",
            neutral_index=neutral_index,
            cnv_index=cnv_index,
            GMM_detection=GMM_detect,
            gmm_comp=BAF_denoise_GMM_comp,
            cell_prop_cutoff=BAF_denoise_cellprop_cutoff,
            out_layer="denoised_posterior_mtx",
        )
    if not develop_mode:
        layers_to_keep = [
            'BAF', 'BAF_phased', 'fill_BAF_phased', 'BAF_phased_KNN', 'BAF_phased_KNN_WMA',
            'BAF_phased_WMA', 'posterior_mtx', 'posterior_mtx_log', 'add_posterior_mtx',
            'denoised_add_posterior_mtx', 'denoised_posterior_mtx'
        ]
        merge_Xdata = xclone.pp.keep_layers(merge_Xdata, layers_to_keep)
    write_adata_safe(merge_Xdata, BAF_final_file, "BAF_final_file")

    log_duration(main_logger, start_time, "XClone BAF module")

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
    _, _, out_plot_dir = resolve_output_dirs(out_dir)

    fig_title = ""
    baf_smooth_fig = os.path.join(out_plot_dir, f"{dataset_name}_BAF_smooth.png")
    baf_final_fig1 = os.path.join(out_plot_dir, f"{dataset_name}_BAF_CNV.png")
    baf_final_fig2 = os.path.join(out_plot_dir, f"{dataset_name}_BAF_CNV_denoise.png")
    baf_final_fig3 = os.path.join(out_plot_dir, f"{dataset_name}_BAF_CNV_3states.png")
    baf_final_fig4 = os.path.join(out_plot_dir, f"{dataset_name}_BAF_CNV_3states_denoise.png")

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

    
    log_duration(sub_logger, start_time, "BAF plot module")
    return None


def run_BAF_complex_plot(merge_Xdata,
                 dataset_name,
                 plot_cell_anno_key,
                 set_figtitle = True,
                 out_dir = None,
                 **kwargs):
    """
    """
    _, _, out_plot_dir = resolve_output_dirs(out_dir)

    fig_title = ""
    baf_smooth_fig = os.path.join(out_plot_dir, f"{dataset_name}_BAF_smooth.png")
    baf_final_fig1 = os.path.join(out_plot_dir, f"{dataset_name}_BAF_CNV.png")
    baf_final_fig2 = os.path.join(out_plot_dir, f"{dataset_name}_BAF_CNV_denoise.png")
    baf_final_fig3 = os.path.join(out_plot_dir, f"{dataset_name}_BAF_CNV_3states.png")
    baf_final_fig4 = os.path.join(out_plot_dir, f"{dataset_name}_BAF_CNV_3states_denoise.png")

    sub_logger = get_logger("BAF plot module")
    sub_logger.info("BAF plot module started")
    start_time = datetime.now(timezone.utc)

    if set_figtitle:
        fig_title = dataset_name + " BAF smooth"
    xclone.pl.BAF_smooth_complex_visualization(merge_Xdata, Xlayer = "BAF_phased_KNN_WMA", 
                                       cell_anno_key = plot_cell_anno_key, 
                                       clusters_display_name = plot_cell_anno_key,
                                       change_colorbar = False, 
                                       colorbar_name = "BAF",
                                       title = fig_title,
                                       save_file = True, 
                                       out_file = baf_smooth_fig,
                                       **kwargs)
    if set_figtitle:
        fig_title = dataset_name + " BAF module"
    xclone.pl.Complex_BAF_CNV_visualization(merge_Xdata, weights = False, 
                                    cell_anno_key = plot_cell_anno_key, 
                                    clusters_display_name = plot_cell_anno_key, 
                                    title = fig_title,
                                    save_file = True, 
                                    out_file = baf_final_fig1,
                                    **kwargs)
                                    
    
    if "denoised_posterior_mtx" in merge_Xdata.layers:
        xclone.pl.Complex_BAF_CNV_visualization(merge_Xdata, Xlayer = "denoised_posterior_mtx",
                                        weights = False, 
                                        cell_anno_key = plot_cell_anno_key, 
                                        clusters_display_name = plot_cell_anno_key,
                                        title = fig_title + " (denoise)",
                                        save_file = True, 
                                        out_file = baf_final_fig2,
                                        **kwargs)
    
    if "add_posterior_mtx" in merge_Xdata.layers:
        xclone.pl.Complex_BAF_CNV_visualization(merge_Xdata, Xlayer = "add_posterior_mtx",
                                        weights = False, 
                                        cell_anno_key = plot_cell_anno_key, 
                                        clusters_display_name = plot_cell_anno_key,
                                        title = fig_title,
                                        save_file = True, 
                                        out_file = baf_final_fig3,
                                        **kwargs)
    if "denoised_add_posterior_mtx" in merge_Xdata.layers:
        xclone.pl.Complex_BAF_CNV_visualization(merge_Xdata, Xlayer = "denoised_add_posterior_mtx",
                                        weights = False, 
                                        cell_anno_key = plot_cell_anno_key, 
                                        clusters_display_name = plot_cell_anno_key,
                                        title = fig_title + " (denoise)",
                                        save_file = True, 
                                        out_file = baf_final_fig4,
                                        **kwargs)

    
    log_duration(sub_logger, start_time, "BAF plot module")
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
    _, _, out_plot_dir = resolve_output_dirs(out_dir)

    fig_title = ""
    
    baf_fig1 = os.path.join(out_plot_dir, f"{dataset_name}_BAF.png")
    baf_fig2 = os.path.join(out_plot_dir, f"{dataset_name}_BAF_phased.png")
    baf_fig3 = os.path.join(out_plot_dir, f"{dataset_name}_BAF_phased_KNN.png")
    baf_fig4 = os.path.join(out_plot_dir, f"{dataset_name}_BAF_phased_KNN_WMA.png")
    
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
