"""Base pipeline for XClone combination module."""

# Author: Rongting Huang
# Date: 2022-11-17
# Update: 2026-01-12 by Jiamu James Qiao

import gc
import os
from datetime import datetime, timezone

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
from xclone.model.clustering import (
    refine_clones_bayesian,
    tumor_classify,
    xclone_subclonal_analysis,
)


def _select_baf_layer(BAF_merge_Xdata, baf_denoise: bool) -> str:
    """Choose the appropriate BAF layer based on availability and config."""
    if baf_denoise:
        if "denoised_posterior_mtx" in BAF_merge_Xdata.layers:
            return "denoised_posterior_mtx"
        print("[XClone] warning: please provide denoised BAF posterior_mtx layer")
    return "posterior_mtx"


def _combine_probability_pipeline(
    combine_Xdata,
    rdr_denoise: bool,
    copyloss_correct: bool,
    copyloss_correct_mode: int,
    copygain_correct: bool,
    copygain_correct_mode,
    RDR_prior: bool,
):
    """Run RDR denoising (optional) then combine RDR/BAF probabilities."""
    rdr_layer = "posterior_mtx"
    if rdr_denoise:
        combine_Xdata = xclone.model.denoise_rdr_by_baf(
            combine_Xdata,
            RDR_layer="posterior_mtx",
            BAF_layer="BAF_extend_post_prob",
            out_RDR_layer="rdr_posterior_mtx_denoised",
        )
        rdr_layer = "rdr_posterior_mtx_denoised"

    return xclone.model.CNV_prob_combination(
        combine_Xdata,
        RDR_layer=rdr_layer,
        BAF_layer="BAF_extend_post_prob",
        copyloss_correct=copyloss_correct,
        copyloss_correct_mode=copyloss_correct_mode,
        copygain_correct=copygain_correct,
        copygain_correct_mode=copygain_correct_mode,
        RDR_prior=RDR_prior,
    )


def _run_wgd_detection(combine_Xdata, config) -> None:
    """Trigger WGD detection with config thresholds."""
    print("[XClone WGD detection performing]")
    xclone.model.WGD_warning(
        combine_Xdata,
        Xlayer="combine_base_prob",
        genome_level=config.WGD_detect_genome_level,
        prop_value_threshold=config.WGD_prop_value_threshold,
        cell_prop_threshold=config.WGD_cell_prop_threshold,
    )



def run_combine(RDR_Xdata,
                BAF_merge_Xdata,
                verbose = True,
                run_verbose = True,
                config_file = None):
    """
    Run the Combine (RDR & BAF) analysis on the provided annotated data.

    This function performs the combine analysis using the provided annotated data (`RDR_adata` & `BAF_merge_Xdata`).
    It allows for verbose output and the use of a custom configuration object. If no configuration
    object is specified, default settings from XClone's Combine module will be used.

    Parameters
    ----------

        RDR_Xdata : anndata.AnnData
            The annotated data matrix output by the RDR moudle analysis.
        BAF_merge_Xdata : anndata.AnnData
            The annotated data matrix output by the BAF moudle analysis.
        verbose : bool, optional
            If True, prints detailed information about the process. Default is True.
        run_verbose : bool, optional
            If True, provides verbose output during the run. Default is True.
        config_file : xclone.XCloneConfig or None, optional
            The XClone configuration object. If None, the default settings in XClone-RDR will be used.
            Default is None.
            The configuration can be created as follows:
            
            .. code-block:: python

                config_file = xclone.XCloneConfig(dataset_name="your_dataset_name", module="Combine")

    Returns
    -------

        combine_Xdata : anndata.AnnData
            The finalized output with multiple layers of information in the `anndata.AnnData` from Combine module.


    Example
    -------

        .. code-block:: python

            import xclone
            
            # Run Combine analysis with default settings
            combine_Xdata = xclone.model.run_combine(RDR_Xdata, BAF_merge_Xdata, verbose=True, run_verbose=True)
            
            # Run Combine analysis with a custom configuration object
            xconfig = xclone.XCloneConfig(dataset_name = dataset_name, module = "Combine")
            xconfig.outdir = out_dir
            #... other specified parameters in `xconfig`
            xconfig.display()
            combine_Xdata = xclone.model.run_combine(RDR_Xdata, BAF_merge_Xdata, verbose=True, run_verbose=True, config_file=xconfig)
    
    """
    config = load_config("Combine", config_file)
    configure_warnings(config.warninig_ignore)

    out_dir, out_data_dir, _ = resolve_output_dirs(config.outdir)
    dataset_name = config.dataset_name
    main_logger = get_logger("Main combine module")
    main_logger.info("XClone combine module Started")
    start_time = datetime.now(timezone.utc)

    if run_verbose:
        print("[XClone Combination module running]************************")

    if config.exclude_XY:
        RDR_Xdata = xclone.pp.exclude_XY_adata(RDR_Xdata)
        BAF_merge_Xdata = xclone.pp.exclude_XY_adata(BAF_merge_Xdata)
        print("[XClone warning] Combine module exclude chr XY analysis.")

    BAF_use_Xlayer = _select_baf_layer(BAF_merge_Xdata, config.BAF_denoise)

    RDR_Xdata, BAF_merge_Xdata = xclone.pp.check_RDR_BAF_samecellnumber(
        RDR_Xdata, BAF_merge_Xdata
    )
    combine_Xdata = xclone.model.bin_to_gene_mapping(
        BAF_merge_Xdata,
        RDR_Xdata,
        Xlayer=BAF_use_Xlayer,
        extend_layer="BAF_extend_post_prob",
        return_prob=False,
    )

    del RDR_Xdata
    if not config.clustering:
        del BAF_merge_Xdata
    gc.collect()

    combine_Xdata = _combine_probability_pipeline(
        combine_Xdata,
        rdr_denoise=config.RDR_denoise,
        copyloss_correct=config.copyloss_correct,
        copyloss_correct_mode=config.copyloss_correct_mode,
        copygain_correct=config.copygain_correct,
        copygain_correct_mode=config.copygain_correct_mode,
        RDR_prior=config.RDR_prior,
    )

    if config.WGD_detection:
        _run_wgd_detection(combine_Xdata, config)

    if config.develop_mode:
        corrected_path = os.path.join(out_data_dir, "combine_adata_corrected.h5ad")
        write_adata_safe(combine_Xdata, corrected_path, "combine_corrected_file")

    combine_Xdata = xclone.model.CNV_prob_merge_for_plot(
        combine_Xdata, Xlayer="corrected_prob"
    )

    if config.tumor_classification:
        print("[XClone tumor classification performing]")
        combine_Xdata = tumor_classify(
            combine_Xdata, config.tumor_classification_layer, out_data_dir
        )

    if config.clustering:
        print("[XClone clustering performing]")
        combine_Xdata = xclone_subclonal_analysis(
            combined_adata=combine_Xdata,
            baf_adata=BAF_merge_Xdata,
            method="combined",
            n_clones=config.n_clones,
            out_dir=out_data_dir,
            sample_name=dataset_name,
        )
        combine_Xdata = refine_clones_bayesian(
            adata=combine_Xdata,
            initial_col="clone_id",
            prob_layer="prob1_merge",
            n_iter=15,
            alpha=20.0,  # higher = smoother, less overfitting
            min_cells=50,
            n_clones=config.n_clones,
            out_dir=out_data_dir,
            sample_name=dataset_name,
        )

    del BAF_merge_Xdata
    gc.collect()

    final_file = os.path.join(out_data_dir, "combined_final.h5ad")
    write_adata_safe(combine_Xdata, final_file, "combine_final_file")

    log_duration(main_logger, start_time, "XClone combine module")

    if config.xclone_plot:
        plot_cell_anno_key = config.plot_cell_anno_key or config.cell_anno_key
        if run_verbose:
            print("[XClone plotting]")
        run_combine_plot(
            combine_Xdata,
            dataset_name,
            plot_cell_anno_key,
            config.merge_loss,
            config.merge_loh,
            config.set_figtitle,
            out_dir,
            config.customizedplotting,
        )
        if config.develop_mode and config.RDR_denoise:
            rdr_denoise_plot(
                combine_Xdata,
                dataset_name,
                plot_cell_anno_key,
                config.set_figtitle,
                out_dir,
            )

    return combine_Xdata

def rdr_denoise_plot(combine_Xdata,
                     dataset_name,
                     plot_cell_anno_key,
                     set_figtitle = True,
                     out_dir = None,
                     **kwargs):
    """
    plotting RDR_Denoised.
    """
    _, _, out_plot_dir = resolve_output_dirs(out_dir)

    rdr_final_denoise_fig = os.path.join(out_plot_dir, f"{dataset_name}_RDR_CNV_denoise.png")
    fig_title = f"{dataset_name} RDR_CNV_denoise" if set_figtitle else ""

    xclone.pl.CNV_visualization(
        combine_Xdata,
        Xlayer="rdr_posterior_mtx_denoised",
        states_weight=np.array([1, 2, 3]),
        weights=True,
        cell_anno_key=plot_cell_anno_key,
        title=fig_title,
        save_file=True,
        out_file=rdr_final_denoise_fig,
        **kwargs,
    )

def run_combine_plot(combine_Xdata,
            dataset_name,
            plot_cell_anno_key,
            merge_loss = True,
            merge_loh = True,
            set_figtitle = True,
            out_dir = None,
            customizedplotting = False,
            **kwargs):
    """
    """
    _, _, out_plot_dir = resolve_output_dirs(out_dir)

    fig_title = ""
    combine_res_base_fig = os.path.join(out_plot_dir, f"{dataset_name}_combine_base.png")
    combine_res_refined_fig = os.path.join(out_plot_dir, f"{dataset_name}_combine_refined.png")
    combine_res_select_fig = os.path.join(out_plot_dir, f"{dataset_name}_combine_select.png")

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

    # cluster refined plot
    if "prob1_merge_refined" in combine_Xdata.layers:
        xclone.pl.Combine_CNV_visualization(combine_Xdata, Xlayer = "prob1_merge_refined", 
            cell_anno_key = plot_cell_anno_key,  
            title = fig_title,
            save_file = True, out_file = combine_res_refined_fig,
            **kwargs)
        
    if customizedplotting:
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