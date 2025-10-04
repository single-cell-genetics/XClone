"""
XClone
Base Configurations class.

Author: Rongting Huang
"""

#! Base Configuration Class
#! Don't use this class directly. 
#! Instead, sub-class it and override the configurations you need to change.

import inspect
from pathlib import Path
from time import time
from typing import Any, Union, Optional
from typing import Tuple

import numpy as np

def _type_check(var: Any, varname: str, types: Union[type, Tuple[type, ...]]):
    if isinstance(var, types):
        return
    if isinstance(types, type):
        possible_types_str = types.__name__
    else:
        type_names = [t.__name__ for t in types]
        possible_types_str = "{} or {}".format(
            ", ".join(type_names[:-1]), type_names[-1]
        )
    raise TypeError(f"{varname} must be of type {possible_types_str}")

## todo: try save dataset_name in uns
class Base_settings():
    def __init__(self):
        self.warninig_ignore = True

class PreprocessingConfig():
    """
    Configuration manager for preprocessing/loading XClone data.

    This class handles the configuration settings required for loading and preprocessing
    XClone data, including dataset details, module-specific configurations, file paths,
    and spatial settings.

    Attributes
    ----------

        module : str
            The module to configure ('pre_check', 'RDR', 'BAF', 'Combine').
        set_spatial : bool
            Flag to indicate if spatial data settings should be applied.
        rdr_data_dir : str or None
            Directory path for RDR data files.
        baf_data_dir : str or None
            Directory path for BAF data files.
        dataset_name : str
            The name of the dataset.
        RDR_file : str or None
            Path to the RDR matrix file.
        mtx_barcodes_file : str or None
            Path to the barcodes file for the matrix.
        regions_anno_file : str or None
            Path to the regions annotation file.
        cell_anno_file : str or None
            Path to the cell annotation file.
        barcodes_key : str
            Key for barcodes in the cell annotation file.
        anno_file_sep : str
            Separator used in the annotation file.
        cell_anno_key : str
            Key for cell type annotations.
        genome_mode : str
            Genome mode setting. One of 'hg38_genes', 'hg38_blocks', 'hg19_genes', 
            'hg19_blocks', or 'mm10_genes'. Default is 'hg38_genes'.
        spot_position_file : str or None
            Path to the spot position file for spatial data.
        AD_file : str or None
            Path to the AD matrix file for BAF data.
        DP_file : str or None
            Path to the DP matrix file for BAF data.
        RDR_adata_file : str or None
            Path to the processed RDR AnnData file. Specific for 'Combine' module.
        BAF_adata_file : str or None
            Path to the processed BAF AnnData file. Specific for 'Combine' module.
    """
    def __init__(
        self,
        dataset_name = "XClone_scDATA",
        module = "RDR",
        cell_anno_file = None,
        barcodes_key = "barcodes",
        anno_file_sep = ",", 
        set_spatial: bool = False,
        spot_position_file = None,
        rdr_data_dir = None,
        baf_data_dir = None):
        """
        Initialize the PreprocessingConfig class with various configuration parameters.

        Parameters
        ----------

            dataset_name : str, optional
                The name of the dataset (default is "XClone_scDATA").
            module : str, optional
                The module to configure, must be one of "pre_check", "RDR", "BAF", or "Combine" 
                (default is "RDR").
            cell_anno_file : str or None, optional
                Path to the cell annotation file (default is None).
            barcodes_key : str, optional
                Key for barcodes in the cell annotation file (default is "barcodes").
            anno_file_sep : str, optional
                Separator used in the annotation file (default is ",").
            set_spatial : bool, optional
                Flag to set spatial specific configurations (default is False).
            spot_position_file : str or None, optional
                Path to the spot position file for spatial data (default is None).
            rdr_data_dir : str or None, optional
                Directory path for RDR data files (default is None).
            baf_data_dir : str or None, optional
                Directory path for BAF data files (default is None).
        """

        
        self.module = module
        self.set_spatial = set_spatial
        self.rdr_data_dir = rdr_data_dir
        self.baf_data_dir = baf_data_dir
        self.dataset_name = dataset_name

        if self.module == "pre_check":
            self.RDR_barcodes_file = self.rdr_data_dir + "barcodes.tsv"
            self.BAF_barcodes_file = self.baf_data_dir + "xcltk.samples.tsv"
        
        if self.module == "RDR":
            self.RDR_file = self.rdr_data_dir + "matrix.mtx"
            self.mtx_barcodes_file = self.rdr_data_dir + "barcodes.tsv"
            self.regions_anno_file = None
            self.cell_anno_file = cell_anno_file
            self.barcodes_key = barcodes_key
            self.anno_file_sep = anno_file_sep
            self.cell_anno_key = "cell_type"
            self.genome_mode = "hg38_genes"
            
            if self.set_spatial:
                self.spot_position_file = spot_position_file

        if self.module == "BAF":
            self.AD_file = self.baf_data_dir + "xcltk.AD.mtx"
            self.DP_file = self.baf_data_dir + "xcltk.DP.mtx"
            self.mtx_barcodes_file = self.baf_data_dir + "xcltk.samples.tsv"
            self.regions_anno_file = None
            self.cell_anno_file = cell_anno_file
            self.barcodes_key = barcodes_key
            self.anno_file_sep = anno_file_sep
            self.cell_anno_key = "cell_type"
            self.genome_mode = "hg38_genes"
            
            if self.set_spatial:
                self.spot_position_file = spot_position_file

        if self.module == "Combine":
            self.RDR_adata_file = self.rdr_data_dir + "RDR_adata_KNN_HMM_post.h5ad"
            self.BAF_adata_file = self.baf_data_dir + "BAF_merge_Xdata_KNN_HMM_post.h5ad"
    
    def display(self):
        """Display Configuration values."""
        print(self.module, "\nConfigurations:")
        for a in dir(self):
            if not a.startswith("__") and not callable(getattr(self, a)):
                print("{:30} {}".format(a, getattr(self, a)))
        print("\n")    

class XCloneGeneral_config():
    """
    General configuration settings for XClone.

    This class manages general configuration settings for the XClone tool,
    including cell annotation keys, reference cell types, data exclusion options,
    KNN smoothing parameters, and plotting preferences.

    Attributes
    ----------

        cell_anno_key : str
            The key for cell type annotations.
        ref_celltype : str
            The reference cell type.
        exclude_XY : bool
            Flag to exclude XY chromosomes from analysis.
        remove_guide_XY : bool
            Flag to remove guide XY chromosomes.
        KNN_neighbors : int
            Number of neighbors for KNN smoothing.
        plot_remove_immune : bool
            Flag to remove immune cells during plotting.
        plot_immune_celltype : str or None
            The cell type to be considered as immune cells for plotting.
        plot_remove_reference : bool
            Flag to remove reference cells during plotting.
        plot_ref_celltype : str or None
            The cell type to be considered as reference cells for plotting.
        develop_mode : bool
            Flag to set develop mode or not (main difference in saved outputs).

    """

    def __init__(self):
        """
        Initialize the XCloneGeneral_config class with default configuration parameters.

        Parameters
        ----------

            None
        """
        self.cell_anno_key = "cell_type"
        self.ref_celltype = "N"
        self.exclude_XY = False
        self.remove_guide_XY = False
        # KNN smoothing
        self.KNN_neighbors = 10
        self.KNN_npcs = 40
        # Plotting
        # self.plot_cell_anno_key =  None
        self.plot_remove_immune = True
        self.plot_immune_celltype = None
        self.plot_remove_reference = True
        self.plot_ref_celltype = None
        # develop 
        self.develop_mode = False


class RDR_prev_General_config():
    """
    Previous version of RDR

    General configuration settings for RDR (Read Depth Ratio) analysis.

    This class manages the configuration settings for RDR analysis, particularly
    for 10X scRNA-seq data, including transformation options, filtering criteria,
    marker gene settings, GLM fitting settings, smoothing parameters, and plotting preferences.

    Attributes
    ----------

        smart_transform : bool
            Flag to apply smart transformation (default is False).
        filter_ref_ave : float
            Threshold for filtering based on reference average (default is 0.5).
        min_gene_keep_num : int
            Minimum number of genes to keep (default is 3000).
        multi_refcelltype : bool
            Flag to use multiple reference cell types.
        marker_group_anno_key : str or None
            Annotation key for marker groups.
        get_marker_genes : bool
            Flag to retrieve marker genes (default is True).
        top_n_marker : int
            Number of top marker genes to select (default is 15).
        remove_marker : bool
            Flag to remove marker genes (default is True).
        fit_GLM_libratio : bool
            Flag to fit GLM using library ratio (default is False, to use counts ratio).
        select_normal_chr_num : int
            Number of normal chromosomes to select (default is 4).
        dispersion_celltype : str or None
            Cell type for dispersion calculation.
        gene_exp_group : int
            Expression group for gene analysis (default is 1).
        gene_exp_ref_log : bool
            Flag to use log transformation for reference expression (should active when exp_group > 1).
        guide_cnv_ratio : float or None
            Ratio for guiding CNV.
        guide_chr_anno_key : str
            Annotation key for chromosome guidance (default is "chr_arm").
        guide_qt_lst : list of float
            List of quantiles for guidance (default is [1e-04, 0.96, 0.99] for "chr_arm",
            We recommend try [0.00001, 0.96, 0.999] for "chr").
        WMA_window_size : int
            Window size for Weighted Moving Average (WMA) smoothing (default is 40).
        WMA_smooth_key : str
            Key for WMA smoothing (default is "chr_arm").
        xclone_plot : bool
            Flag to enable XClone plotting (default is True).
        plot_cell_anno_key : str or None
            Annotation key for plotting cell annotations.
        rdr_plot_vmin : float
            Minimum value for RDR plot color scale.
        rdr_plot_vmax : float
            Maximum value for RDR plot color scale.
        set_figtitle : bool
            Flag to set figure titles in plots (default is True).
        """
    def __init__(self):
        """
        Initialize the class RDR_General_config() class with default configuration parameters.

        Parameters
        ----------

            None
        """
        self.smart_transform = False
        self.filter_ref_ave = 0.5
        self.min_gene_keep_num = 3000
        self.multi_refcelltype = False
        self.marker_group_anno_key = None
        self.get_marker_genes = True
        self.top_n_marker = 15
        self.remove_marker = True
        self.fit_GLM_libratio = False # default use counts ratio
        self.select_normal_chr_num = 4
        self.dispersion_celltype = None
        self.gene_exp_group = 1
        # active when exp_group larger than 2
        self.gene_exp_ref_log = True
        self.guide_cnv_ratio = None
        self.guide_chr_anno_key = "chr_arm"
        self.guide_qt_lst = [1e-04, 0.96, 0.99]
        ## smoothing
        self.WMA_window_size = 40
        self.WMA_smooth_key = "chr_arm"
        # Notes: WMA_smooth_key may update to predefined segment for simple clones
        
        ## RDR plotting
        self.xclone_plot = True
        self.plot_cell_anno_key =  None
        self.rdr_plot_vmin = -0.7
        self.rdr_plot_vmax = 0.7
        self.set_figtitle = True

class RDR_General_config():
    """
    General configuration settings for new version of RDR (RDR_gaussian)
    Adapted from RDR_prev_General_config for Gaussian-based RDR analysis.

    This class manages the configuration settings for RDR analysis, particularly
    for 10X scRNA-seq data, including transformation options, filtering criteria,
    marker gene settings, GLM fitting settings, smoothing parameters, and plotting preferences.

    Attributes
    ----------

        smart_transform : bool
            Flag to apply smart transformation (default is False).
        filter_ref_ave : float
            Threshold for filtering based on reference average (default is 0.5).
        min_gene_keep_num : int
            Minimum number of genes to keep (default is 3000).
        multi_refcelltype : bool
            Flag to use multiple reference cell types.
        marker_group_anno_key : str or None
            Annotation key for marker groups.
        get_marker_genes : bool
            Flag to retrieve marker genes (default is True).
        top_n_marker : int
            Number of top marker genes to select (default is 15).
        remove_marker : bool
            Flag to remove marker genes (default is True).
        fit_GLM_libratio : bool
            Flag to fit GLM using library ratio (default is False, to use counts ratio).
        select_normal_chr_num : int
            Number of normal chromosomes to select (default is 4).
        dispersion_celltype : str or None
            Cell type for dispersion calculation.
        gene_exp_group : int
            Expression group for gene analysis (default is 1).
        gene_exp_ref_log : bool
            Flag to use log transformation for reference expression (should active when exp_group > 1).
        guide_cnv_ratio : float or None
            Ratio for guiding CNV.
        guide_chr_anno_key : str
            Annotation key for chromosome guidance (default is "chr_arm").
        guide_qt_lst : list of float
            List of quantiles for guidance (default is [1e-04, 0.96, 0.99] for "chr_arm",
            We recommend try [0.00001, 0.96, 0.999] for "chr").
        WMA_window_size : int
            Window size for Weighted Moving Average (WMA) smoothing (default is 40).
        WMA_smooth_key : str
            Key for WMA smoothing (default is "chr_arm").
        xclone_plot : bool
            Flag to enable XClone plotting (default is True).
        plot_cell_anno_key : str or None
            Annotation key for plotting cell annotations.
        rdr_plot_vmin : float
            Minimum value for RDR plot color scale.
        rdr_plot_vmax : float
            Maximum value for RDR plot color scale.
        set_figtitle : bool
            Flag to set figure titles in plots (default is True).
        """
    def __init__(self):
        """
        Initialize the class RDR_General_config() class with default configuration parameters.

        Parameters
        ----------

            None
        """
        self.smart_transform = False
        self.filter_ref_ave = 0.5
        self.min_gene_keep_num = 3000
        self.multi_refcelltype = False
        self.marker_group_anno_key = None
        self.get_marker_genes = True
        self.top_n_marker = 15
        self.remove_marker = True
        self.fit_GLM_libratio = False # default use counts ratio
        self.select_normal_chr_num = 4
        self.dispersion_celltype = None
        self.gene_exp_group = 1
        # active when exp_group larger than 2
        self.gene_exp_ref_log = True
        self.guide_cnv_ratio = None
        self.guide_chr_anno_key = "chr_arm"
        self.guide_qt_lst = [1e-04, 0.96, 0.99]
        ## smoothing
        self.WMA_window_size = 40
        self.WMA_smooth_key = "chr_arm"
        # Notes: WMA_smooth_key may update to predefined segment for simple clones
        
        ## RDR plotting
        self.xclone_plot = True
        self.plot_cell_anno_key =  None
        self.rdr_plot_vmin = -0.7
        self.rdr_plot_vmax = 0.7
        self.set_figtitle = True

        # adaptive baseline
        self.ab_k_neighbors = 5
        self.ab_pseudo_count = 1e-6

        # denoise
        self.denoise_sd_amplifier = 1.5

        # GMM
        self.c_k = np.array([0.5, 1, 1.5])

        # low rank
        self.low_rank = True
        self.low_rank_n_components = 10


class BAF_General_config():
    """
    General configuration settings for BAF (B Allele Frequency) analysis.

    This class manages the configuration settings for BAF analysis, particularly
    for 10X scRNA-seq data, including bias mode settings, related RDR settings,
    KNN connectivity options, theoretical CNV states, phasing parameters, smoothing,
    postprocessing options, and plotting preferences.

    Attributes
    ----------

        baf_bias_mode : int
            Mode for BAF bias setting.
        CNV_N_components : int
            Number of components for CNV analysis, dependent on baf_bias_mode.
        BAF_add : bool or None
            Additional BAF configuration based on baf_bias_mode.
        update_info_from_rdr : bool
            Flag to update information from RDR.
        RDR_file : str or None
            Path to the RDR file.
        remove_marker_genes : bool
            Flag to remove marker genes.
        KNN_connect_use_key : str
            Key for KNN connectivity usage.
        get_BAF_KNN_connectivities : bool
            Flag to get BAF KNN connectivities.
        KNN_Xlayer : str
            X layer for KNN smoothing.
        guide_theo_CNV_states : str or None
            Theoretical CNV states for guidance.
        theo_neutral_BAF : str or None
            Theoretical neutral BAF value.
        ref_BAF_clip : bool
            Flag to clip reference BAF.
        concentration : int
            Concentration value for analysis.
        extreme_count_cap : bool
            Flag to cap extreme counts.
        gene_specific_concentration : bool
            Flag for gene-specific concentration.
        concentration_lower : int
            Lower bound for concentration.
        concentration_upper : int
            Upper bound for concentration.
        feature_mode : str
            Mode for feature selection.
        phasing_region_key : str
            Key for phasing regions.
        phasing_len : int
            Length of phasing regions.
        bin_nproc : int
            Number of processes for binning.
        WMA_window_size : int
            Window size for Weighted Moving Average (WMA) smoothing.
        WMA_smooth_key : str
            Key for WMA smoothing.
        BAF_denoise : bool
            Flag to enable BAF denoising.
        BAF_denoise_GMM_detection : bool
            Flag to enable GMM detection in BAF denoising.
        BAF_denoise_GMM_comp : int
            Number of components for GMM in BAF denoising.
        BAF_denoise_cellprop_cutoff : float
            Cutoff for cell proportion in BAF denoising.
        xclone_plot : bool
            Flag to enable XClone plotting.
        plot_cell_anno_key : str or None
            Annotation key for plotting cell annotations.
        set_figtitle : bool
            Flag to set figure titles in plots.
    """
    def __init__(self, baf_bias_mode):
        """
        Initialize the BAF_General_config class with configuration parameters based on the BAF bias mode.

        This method sets default values for BAF analysis settings, optimized for 10X scRNA-seq data.

        Parameters
        ----------

            baf_bias_mode : int
                Mode for BAF bias setting, determines specific configurations.
                Default 1 (5 allele bias states), otherwise 0 (3 allele bias states).
        """
        self.baf_bias_mode = baf_bias_mode
        if self.baf_bias_mode == 0:
            self.CNV_N_components = 3
            self.BAF_add = False
        elif self.baf_bias_mode == 1:
            self.CNV_N_components = 5
            self.BAF_add = None
        
        ## related to RDR
        self.update_info_from_rdr = True
        self.RDR_file = None
        self.remove_marker_genes = True
        self.KNN_connect_use_key = "connectivities_expr"
        
        self.get_BAF_KNN_connectivities = False
        self.KNN_Xlayer = "fill_BAF_phased"
        
        self.guide_theo_CNV_states = None
        self.theo_neutral_BAF = None
        self.ref_BAF_clip = False
        self.concentration = 100
        self.extreme_count_cap = False
        self.gene_specific_concentration = False
        self.concentration_lower = 20
        self.concentration_upper = 100
        
        ## loacal phasing
        self.feature_mode = "GENE"
        # self.init_mode = "warm"
        self.phasing_region_key = "chr"
        self.phasing_len = 100
        self.bin_nproc = 20
        ## smoothing
        self.WMA_window_size = 101
        self.WMA_smooth_key = "chr_arm"
        self.HMM_nproc = 40
        ## postprocessing
        self.BAF_denoise = True
        self.BAF_denoise_GMM_detection = True
        self.BAF_denoise_GMM_comp = 2
        self.BAF_denoise_cellprop_cutoff = 0.05
        
        ## BAF plotting
        self.xclone_plot = True
        self.plot_cell_anno_key =  None
        self.set_figtitle = True



class Combine_General_config():
    """
    General configuration settings for combining BAF and RDR analyses.

    This class manages the configuration settings for combining BAF and RDR analyses,
    including denoising options, copy number correction settings, plotting preferences,
    and whole-genome duplication (WGD) detection parameters.

    Attributes
    ----------

        BAF_denoise : bool
            Flag to enable BAF denoising.
        RDR_denoise : bool
            Flag to enable RDR denoising.
        copyloss_correct : bool
            Flag to apply copy loss correction (default is True).
        copyloss_correct_mode : int
            Mode for copy loss correction.
        copygain_correct : bool
            Flag to apply copy gain correction (default is False).
        copygain_correct_mode : int or None
            Mode for copy gain correction, None if copygain_correct is False.
        RDR_prior : bool
            Flag to prioritize RDR in combination analysis.
        xclone_plot : bool
            Flag to enable XClone plotting.
        plot_cell_anno_key : str or None
            Annotation key for plotting cell annotations.
        merge_loss : bool
            Flag to merge loss segments in plots.
        merge_loh : bool
            Flag to merge loss of heterozygosity (LOH) segments in plots.
        set_figtitle : bool
            Flag to set figure titles in plots.
        customizedplotting : bool
            Flag to get customized plot, related to the setting in `merge_loss` and `merge_loh`.
        WGD_detection : bool
            Flag to enable whole-genome duplication (WGD) detection.
        WGD_detect_genome_level : str
            Genome level for WGD detection (e.g., "chr_arm").
        WGD_prop_value_threshold : float
            Proportion value threshold for WGD detection.
        WGD_cell_prop_threshold : int
            Cell proportion threshold for WGD detection.
    """

    def __init__(self):
        """
        Initialize the Combine_General_config class with default combination parameters.

        This method sets default values for combining BAF and RDR analyses, optimized for default settings.

        Parameters
        ----------

            None
        """
        ## combine performing
        self.BAF_denoise = False
        self.RDR_denoise = False
        self.copyloss_correct = True # default
        self.copyloss_correct_mode = 1
        self.copygain_correct= False # default
        if  self.copygain_correct== False:
            self.copygain_correct_mode = None
        self.RDR_prior = True

        ## combine plotting
        self.xclone_plot = True
        self.plot_cell_anno_key =  None
        self.merge_loss = True
        self.merge_loh = True
        self.set_figtitle = True
        self.customizedplotting = False
        
        ## function in combine module
        self.WGD_detection = True
        self.WGD_detect_genome_level = "chr_arm"
        self.WGD_prop_value_threshold = 0.9
        self.WGD_cell_prop_threshold = 50

class HMM_Configs():
    """
    Configuration settings for Hidden Markov Model (HMM) smoothing.

    This class manages the configuration settings for HMM smoothing, 
    including transition probabilities, HMM breaks, and module-specific settings for RDR and BAF.

    Attributes
    ----------

        trans_t : float
            base setting for transition probability.
        HMM_brk : str
            Break point for HMM analysis (e.g., "chr_arm").
        start_prob : numpy.ndarray
            Starting probabilities for the HMM states.
        trans_prob : numpy.ndarray
            Transition probabilities for the HMM states.
        max_iter : int
            Maximum number of iterations for HMM.
        min_iter : int
            Minimum number of iterations for HMM.
        module : str
            Module type for HMM configuration ("RDR" or "BAF").
        CNV_N_components : int
            Number of CNV components, relevant for BAF module.
    """
    def __init__(self):
        """
        Initialize the HMM_Configs class with default HMM smoothing parameters.

        This method sets default values for HMM smoothing parameters, optimized for different modules (RDR and BAF).

        Parameters
        ----------

            None
        """
        ## base setting
        self.trans_t = 1e-6
        # self.trans_prob = None
        t = self.trans_t
        
        self.HMM_brk = "chr_arm"
        
        if self.module == "RDR" or self.module == "RDR_prev":
            self.start_prob = np.array([0.1, 0.8, 0.1])
            self.trans_prob = np.array([[1-2*t, t, t],[t, 1-2*t, t],[t, t, 1-2*t]])

            ## HMM iteration
            self.max_iter = 2
            self.min_iter = 1

        if self.module == "BAF":
            if self.CNV_N_components == 3:
                self.start_prob = np.array([0.3, 0.4, 0.3])
                self.trans_prob = np.array([[1-2*t, t, t],[t, 1-2*t, t],[t, t, 1-2*t]])
            elif self.CNV_N_components == 5:
                self.start_prob = np.array([0.2, 0.15, 0.3, 0.15, 0.2])
                self.trans_prob = np.array([[1-4*t, t, t, t,t],[t, 1-4*t, t, t,t],[t, t, 1-4*t, t,t], [t, t, t, 1-4*t, t], [t, t, t, t, 1-4*t]])

        

class Smartseq_Config():
    """
    Configuration settings specific to Smart-seq data analysis.

    This class manages the configuration settings for Smart-seq data, including
    module-specific settings for RDR, BAF, and their combination, as well as general settings.

    Attributes
    ----------

        module : str
            Module type for configuration ("RDR", "BAF", or "Combine").
        smart_transform : bool
            Flag to enable smart transformation (RDR module).
        filter_ref_ave : float
            Reference average filter threshold (RDR module).
        start_prob : numpy.ndarray
            Starting probabilities for the HMM states (RDR module).
        extreme_count_cap : bool
            Flag to cap extreme counts (BAF module).
        gene_specific_concentration : bool
            Flag to enable gene-specific concentration (BAF module).
        concentration : None
            Concentration value, None if gene_specific_concentration is True (BAF module).
        WMA_window_size : int
            Window size for Weighted Moving Average (WMA) smoothing (BAF module).
        exclude_XY : bool
            Flag to exclude sex chromosomes (general setting).
        remove_guide_XY : bool
            Flag to remove guide for sex chromosomes (general setting).
    """
    def __init__(self):
        """
        Initialize the Smartseq_Config class with default Smart-seq specific parameters.

        This method sets default values for Smart-seq specific parameters, optimized for different modules (RDR, BAF, Combine).

        Parameters
        ----------

            None
        """

        if self.module == "RDR" or self.module == "RDR_prev":
            self.smart_transform = True
            self.filter_ref_ave = 1.8
            ## HMM related
            self.start_prob = np.array([0.3, 0.4, 0.3])

        if self.module == "BAF":
            self.extreme_count_cap = False
            self.gene_specific_concentration = True
            if self.gene_specific_concentration == True:
                self.concentration = None

            ## smoothing related
            self.WMA_window_size = 6
        
        if self.module == "Combine":
            pass
        
        ## general settings
        self.exclude_XY = True
        self.remove_guide_XY = True
        # self.KNN_neighbors = 5


class Spatial_Config():
    """
    Configuration settings specific to spatial transcriptomics analysis.

    This class manages the configuration settings for spatial transcriptomics data, 
    including module-specific settings for RDR, BAF, and their combination, as well as general settings.

    Attributes
    ----------

        module : str
            Module type for configuration ("RDR", "BAF", or "Combine").
        smart_transform : bool
            Flag to enable smart transformation (RDR module).
        spatial_transform : bool
            Flag to enable spatial transformation (RDR module).
        filter_ref_ave : float
            Reference average filter threshold (RDR module).
        min_gene_keep_num : int
            Minimum number of genes to keep (RDR module).
        start_prob : numpy.ndarray
            Starting probabilities for the HMM states (RDR module).
        extreme_count_cap : bool
            Flag to cap extreme counts (BAF module).
        gene_specific_concentration : bool
            Flag to enable gene-specific concentration (BAF module).
        concentration : None
            Concentration value, None if gene_specific_concentration is True (BAF module).
        WMA_window_size : int
            Window size for Weighted Moving Average (WMA) smoothing (BAF module).
        exclude_XY : bool
            Flag to exclude sex chromosomes (general setting).
        remove_guide_XY : bool
            Flag to remove guide for sex chromosomes (general setting).
    """
    def __init__(self):
        """
        Initialize the Spatial_Config class with default spatial transcriptomics specific parameters.

        This method sets default values for spatial transcriptomics specific parameters, 
        optimized for different modules (RDR, BAF, Combine).

        Parameters
        ----------

            None
        """

        if self.module == "RDR" or self.module == "RDR_prev":
            # self.spatail_imputation = True
            self.smart_transform = False
            self.spatial_transform = True
            self.filter_ref_ave = 0.5
            self.min_gene_keep_num = 1500
            ## HMM related
            # self.start_prob = np.array([0.1, 0.8, 0.1])
            self.start_prob = np.array([0.3, 0.4, 0.3])

        if self.module == "BAF":
            self.extreme_count_cap = False
            self.gene_specific_concentration = True
            if self.gene_specific_concentration == True:
                self.concentration = None

            ## smoothing related
            self.WMA_window_size = 6
        
        if self.module == "Combine":
            pass
        
        ## general settings
        self.exclude_XY = True
        self.remove_guide_XY = True
        # self.KNN_neighbors = 5


## todo: denoise part

class XCloneConfig():
    """
    Config manager for XClone.

    This class manages the configuration settings for the XClone tool,
    including dataset details, module-specific configurations, file formats,
    and output directories.
    """

    def __init__(
        self,
        dataset_name: str = "XClone_scDATA",
        set_smartseq: bool = False,
        set_spatial: bool = False,
        module: str = "RDR",
        baf_bias_mode = 1,
        plot_suffix: str = "",
        file_format_data: str = "h5ad",
        file_format_figs: str = "pdf",
        outdir: Union[str, Path] = "./XCLONE_OUT/",
        _frameon: bool = True,
        _vector_friendly: bool = False
    ):
        """\
        Initialize the XCloneConfig class with various configuration parameters.

        Parameters
        ----------

            dataset_name : str, optional
                The name of the dataset (default is "XClone_scDATA").
            set_smartseq : bool, optional
                Flag to set Smart-seq specific configurations (default is False).
            set_spatial : bool, optional
                Flag to set spatial specific configurations (default is False).
            module : str, optional
                The module to configure, must be one of "RDR", "BAF", or "Combine" 
                (default is "RDR").
            baf_bias_mode : int, optional
                BAF bias mode (default is 1).
            plot_suffix : str, optional
                Suffix to append to plot filenames (default is "").
            file_format_data : str, optional
                Format for saving data files (default is "h5ad").
            file_format_figs : str, optional
                Format for saving figure files (default is "pdf").
            outdir : Union[str, Path], optional
                Directory to save output files (default is "./XCLONE_OUT/").
            _frameon : bool, optional
                Flag to add frames and axes labels to scatter plots (default is True).
            _vector_friendly : bool, optional
                Flag for vector-friendly plotting (default is False).
        """
        self.dataset_name = dataset_name
        self.set_smartseq = set_smartseq
        self.set_spatial = set_spatial
        if module in ["RDR", "BAF", "Combine", "RDR_prev"]:
            pass
        else:
            print(module)
            raise ValueError("module is NOT supported in XClone.")
        
        # module specific
        self.module = module
        if self.module == "RDR":
            RDR_General_config.__init__(self)
            HMM_Configs.__init__(self)
        if self.module == "RDR_prev":
            RDR_prev_General_config.__init__(self)
            HMM_Configs.__init__(self)
        if self.module == "BAF":
            BAF_General_config.__init__(self, baf_bias_mode)
            HMM_Configs.__init__(self)
        if self.module == "Combine":
            Combine_General_config.__init__(self)
        
        # xclone general config
        Base_settings.__init__(self)
        XCloneGeneral_config.__init__(self)
        
        if self.set_smartseq:
            Smartseq_Config.__init__(self)
        if self.set_spatial:
            Spatial_Config.__init__(self)


        # other general config
        self.plot_suffix = plot_suffix
        self.file_format_data = file_format_data
        self.file_format_figs = file_format_figs
        self.outdir = outdir
        self._frameon = _frameon
        """bool: See set_figure_params."""

        self._vector_friendly = _vector_friendly
        """Set to true if you want to include pngs in svgs and pdfs."""

        self._start = time()
        """Time when the settings module is first imported."""

    @property
    def plot_suffix(self) -> str:
        """Global suffix that is appended to figure filenames."""
        return self._plot_suffix

    @plot_suffix.setter
    def plot_suffix(self, plot_suffix: str):
        _type_check(plot_suffix, "plot_suffix", str)
        self._plot_suffix = plot_suffix

    @property
    def file_format_data(self) -> str:
        """File format for saving AnnData objects.
        Allowed are 'txt', 'csv' (comma separated value file) for exporting and 'h5ad'
        (hdf5) for lossless saving.
        """
        return self._file_format_data

    @file_format_data.setter
    def file_format_data(self, file_format: str):
        _type_check(file_format, "file_format_data", str)
        file_format_options = {"txt", "csv", "h5ad"}
        if file_format not in file_format_options:
            raise ValueError(
                f"Cannot set file_format_data to {file_format}. "
                f"Must be one of {file_format_options}"
            )
        self._file_format_data = file_format

    @property
    def file_format_figs(self) -> str:
        """File format for saving figures.
        For example 'png', 'pdf' or 'svg'. Many other formats work as well (see
        `matplotlib.pyplot.savefig`).
        """
        return self._file_format_figs

    @file_format_figs.setter
    def file_format_figs(self, figure_format: str):
        _type_check(figure_format, "figure_format_data", str)
        self._file_format_figs = figure_format

    @property
    def outdir(self) -> Path:
        """\
        Directory where the function xclone.write writes to by default.
        """
        return self._outdir

    @outdir.setter
    def outdir(self, outdir: Union[str, Path]):
        _type_check(outdir, "outdir", (str, Path))
        self._outdir = Path(outdir)

    def display(self):
        """Display Configuration values."""
        print(self.module, "\nConfigurations:")
        for a in dir(self):
            if not a.startswith("__") and not callable(getattr(self, a)):
                print("{:30} {}".format(a, getattr(self, a)))
        print("\n")

    # --------------------------------------------------------------------------------
    # Plotting config settings
    # --------------------------------------------------------------------------------

    # Collected from the print_* functions in matplotlib.backends
    # fmt: off
    try:
        from typing import Literal
    except ImportError:
        try:
            from typing_extensions import Literal
        except ImportError:

            class LiteralMeta(type):
                def __getitem__(cls, values):
                    if not isinstance(values, tuple):
                        values = (values,)
                    return type('Literal_', (Literal,), dict(__args__=values))

            class Literal(metaclass=LiteralMeta):
                pass

    _Format = Literal[
        'png', 'jpg', 'tif', 'tiff',
        'pdf', 'ps', 'eps', 'svg', 'svgz', 'pgf',
        'raw', 'rgba',
    ]
    # fmt: on

    def set_figure_params(
        self,
        xclone: bool = True,
        dpi: int = 80,
        dpi_save: int = 150,
        frameon: bool = True,
        vector_friendly: bool = True,
        fontsize: int = 14,
        figsize: Optional[int] = None,
        color_map: Optional[str] = None,
        grid_alpha: Optional[float] = None,
        format: _Format = "pdf",
        facecolor: Optional[str] = None,
        transparent: bool = False,
        ipython_format: str = "png2x",
    ):
        """\
        Set resolution/size, styling and format of figures.

        Parameters
        ----------

            xclone
                Init default values for :obj:`matplotlib.rcParams` suited for XClone.
            dpi
                Resolution of rendered figures, this influences the size of figures in notebooks.
            dpi_save
                Resolution of saved figures. This should typically be higher to achieve publication quality.
            frameon
                Add frames and axes labels to scatter plots.
            vector_friendly
                Plot scatter plots using `png` backend even when exporting as `pdf` or `svg`.
            fontsize
                Set the fontsize for several `rcParams` entries. Ignored if `xclone=False`.
            figsize
                Set plt.rcParams['figure.figsize'].
            color_map
                Convenience method for setting the default color map. Ignored if `xclone=False`.
            format
                This sets the default format for saving figures: `file_format_figs`.
            facecolor
                Sets backgrounds via `rcParams['figure.facecolor'] = facecolor` and
                `rcParams['axes.facecolor'] = facecolor`.
            transparent
                Save figures with transparent back ground. Sets
                `rcParams['savefig.transparent']`.
            ipython_format
                Only concerns the notebook/IPython environment; see
                :func:`~IPython.display.set_matplotlib_formats` for details.
        """
        if self._is_run_from_ipython():
            import IPython

            if isinstance(ipython_format, str):
                ipython_format = [ipython_format]
            IPython.display.set_matplotlib_formats(*ipython_format)

        from matplotlib import rcParams

        self._vector_friendly = vector_friendly
        self.file_format_figs = format
        if dpi is not None:
            rcParams["figure.dpi"] = dpi
        if dpi_save is not None:
            rcParams["savefig.dpi"] = dpi_save
        if transparent is not None:
            rcParams["savefig.transparent"] = transparent
        if facecolor is not None:
            rcParams['figure.facecolor'] = facecolor
            rcParams['axes.facecolor'] = facecolor
        if xclone:
            from .plot._plot_style import set_style_xclone
            set_style_xclone(fontsize=fontsize, color_map=color_map, grid_alpha = grid_alpha)
        if figsize is not None:
            rcParams['figure.figsize'] = figsize
        self._frameon = frameon

    @staticmethod
    def _is_run_from_ipython():
        """Determines whether we're currently in IPython."""
        import builtins

        return getattr(builtins, "__IPYTHON__", False)

    def __str__(self) -> str:
        return '\n'.join(
            f'{k} = {v!r}'
            for k, v in inspect.getmembers(self)
            if not k.startswith("_") and not k == 'getdoc'
        )


settings = XCloneConfig()