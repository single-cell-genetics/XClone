"""XClone modelling.
"""

from .phasing import Local_Phasing, Global_Phasing
from .analysis_utils import filter_data, filter_2data, sub_chr


from .base_utils import normalize
from .data_base_utils import mtx_describe
from .data_base_utils import mtx_bound

from .data_base_utils import Xdata_mapping

## RDR module
from ._RDR_libratio import total_count_ratio
from ._RDR_libratio import optmize_libratio2, get_libsize2, libratio_elbo_func, libratio_elbo_fun

from ._RDR_libratio import select_normal_CHR, extend_chr_index, select_normal_CHR1
from ._RDR_libratio import extend_chr_index_to_singlecell
from ._RDR_libratio import cell_filter_by_celltype
from ._RDR_libratio import get_libsize

from ._RDR_libratio import fit_lib_ratio, fit_lib_ratio_accelerate
from ._RDR_libratio import remove_cells
from ._RDR_libratio import check_libratio
from ._RDR_libratio import libsize_clip, libsize_select


from ._RDR_dispersion import fit_Dispersions 
from ._RDR_dispersion import view_celltype, select_celltype
from ._RDR_dispersion import check_dispersion
from ._RDR_dispersion import dispersion_select, dispersion_clip

from ._RDR_dispersion import Estimate_Dispersions
from ._RDR_dispersion import map_var_info
from ._RDR_dispersion import remove_genes

from ._RDR_CNVratio  import guide_CNV_chrs, guide_CNV_states
from ._RDR_CNVratio  import gene_exp_group
from ._RDR_CNVratio import gene_specific_states
from ._RDR_CNVratio  import fit_CNV_ratio

from ._RDR_CNVratio  import generate_emm_GT
from ._RDR_CNVratio  import hard_assign_prob

from ._RDR_CNVratio  import guide_RDR_CNV_ratio

## iteration part
from ._RDR_CNVratio import combine_BAF_emmprob
from ._RDR_CNVratio import CNV_optimazation
from .base_utils import cal_log_lik

##debug
from ._RDR_cnv import fit_CNV_ratio_update

## RDR remove ref bio
from ._RDR_process import NMF_confounder
## RDR analysis
from ._RDR_process import extra_preprocess, extra_preprocess2
from ._RDR_smoothing import RDR_smoothing_base

## HMM

from .HMM_base import XHMM_smoothing
from .HMM_base import XC_HMM_base

from .HMM_NB import calculate_Xemm_prob
from .HMM_NB import calculate_Xemm_prob2
from .HMM_NB import calculate_Xemm_probTry


## BAF module
from ._BAF_base import BAF_Local_phasing
from ._BAF_base import BAF_Global_phasing, BAF_Global_phasing_rev
from ._BAF_base import BAF_fillna
from ._BAF_base import BAF_remove_ref

from .HMM_BB import calculate_Xemm_prob_bb
from .HMM_BB import generate_bb_logprob
from .HMM_BB import get_BAF_ref, get_BAF_ref_limited, gene_specific_BAF, specific_BAF

from ._BAF import extrme_count_capping
from ._BAF import concentration_mapping
## optimize BAF theoretical values
from ._BAF import BAF_theoretical_value_optimization

## smoothing
from ._BAF import BAF_smoothing
from .smoothing import WMA_smooth, KNN_smooth
from ._BAF_process import get_KNN_connectivities_from_expr
from ._BAF_process import extra_preprocess_BAF

from ._Gaussian_mixture import get_CNV_states, guide_BAF_theo_states

## combination
from .XClone_combine import gene_to_bin_mapping
# from .XClone_combine import BAF_LOH_merge, BAF_LOH_corrected
# from .XClone_combine import CNV_states_combination, copyloss_adjust
# from .XClone_combine import visualize_states_process

from .XClone_combine import bin_to_gene_mapping
from .XClone_combine import CNV_prob_combination
from .XClone_combine import CNV_prob_merge_for_plot

from .XClone_combine import WGD_warning

## post-step[analysis for specific dataset]
from ._denoise import denoise_gene_scale, merge_denoise_prob
from ._denoise import denoise_by_cluster
from ._denoise import denoise_rdr_by_baf

## wrap
from .xclone_rdr_wrap import preview_RDR
from .xclone_rdr_wrap import run_RDR_prev, run_RDR_plot, run_RDR_complex_plot
from .xclone_rdr_wrap import plot_processed_RDR

from .xclone_baf_wrap import preview_BAF
from .xclone_baf_wrap import run_BAF, run_BAF_plot, run_BAF_complex_plot
from .xclone_baf_wrap import plot_processed_BAF

from .xclone_combine_wrap import run_combine, run_combine_plot

## debug-develop-improvement
from ._HMM import Model_NB, HMM_Frame


## RDR_gaussian

### HMM
from ._RDR_gaussian_HMM import calculate_z_post_parallel

### preprocessing and smoothing
from ._RDR_gaussian_base import preprocess_adaptive_baseline, WMA_smoothing, denoise

### GMM
from ._RDR_gaussian_base import compute_gaussian_probabilities

### probability smoothing
from ._RDR_gaussian_base import smooth_and_normalize, smooth_anndata_layer, low_rank_approximation, knn_smooth_prob_layer

### plot
from ._RDR_gaussian_base import run_RDR_gaussian_plot

### RDR Gaussian
from .xclone_rdr_gaussian import run_RDR