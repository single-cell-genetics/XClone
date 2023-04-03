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
    Xdata loading params.
    """
    def __init__(
        self,
        dataset_name = "XClone_scDATA",
        module = "RDR",
        cell_anno_file = None,
        rdr_data_dir = None,
        baf_data_dir = None):
        
        self.module = module
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
            # self.regions_anno_file = self.rdr_data_dir + "features.tsv"
            self.cell_anno_file = cell_anno_file
            self.cell_anno_key = "cell_type"
            self.genome_mode = "hg38_genes"

        if self.module == "BAF":
            self.AD_file = self.baf_data_dir + "xcltk.AD.mtx"
            self.DP_file = self.baf_data_dir + "xcltk.DP.mtx"
            self.mtx_barcodes_file = self.baf_data_dir + "xcltk.samples.tsv"
            self.regions_anno_file = None
            # self.regions_anno_file = self.baf_data_dir + "xcltk.region.tsv"
            self.cell_anno_file = cell_anno_file
            self.cell_anno_key = "cell_type"
            self.genome_mode = "hg38_genes"

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
    def __init__(self):
        self.cell_anno_key = "cell_type"
        self.ref_celltype = "N"
        self.exclude_XY = False
        self.remove_guide_XY = False
        # KNN smoothing
        self.KNN_neighbors = 10

class RDR_General_config():
    def __init__(self):
        """
        RDR params init settings. 
        default for 10X scRNA-seq data.
        """
        self.smart_transform = False
        self.filter_ref_ave = 0.5
        self.marker_group_anno_key = None
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
        # [0.00001, 0.96, 0.99999]
        # self.guide_chr_anno_key = "chr"
        # self.guide_qt_lst = [0.00001, 0.96, 0.999]
        ## smoothing
        self.WMA_window_size = 40
        self.WMA_smooth_key = "chr_arm"
        # WMA_smooth_key may update to predefined segment for simple clones
        ## RDR plotting
        self.xclone_plot = True
        self.plot_cell_anno_key =  None
        self.rdr_plot_vmin = -0.7
        self.rdr_plot_vmax = 0.7
        self.set_figtitle = True


class BAF_General_config():
    def __init__(self, baf_bias_mode):
        """
        BAF params init settings.
        default for 10X scRNA-seq data.
        """
        self.baf_bias_mode = baf_bias_mode
        if self.baf_bias_mode == 0:
            self.CNV_N_components = 3
        elif self.baf_bias_mode == 1:
            self.CNV_N_components = 5
            self.BAF_add = None

        self.RDR_file = None
        self.guide_theo_CNV_states = None
        self.theo_neutral_BAF = None
        self.ref_BAF_clip = False
        self.concentration = 100
        self.extreme_count_cap = False
        self.gene_specific_concentration = False
        self.concentration_lower = 20
        self.concentration_upper = 100
        ## related to RDR
        self.remove_marker_genes = True
        ## loacal phasing
        self.feature_mode = "GENE"
        # self.init_mode = "warm"
        self.phasing_region_key = "chr"
        self.phasing_len = 100
        self.bin_nproc = 20
        ## smoothing
        self.WMA_window_size = 101
        self.WMA_smooth_key = "chr_arm"
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
    def __init__(self):
        """
        Combination parmas init settings.
        default settings.
        """
        ## combine performing
        self.BAF_denoise = False
        self.copyloss_correct = True # default
        self.copyloss_correct_mode = 1
        self.copygain_correct= False # default
        if  self.copygain_correct== False:
            self.copygain_correct_mode = None

        ## combine plotting
        self.merge_loss = True
        self.merge_loh = True
        self.set_figtitle = True
        self.xclone_plot = True
        self.plot_cell_anno_key =  None
        

class HMM_Configs():
    def __init__(self):
        """
        HMM smoothing params init settings.
        default settings.
        """
        ## base setting
        self.trans_t = 1e-6
        # self.trans_prob = None
        t = self.trans_t
        
        self.HMM_brk = "chr_arm"
        
        if self.module == "RDR":
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
    def __init__(self):
        """
        Smartseq specific config settings.
        default settings.
        """
        if self.module == "RDR":
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

## todo: denoise part

class XCloneConfig():
    """\
    Config manager for xclone.
    """

    def __init__(
        self,
        dataset_name: str = "XClone_scDATA",
        set_smartseq: bool = False,
        module: str = "RDR",
        baf_bias_mode = 1,
        plot_suffix: str = "",
        file_format_data: str = "h5ad",
        file_format_figs: str = "pdf",
        outdir: Union[str, Path] = "./XCLONE_OUT/",
        _frameon: bool = True,
        _vector_friendly: bool = False
    ):
        """
        """
        self.dataset_name = dataset_name
        self.set_smartseq = set_smartseq
        if module in ["RDR", "BAF", "Combine"]:
            pass
        else:
            print(module)
            raise ValueError("module is NOT supported in XClone.")
        
        # module specific
        self.module = module
        if self.module == "RDR":
            RDR_General_config.__init__(self)
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
            Resolution of rendered figures – this influences the size of figures in notebooks.
        dpi_save
            Resolution of saved figures. This should typically be higher to achieve
            publication quality.
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