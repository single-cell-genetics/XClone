.. _getting started:
===============
Getting Started
===============
XClone is a software tool designed for analyzing single-cell RNA sequencing data to identify copy number variations within individual cells.  
CNVs are important in a variety of biological processes, especially in cancer development.

XClone uses probabilistic modelling to estimate the copy number variation status of each gene within each cell in a sample. 
By analyzing the expression levels of genes and B allele frequency of gene_bins, 
XClone can determine which cells are likely to have copy number variations and which genes are affected.

In this tutorial, we'll walk through the steps of using XClone to analyze single-cell RNA sequencing data for CNV analysis. 
We'll cover everything from preparing the input data to interpreting the output of the analysis. 
Whether you're new to single-cell RNA sequencing or an experienced user, this tutorial will provide a comprehensive introduction to using XClone for CNV analysis.

Public Datasets
===============

`TNBC1 data`_: a triple-negative breast cancer (TNBC) sample that was assayed by droplet-based scRNA-seq (10x Genomics), from `Gao et al. (2021)`_.

`BCH869 data`_: a glioma sample BCH869 with histone H3 lysine27-to-methionine mutations (H3K27M-glioma), where 489 malignant cells and 3 non-tumour cells were probed by smart-seq2, from `Filbin et al. (2018)`_.

Raw Data matrix (in anndata format) for RDR module and BAF moudle can be downloaded and imported directly (see notebook tutorials). Data download at `demo_data`_.

Examples of XClone and steps for reproducible results are provided in **Jupyter Notebook** under `examples`_ folder. 

For start, please refer to tutorials analyzing `TNBC1`_ and `BCH869`_ datasets, or `TNBC1_tutorial`_ and `BCH869_tutorial`_ records step by step.

.. _examples: https://github.com/Rongtingting/xclone-data/tree/main/examples
.. _TNBC1: ./TNBC1_XClone_tutorials.html
.. _BCH869: ./BCH869_XClone_tutorials.html
.. _TNBC1_tutorial: ./TNBC1_XClone_demo_v2.html
.. _BCH869_tutorial: ./BCH869_XClone_demo_v2.html
.. _demo_data: https://connecthkuhk-my.sharepoint.com/:f:/g/personal/rthuang_connect_hku_hk/EnKri0rS-ZpHl0VGVHUp4k0B_3iZ_gpD-obVuDwEMQUieQ?e=k0eR4T
.. _TNBC1 data: https://connecthkuhk-my.sharepoint.com/:f:/g/personal/rthuang_connect_hku_hk/Etlhi3gMu_VJuhmtrQiQRO4BRu4VVxIE_yL3Mt6iQ10kkA?e=zV0qbe
.. _BCH869 data: https://connecthkuhk-my.sharepoint.com/:f:/g/personal/rthuang_connect_hku_hk/EhnxMmkOFsNOto8XN0OYNr0BNVAvZOem3SKFcpjBKMTJFw?e=0e73Rg
.. _Gao et al. (2021): https://www.nature.com/articles/s41587-020-00795-2
.. _Filbin et al. (2018): DOI: 10.1126/science.aao4750

XClone on new datasets
======================
- import package

::

    import xclone

XClone provides integrated functions (RDR BAF and Combination) for CNV analysis by default 
whilst specific configurations might need to be adjusted accordingly. For each XClone module, we provide
independent Class for configuration. You can specify which module to use by set `module = "RDR"`, `module = "BAF"`
or `module = "Combine"` in `xclone.XCloneConfig(dataset_name = dataset_name, module = "RDR")`.
We provide two sets of default base configurations for 10X Genomics scRNA-seq and SMART-seq, default is 10X scRNA-seq.
If you want to get default base configurations for SMART-seq dataset, you can set `set_smartseq = True` in `xclone.XCloneConfig()`.

If you want to change params setting, Please Sub-class and override base configuration file (here lists a few frequently params used), 
please refer config.py for detailed arguments. After overriding, you can print `xconfig.display()` to show the updated configurations 
before you run the module, e.g., `RDR_Xdata = xclone.model.run_RDR(RDR_adata, config_file = xconfig)`.

- RDR module

::

    xconfig = xclone.XCloneConfig(dataset_name = dataset_name, module = "RDR", set_smartseq = True)
    xconfig.set_figure_params(xclone= True, fontsize = 18)
    xconfig.outdir = out_dir
    xconfig.cell_anno_key = "cell_type"
    xconfig.ref_celltype = "N"
    xconfig.top_n_marker = 25
    xconfig.marker_group_anno_key = "cell_type"
    xconfig.xclone_plot= True
    xconfig.plot_cell_anno_key = "Clone_ID"
    xconfig.trans_t = 1e-6
    xconfig.start_prob = np.array([0.3, 0.4, 0.3])


    xconfig.display()

    RDR_Xdata = xclone.model.run_RDR(RDR_adata,
                config_file = xconfig)


- BAF module

::

    xconfig = xclone.XCloneConfig(dataset_name = dataset_name, module = "BAF", set_smartseq = True)
    xconfig.set_figure_params(xclone= True, fontsize = 18)
    xconfig.outdir = out_dir
    xconfig.cell_anno_key = "cell_type"
    xconfig.ref_celltype = "N"
    xconfig.concentration = 35.5
    xconfig.concentration_lower = 20
    xconfig.concentration_upper = 100
    xconfig.theo_neutral_BAF = 0.5

    xconfig.xclone_plot= True
    xconfig.plot_cell_anno_key = "Clone_ID"
    xconfig.phasing_region_key = "chr"
    xconfig.phasing_len = 100

    xconfig.WMA_window_size = 6

    xconfig.trans_t = 1e-6
    xconfig.start_prob = np.array([0.2, 0.15,  0.3, 0.15, 0.2])

    t = xconfig.trans_t
    xconfig.trans_prob = np.array([[1-4*t, t, t, t,t],[t, 1-4*t, t, t,t],[t, t, 1-4*t, t,t], [t, t, t, 1-4*t, t], [t, t, t, t, 1-4*t]])
    xconfig.CNV_N_components = 5

    ## update
    xconfig.BAF_denoise = True
    xconfig.display()

    BAF_merge_Xdata = xclone.model.run_BAF(BAF_adata,
                config_file = xconfig)


- Combine module

::

    xconfig = xclone.XCloneConfig(dataset_name = dataset_name, module = "Combine")
    xconfig.set_figure_params(xclone= True, fontsize = 18)
    xconfig.outdir = out_dir

    xconfig.cell_anno_key = "cell_type"
    xconfig.ref_celltype = "N"


    xconfig.copygain_correct= False

    xconfig.xclone_plot= True
    xconfig.plot_cell_anno_key = "Clone_ID"
    xconfig.merge_loss = False
    xconfig.merge_loh = True

    xconfig.BAF_denoise = True
    xconfig.display()

    combine_Xdata = xclone.model.run_combine(RDR_Xdata,
                    BAF_merge_Xdata,
                    verbose = True,
                    run_verbose = True,
                    config_file = xconfig)



XClone on GX109-T1c
===================

