Release History
===============

Version 0.3.5
-------------
- [BAF] add allele flip status for each gene in Local phasing and Global phasing
- [minor debug] baf_bias_mode configs update in BAF module
- [minor debug] add check cell quality in RDR `extra_preprocess`.[DEVELOP FOR ST](need test more)
- [Analysis] initialize reference cells detection function (from BAF) (todo)
- [BAF] update BAF theo ref strategy in wrap BAF module (if refernce cells is limited)
- [BAF] test 2 BAF theo ref strategies in wrap BAF module (if refernce cells is limited)
- [RDR] update RDR guide ratio correct strategy
- [Main update] default params setting for low coverage datasets
- [mm10 extension init] update mm10 genes and extend to mm10 application
- [Spatial extension init] add some params for spatial transcriptomics
- [Spatial filter Genes] set lower `min_gene_keep_num` in RDR module
- [RDR denoise] add function to use BAF information denoisde RDR layer before combination[test now]
- [BAF] minor update to get its own KNN `connectivities` instead of `connectivities_expr`
- [RDR] Reference cell usage strategy update (for multiple celltyps as reference)[init]
- [Combine] Special WGD warning function in the combine module.
- [Spatial] add spatial coordinates in adata.
- [Spatial] debug some preprocessins issues to adapt to spatial cases.






Version 0.3.4
-------------
- remove deprecated funcs
- update plotting tool, mainly in color settings
- update xclone tutorials in version 0.3.4
- BAF local phasing improvement
- BAF module support 3 and 5 allele bias states
- Combine module support prob correct mode for copy gain
- Combine module support prob correct mode for copy loss (default)
- update libratio capping and add cell filter limits
- use counts ratio as default
- libratio fitting update
- provide new denoise function (denoise by cluster)
- add ROC evaluation function in analysis module
- add functions for XClone clustering
- wrap independent plotting modules
- potential init to add predefined segmentation
- minor debug on BAF guide allele bias ratio
- add tutorials in xclone plotting #
- add tutorials in xclone data analysis #


Version 0.3.3
-------------
- Add tutorials for data loadding
- Wrap data loading for 3 main modules
- try supporting h5py higher version
- provide func to get results of different resolution (for benchmarking)
- updates in config: add default Smartseq_Config
- minor updates in other modules' configs and general configs
- BAF module update: feature specific concentration; extreme DP detection and rescaling(test)
- RDR module update: guide chr exclude chr X&Y
- more protecting steps in wrapped main module

Version 0.3.2
-------------
- fix bugs in loadding demo datasets

Version 0.3.1
-------------
- Wrap 3 main modules: RDR, BAF and Combine
- Add tutorials and documentations
- Provide demo data

Version 0.3.0
-------------
- Alpha version of XClone released