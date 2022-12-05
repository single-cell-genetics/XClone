"""Base pipeline for XClone data prerpocessing and loadding"""

# Author: Rongting Huang
# Date: 2022-12-03
# update: 2022-12-03

import xclone

def load_Xdata(module = "RDR", 
               rdr_data_dir = None,
               baf_data_dir = None,
               config_file = None):
    """
    """
    ## settings
    from .._config import PreprocessingConfig
    
    if config_file == None:
        print (
            f'Model configuration file not specified.\n'
            f'Default settings in preprocessing will be used.'
        )
        config = PreprocessingConfig(module = module, 
        rdr_data_dir = rdr_data_dir,
        baf_data_dir = baf_data_dir)

    else:
        config = config_file
    ## general settings
    dataset_name = config.dataset_name
    
    ## moudle specific settings
    if module == "pre_check":
        BAF_barcodes_file = config.BAF_barcodes_file
        RDR_barcodes_file = config.RDR_barcodes_file
        success_flag = xclone.pp.check_RDR_BAF_anno(BAF_barcodes_file, RDR_barcodes_file)
        print("""[XClone data hint] Pls also check data chr order, if not in genome order, 
                  pls use FUNC 'xclone.pp.resort_mtx_bychr() to resort the mtx.'""")
        return success_flag
    
    if module == "RDR":
        print("[XClone RDR data loading]************************")
        RDR_file = config.RDR_file
        mtx_barcodes_file = config.mtx_barcodes_file
        regions_anno_file = config.regions_anno_file
        anno_file = config.cell_anno_file
        cell_anno_key = config.cell_anno_key

        genome_mode = config.genome_mode

        RDR_adata = xclone.pp.xclonedata(RDR_file, 'RDR', mtx_barcodes_file, 
                                         regions_anno_file,
                                         genome_mode = genome_mode,
                                         data_notes = dataset_name)
        RDR_adata = xclone.pp.extra_anno(RDR_adata, anno_file, barcodes_key = "cell", 
                                         cell_anno_key = cell_anno_key, sep =",") 
        return RDR_adata

    if module == "BAF":
        print("[XClone BAF data loading]************************")
        AD_file = config.AD_file
        DP_file = config.DP_file
        mtx_barcodes_file = config.mtx_barcodes_file
        regions_anno_file = config.regions_anno_file
        anno_file = config.cell_anno_file
        cell_anno_key = config.cell_anno_key
        genome_mode = config.genome_mode

        BAF_adata = xclone.pp.xclonedata([AD_file, DP_file], 'BAF', mtx_barcodes_file, 
                                         regions_anno_file,
                                         genome_mode = genome_mode,
                                         data_notes = dataset_name)
        BAF_adata = xclone.pp.extra_anno(BAF_adata, anno_file, barcodes_key = "cell", 
                                         cell_anno_key = cell_anno_key, sep =",")
        
        return BAF_adata

    if module == "Combine":
        print("[XClone Combine data loading]************************")
        import anndata as an
        RDR_Xdata = an.read_h5ad(config.RDR_adata_file)
        BAF_merge_Xdata = an.read_h5ad(config.BAF_adata_file)
        ## check
        xclone.pp.check_data_combine(RDR_Xdata, BAF_merge_Xdata)
        return RDR_Xdata, BAF_merge_Xdata