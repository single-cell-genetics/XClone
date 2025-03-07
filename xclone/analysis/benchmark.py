"""Base functions for Benchmarking task: performance evaluation.
"""

## Step1: Extract results from all CNA Detection methods
import os
import pandas as pd
from .utils import reg2gene

def extract_numbat(sid, dat_dir, cnv_type, method_sub="numbat", mtx_type="expr", cnv_scale="gene", gene_anno=None):
    obj_fn = f"{dat_dir}/joint_post_2.tsv"
    mtx = pd.read_csv(obj_fn, sep='\t')

    mtx['p_amp'] += mtx['p_bamp']
    mtx['p_del'] += mtx['p_bdel']
    mtx['chrom'] = mtx['CHROM'].str.replace('chr', '')
    mtx['reg_id'] = mtx.apply(lambda x: f"{x['chrom']}:{x['seg_start']}-{x['seg_end']}", axis=1)
    
    if cnv_type == "copy_gain":
        mtx = mtx.pivot(index='cell', columns='reg_id', values='p_amp')
    elif cnv_type == "copy_loss":
        mtx = mtx.pivot(index='cell', columns='reg_id', values='p_del')
    elif cnv_type == "loh":
        mtx = mtx.pivot(index='cell', columns='reg_id', values='p_loh')
    else:
        raise ValueError(f"Error: unknown cnv type '{cnv_type}'.")

    mtx_matrix = mtx.to_numpy()
    mtx.index = mtx.index.tolist()  # Set rownames to cells

    dat = reg2gene(mtx_matrix, gene_anno)  # cell x gene matrix
    return dat

def extract_numbat(sid, dat_dir, cnv_type, method_sub="numbat", mtx_type="expr", cnv_scale="gene", gene_anno=None):
    obj_fn = f"{dat_dir}/joint_post_2.tsv"
    mtx = pd.read_csv(obj_fn, sep='\t')

    mtx['p_amp'] += mtx['p_bamp']
    mtx['p_del'] += mtx['p_bdel']
    mtx['chrom'] = mtx['CHROM']
    mtx['reg_id'] = mtx.apply(lambda x: f"{x['chrom']}:{x['seg_start']}-{x['seg_end']}", axis=1)

    # mtx['reg_id'].unique()
    # array(['3:3057678-197960200', '9:14475-137619005',
    #    '15:22786225-101654085'], dtype=object)

    
    if cnv_type == "copy_gain":
        mtx = mtx.pivot(index='cell', columns='reg_id', values='p_amp')
    elif cnv_type == "copy_loss":
        mtx = mtx.pivot(index='cell', columns='reg_id', values='p_del')
    elif cnv_type == "loh":
        mtx = mtx.pivot(index='cell', columns='reg_id', values='p_loh')
    else:
        raise ValueError(f"Error: unknown cnv type '{cnv_type}'.")

    mtx_matrix = mtx.to_numpy()
    mtx.index = mtx.index.tolist()  # Set rownames to cells

    dat = reg2gene(mtx_matrix, gene_anno)  # cell x gene matrix
    return dat


def extract_copykat(sid, dat_dir, cnv_type, method_sub="copykat", mtx_type="expr", cnv_scale="gene", gene_anno=None):
    """# Example usage
# dat = extract_copykat('sample_id', 'data_directory', 'cnv_type')

    Args:
        sid (_type_): _description_
        dat_dir (_type_): _description_
        cnv_type (_type_): _description_
        method_sub (str, optional): _description_. Defaults to "copykat".
        mtx_type (str, optional): _description_. Defaults to "expr".
        cnv_scale (str, optional): _description_. Defaults to "gene".
        gene_anno (_type_, optional): _description_. Defaults to None.

    Returns:
        _type_: _description_
    """
    obj_fn = f"{dat_dir}/{sid}_copykat_CNA_raw_results_gene_by_cell.txt"
    mtx = pd.read_csv(obj_fn, sep='\t')
    
    # Set the row names to 'hgnc_symbol'
    mtx.set_index('hgnc_symbol', inplace=True)
    
    # Remove the first 7 columns and transpose
    mtx = mtx.iloc[:, 7:].transpose()
    
    dat = {'mtx': mtx, 'overlap': None}
    return dat

# extract_infercnv <- function(
#   sid, dat_dir, cnv_type, method_sub = "infercnv", mtx_type = "expr",
#   cnv_scale = "gene", gene_anno = NULL) 
# {
#   obj_fn <- sprintf(
#     "%s/BayesNetOutput.HMMi6.hmm_mode-samples/MCMC_inferCNV_obj.rds", dat_dir)
#   obj <- readRDS(obj_fn)
#   mtx <- obj@expr.data
#   mtx <- t(mtx)    # cell x gene matrix
#   dat <- list(mtx = mtx, overlap = NULL)
#   return(dat)
# }


def read_rds(file_path):
    # Placeholder function to read RDS files
    # You can use rpy2 or other methods to read RDS files in Python
    import rpy2.robjects as ro
    from rpy2.robjects import pandas2ri
    pandas2ri.activate()
    
    readRDS = ro.r['readRDS']
    rds_data = readRDS(file_path)
    return pandas2ri.ri2py(rds_data)

def extract_infercnv(sid, dat_dir, cnv_type, method_sub="infercnv", mtx_type="expr", cnv_scale="gene", gene_anno=None):
    obj_fn = os.path.join(dat_dir, "BayesNetOutput.HMMi6.hmm_mode-samples", "MCMC_inferCNV_obj.rds")
    
    # Assuming you have a function to read RDS files in Python
    obj = read_rds(obj_fn)
    
    mtx = obj['expr.data']
    mtx = mtx.T  # Transpose to get cell x gene matrix
    
    dat = {'mtx': mtx, 'overlap': None}
    return dat