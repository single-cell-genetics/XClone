"""Base functions for Benchmarking task: performance evaluation.
"""

import pandas as pd
import numpy as np
import re


# parse_region1 <- function(region) {
#   chrom <- stringr::str_extract(region, "(.*?)(?=:|$)")

#   start <- stringr::str_remove_all(
#     stringr::str_extract(region, "(?<=:)(.*?)(?=-|$)"), ",")
#   if (is.na(start) || is.na(as.numeric(start)))
#     start <- -Inf

#   end <- stringr::str_remove_all(
#     stringr::str_extract(region, "(?<=-)(.*)"), ",")
#   if (is.na(end) || is.na(as.numeric(end)))
#     end <- Inf

#   return(c(chrom = chrom, start = start, end = end))
# }

import pandas as pd

def parse_region1(region):
    chrom, *rest = region.split(':')
    start, end = (x.strip() for x in rest[0].split('-') if '-' in rest[0]) or (None, None)
    start = float('-inf') if start is None else float(start)
    end = float('inf') if end is None else float(end)
    return {'chrom': chrom, 'start': start, 'end': end}

def parse_regions(regions):
    parsed = [parse_region1(region) for region in regions]
    res = pd.DataFrame(parsed)
    res['start'] = res['start'].astype(float)
    res['end'] = res['end'].astype(float)
    return res


def parse_region1(region):
    """Parse Region String
    #' 
    #' Parse a region string '[chr]xxx[:xxx-xxx]' to a vector of 
    #' (chrom, start, end).
    #' 
    #' @param region A string representing a genomic region in the format of 
    #'   '[chr]xxx[:xxx-xxx]'.
    #' @return A string vector of 3 extracted elements: (chrom, start, end).
    #'   Note that start and end would be -Inf or Inf if not specified in the 
    #'   string.
    #' @examples
    #' parse_region1("chr1")         # c("chr1", "-Inf", "Inf")
    #' parse_region1("22")           # c("22", "-Inf", "Inf")
    #' parse_region1("chr22:-")      # c("chr22", "-Inf", "Inf")
    #' parse_region1("chrX:1230")          # c("chrX", "1230", "Inf")
    #' parse_region1("chrX:1230-")         # c("chrX", "1230", "Inf")
    #' parse_region1("X:-20221230")        # c("X", "-Inf", "20221230")
    #' parse_region1("X:1230-20221230")    # c("X", "1230", "20221230")

    Args:
        region (str): A string representing a genomic region in the format of 
        '[chr]xxx[:xxx-xxx]'.

    Returns:
        _type_: _description_
    """
    chrom_match = re.search(r"(.*?)(?=:|$)", region)
    chrom = chrom_match.group(0) if chrom_match else None

    start_match = re.search(r"(?<=:)(.*?)(?=-|$)", region)
    start = start_match.group(0).replace(",", "") if start_match else None
    start = float(start) if start and start.isdigit() else -np.inf

    end_match = re.search(r"(?<=-)(.*)", region)
    end = end_match.group(0).replace(",", "") if end_match else None
    end = float(end) if end and end.isdigit() else np.inf

    return {"chrom": chrom, "start": start, "end": end}


#' Parse Multiple Regions
#' 
#' Parse a set of region strings into a dataframe of 3 columns.
#' 
#' @param regions A string vector of genomic regions, see `parse_region1`
#'   for details of the format of regions.
#' @return A data frame containing 3 columns: 1) chrom <chr>; 2) start <num>;
#'   3) end <num>. 
#' @examples
#' regions <- c("chr1", "22:1556", "chrX:1230-20221230")
#' parse_regions(regions)
# parse_regions <- function(regions) {
#   res <- sapply(regions, parse_region1)
#   res <- as.data.frame(t(res))
#   res$start <- as.numeric(res$start)
#   res$end <- as.numeric(res$end)
#   return(res)
# }


def parse_regions(regions):
    """Parse Multiple Regions
    #' Parse a set of region strings into a dataframe of 3 columns.
    #' @param regions A string vector of genomic regions, see `parse_region1`
    #'   for details of the format of regions.
    #' @return A data frame containing 3 columns: 1) chrom <chr>; 2) start <num>;
    #'   3) end <num>. 
    #' @examples
    #' regions <- c("chr1", "22:1556", "chrX:1230-20221230")
    #' parse_regions(regions)

    Args:
        regions (_type_): _description_

    Returns:
        _type_: _description_
    """
    parsed = [parse_region1(region) for region in regions]
    res = pd.DataFrame(parsed)

    return res

def reg2gene(mtx, gene_anno, verbose=True):
    """_summary_

    Args:
        mtx (_type_): _description_
        gene_anno (_type_): _description_
        verbose (bool, optional): _description_. Defaults to True.
    """
    def parse_regions(region_ids):
        # Dummy implementation: Replace with actual parsing logic
        return pd.DataFrame([r.split(':') + r.split('-') for r in region_ids], columns=['chrom', 'start', 'end'])

    # Parse regions
    regions = parse_regions(mtx.columns)
    regions['reg_id'] = mtx.columns
    regions['start'] = pd.to_numeric(regions['start'])
    regions['end'] = pd.to_numeric(regions['end'])

    # Overlap calculation (Placeholder)
    res = overlap_gene_anno(regions, gene_anno)
    gene_overlap = res['gene_overlap']
    n_dup = res['n_dup']
    gene_uniq = res['gene_uniq']

    if verbose:
        print(f"[I::reg2gene] gene_overlap:")
        print(gene_overlap)

        if n_dup > 0:
            print(f"[W::reg2gene] there are {n_dup} genes overlap with >1 regions!")
        
        print(f"[I::reg2gene] {len(gene_uniq)} genes overlap with 1 region.")

    mtx_gene = mtx.loc[:, gene_uniq['reg_id']].copy()
    mtx_gene.columns = gene_uniq['Gene']

    return {'mtx': mtx_gene, 'overlap': gene_overlap}




def overlap_gene_anno(regions, gene_anno):
    """
    #' Overlap Genes with Regions
    #' @param regions A dataframe containing 4 columns: reg_id <str> Region ID,
    #'   chrom <str>, start <int>, end <int>.
    #' @param gene_anno A data frame. It should contain at least 4 columns:
    #'   1) Gene <chr> gene name; 2) Chr <chr> chrom; 3) start <num> 
    #'   1-based start pos, inclusive; 4) end <num> 1-based end pos, inclusive.
    #' @return A list of 3 elements: 
    #'   gene_overlap, a dataframe containing 2 columns: Gene, <str>, gene name;
    #'     reg_id, <str> ID of overlapping region(s).
    #'   n_dup, an integer, number of genes overlapping more than one region.
    #'   gene_uniq, a dataframe containing 2 columns: Gene and reg_id, only for
    #'     genes overlapping with one region.

    Args:
        regions (_type_): _description_
        gene_anno (_type_): _description_

    Returns:
        _type_: _description_
    """
    # Remove 'chr' from chromosome identifiers
    regions['chrom'] = regions['chrom'].str.replace('chr', '')
    gene_anno['Chr'] = gene_anno['Chr'].str.replace('chr', '')
    
    # Find overlapping genes
    def find_overlaps(gene_row):
        overlaps = regions[
            (regions['chrom'] == gene_row['Chr']) &
            (regions['start'] <= gene_row['end']) &
            (regions['end'] >= gene_row['start'])
        ]
        return pd.DataFrame({'Gene': gene_row['Gene'], 'reg_id': overlaps['reg_id']})
    
    # def find_overlapping_genes_for_region(region, gene_anno):
    #     chrom, start, end = region['chrom'], region['start'], region['end']
        
    #     # Filter genes that overlap with the given region
    #     overlapping_genes = gene_anno[
    #         (gene_anno['Chr'] == chrom) &
    #         (gene_anno['start'] <= end) &
    #         (gene_anno['end'] >= start)
    #     ]
        
    #     return overlapping_genes['Gene'].tolist()
    # # Apply the function to each region
    # regions['overlapping_genes'] = regions.apply(
    #     find_overlapping_genes_for_region, gene_anno=gene_anno, axis=1
    # )

    overlaps_list = gene_anno.apply(find_overlaps, axis=1).tolist()
    gene_overlap = pd.concat(overlaps_list).reset_index(drop=True)

    # Calculate statistics
    gene_stat = gene_overlap.groupby('Gene').size().reset_index(name='n')
    gene_dup = gene_stat[gene_stat['n'] > 1]
    gene_uniq = gene_stat[gene_stat['n'] == 1].merge(gene_overlap, on='Gene')[['Gene', 'reg_id']]

    return {
        'gene_overlap': gene_overlap,
        'n_dup': len(gene_dup),
        'gene_uniq': gene_uniq
    }