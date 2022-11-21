"""Base functions for XClone data manipulation
"""
# Author: Rongting Huang
# Date: 2022/06/29
# update: 2022/06/29

## Part I
def Xdata_region_selection(Xdata,
                           select = False,
                           chr_anno_key = None,
                           chr_lst = None,
                           update_uns = True,
                           uns_anno_key = None):
    """
    Function:
    select/remove chr regions.
    if the chr_anno_key is None, then use the chr.
    chr_anno_key could be "bin_idx_cum", "chr", "chr_arm"...
    
    uns part is taken into consideration.
    if the uns_anno_key is None, then use the default important lst,
    mainly include the calculated probability.
    """
    if chr_anno_key is None:
        chr_anno_key = "chr"
    if chr_lst is None:
        raise ValueError("[XClone data selection]-pls check chr_lst.")
    
    gene_total_num = Xdata.var.shape[0]
    region_flag = Xdata.var[chr_anno_key].isin(chr_lst)
    if select == True:
        print("[XClone-data selection]:")
        print("select %d regions / %d total regions" 
               %(region_flag.sum(), gene_total_num))

    else:
        print("[XClone-data removing]:")
        print("Filter out %d genes / %d total genes, remain %d regions" 
               %(region_flag.sum(), gene_total_num, gene_total_num-region_flag.sum()))
        region_flag = ~region_flag
    
    Xdata_subset = Xdata[:, region_flag].copy()

    ## uns selection:mainly for probability in shape[cells, regions, states]
    if update_uns == True:
        if uns_anno_key is None:
            ## dic key search include _prob_ item.
            for key_ in Xdata_subset.uns:
                if "_prob_" in key_:
                    Xdata_subset.uns.update({key_: Xdata.uns[key_][:, region_flag,:].copy()})
        else:
            for key_ in uns_anno_key:
                Xdata_subset.uns.update({key_: Xdata.uns[key_][:, region_flag,:].copy()})

    return Xdata_subset

def Xdata_cell_selection(Xdata,
                         select = False,
                         cell_anno_key = None,
                         cell_lst = None,
                         update_uns = True,
                         uns_anno_key = None):
    """
    Function:
    select/remove cells.
    if cell_anno_key is None, then use cell_type;

    uns part is taken into consideration.
    if the uns_anno_key is None, then use the default important lst,
    mainly include the calculated probability.
    """

    if cell_anno_key is None:
        cell_anno_key = "cell_type"
    if cell_lst is None:
        raise ValueError("[XClone data selection]-pls check cell_lst.")
    cell_total_num = Xdata.obs.shape[0]
    cell_flag = Xdata.obs[cell_anno_key].isin(cell_lst)
    if select == True:
        print("[XClone-data selection]:")
        print("select %d cells / %d total cells" 
               %(cell_flag.sum(), cell_total_num))
    else:
        print("[XClone-data removing]:")
        print("Filter out %d cells / %d total cells, remain %d cells" 
               %(cell_flag.sum(), cell_total_num, cell_total_num-cell_flag.sum()))
        cell_flag = ~cell_flag
    
    Xdata_subset = Xdata[cell_flag, :].copy()

    ## uns selection:mainly for probability in shape[cells, regions, states]
    if update_uns == True:
        if uns_anno_key is None:
            ## dic key search include _prob_ item.
            for key_ in Xdata_subset.uns:
                if "_prob_" in key_:
                    Xdata_subset.uns.update({key_: Xdata.uns[key_][cell_flag,:,:].copy()})
        else:
            for key_ in uns_anno_key:
                Xdata_subset.uns.update({key_: Xdata.uns[key_][cell_flag ,:,:].copy()})
    return Xdata_subset


