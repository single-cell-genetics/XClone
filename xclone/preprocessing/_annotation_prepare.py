"""Base functions for XClone annotation data 
preprocessing/preparation.
"""

# Author: Rongting Huang
# Date: 2021/05/04
# update: 2021/07/21


import pandas as pd
import numpy as np

## Part I: For hg19/38 blocks/genes chr_region annotation preprocessing (2021/05/05)
def resort_chr_df(chr_df):
    """
    Function:
    sort the dataframe as the assignned chromosome index.
    """
    assign_sort_index=['1', '2','3', '4', '5', '6', '7', '8', '9','10', '11', '12', '13', '14', '15', '16', '17', '18', '19','20', '21', '22', 'X', 'Y']
    cnt = 0
    for i in assign_sort_index:
        tmp_df = chr_df[chr_df["chr"] == i]
        if cnt==0:
            resort_df = tmp_df
        else:
            resort_df = pd.concat([resort_df,tmp_df],ignore_index=True)
        cnt+=1
    if cnt == len(assign_sort_index):
        print("sorted %d chromosomes" %cnt)
    return resort_df


def get_chr_arm_df(chr_breakpoint_df,regions,sort=True):
    """
    Function:
    annotate p/q info for regions according to the chr arm breakpoints;
    and sort with normal chromosome index.
    
    So, actually just processing 24 chromosomes and filter out the unknown genes.
    
    Example:
    ========================
    input:
    1. hg19_chr_breakpoint with columns=["chr","p_start","p_stop"]
    2. blocks_50kb_region with columns=["chr","start","stop"]
      or genes_region with columns=["chr","start","stop"]
    usage:
    chr_arm_blocks = get_chr_arm_df(hg19_chr_breakpoint, blocks_50kb_region)
    """
    print("build chr_list and arm_list------------")
    arm_list = []
    chr_list = []
    for chr_ in np.unique(chr_breakpoint_df["chr"]):
        tmp_breakpoint = int(chr_breakpoint_df[chr_breakpoint_df["chr"]==chr_]["p_stop"])
        print("chr:", chr_, "  breakpoint: ", tmp_breakpoint)
        flag_ = regions["chr"]==chr_
        flag1 = regions[flag_]["stop"] <= tmp_breakpoint
        print("chr records: ",flag_.sum(), "p_arm records: ", flag1.sum())
        arm_list.append(['p']*flag1.sum())
        arm_list.append(['q']*(flag_.sum() - flag1.sum()))
        chr_list.append([chr_]*flag_.sum())
    print("processing chr--------------------")
    chr_list_all = []
    for chrlst_item in chr_list:
        for chr_item in chrlst_item:
            chr_list_all.append(chr_item)
    print(len(chr_list_all))
    print("processing arm---------------------")
    arm_list_all = []
    for pqlst_item in arm_list:
        for pq_item in pqlst_item:
            arm_list_all.append(pq_item)
    print(len(arm_list_all))
    
    chr_arm_df = pd.DataFrame({"chr":chr_list_all,"arm":arm_list_all})
    if sort:
        resort_chr_arm_df = resort_chr_df(chr_arm_df)
        return resort_chr_arm_df
    else:
        return chr_arm_df


def concat_df(sort_regions, sort_chr_arm_df):
    """
    concat the region info and chr_arm info
    These two dfs should be in the same chromosome order.
    ========================
    Example:
    annotate_blocks = concat_df(blocks_50kb_region, chr_arm_blocks)
    annotate_genes = concat_df(sort_genes_region, chr_arm_genes)
    """
    concat_df = pd.concat([sort_regions, sort_chr_arm_df],axis=1)
    concat_df.columns = ["chr","start","stop","another_chr","arm"]
    ## check the concat status
    check_num = (concat_df["chr"] == concat_df["another_chr"]).sum()
    if len(concat_df) == check_num:
        print("concat success!!")
    else:
        print("Error! Pls check chromosome order!")
    out_put_df = concat_df.drop(columns=['another_chr'])
    return out_put_df



## chr info prepare
