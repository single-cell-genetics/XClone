import random
import pandas as pd
import numpy as np
import anndata as ad

from ..preprocessing import load_anno
from .CNV_plot import Combine_CNV_visualization



def calculate_count(total_cells, clone_meta):
    clone_meta['generated_count'] = (clone_meta['count'] / clone_meta['count'].sum()) * total_cells

    # Round to nearest integer
    clone_meta['generated_count'] = clone_meta['generated_count'].round()

    # Adjust the total to match exactly the total_cells
    difference = int(total_cells - clone_meta['generated_count'].sum())

    # Distribute the difference
    indices = clone_meta.index.to_list()
    adjustments = np.random.choice(indices, abs(difference), replace=True)

    for i in adjustments:
        clone_meta.at[i, 'generated_count'] += np.sign(difference)

    # Ensure the sum matches the total_cells
    assert clone_meta['generated_count'].sum() == total_cells
    
    return clone_meta



def generate_barcode():
    return ''.join(random.choices('ATCG', k=16))


def generate_barcode_lst(barcodes_num):
    # Set to store unique barcodes
    unique_barcodes = set()

    # Generate 400 unique barcodes
    while len(unique_barcodes) < barcodes_num:
        unique_barcodes.add(generate_barcode())

    # Convert to a list for easy handling
    barcode_list = list(unique_barcodes)
    return barcode_list  
        
    

def Generate_adata(total_cells, clone_meta):
    barcode_list = generate_barcode_lst(total_cells)
    clones_meta = calculate_count(total_cells, clone_meta)
    
    # Generate the list of annotation labels
    label_list = []
    # for _, row in clones_meta.iterrows():
    #     label_list.extend([row['simulated_label']] * row['generated_count'])
    for _, row in clones_meta.iterrows():
        label_list.extend([row['simulated_label']] * int(row['generated_count']))
        
        
    adata_obs = pd.DataFrame(label_list, index=barcode_list)
    adata_obs.index.name = None
    adata_obs.columns = ["annotation"]
    
    adata_var = load_anno(genome_mode = "hg38_genes_select")
    

    # Create AnnData object
    adata = ad.AnnData(X=np.zeros((len(adata_obs), len(adata_var))),
                       obs=adata_obs,
                       var=adata_var)
    adata.var["chr"] = adata.var["chr"].astype(str)
    return adata
       
    
    
# Function to determine label
def assign_label(row, mode=1):
    total = row['allele A'] + row['allele B']
    if mode == 1:
        if total > 2:
            return 'copy gain'
        elif total < 2:
            return 'copy loss'
        elif total == 2 and (row['allele A'] == 0 or row['allele B'] == 0):
            return 'loh'
        else:
            return 'copy neutral'
    elif mode == 2:
        if total > 2:
            return 'copy gain'
        elif total < 2 and row['allele A'] == 0:
            return 'allele A loss'
        elif total < 2 and row['allele B'] == 0:
            return 'allele B loss'
        elif total == 2 and (row['allele A'] == 0 or row['allele B'] == 0):
            return 'loh'
        else:
            return 'copy neutral'
        


def generate_illustrate_GT(adata, cnv_profile, mode=1):
    ## init layer for plotting
    Xlayer = "illustrate_GT"
    
    if mode == 1:
        adata.layers[Xlayer] = np.zeros((len(adata.obs), len(adata.var), 4))
        adata.layers[Xlayer][:,:, 2] = 1
    if mode == 2:
        adata.layers[Xlayer] = np.zeros((len(adata.obs), len(adata.var), 5))
        adata.layers[Xlayer][:,:, 3] = 1
        
    
    # cnv_profile['cnvlabel'] = cnv_profile.apply(assign_label, axis=1)
    cnv_profile['cnvlabel'] = cnv_profile.apply(lambda row: assign_label(row, mode=mode), axis=1)
    
    
    anndata_lst = []
    adata_normal = adata[adata.obs["annotation"] == "normal", :].copy()
    anndata_lst.append(adata_normal)
    
   
    clone_lst = cnv_profile["clone"].unique()
    for clone_label in clone_lst:
        adata_tmp = adata[adata.obs["annotation"] == clone_label, :].copy()
        cnv_profile_tmp = cnv_profile[cnv_profile["clone"] == clone_label]
        for idx_ in cnv_profile_tmp.index:
            if cnv_profile_tmp.loc[idx_]["start"] == 1 and cnv_profile_tmp.loc[idx_]["end"] == np.inf:
                flag = adata_tmp.var["chr"].isin([cnv_profile_tmp.loc[idx_]["chr"].strip("chr"), ])
            else:
                flag = np.logical_and(
                    np.logical_and(
                        adata_tmp.var["chr"] == cnv_profile_tmp.loc[idx_]["chr"].strip("chr"),
                        adata_tmp.var["start"] <= cnv_profile_tmp.loc[idx_]["end"]
                    ),
                    adata_tmp.var["stop"] >= cnv_profile_tmp.loc[idx_]["start"]
                )
            
            if mode == 1:
                if cnv_profile_tmp.loc[idx_]["cnvlabel"] == "copy gain":
                    adata_tmp.layers[Xlayer][:, flag, 3] = 1
                    adata_tmp.layers[Xlayer][:, flag, 2] = 0
                    adata_tmp.layers[Xlayer][:, flag, 1] = 0
                    adata_tmp.layers[Xlayer][:, flag, 0] = 0
                if cnv_profile_tmp.loc[idx_]["cnvlabel"] == "copy loss":
                    adata_tmp.layers[Xlayer][:, flag, 3] = 0
                    adata_tmp.layers[Xlayer][:, flag, 2] = 0
                    adata_tmp.layers[Xlayer][:, flag, 1] = 0
                    adata_tmp.layers[Xlayer][:, flag, 0] = 1
                if cnv_profile_tmp.loc[idx_]["cnvlabel"] == "loh":
                    adata_tmp.layers[Xlayer][:, flag, 3] = 0
                    adata_tmp.layers[Xlayer][:, flag, 2] = 0
                    adata_tmp.layers[Xlayer][:, flag, 1] = 1
                    adata_tmp.layers[Xlayer][:, flag, 0] = 0
            elif mode == 2:
                if cnv_profile_tmp.loc[idx_]["cnvlabel"] == "copy gain":
                    adata_tmp.layers[Xlayer][:, flag, 4] = 1
                    adata_tmp.layers[Xlayer][:, flag, 3] = 0
                    adata_tmp.layers[Xlayer][:, flag, 2] = 0
                    adata_tmp.layers[Xlayer][:, flag, 1] = 0
                    adata_tmp.layers[Xlayer][:, flag, 0] = 0
                if cnv_profile_tmp.loc[idx_]["cnvlabel"] == "allele A loss":
                    adata_tmp.layers[Xlayer][:, flag, 4] = 0
                    adata_tmp.layers[Xlayer][:, flag, 3] = 0
                    adata_tmp.layers[Xlayer][:, flag, 2] = 0
                    adata_tmp.layers[Xlayer][:, flag, 1] = 0
                    adata_tmp.layers[Xlayer][:, flag, 0] = 1
                if cnv_profile_tmp.loc[idx_]["cnvlabel"] == "allele B loss":
                    adata_tmp.layers[Xlayer][:, flag, 4] = 0
                    adata_tmp.layers[Xlayer][:, flag, 3] = 0
                    adata_tmp.layers[Xlayer][:, flag, 2] = 0
                    adata_tmp.layers[Xlayer][:, flag, 1] = 1
                    adata_tmp.layers[Xlayer][:, flag, 0] = 0
                if cnv_profile_tmp.loc[idx_]["cnvlabel"] == "loh":
                    adata_tmp.layers[Xlayer][:, flag, 4] = 0
                    adata_tmp.layers[Xlayer][:, flag, 3] = 0
                    adata_tmp.layers[Xlayer][:, flag, 2] = 1
                    adata_tmp.layers[Xlayer][:, flag, 1] = 0
                    adata_tmp.layers[Xlayer][:, flag, 0] = 0
        anndata_lst.append(adata_tmp)
        
    adata_concatenated = ad.concat(anndata_lst, axis=0)
    adata_concatenated.var = adata_normal.var
    
    return adata_concatenated 



def Plot_Simulation_GT(clone_meta_path, cnv_profile_path, out_fig_path = "auto", total_cells = 400, normal = True, mode = 1, plot = True):
    """
    Plotting for simulation ground truth.
    
    mode = 1: copy gain, copy loss, loh and copy neutral.
    mode = 2: copy gain, allele-specific copy loss, loh and copy neutral.
    """
    column_names = ["simulated_label", 'seed_label', 'count']
    clone_meta = pd.read_table(clone_meta_path, names=column_names, header=None)
    
    cnv_profile = pd.read_table(cnv_profile_path, header=None)
    if cnv_profile.shape[1] == 6:
        column_names = ["chr", 'start', 'end', "clone", "allele A", "allele B"]
        cnv_profile = pd.read_table(cnv_profile_path, names=column_names, header=None)
    elif cnv_profile.shape[1] == 7:
        column_names = ["chr", 'start', 'end', "label", "clone", "allele A", "allele B"]
        cnv_profile = pd.read_table(cnv_profile_path, names=column_names, header=None)
    else:
        raise ValueError("invalid CNA profile, expect 6 or 7 columns.")
    
    adata = Generate_adata(total_cells, clone_meta)
    
    adata_concatenated = generate_illustrate_GT(adata, cnv_profile, mode=mode)
    if normal:
        pass
    else:
        adata_concatenated = adata_concatenated[adata_concatenated.obs["annotation"] != "normal"]
    if plot:
        if mode == 1:
            #combine_res_base_fig = "GT_Simulation_test.png"
            if out_fig_path == "auto":
                out_fig_path = "GT_Simulation_test.png"
            plot_cell_anno_key = "annotation"
            Combine_CNV_visualization(adata_concatenated, Xlayer = "illustrate_GT", 
                    cell_anno_key = plot_cell_anno_key, save_file = True, out_file = out_fig_path, title = None)
        elif mode == 2:
            #combine_res_base_fig = "GT_Simulation_allele_test.png"
            if out_fig_path == "auto":
                out_fig_path = "GT_Simulation_allele_test.png"
            plot_cell_anno_key = "annotation"
            colorbar_ticks = [0,1,2,3,4]
            colorbar_label = ["copy lossA", "copy lossB", "LOH", "copy neutral", "copy gain"]
            Combine_CNV_visualization(adata_concatenated, Xlayer = "illustrate_GT", 
                                                    cell_anno_key = plot_cell_anno_key, 
                                                    color_map_name = "combine_cmap2", 
                                                    states_num = 5,
                                                    colorbar_ticks = colorbar_ticks,
                                                    colorbar_label = colorbar_label,
                                                    save_file = True, out_file = out_fig_path, title = None)

    return adata_concatenated