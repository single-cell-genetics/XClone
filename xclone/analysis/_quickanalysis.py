
import matplotlib.pylab as plt
import matplotlib.patches as mpatches
from ._clustering import XClustering, Cluster_mapping
from .evaluation import get_confusion, get_confuse_mat_df
from ..plot.base_plot import confuse_heatmap
from ..plot.CNV_plot import Complex_BAF_CNV_visualization, BAF_CNV_visualization

def exploreClustering(adata, ref_anno_key = "spot_anno", Xlayer = "posterior_mtx", 
                      min_clusters = 2, max_clusters=5):
        """
        Args:
            adata (_type_): _description_
            Xlayer (str, optional): _description_. Defaults to "posterior_mtx".
            min_clusters (int, optional): _description_. Defaults to 2.
            max_clusters (int, optional): _description_. Defaults to 5.
        """
        for n_clusters in range(min_clusters, max_clusters+1):
             Z, embedding = XClustering(adata, Xlayer, n_clusters = n_clusters)
             plt.scatter(embedding[:, 0], embedding[:, 1], c = Z)
             plt.title(str(n_clusters) + " clsuters")
             plt.show()
        data = adata.obs[ref_anno_key]
        category_mapping = {category: idx for idx, category in enumerate(data.cat.categories)}
        numerical_data = data.map(category_mapping)

        # plt.scatter(embedding[:, 0], embedding[:, 1], c = numerical_data)
        scatter = plt.scatter(embedding[:, 0], embedding[:, 1], c=numerical_data, alpha=0.2, s=12, cmap = "summer") #, cmap='viridis'
        plt.title(ref_anno_key)
        # Create a custom legend
        handles = [
        mpatches.Patch(color=scatter.cmap(scatter.norm(category_mapping[cat])), label=cat)
        for cat in data.cat.categories
        ]
        plt.legend(handles=handles, title="ref_anno_key")
        plt.show()
        return None

    

def OnestopBAFClustering(adata, Xlayer = "posterior_mtx", 
                         n_clusters = 2, ref_anno_key = "spot_anno", clone_anno_key = "clone(2)", 
                         plot_title = "HCC-3 (XClone)",
                         file_save_path = "./", file_save_name =  "Liver_HCC-3_2clones",
                         complex_plot = False):
    """
    OnestopBAFClustering for allele bias states probability analysis.
    """
    # clustering
    Z, embedding = XClustering(adata, Xlayer, n_clusters = n_clusters)

    plt.scatter(embedding[:, 0], embedding[:, 1], c = Z)
    plt.show()

    adata.obs[clone_anno_key] = Z
#     adata.obs[clone_anno_key] = adata.obs[clone_anno_key].astype(str)
#     adata.obs[clone_anno_key] = adata.obs[clone_anno_key].astype('category')
    adata = Cluster_mapping(adata,
                    pre_anno = clone_anno_key,
                    after_anno = clone_anno_key,
                    mapping_mode = 1)

    # confuse analysis
    confuse_mat, ids1_uniq, ids2_uniq = get_confusion(adata.obs[ref_anno_key], adata.obs[clone_anno_key])
    ids1_uniq
    ids2_uniq

    confuse_mat_df = get_confuse_mat_df(confuse_mat,
            index_names=ids1_uniq,
            columns_names=ids2_uniq)
    
    confuse_heatmap_savefile = file_save_path + file_save_name + "confusemap" + str(n_clusters) + "cluster.pdf"

    confuse_heatmap(confuse_mat_df,  plot_xlabel = "CNA_Clones",
                    plot_ylabel = "ref_annos",
                    plt_title = "Concordance in subclone identification",
                    save_file_name = confuse_heatmap_savefile)
    
    # heatmap plotting
    out_path = file_save_path + file_save_name
    if complex_plot:
        cell_anno_key_plot = [clone_anno_key, ref_anno_key]
        clusters_display_name = [clone_anno_key, ref_anno_key]
        Complex_BAF_CNV_visualization(adata, Xlayer = Xlayer,
                                  cell_anno_key = cell_anno_key_plot, clusters_display_name = clusters_display_name, 
                                  title = plot_title, save_file = True, out_file = out_path)
    else:
        BAF_CNV_visualization(adata, Xlayer = Xlayer,
                              cell_anno_key = clone_anno_key, clusters_display_name = clone_anno_key,
                                  title = plot_title, save_file = True, out_file = out_path)
    
    
    anno_save_file = file_save_path + file_save_name + str(n_clusters) + "clone_anno.csv"
    adata.obs.to_csv(anno_save_file)
    
    return adata.obs
    
    
    