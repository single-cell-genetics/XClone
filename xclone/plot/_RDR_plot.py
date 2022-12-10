"""Base functions for XClone RDR plotting.
"""

# Author: Rongting Huang
# Date: 2021/07/22
# update: 2021/12/16

##==============================
## color
# https://cran.r-project.org/web/packages/viridis/vignettes/intro-to-viridis.html

## develop stage - check lib size work or not
import matplotlib.pyplot as plt
import matplotlib
from matplotlib import colors
# from matplotlib.ticker import PercentFormatter

def X_hist_orig(X, bins_num = 50, cmap_name = 'viridis', plot_title = None, legend = "distribution",
           xlabel = 'xaxis', ylabel = 'yaxis', add_text = None, **kwargs):
    """
    For testing---
    X:dataset
    bins_num
    
    cmap_name: 'viridis', 'magma', 'magma_r'
    
    **kwargs  for plt.hist func
    
    # https://www.geeksforgeeks.org/plotting-histogram-in-python-using-matplotlib/
    
    """
    fig, axs = plt.subplots(1, 1,
                        figsize =(10, 7),
                        tight_layout = True)
 
    # Remove axes splines
    for s in ['top', 'bottom', 'left', 'right']:
        axs.spines[s].set_visible(False)
        
    # Remove x, y ticks
    axs.xaxis.set_ticks_position('none')
    axs.yaxis.set_ticks_position('none')
    
    # Add padding between axes and labels
    axs.xaxis.set_tick_params(pad = 5)
    axs.yaxis.set_tick_params(pad = 10)
    
    # Add x, y gridlines
    axs.grid(b = True, color ='grey',
        linestyle ='-.', linewidth = 0.5,
        alpha = 0.6)
    
     
    # Add Text mark
    axs.text(0.7, 0.3, add_text,
         fontsize = 12,
         color ='red',
         ha ='left',
         va ='bottom',
         alpha = 0.7)
    
    ## hist plot params setting
    paras = {}
    paras['x'] = X
    paras['bins'] = bins_num
    paras['density'] = False
    
#     ### update the default paras
#     args, legend_paras = set_XXheatmap_args(legend_pos,legend_mode)
#     paras.update(**args)
    ### update the user assign paras
    paras.update(**kwargs)
    
    # Creating histogram
#     N_count, bins, patches = axs.hist(X, bins = bins_num, density=False)
    N_count, bins, patches = axs.hist(**paras)
    
    # Setting color
    ## todo maybe log N_count?
    fracs = ((N_count**(1 / 5)) / N_count.max())
    norm = colors.Normalize(fracs.min(), fracs.max())
    
    hist_cmap = matplotlib.cm.get_cmap(cmap_name)
    for thisfrac, thispatch in zip(fracs, patches):
        color_ = hist_cmap(norm(thisfrac))
        thispatch.set_facecolor(color_)
    
#     # Adding extra features   
#     axs.set_xlabel(xlabel)
#     axs.set_ylabel(ylabel)
#     axs.legend(legend)
#     axs.set_title(plot_title)
    
    # Adding extra features   
    plt.xlabel(xlabel, fontsize=16)
    plt.ylabel(ylabel, fontsize=16)
    plt.legend(legend)
    plt.title(plot_title, fontsize=20)
 
    # Show plot
    plt.show()
    

def X_hist_base(X, axs, bins_num = 50, cmap_name = 'viridis', plot_title = None, label = "distribution",
                xlabel = 'xaxis', ylabel = 'count', add_text = None, **kwargs):
    """
    Base hist function in XClone
    color setting for better reading
    params:
    -----
    X:dataset
    bins_num

    ----
    

    Notes:
    -----
    take care for plt.hist  data is x
    for my X_hist_base data is X
    """
    # Setting for the base hist
    ## Remove axes splines
    for s in ['top', 'bottom', 'left', 'right']:
        axs.spines[s].set_visible(False)
        
    ## Remove x, y ticks
    axs.xaxis.set_ticks_position('none')
    axs.yaxis.set_ticks_position('none')
    
    ## Add padding between axes and labels
    axs.xaxis.set_tick_params(pad = 5)
    axs.yaxis.set_tick_params(pad = 10)
    
    ## Add x, y gridlines
    axs.grid(b = True, color ='grey', linestyle ='-.', linewidth = 0.5, alpha = 0.5)

    # Creating histogram
    ## params setting for plt.hist
    paras = {}
    paras['x'] = X
    paras['bins'] = bins_num
#     paras['density'] = False
    
#     ### update the default paras
#     args, legend_paras = set_XXheatmap_args(legend_pos,legend_mode)
#     paras.update(**args)
    ### update the user assign paras
    paras.update(**kwargs)
    
    # plot histogram
#     N_count, bins, patches = axs.hist(X, bins = bins_num, density=False)
    N_count, bins, patches = axs.hist(**paras)
    

    ##  Setting color
    fracs = ((N_count**(1 / 5)) / N_count.max())
    norm = colors.Normalize(fracs.min(), fracs.max())
    
    hist_cmap = matplotlib.cm.get_cmap(cmap_name)
    for thisfrac, thispatch in zip(fracs, patches):
        color_ = hist_cmap(norm(thisfrac))
        thispatch.set_facecolor(color_)
    
    # Adding extra features
    ## Add Text mark
    m_brk = round(bins_num/2)
    if add_text != None:
        axs.text(bins[m_brk], N_count[m_brk], add_text, fontsize = 12,
                 color ='black', ha ='right', va ='bottom', alpha = 0.7)
    
    ## label and title
    axs.set_xlabel(xlabel, fontsize=16)
    axs.set_ylabel(ylabel, fontsize=16)
    # axs.set_title(plot_title)


def X_hist_multiplots(X_data, bins_num = 50, cmap_name = 'viridis', plot_title = 'customised histogram', sharex=True, figsize = (10,22), save_fig = False, **kwargs):
    """
    Function:
    -----
    genome based
    
    Example:
    xclone.pl.X_hist_multiplots(rr_ad, plot_title = 'ratio histogram genome') #, save_fig = "ratio_hist_genome_based.pdf"
    """
    # Creating fig and subplots
    n_plots = X_data.X.shape[0] 
    fig, axes = plt.subplots(n_plots, 1, figsize = figsize, tight_layout = True, sharex=sharex)
    
    # common params setting for X_hist_base
    paras = {}
    paras['bins_num'] = bins_num
    paras['cmap_name'] = cmap_name
    paras['density'] = False
    
#     ### update the default paras
#     args, legend_paras = set_XXheatmap_args(legend_pos,legend_mode)
#     paras.update(**args)
    ### update the user assign paras
    paras.update(**kwargs)
    
    for i in range(n_plots):
        paras['X'] = X_data.X[i]
        paras['axs'] = axes[i]
#         paras['bins_num'] = bins_num
        paras['xlabel'] = X_data.obs["cell_type"][i]
#         X_hist_base(X_data.X[i], axes[i], bins_num = bins_num, xlabel = X_data.obs["cell_type"][i])
        X_hist_base(**paras)

    # Adding extra features for whole plot
    fig.suptitle(plot_title, fontsize=20)
    
    # Show plot
    if save_fig != False:
        plt.savefig(save_fig)
    else:
        plt.show()
        
def X_hist_multiplots_CHR(X_data, bins_num = 50, cmap_name = 'viridis', 
                          plot_title = 'customised histogram', sharex=True, sharey=True, figsize = (50,22), save_fig = False, 
                          chr_show_num = 22, chr_lst = None, **kwargs):
    """
    Function:
    chr based

    --
    if chr_lst is specified, then chr_show_num here is useless
    
    -----
    Example:
    ## GX109 Raw ratio
    ### show first 4 col here, default show chr1-22
    xclone.pl.X_hist_multiplots_CHR(rr_ad, plot_title = 'ratio histogram', chr_show_num = 4, figsize = (60,22))
    ### select chr to shows
    xclone.pl.X_hist_multiplots_CHR(rr_ad, plot_title = 'ratio histogram', chr_lst = ["8", "18", "20" ], figsize = (60,22))

    """
    
    # Creating fig and subplots
    plots_row = X_data.X.shape[0]
    if chr_lst == None:
        plots_col = chr_show_num
    else:
        plots_col = len(chr_lst)
        
    fig, axs = plt.subplots(plots_row, plots_col, figsize = figsize, tight_layout = True, sharex=sharex, sharey=sharey)
    
    # common params setting for X_hist_base
    paras = {}
    paras['bins_num'] = bins_num
    paras['cmap_name'] = cmap_name
    paras['density'] = False
    
#     ### update the default paras
#     args, legend_paras = set_XXheatmap_args(legend_pos,legend_mode)
#     paras.update(**args)
    ### update the user assign paras
    paras.update(**kwargs)
    
    if chr_lst == None:
        for i in range(plots_row):
            for j in range(chr_show_num):
                tmp_flag = X_data.var["chr"] == str(j+1)
                chr_num = tmp_flag.sum()
                tmp_data = X_data.X[i][tmp_flag]
            
                paras['X'] = tmp_data
                paras['axs'] = axs[i, j]
                paras['xlabel'] = X_data.obs["cell_type"][i]
                paras['ylabel'] = "chr" + str(j+1)
                paras['add_text'] = str(chr_num)
                # X_hist_base(tmp_data, axs[i, j], bins_num = bins_num, xlabel = X_data.obs["cell_type"][i], ylabel = "chr" + str(j+1), add_text = str(chr_num))
                X_hist_base(**paras)
                
    else:
        for i in range(plots_row):
            for j in range(len(chr_lst)):
                chr_ = chr_lst[j]
                tmp_flag = X_data.var["chr"] == chr_
                chr_num = tmp_flag.sum()
                tmp_data = X_data.X[i][tmp_flag]
            
                paras['X'] = tmp_data
                paras['axs'] = axs[i, j]
                paras['xlabel'] = X_data.obs["cell_type"][i]
                paras['ylabel'] = "chr" + chr_
                paras['add_text'] = str(chr_num)
                X_hist_base(**paras)
        
    # Adding extra features
    fig.suptitle(plot_title, fontsize=20)
    
    # Show plot
    if save_fig != False:
        plt.savefig(save_fig)
    else:
        plt.show()