import numpy as np
import seaborn as sns

WZcolors = np.array(['#4796d7', '#f79e54', '#79a702', '#df5858', '#556cab', 
                     '#de7a1f', '#ffda5c', '#4b595c', '#6ab186', '#bddbcf', 
                     '#daad58', '#488a99', '#f79b78', '#ffba00'])

def annoHeatMap(X, anno, cax_pos=[.99, .2, .03, .45], cmap="GnBu", **kwargs):
    """Heatmap with annotations"""
    idx = np.argsort(np.dot(X, 2**np.arange(X.shape[1])) + 
                     anno * 2**X.shape[1])

    g = sns.clustermap(X[idx], cmap=cmap, yticklabels=False,
                       col_cluster=False, row_cluster=False,
                       row_colors=WZcolors[anno][idx], **kwargs)
    
    for idx, val in enumerate(np.unique(anno)):
        g.ax_col_dendrogram.bar(0, 0, color=WZcolors[idx],
                                label=val, linewidth=0)
    
    g.ax_col_dendrogram.legend(loc="center", ncol=6, title="True clone")
    g.cax.set_position(cax_pos)

    return g

import matplotlib.pyplot as plt
def scatter_anno(x, y, anno_text = None, alpha = 0.6, 
                 xlabel="sample_chr_total", 
                 ylabel="sample_libsize",
                 plot_title = "check libsize", 
                 **kwargs):
    """
    Parameters
    ----------

    https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.scatter.html
    
    
    Example
    -------
    """
    fig, axs = plt.subplots(1, 1,
                        figsize =(10, 7),
                        tight_layout = True)
    
    # Remove axes splines
    for s in ['top',  'right']:
        axs.spines[s].set_visible(False)
    # 'bottom', 'left',
        
    # Remove x, y ticks
    axs.xaxis.set_ticks_position('none')
    axs.yaxis.set_ticks_position('none')
    
    # Add x, y gridlines
    axs.grid(b = True, color ='grey',
        linestyle ='-.', linewidth = 0.5,
        alpha = 0.6)
    
    ## scatter params
    paras = {}
    paras['x'] = x
    paras['y'] = y
    paras['alpha'] = alpha
    ### update the user assign paras
    paras.update(**kwargs)
    
    axs.scatter(**paras)
    
    if anno_text == None:
        pass
    else:
        for i, txt_ in enumerate(anno_text):
#             axs.annotate(txt_, (x[i], y[i]), fontsize=14)
            axs.annotate(txt_, xy = (x[i], y[i]), xycoords = 'data', xytext=(+20,-20),  textcoords = 'offset points', fontsize=14)
            
    # Adding extra features   
    plt.xlabel(xlabel, fontsize=16)
    plt.ylabel(ylabel, fontsize=16)
#     plt.legend(legend)
    plt.title(plot_title, fontsize=20)
 
    # Show plot
    plt.show()

import scipy.stats as st
from sklearn import linear_model

## from yuanhua hilearn
def corr_plot(x, y, max_num=10000, outlier=0.01, line_on=True,
    legend_on=True, size=30, dot_color=None, outlier_color="r",
    alpha=0.8, color_rate=10, corr_on=None):
    """Correlation plot for large number of dots by showing the density and 
    Pearson's correlatin coefficient
    Parameters
    ----------
    x: `array_like`, (1, ) 
        Values on x-axis
    y: `array_like`, (1, ) 
        Values on y-axis
    max_num: int
        Maximum number of dots to plotting by subsampling
    outlier: float
        The proportion of dots as outliers in different color
    line_on : bool
        If True, show the regression line
    legend_on: bool
        If True, show the Pearson's correlatin coefficient in legend. Replace
        of *corr_on*
    size: float
        The dot size
    dot_color: string
        The dot color. If None (by default), density color will be use
    outlier_color: string
        The color for outlier dot
    alpha : float
        The transparency: 0 (fully transparent) to 1
    color_rate: float
        Color rate for density
    Returns
    -------
    ax: matplotlib Axes
        The Axes object containing the plot.
    Examples
    --------
    .. plot::
        >>> import numpy as np
        >>> from hilearn.plot import corr_plot
        >>> np.random.seed(1)
        >>> x = np.append(np.random.normal(size=200), 5 + np.random.normal(size=500))
        >>> y = 2 * x + 3 * np.random.normal(size=700)
        >>> corr_plot(x, y)
    """
    
    score = st.pearsonr(x, y)
    np.random.seed(0)
    if len(x) > max_num:
        idx = np.random.permutation(len(x))[:max_num]
        x, y = x[idx], y[idx]
    outlier = int(len(x) * outlier)
    
    xy = np.vstack([x,y])
    z = st.gaussian_kde(xy)(xy)
    idx = z.argsort()
    idx1, idx2 = idx[outlier:], idx[:outlier]
    
    if dot_color is None: 
        #c_score = np.log2(z[idx]+100)
        c_score = np.log2(z[idx] + color_rate*np.min(z[idx]))
    else:
        #idx2 = []
        c_score = dot_color
    
    plt.set_cmap("Blues")
    plt.scatter(x[idx], y[idx], c=c_score, edgecolor=None, s=size, alpha=alpha)
    plt.scatter(x[idx2], y[idx2], c=outlier_color, edgecolor=None, s=size/5, 
                alpha=alpha/3.0)#/5

    # plt.grid(alpha=0.4)

    if line_on:
        clf = linear_model.LinearRegression()
        clf.fit(x.reshape(-1,1), y)
        xx = np.linspace(x.min(), x.max(), 1000).reshape(-1,1)
        yy = clf.predict(xx)
        plt.plot(xx, yy, "k--", label="R=%.3f" %score[0])
        # plt.plot(xx, yy, "k--")

    if legend_on or corr_on:
        plt.legend(loc="best", fancybox=True, ncol=1)
        # plt.annotate("R=%.3f\np=%.1e" %score, fontsize='x-large', 
        #             xy=(0.97, 0.05), xycoords='axes fraction',
        #             textcoords='offset points', ha='right', va='bottom')


## confusion matrix heatmap
def confuse_heatmap(confuse_mat_df, cmap = "Blues", version1 = True, version2 =True, plot_xlabel = None,
                    plot_ylabel = None, plt_title = None, save_file_name = None):
    """
    Function:
    confusion matrix heatmap
    
    confuse_mat_df: confusion matrix in pd.DataFrame
    
    Example:
    --------
    >>> confuse_mat, ids1_uniq, ids2_uniq = get_confusion(CopyKAT_lst, expression_lst)
    >>> confuse_mat_df = get_confuse_mat_df(confuse_mat,
            index_names=["cloneA", "cloneB","Normal"],
            clolums_names=["cloneA", "cloneB","Normal"])
    >>> confuse_heatmap(confuse_mat_df,  plot_xlabel = "CopyKAT",
                    plot_ylabel = "Ground Truth", 
                    plt_title = "Concordance in subclone identification", 
                    save_file_name = 'CopyKAT_vs_Groudtruth1.pdf')
    """
    import matplotlib.pyplot as plt
    import seaborn as sns
    ## setting
    font = {'family' : 'DejaVu Sans',
            'size'   : 15}
    plt.rc('font', **font)
    if version1:
        sns.heatmap(confuse_mat_df, annot=True, fmt="d", cmap=cmap)
    if version2:
        ## data normalization
        plot_df = (confuse_mat_df.T/confuse_mat_df.sum(axis=1).values).T
        ## plot percentage
        plt.figure(figsize=[8,6])
        ax = sns.heatmap(plot_df, cmap=cmap) # annot=True,
        ## annot original count
        height, width = np.array(confuse_mat_df).shape
        text_colors = ['black', 'white']
        for x in range(width):
            for y in range(height):
                ax.annotate(str(np.array(confuse_mat_df)[y][x]), xy=(x+0.5, y+0.5), 
                            horizontalalignment='center',
                            verticalalignment='center', color=text_colors[int(np.array(plot_df)[y][x] > 0.5)], fontsize = 15)
        if plot_xlabel is None:
            pass
        else:
            plt.xlabel(plot_xlabel)
        if plot_ylabel is None:
            pass
        else:
            plt.ylabel(plot_ylabel)
        # plt.yticks(rotation=45)
        # plt.xticks(rotation=45)
        # plt.xticks(range(3), confusion_matrix.columns) #, rotation=315)
        # plt.yticks(range(len(norm_conf)), set(clone_id))
        plt.tight_layout()
        if plt.title is None:
            pass
        else:
            plt.title(plt_title, fontsize = 18)          
            # plt.title('Concordance in subclone identification', fontsize = 18)
        if save_file_name is None:
            plt.show()
        else:
            plt.savefig(save_file_name, dpi=300, bbox_inches='tight') 
            plt.show()
            #plt.savefig('CopyKAT_vs_Groudtruth.pdf', dpi=300, bbox_inches='tight') 