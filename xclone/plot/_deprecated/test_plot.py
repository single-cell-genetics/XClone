def test_plot():
    """
    From weizhong.
    """
    import matplotlib.patches as mpatches
    sns.set(rc={"figure.dpi":600, 'savefig.dpi':600})
    f=sns.clustermap(pr_ay,col_cluster=True, cmap='RdBu_r', 
             colors_ratio=(0.13,0),
             dendrogram_ratio=(0.42,0.15), figsize=(4,6), row_colors=row_color.to_numpy(),
             cbar_pos=(0.05, 0.88, 0.4, 0.03), center=0, 
             cbar_kws={'orientation':'horizontal'},
                )

    ax2 = f.fig.add_axes((0.05, 0.91, 0.4, 0.05))
    # ax2.set_aspect(aspect=0.05)
    # ax2.set_ylim(-0.1,3)
    ax2.set_axis_off()

    f.ax_heatmap.set_yticks([])
    f.ax_heatmap.set_xticklabels(['MPA-H', 'RP-293T','RP-pc3'], fontsize=8)

    for i in range(3):
        sns.kdeplot(pr_ay[:,i], bw_adjust=0.3, ax = ax2, alpha=0.8, linewidth=1.2)

    f.ax_row_dendrogram.legend(handles=handle_patches,
                           frameon=False, fontsize=8.2, 
                           bbox_to_anchor=(0.87, 0.5),
                          )

    handle_patches = [mpatches.Patch(label=lc, color=c) for lc,c in zip(label, motifs_color)]