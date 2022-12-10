## From BAF_plot
def visualize_cell_BAF1(Xdata,
                       Xlayer = "BAF", 
                       chr_lst = None, 
                       plot_title = "visualize_cell_BAF test",
                       save_fig = None, 
                       **kwargs):
    """
    ## no chr annotation
    ## plt.imshow
    ## will be deprecated in future
    deprecated 2022-07-27
    """
    import matplotlib.pyplot as plt
    # fig = plt.figure(figsize=(6, 5))

    if chr_lst is None:
        plt.imshow(Xdata.layers[Xlayer], interpolation='none', 
                aspect='auto', cmap='coolwarm') #seismic, bwr

    else:
        FLAG_ = Xdata.var["chr"].isin(chr_lst)
        plt.imshow(Xdata[:,FLAG_].layers[Xlayer], interpolation='none', 
                aspect='auto', cmap='coolwarm') #seismic, bwr
        plt.xlabel("Chr" + str(chr_lst))

    plt.colorbar()
    plt.tight_layout()
    plt.ylabel("cells")
    plt.title(plot_title)
    plt.show()

## from _visualize.py
def convert_res_ann(res_prob, obs_Xdata, states_weight = np.array([1,2,3])):
    """
    deprecated 2022-06-27
    Function: 
    Convert cell based res to anndata for visualization.
    Check nan value.
    todo: 1) why -1 for visualization (center?)
          2）change legend for copy states
    Parameters:
    ----------
    res: np.array. same order with Xdata.X
    Xdata: anndata. Provide the annotation.
    Return:
    ------
    
    Example:
    -------
    """
    res_cnv = np.argmax(res_prob, axis=-1)
    res_cnv_ad = obs_Xdata.copy()
    # res_cnv_ad.X = res_cnv.copy() - 1
    res_cnv_ad.X = res_cnv.copy()

    res_cnv_weights = (res_prob * states_weight).sum(axis=-1)

    res_cnv_weights_ad = obs_Xdata.copy()
    res_cnv_weights_ad.X = res_cnv_weights.copy()

    res_cnv_weights_ad1 = obs_Xdata.copy()
    # res_cnv_weights_ad1.X = res_cnv_weights.copy() - 1
    res_cnv_weights_ad1.X = res_cnv_weights.copy()

    # return res_cnv_weights_ad, res_cnv_weights_ad1
    return res_cnv_ad, res_cnv_weights_ad1