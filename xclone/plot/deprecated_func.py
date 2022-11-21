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