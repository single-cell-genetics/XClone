# init version
import matplotlib.colors as clr
cnt_color = clr.LinearSegmentedColormap.from_list('magma', ["#000003",  "#3B0F6F",  "#8C2980",   "#F66E5B", "#FD9F6C", "#FBFCBF"], N=256)
def simple_spatial_scatter(spatial_loc, data, title, ax):
    """
    

    Args:
        spatial_loc (_type_): _description_
        data (_type_): _description_
        title (_type_): _description_
        ax (_type_): _description_

    Returns:
        _type_: _description_
    # Notes: add functions for Discrete variables
    """
    scatter = ax.scatter(spatial_loc[:,0], spatial_loc[:,1], c=data, cmap=cnt_color)    # OrRd   Greens   YlOrRd  autumn   autumn_r
    ax.invert_yaxis()
    ax.set_title(title)
    return scatter