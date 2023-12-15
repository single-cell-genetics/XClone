"""Base functions for XClone probability plotting
"""

# Author: Rongting Huang
# Date: 2023/12/13
# update: 2023/12/13


from ._data import reorder_data_by_cellanno
from ._base_xanndata import Xheatmap, XXheatmap

def prob_visualization(Xdata, cell_anno_key = "cell_type", 
                         Xlayer = "Prob", vmin=0, vmax=1,
                         colorbar_name = "LOH probability", 
                         **kwargs):
    """
    RENAME TO `RDR_smooth_visualization`.
    """
    re_Xdata = reorder_data_by_cellanno(Xdata, cell_anno_key = cell_anno_key)
    Xheatmap(re_Xdata, Xlayer = Xlayer, cell_anno_key = cell_anno_key, 
            vmin=vmin, vmax=vmax,
            colorbar_ticks = [vmin, 0.5, vmax], 
            colorbar_label = ["0",  "0.5",  "1"], 
            colorbar_name = colorbar_name, **kwargs)

