"""XClone plotting.
"""

## develop
from ._test import random_assign_celltype

## data
from ._data import reorder_data_by_cellanno

## base plotting
from .base_plot import annoHeatMap
from .base_plot import scatter_anno
from .base_plot import corr_plot
from .base_plot import confuse_heatmap

from ._base_xanndata import Xheatmap, Xheatmap_addref
from ._base_xanndata import XXheatmap, XXheatmap_addref

## RDR plotting
from ._RDR_plot import X_hist_orig, X_hist_base, X_hist_multiplots, X_hist_multiplots_CHR


## BAF plotting

from ._BAF_plot import prepare_manhattan_df, visualize_bulk_BAF
from ._BAF_plot import compare_visualize_bulk_BAF

from ._BAF_plot import calculate_cell_BAF, visualize_cell_BAF

### debug
from ._BAF_debug import state_specific_df, transfer_prob_to_df, visualize_emm_prob
from ._BAF_debug import visualize_BAF_distribution

from ._BAF_debug import get_count_df, visualize_DP
from ._BAF_debug import bulk_AF_chr

from ._explore import visualize_prob, compare_prob
from ._explore import visualize_count

## smoothing plotting
from .smooth_plot import raw_ratio_visualization
from .smooth_plot import smooth_visualization, RDR_smooth_visualization, BAF_smooth_visualization, BAF_smooth_complex_visualization
from .smooth_plot import smooth_visualization2, smooth_visualization3

## CNV plotting
### update
from .CNV_plot import CNV_visualization, Complex_CNV_visualization
from .CNV_plot import BAF_CNV_visualization, Complex_BAF_CNV_visualization
from .CNV_plot import Combine_CNV_visualization, Complex_Combine_CNV_visualization

# from .CNV_plot import CNV_visualization, CNV_visualization_complex
# from .CNV_plot import CNV_visualization2, CNV_visualization2_complex
# from .CNV_plot import CNV_visualization_combine, complex_CNV_visualization_combine

# from .CNV_plot import CNV_combine_visualization, CNV_combine_visualization_complex
# from .CNV_plot import CNV_LOH_visualization
