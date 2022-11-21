"""Set the default matplotlib.rcParams for XClone.
"""

import matplotlib
from matplotlib import rcParams
# from cycler import cycler


def set_style_xclone(fontsize=14, color_map=None, grid_alpha = None):
    """Set matplotlib.rcParams to XClone defaults.
    Call this through `settings.set_figure_params`.
    """

    # figure
    rcParams['figure.figsize'] = (12, 8)
    rcParams['figure.subplot.left'] = 0.18
    rcParams['figure.subplot.right'] = 0.96
    rcParams['figure.subplot.bottom'] = 0.15
    rcParams['figure.subplot.top'] = 0.91

    rcParams['lines.linewidth'] = 1.5  # the line width of the frame
    rcParams['lines.markersize'] = 6
    rcParams['lines.markeredgewidth'] = 1

    # font
    rcParams['font.sans-serif'] = [
        'Arial',
        'Helvetica',
        'DejaVu Sans',
        'Bitstream Vera Sans',
        'sans-serif',
    ]
    fontsize = fontsize
    rcParams['font.size'] = fontsize
    rcParams['legend.fontsize'] = 1.1 * fontsize
    rcParams['axes.labelsize'] = 1.1* fontsize
    rcParams['axes.titlesize'] = 1.2 * fontsize
    rcParams['axes.titleweight'] = 'bold'


    # legend
    rcParams['legend.numpoints'] = 1
    rcParams['legend.scatterpoints'] = 1
    rcParams['legend.handlelength'] = 0.5
    rcParams['legend.handletextpad'] = 0.4

    # lines
    rcParams['axes.linewidth'] = 0.8
    rcParams['axes.edgecolor'] = 'black'
    rcParams['axes.facecolor'] = 'white'

    # ticks
    rcParams['xtick.color'] = 'k'
    rcParams['ytick.color'] = 'k'
    rcParams['xtick.labelsize'] = fontsize
    rcParams['ytick.labelsize'] = fontsize

    # axes grid
    rcParams['axes.grid'] = False
    if grid_alpha is not None and grid_alpha >= 0 and grid_alpha <= 1:
        rcParams.update({"axes.grid" : True, "grid.alpha": grid_alpha})

    # color map
    rcParams['image.cmap'] = rcParams['image.cmap'] if color_map is None else color_map


def set_style_defaults():
    """Reset `matplotlib.rcParams` to defaults."""
    rcParams.update(matplotlib.rcParamsDefault)