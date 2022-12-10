"""[XClone] 
Base manhattan plotting functions for XClone BAF analysis.
"""

# Author: Rongting Huang
# rthuang@connect.hku.hk

### update from https://github.com/Rongtingting/limix-plot
### [orogin: https://github.com/limix/limix-plot]
###======================================== manhattan plot for BAF[start]
def manhattan(data, colora="#5689AC", colorb="#21334F", 
              anno_pv_max=None, pts_kws=None, ax=None, figsize = None):
    """
    Produce a manhattan plot.
    Parameters
    ----------
    data : DataFrame, dict
        DataFrame containing the chromosome, base-pair positions, and
        BAF-values.
    colora : matplotlib color
        Points color of the first group.
    colorb : matplotlib color
        Points color of the second group.
    anno_pv_max : float
        Threshold of maximum p value for annotating data['id'] in the plot.
    pts_kws : dict, optional
        Keyword arguments forwarded to the matplotlib function used for
        plotting the points.
    ax : matplotlib Axes, optional
        The target handle for this figure. If ``None``, the current axes is
        set.
    Example
    -------
    .. plot::
        >>> import limix_plot as lp
        >>> from numpy import log10
        >>>
        >>> df = lp.load_dataset('gwas')
        >>> df = df.rename(columns={"chr": "chrom"})
        >>> df = df.rename(columns={"pv": "BAF"})
        >>> print(df.head())
            chrom     pos       BAF
        234    10  224239  0.00887
        239    10  229681  0.00848
        253    10  240788  0.00721
        258    10  246933  0.00568
        266    10  255222  0.00593
        >>> lp.manhattan(df)
        >>> plt = lp.get_pyplot()
        >>> _ = plt.axhline(-log10(1e-7), color='red')
        >>> _ = plt.ylim(2, plt.ylim()[1])
    """
    from numpy import log10, unique, where
    from xarray import DataArray
    import pandas as pd
    
    import matplotlib.pyplot as plt

    # plt = get_pyplot()
    if figsize == None:
        plt.figure(figsize=(24, 6))
    else:
        plt.figure(figsize=figsize)


    if isinstance(data, pd.DataFrame):
        data = DataArray(
            data=data["BAF"],
            dims=["candidate"],
            coords={k: ("candidate", data[k]) for k in data.columns},
        )
    else:
        data = DataArray(data=data)

    if len(data) == 0:
        raise ValueError("DataFrame is empty.")

    if pts_kws is None:
        pts_kws = dict()

    ax = plt.gca() if ax is None else ax

    data["chrom"] = data["chrom"].astype(str)
    data["pos"] = data["pos"].astype(int)
    chr_order = _chr_precedence(data)
    data["order"] = ("candidate", [chr_order[i] for i in data["chrom"].values])

    data = data.sortby(["order", "pos"])

    data = _abs_pos(data)

    if "markersize" not in pts_kws:
        pts_kws["markersize"] = 2
    if "marker" not in pts_kws:
        pts_kws["marker"] = "."
    if "linestyle" not in pts_kws:
        pts_kws["linestyle"] = ""

    colors = {0: colora, 1: colorb}

    for i, c in enumerate(unique(data["order"])):
        ok = data["order"] == c
        pts_kws["color"] = colors[i % 2]
        x = data.loc[ok]["abs_pos"]
        y = data.loc[ok].values
        ax.plot(x, y, **pts_kws)
        
        if anno_pv_max is not None:
            _idx = where(y > anno_pv_max)[0]
            for _ii in _idx:
                if 'id' in data.coords:
                    _txt = data['id'].loc[ok].loc[_ii].values
                else:
                    _txt = (str(data['chrom'].loc[ok].loc[_ii].values) + "_" + 
                            str(data['pos'].loc[ok].loc[_ii].values))
                ax.annotate(_txt, (x[_ii], y[_ii]), ha='center')

    ax.set_xlim(data["abs_pos"].min(), data["abs_pos"].max())
#     ax.set_ylim(0, ax.get_ylim()[1])
    ax.set_ylim(-0.01, 1.01)

    ax.set_ylabel("BAF", fontsize=15)
    ax.set_xlabel("chromosome", fontsize=15)

    u = unique(data["chrom"].values)
    chrom_labels = sorted(u, key=lambda x: chr_order[x])
    _set_ticks(ax, _chrom_bounds(data), chrom_labels)


def _plot_points(ax, data, alpha, null_style, alt_style):
    from numpy import log10

    null = data.loc[data.values >= alpha, :]
    alt = data.loc[data.values < alpha, :]

    ax.plot(null["abs_pos"], null.values, ".", ms=7, **null_style)
    ax.plot(alt["abs_pos"], alt.values, ".", ms=7, **alt_style)


def _set_ticks(ax, chrom_bounds, chrom_labels):
    from numpy import asarray, mean

    n = len(chrom_bounds)
    xticks = asarray([mean(chrom_bounds[i]) for i in range(n)])
    ax.set_xticks(xticks)
    ax.set_xticklabels(chrom_labels, fontsize=15)
    
    ax.tick_params(axis='both', which='major', labelsize=15)
#     ax.tick_params(axis='both', which='minor', labelsize=8)


def _abs_pos(data):
    from numpy import cumsum, flipud, unique

    order = unique(data["order"].values)
    chrom_ends = [data["pos"].values[data["order"].values == c].max() for c in order]

    offset = flipud(cumsum(chrom_ends)[:-1])

    data["abs_pos"] = data["pos"].copy()

    order = list(reversed(order))
    for i, oi in enumerate(offset):
        ix = data["order"] == order[i]
        data["abs_pos"].values[ix] = data.loc[ix]["abs_pos"] + oi

    return data


def _chrom_bounds(data):
    from numpy import unique

    order = unique(data["order"])
    v = []
    for c in order:
        vals = data["abs_pos"][data["order"] == c]
        v += [(vals.min(), vals.max())]
    return v


def _isint(i):
    try:
        int(i)
    except ValueError:
        return False
    else:
        return True


def _chr_precedence(data):
    from numpy import unique

    uchr = unique(data["chrom"].values)
    nchr = [int(i) for i in uchr if _isint(i)]
    if len(nchr) > 0:
        offset = max(nchr)
    else:
        offset = -1
    precedence = {str(i): i for i in nchr}
    schr = sorted([i for i in uchr if not _isint(i)])
    for i, s in enumerate(schr):
        precedence[s] = offset + i + 1
    return precedence

###======================================== manhattan plot for BAF[end]


def get_pyplot():
    """
    Get :mod:`matplotlib.pyplot`.
    Returns
    -------
    pyplot : :mod:`matplotlib.pyplot`
        MATLAB-like interface.
    """
    from sys import platform as sys_pf
    from matplotlib import rcParams

    if get_pyplot.pyplot is not None:
        return get_pyplot.pyplot

    # Bug with matplotlib on macos:
    # https://github.com/matplotlib/matplotlib/issues/10239
    if "backend" not in rcParams and sys_pf == "darwin":
        from matplotlib import use as _backend_use

        _backend_use("TkAgg")

    from matplotlib import pyplot

    get_pyplot.pyplot = pyplot
    return pyplot


get_pyplot.pyplot = None

###======================================== function from limix-plot for manhattan