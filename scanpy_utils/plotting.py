from typing import Optional, Dict, Union, List, Tuple

from matplotlib.axes import Axes
from matplotlib import gridspec
import matplotlib.lines as mlines
import matplotlib.transforms as mtransforms
import matplotlib.pyplot as plt

from anndata import AnnData

from scanpy.plotting._tools.scatterplots import plot_scatter, _wraps_plot_scatter
from scanpy import logging as logg
from scanpy.plotting._utils import savefig_or_show
from scanpy.utils import doc_params
from scanpy.plotting._docs import doc_show_save_ax, doc_common_plot_args
from scanpy.plotting._anndata import _check_var_names_type, _prepare_dataframe

from scvelo.plotting.utils import (
    is_categorical,
    make_dense,
    interpret_colorkey,
    default_basis,
    default_arrow,
    default_size,
    get_components,
    set_colorbar,
)
from scvelo.plotting.velocity_embedding_grid import compute_velocity_on_grid

import numpy as np
import pandas as pd


@_wraps_plot_scatter
def dimplot(adata, reduction="umap", **kwargs) -> Union[Axes, List[Axes], None]:
    """\
    A generalized interface for plotting dimensional reductions.
    Parameters
    ----------
    {adata_color_etc}
    {edges_arrows}
    {scatter_bulk}
    {show_save_ax}
    Returns
    -------
    If `show==False` a :class:`~matplotlib.axes.Axes` or a list of it.
    """
    return plot_scatter(adata, reduction, **kwargs)


def velocity_compare(
    adatas, feature_list: List[str], xkey: str = "Ms", ykey: str = "Mu"
):

    # TODO: add code here to generate the color list!

    if not isinstance(adatas, list):
        adatas = [adatas]
    for _ in feature_list:
        fig, axs = plt.subplots(
            1, len(adatas), sharey=True, sharex=True, figsize=(12, 4)
        )
        basis = _

        for index, item in enumerate(adatas):
            x = make_dense(adatas[item][:, basis].layers[xkey]).flatten()
            y = make_dense(adatas[item][:, basis].layers[ykey]).flatten()

            axs[index].scatter(x, y, alpha=0.5, color=adatas[item], s=1)

            xnew = np.linspace(
                0, np.percentile(make_dense(adatas[item][:, basis].layers[xkey]), 98)
            )
            gamma = adatas[item][:, basis].var["velocity_gamma"].values
            beta = 1
            offset = 0
            line = mlines.Line2D(
                xnew, gamma / beta * xnew + offset / beta, color="black"
            )

            axs[index].set_title(f"{item} {basis}")
            axs[index].set_xlabel("spliced")
            axs[index].set_ylabel("unspliced")
            axs[index].add_line(line)
        plt.show()


def plot_grouped_velo_grid(
    adata: AnnData,
    group: str,
    basis: Optional[str] = None,
    color: str = "YES1",
    layer: str = "velocity",
    vkey: str = "velocity",
    gridlines: bool = False,
    show_quiver: bool = True,
    arrow_length: Optional[float] = None,
    arrow_color: Optional[str] = None,
    arrow_size: Optional[int] = None,
    colormap: str = None,
    components: Optional[str] = None,
    size: Optional[int] = None,
    show: bool = True,
) -> plt.figure:

    basis = default_basis(adata) if basis is None else basis
    colormap = "viridis" if colormap is None else colormap
    unique_group_vals = np.unique(adata.obs[group])

    obj_dict = {val: adata[adata.obs[group] == val] for val in unique_group_vals}

    scale = 1 / arrow_length if arrow_length is not None else 1

    hl, hw, hal = default_arrow(arrow_size)
    quiver_kwargs = {
        "scale": scale,
        "angles": "xy",
        "scale_units": "xy",
        "width": 0.001,
        "color": "grey" if arrow_color is None else arrow_color,
        "edgecolors": "k",
        "headlength": hl / 2,
        "headwidth": hw / 2,
        "headaxislength": hal / 2,
        "linewidth": 0.2,
    }

    quiver_dict = {
        item: {
            k: v
            for k, v in zip(
                ["XY_grid", "V_grid"],
                compute_velocity_on_grid(
                    X_emb=np.array(
                        obj_dict[item].obsm[f"X_{basis}"][
                            :, get_components(components, basis)
                        ]
                    ),
                    V_emb=np.array(
                        obj_dict[item].obsm[f"{vkey}_{basis}"][
                            :, get_components(components, basis)
                        ]
                    ),
                    autoscale=True,
                ),
            )
        }
        for item in obj_dict
    }

    size = 4 * default_size(adata) if size is None else size
    fig, axs = plt.subplots(1, len(obj_dict), figsize=(12, 4))

    for index, item in enumerate(obj_dict):
        embed_points = pd.DataFrame(obj_dict[item].obsm[f"X_{basis}"]).rename(
            columns={0: f"{basis}_1", 1: f"{basis}_2"}
        )
        exprs = interpret_colorkey(obj_dict[item], color, layer, 0.95)

        axs[index].scatter(
            x=embed_points[f"{basis}_1"],
            y=embed_points[f"{basis}_2"],
            c=exprs,
            s=1,
            cmap=colormap,
        )

        if show_quiver:
            axs[index].quiver(
                quiver_dict[item]["XY_grid"][:, 0],
                quiver_dict[item]["XY_grid"][:, 1],
                quiver_dict[item]["V_grid"][:, 0],
                quiver_dict[item]["V_grid"][:, 1],
                zorder=3,
                **quiver_kwargs,
            )
        axs[index].set_title(f"{item}: {color}")
        if not gridlines:
            axs[index].grid()
    if show:
        plt.show()
    return fig


def plot_dictionary_velo_grid(
    adata_dict: Dict[str, AnnData],
    basis: Optional[str] = None,
    color: Optional[str] = None,
    layer: str = "velocity",
    vkey: str = "velocity",
    arrow_length: float = None,
    arrow_color: str = None,
    figsize: Tuple[int] = (12, 4),
    gridlines: bool = False,
    show_quiver: bool = True,
    colormap: Union[str, List[str], Dict[str, str]] = None,
    arrow_size: Optional[int] = None,
    components: Optional[str] = None,
    size: Optional[int] = None,
    show: bool = True,
) -> plt.figure:

    if colormap is None:
        colormap = {_: "viridis" for _ in adata_dict}
    elif isinstance(colormap, str):
        colormap = {_: colormap for _ in adata_dict}
    elif isinstance(colormap, list):
        colormap = {item: colormap[index] for index, item in enumerate(adata_dict)}

    scale = 1 / arrow_length if arrow_length is not None else 1

    hl, hw, hal = default_arrow(arrow_size)
    quiver_kwargs = {
        "scale": scale,
        "angles": "xy",
        "scale_units": "xy",
        "width": 0.001,
        "color": "grey" if arrow_color is None else arrow_color,
        "edgecolors": "k",
        "headlength": hl / 2,
        "headwidth": hw / 2,
        "headaxislength": hal / 2,
        "linewidth": 0.2,
    }

    fig = plt.figure(figsize=figsize)
    for index, item in enumerate(adata_dict):
        basis = default_basis(adata_dict[item]) if basis is None else basis
        quiver_data = {
            k: v
            for k, v in zip(
                ["XY_grid", "V_grid"],
                compute_velocity_on_grid(
                    X_emb=np.array(
                        adata_dict[item].obsm[f"X_{basis}"][
                            :, get_components(components, basis)
                        ]
                    ),
                    V_emb=np.array(
                        adata_dict[item].obsm[f"{vkey}_{basis}"][
                            :, get_components(components, basis)
                        ]
                    ),
                    autoscale=True,
                ),
            )
        }

        size = 4 * default_size(adata_dict[item]) if size is None else size

        embed_points = pd.DataFrame(adata_dict[item].obsm[f"X_{basis}"]).rename(
            columns={0: f"{basis}_1", 1: f"{basis}_2"}
        )
        
        if color is not None:
            exprs = interpret_colorkey(adata_dict[item], color, layer, 0.95)
        else:
            exprs = np.zeros(shape=(1,adata.shape[1]))

        ax = fig.add_subplot(1, len(adata_dict), index + 1)
        ax.scatter(
            x=embed_points[f"{basis}_1"],
            y=embed_points[f"{basis}_2"],
            c=exprs,
            s=1,
            cmap=colormap[item],
        )

        if show_quiver:
            ax.quiver(
                quiver_data["XY_grid"][:, 0],
                quiver_data["XY_grid"][:, 1],
                quiver_data["V_grid"][:, 0],
                quiver_data["V_grid"][:, 1],
                zorder=3,
                **quiver_kwargs,
            )
        ax.set_title(f"{item}: {color}")
        if not gridlines:
            ax.grid()
    if show:
        plt.show()
    return fig


def dotplot(adata,
            var_names,
            groupby=None,
            use_raw=None,
            log=False,
            num_categories=7,
            expression_cutoff=0.,
            mean_only_expressed=False,
            color_map='Reds',
            dot_max=None,
            dot_min=None,
            show_grid=False,
            grid_color='#CCCCCC',
            grid_linewidth=1,
            grid_linestyle="-",
            grid_alpha=0.25,
            figsize=None,
            dendrogram=False,
            gene_symbols=None,
            var_group_positions=None,
            standard_scale=None,
            smallest_dot=0.,
            var_group_labels=None,
            var_group_rotation=None,
            layer=None,
            show=None,
            save=None,
            **kwds):
    """\
    Makes a *dot plot* of the expression values of `var_names`.
    For each var_name and each `groupby` category a dot is plotted. Each dot
    represents two values: mean expression within each category (visualized by
    color) and fraction of cells expressing the var_name in the
    category (visualized by the size of the dot).  If groupby is not given, the
    dotplot assumes that all data belongs to a single category.
    **Note**: A gene is considered expressed if the expression value in the adata
    (or adata.raw) is above the specified threshold which is zero by default.
    An example of dotplot usage is to visualize, for multiple marker genes,
    the mean value and the percentage of cells expressing the gene accross multiple clusters.
    See also :func:`~scanpy.pl.rank_genes_groups_dotplot` to plot marker genes
    identified using the :func:`~scanpy.tl.rank_genes_groups` function.
    Parameters
    ----------
    {common_plot_args}
    expression_cutoff : `float` (default: `0.`)
        Expression cutoff that is used for binarizing the gene expression and determining the fraction
        of cells expressing given genes. A gene is expressed only if the expression value is greater than
        this threshold.
    mean_only_expressed : `bool` (default: `False`)
        If True, gene expression is averaged only over the cells expressing the given genes.
    color_map : `str`, optional (default: `Reds`)
        String denoting matplotlib color map.
    dot_max : `float` optional (default: `None`)
        If none, the maximum dot size is set to the maximum fraction value found (e.g. 0.6). If given,
        the value should be a number between 0 and 1. All fractions larger than dot_max are clipped to
        this value.
    dot_min : `float` optional (default: `None`)
        If none, the minimum dot size is set to 0. If given,
        the value should be a number between 0 and 1. All fractions smaller than dot_min are clipped to
        this value.
    standard_scale : {{'var', 'group'}}, optional (default: None)
        Whether or not to standardize that dimension between 0 and 1, meaning for each variable or group,
        subtract the minimum and divide each by its maximum.
    smallest_dot : `float` optional (default: 0.)
        If none, the smallest dot has size 0. All expression levels with `dot_min` are potted with
        `smallest_dot` dot size.
    {show_save_ax}
    **kwds : keyword arguments
        Are passed to `matplotlib.pyplot.scatter`.
    Returns
    -------
    List of :class:`~matplotlib.axes.Axes`
    Examples
    -------
    >>> adata = sc.datasets.pbmc68k_reduced()
    >>> sc.pl.dotplot(adata, ['C1QA', 'PSAP', 'CD79A', 'CD79B', 'CST3', 'LYZ'],
    ...               groupby='bulk_labels', dendrogram=True)
    Using var_names as dict:
    >>> markers = {{'T-cell': 'CD3D', 'B-cell': 'CD79A', 'myeloid': 'CST3'}}
    >>> sc.pl.dotplot(adata, markers, groupby='bulk_labels', dendrogram=True)
    """
    if use_raw is None and adata.raw is not None: use_raw = True
    var_names, var_group_labels, var_group_positions = _check_var_names_type(var_names,
                                                                             var_group_labels, var_group_positions)
    categories, obs_tidy = _prepare_dataframe(adata, var_names, groupby, use_raw, log, num_categories,
                                              layer=layer, gene_symbols=gene_symbols)

    # for if category defined by groupby (if any) compute for each var_name
    # 1. the fraction of cells in the category having a value > expression_cutoff
    # 2. the mean value over the category

    # 1. compute fraction of cells having value > expression_cutoff
    # transform obs_tidy into boolean matrix using the expression_cutoff
    obs_bool = obs_tidy > expression_cutoff

    # compute the sum per group which in the boolean matrix this is the number
    # of values > expression_cutoff, and divide the result by the total number of values
    # in the group (given by `count()`)
    fraction_obs = obs_bool.groupby(level=0).sum() / obs_bool.groupby(level=0).count()

    # 2. compute mean value
    if mean_only_expressed:
        mean_obs = obs_tidy.mask(~obs_bool).groupby(level=0).mean().fillna(0)
    else:
        mean_obs = obs_tidy.groupby(level=0).mean()

    if standard_scale == 'group':
        mean_obs = mean_obs.sub(mean_obs.min(1), axis=0)
        mean_obs = mean_obs.div(mean_obs.max(1), axis=0).fillna(0)
    elif standard_scale == 'var':
        mean_obs -= mean_obs.min(0)
        mean_obs = (mean_obs / mean_obs.max(0)).fillna(0)
    elif standard_scale is None:
        pass
    else:
        logg.warning('Unknown type for standard_scale, ignored')

    dendro_width = 0.8 if dendrogram else 0
    colorbar_width = 0.2
    colorbar_width_spacer = 0.5
    size_legend_width = 0.25
    if figsize is None:
        height = len(categories) * 0.3 + 1  # +1 for labels
        # if the number of categories is small (eg 1 or 2) use
        # a larger height
        height = max([1.5, height])
        heatmap_width = len(var_names) * 0.35
        width = heatmap_width + colorbar_width + size_legend_width + dendro_width + colorbar_width_spacer
    else:
        width, height = figsize
        heatmap_width = width - (colorbar_width + size_legend_width + dendro_width + colorbar_width_spacer)

    # colorbar ax width should not change with differences in the width of the image
    # otherwise can become too small

    if var_group_positions is not None and len(var_group_positions) > 0:
        # add some space in case 'brackets' want to be plotted on top of the image
        height_ratios = [0.5, 10]
    else:
        height_ratios = [0, 10.5]

    # define a layout of 2 rows x 5 columns
    # first row is for 'brackets' (if no brackets needed, the height of this row is zero)
    # second row is for main content. This second row
    # is divided into 4 axes:
    #   first ax is for the main figure
    #   second ax is for dendrogram (if present)
    #   third ax is for the color bar legend
    #   fourth ax is for an spacer that avoids the ticks
    #             from the color bar to be hidden beneath the size lengend axis
    #   fifth ax is to plot the size legend
    fig = plt.figure(figsize=(width, height))
    axs = gridspec.GridSpec(nrows=2, ncols=5, wspace=0.02, hspace=0.04,
                            width_ratios=[heatmap_width, dendro_width, colorbar_width, colorbar_width_spacer, size_legend_width],
                            height_ratios=height_ratios)
    if len(categories) < 4:
        # when few categories are shown, the colorbar and size legend
        # need to be larger than the main plot, otherwise they would look
        # compressed. For this, the dotplot ax is split into two:
        axs2 = gridspec.GridSpecFromSubplotSpec(2, 1, subplot_spec=axs[1, 0],
                                                height_ratios=[len(categories) * 0.3, 1])
        dot_ax = fig.add_subplot(axs2[0])
    else:
        dot_ax = fig.add_subplot(axs[1, 0])

    color_legend = fig.add_subplot(axs[1, 2])

    if groupby is None or len(categories) <= 1:
        # dendrogram can only be computed  between groupby categories
        dendrogram = False

    if dendrogram:
        dendro_data = _reorder_categories_after_dendrogram(adata, groupby, dendrogram,
                                                           var_names=var_names,
                                                           var_group_labels=var_group_labels,
                                                           var_group_positions=var_group_positions)

        var_group_labels = dendro_data['var_group_labels']
        var_group_positions = dendro_data['var_group_positions']

        # reorder matrix
        if dendro_data['var_names_idx_ordered'] is not None:
            # reorder columns (usually genes) if needed. This only happens when
            # var_group_positions and var_group_labels is set
            mean_obs = mean_obs.iloc[:,dendro_data['var_names_idx_ordered']]
            fraction_obs = fraction_obs.iloc[:, dendro_data['var_names_idx_ordered']]

        # reorder rows (categories) to match the dendrogram order
        mean_obs = mean_obs.iloc[dendro_data['categories_idx_ordered'], :]
        fraction_obs = fraction_obs.iloc[dendro_data['categories_idx_ordered'], :]

        y_ticks = range(mean_obs.shape[0])
        dendro_ax = fig.add_subplot(axs[1, 1], sharey=dot_ax)
        _plot_dendrogram(dendro_ax, adata, groupby, dendrogram_key=dendrogram, ticks=y_ticks)

    # to keep the size_legen of about the same height, irrespective
    # of the number of categories, the fourth ax is subdivided into two parts
    size_legend_height = min(1.3, height)
    # wspace is proportional to the width but a constant value is
    # needed such that the spacing is the same for thinner or wider images.
    wspace = 10.5 / width
    axs3 = gridspec.GridSpecFromSubplotSpec(2, 1, subplot_spec=axs[1, 4], wspace=wspace,
                                            height_ratios=[size_legend_height / height,
                                                           (height - size_legend_height) / height])
    # make scatter plot in which
    # x = var_names
    # y = groupby category
    # size = fraction
    # color = mean expression

    y, x = np.indices(mean_obs.shape)
    y = y.flatten()
    x = x.flatten()
    frac = fraction_obs.values.flatten()
    mean_flat = mean_obs.values.flatten()
    cmap = plt.get_cmap(color_map)
    if dot_max is None:
        dot_max = np.ceil(max(frac) * 10) / 10
    else:
        if dot_max < 0 or dot_max > 1:
            raise ValueError("`dot_max` value has to be between 0 and 1")
    if dot_min is None:
        dot_min = 0
    else:
        if dot_min < 0 or dot_min > 1:
            raise ValueError("`dot_min` value has to be between 0 and 1")

    if dot_min != 0 or dot_max != 1:
        # clip frac between dot_min and  dot_max
        frac = np.clip(frac, dot_min, dot_max)
        old_range = dot_max - dot_min
        # re-scale frac between 0 and 1
        frac = ((frac - dot_min) / old_range)

    size = (frac * 10) ** 2
    size += smallest_dot
    import matplotlib.colors

    normalize = matplotlib.colors.Normalize(vmin=kwds.get('vmin'), vmax=kwds.get('vmax'))
    colors = cmap(normalize(mean_flat))
    dot_ax.scatter(x, y, color=colors, s=size, cmap=cmap, norm=None, edgecolor='none', **kwds)
    y_ticks = range(mean_obs.shape[0])
    dot_ax.set_yticks(y_ticks)
    dot_ax.set_yticklabels([mean_obs.index[idx] for idx in y_ticks])

    x_ticks = range(mean_obs.shape[1])
    dot_ax.set_xticks(x_ticks)
    dot_ax.set_xticklabels([mean_obs.columns[idx] for idx in x_ticks], rotation=90)
    dot_ax.tick_params(axis='both', labelsize='small')
    if show_grid:
        dot_ax.grid(show_grid, which="both", axis="both", color=grid_color, linestyle=grid_linestyle, linewidth=grid_linewidth, alpha=grid_alpha)
    else:
        dot_ax.grid(False)
    dot_ax.set_xlim(-0.5, len(var_names) + 0.5)
    dot_ax.set_ylabel(groupby)

    # to be consistent with the heatmap plot, is better to
    # invert the order of the y-axis, such that the first group is on
    # top
    ymin, ymax = dot_ax.get_ylim()
    dot_ax.set_ylim(ymax+0.5, ymin - 0.5)

    dot_ax.set_xlim(-1, len(var_names))

    # plot group legends on top of dot_ax (if given)
    if var_group_positions is not None and len(var_group_positions) > 0:
        gene_groups_ax = fig.add_subplot(axs[0, 0], sharex=dot_ax)
        _plot_gene_groups_brackets(gene_groups_ax, group_positions=var_group_positions,
                                   group_labels=var_group_labels,
                                   rotation=var_group_rotation)

    # plot colorbar
    import matplotlib.colorbar
    matplotlib.colorbar.ColorbarBase(color_legend, cmap=cmap, norm=normalize)

    # for the dot size legend, use step between dot_max and dot_min
    # based on how different they are.
    diff = dot_max - dot_min
    if 0.3 < diff <= 0.6:
        step = 0.1
    elif diff <= 0.3:
        step = 0.05
    else:
        step = 0.2
    # a descending range that is afterwards inverted is used
    # to guarantee that dot_max is in the legend.
    fracs_legends = np.arange(dot_max, dot_min, step * -1)[::-1]
    if dot_min != 0 or dot_max != 1:
        fracs_values = ((fracs_legends - dot_min) / old_range)
    else:
        fracs_values = fracs_legends
    size = (fracs_values * 10) ** 2
    size += smallest_dot
    color = [cmap(normalize(value)) for value in np.repeat(max(mean_flat) * 0.7, len(size))]

    # plot size bar
    size_legend = fig.add_subplot(axs3[0])

    size_legend.scatter(np.repeat(0, len(size)), range(len(size)), s=size, color=color)
    size_legend.set_yticks(range(len(size)))
    labels = ["{:.0%}".format(x) for x in fracs_legends]
    if dot_max < 1:
        labels[-1] = ">" + labels[-1]
    size_legend.set_yticklabels(labels)
    size_legend.set_yticklabels(["{:.0%}".format(x) for x in fracs_legends])

    size_legend.tick_params(axis='y', left=False, labelleft=False, labelright=True)

    # remove x ticks and labels
    size_legend.tick_params(axis='x', bottom=False, labelbottom=False)

    # remove surrounding lines
    size_legend.spines['right'].set_visible(False)
    size_legend.spines['top'].set_visible(False)
    size_legend.spines['left'].set_visible(False)
    size_legend.spines['bottom'].set_visible(False)
    #size_legend.grid(False)

    ymin, ymax = size_legend.get_ylim()
    size_legend.set_ylim(ymin, ymax+0.5)

    savefig_or_show('dotplot', show=show, save=save)
    return axs