from typing import Optional, Dict, Union, List, Tuple

from matplotlib.axes import Axes
import matplotlib.lines as mlines
import matplotlib.transforms as mtransforms
import matplotlib.pyplot as plt

from anndata import AnnData

from scanpy.plotting._tools.scatterplots import plot_scatter, _wraps_plot_scatter
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
        fig.show()
    return fig


def plot_dictionary_velo_grid(
    adata_dict: Dict[str, AnnData],
    basis: Optional[str] = None,
    color: str = "YES1",
    layer: str = "velocity",
    vkey: str = "velocity",
    arrow_length: float = None,
    arrow_color: str = None,
    figsize: Tuple[int] = (12, 4),
    gridlines: bool = False,
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
        exprs = interpret_colorkey(adata_dict[item], color, layer, 0.95)

        ax = fig.add_subplot(1, len(adata_dict), index + 1)
        ax.scatter(
            x=embed_points[f"{basis}_1"],
            y=embed_points[f"{basis}_2"],
            c=exprs,
            s=1,
            cmap=colormap[item],
        )

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
        fig.show()
    return fig
