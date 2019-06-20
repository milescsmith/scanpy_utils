from typing import Union, List
from matplotlib.axes import Axes
from scanpy.plotting._tools.scatterplots import plot_scatter


def _wraps_plot_scatter(wrapper):
    annots_orig = {
        k: v for k, v in wrapper.__annotations__.items() if k not in {"adata", "kwargs"}
    }
    annots_scatter = {
        k: v for k, v in plot_scatter.__annotations__.items() if k != "basis"
    }
    wrapper.__annotations__ = {**annots_scatter, **annots_orig}
    wrapper.__wrapped__ = plot_scatter
    return wrapper


@_wraps_plot_scatter
def dimplot(adata, reduction="umap", **kwargs) -> Union[Axes, List[Axes], None]:
    """\
    Scatter plot in UMAP basis.
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
