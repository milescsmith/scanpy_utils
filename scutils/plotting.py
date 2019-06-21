from typing import Union, List
from matplotlib.axes import Axes
from scanpy.plotting._tools.scatterplots import plot_scatter, _wraps_plot_scatter

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
