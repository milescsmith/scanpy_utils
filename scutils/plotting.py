from typing import Union, List
from matplotlib.axes import Axes
from scanpy.plotting._tools.scatterplots import plot_scatter, _wraps_plot_scatter
import matplotlib.lines as mlines
import matplotlib.transforms as mtransforms
from anndata import AnnData


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
    adatas: Union(AnnData, List[AnnData]),
    feature_list: List[str],
    xkey: str = "Ms",
    ykey: str = "Mu",
):

    # TODO: add code here to generate the color list!
    
    if not isinstance(adatas, list):
        adatas = [adatas]
    for _ in feature_list:
        fig, axs = plt.subplots(
            1, len(adatas), sharey=True, sharex=True, figsize=(12, 4)
        )
        basis = _

        for index, item in enumerate(class_dict):
            x = make_dense(class_dict[item][:, basis].layers[xkey]).flatten()
            y = make_dense(class_dict[item][:, basis].layers[ykey]).flatten()

            axs[index].scatter(x, y, alpha=0.5, color=class_colors[item], s=1)

            xnew = np.linspace(
                0,
                np.percentile(make_dense(class_dict[item][:, basis].layers[xkey]), 98),
            )
            gamma = class_dict[item][:, basis].var["velocity_gamma"].values
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
