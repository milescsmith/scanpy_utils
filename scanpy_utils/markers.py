from itertools import chain
import pandas as pd
from anndata import AnnData


def get_marker_genes(
    adata: AnnData, rank_key: str = "rank_genes_groups"
) -> pd.DataFrame:

    assert (
        rank_key in adata.uns
    ), "top_markers() requires you to run `sc.tl.rank_genes_groups` first."
    marker_genes_columns = [
        pd.DataFrame.from_records(adata.uns[rank_key][_]).melt()
        for _ in adata.uns[rank_key]
        if _ != "params"
    ]
    marker_genes_df = pd.concat(marker_genes_columns, axis=1)
    marker_genes_df.columns = [
        _
        for _ in chain.from_iterable(
            zip(["variable"] * 5, [_ for _ in bcells.uns[rank_key] if _ != "params"])
        )
    ]
    marker_genes_df = marker_genes_df.loc[
        :, ~marker_genes_df.columns.duplicated()
    ].rename(index=str, columns={"variable": "cluster"})
    return marker_genes_df


def top_markers(
    adata: AnnData, topn: int = 5, rank_key: str = "rank_genes_groups"
) -> pd.DataFrame:

    assert (
        rank_key in adata.uns
    ), "top_markers() requires you to run `sc.tl.rank_genes_groups` first."
    marker_genes_df = get_marker_genes(adata=adata, rank_key=rank_key)
    top_markers = marker_genes_df.iloc[
        marker_genes_df.groupby("cluster").apply(
            lambda x: x["logfoldchanges"].nlargest(topn).reset_index()
        )["index"]
    ]["names"].values
    return top_markers


def top_markers_for_group(
    adata: AnnData,
    group: str,
    topn: int = 5,
    selection_key: str = "logfoldchanges",
    rank_key: str = "rank_genes_groups",
) -> pd.DataFrame:

    assert (
        rank_key in adata.uns
    ), "top_markers() requires you to run `sc.tl.rank_genes_groups` first."
    marker_genes_df = get_marker_genes(adata=adata, rank_key=rank_key)
    marker_genes_df = marker_genes_df[marker_genes_df["cluster"] == str(group)]
    top_markers = marker_genes_df.nlargest(topn, columns=selection_key)
    return top_markers
