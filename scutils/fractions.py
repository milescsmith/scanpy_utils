from typing import Optional
import pandas as pd
from anndata import AnnData

def cluster_fractions(
    data_obj: AnnData, numerator: str, denominator: str
) -> pd.DataFrame:
    """Calculate the relative fraction of a population within a given set of groups.
    More clearly, calculate how many cells of, say, a treatment group make up a 
    cluster.

    Parameters
    ----------
    data_obj : :class:`~anndata.AnnData
        anndata object with metadata corresponding to the numerator and denominator arguments
    numerator : `str`
        Inner group to calculate relative fraction of
    denominator : `str`
        Outer group to calculate relative fraction of

    Return
    ------
    :class:`~pandas.DataFrame`
    """

    if isinstance(data_obj, "AnnData"):
        df = data_obj.obs
    df["index"] = df.index
    beta = df.groupby([numerator, denominator]).count()
    beta = beta.loc[:, ["index"]].rename(columns={"index": "numerator_count"})
    gamma = df.groupby(denominator).count()
    gamma = gamma.loc[:, ["index"]].rename(columns={"index": "denominator_count"})
    delta = beta.reset_index().merge(gamma.reset_index(), on=denominator)
    delta["fraction"] = delta["numerator_count"] / delta["denominator_count"]
    return delta
