import pandas as pd

def cluster_fractions(
    df: pd.DataFrame, numerator: str, denominator: str
) -> pd.DataFrame:

    df["index"] = df.index
    beta = df.groupby([numerator, denominator]).count()
    beta = beta.loc[:, ["index"]].rename(columns={"index": "numerator_count"})
    gamma = df.groupby(denominator).count()
    gamma = gamma.loc[:, ["index"]].rename(columns={"index": "denominator_count"})
    delta = beta.reset_index().merge(gamma.reset_index(), on=denominator)
    delta["fraction"] = delta["numerator_count"] / delta["denominator_count"]
    return delta
