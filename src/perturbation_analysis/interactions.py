"""Gene interaction inference (Norman 2019 approach)."""
import numpy as np
import pandas as pd


def compute_interaction_score(de_single_a, de_single_b, de_double,
                               score_col='log2FoldChange'):
    """Compare double-perturbation DE to additive expectation.

    Parameters
    ----------
    de_single_a, de_single_b, de_double : pd.DataFrame
        DE results with score_col column and gene index.
    score_col : str

    Returns
    -------
    pd.DataFrame with 'observed', 'expected', 'interaction' columns
    """
    genes = de_single_a.index.intersection(de_single_b.index).intersection(de_double.index)
    result = pd.DataFrame(index=genes)
    result['single_a'] = de_single_a.loc[genes, score_col]
    result['single_b'] = de_single_b.loc[genes, score_col]
    result['expected'] = result['single_a'] + result['single_b']
    result['observed'] = de_double.loc[genes, score_col]
    result['interaction'] = result['observed'] - result['expected']
    return result
