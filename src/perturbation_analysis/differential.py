"""Pseudobulk differential expression wrappers (PyDESeq2)."""
import numpy as np
import pandas as pd
import scanpy as sc
from src.preprocessing.qc import CONTROL_LABELS


def pseudobulk_de(adata, perturbation, control_labels=None, obs_key='perturbation',
                  sample_col=None):
    """Run pseudobulk DE for a single perturbation vs. controls using PyDESeq2.

    Parameters
    ----------
    adata : AnnData  (raw counts required)
    perturbation : str
    control_labels : list of str, optional
    obs_key : str
    sample_col : str, optional  — column for donor/replicate aggregation

    Returns
    -------
    pd.DataFrame with DE results
    """
    from pydeseq2.dds import DeseqDataSet
    from pydeseq2.ds import DeseqStats

    if control_labels is None:
        control_labels = CONTROL_LABELS

    mask = adata.obs[obs_key].isin([perturbation] + control_labels)
    sub = adata[mask].copy()
    sub.obs['group'] = np.where(sub.obs[obs_key] == perturbation, 'perturbed', 'control')

    # Pseudobulk: if sample_col provided, aggregate per sample; else treat each cell as replicate
    counts = pd.DataFrame(
        sub.X.toarray() if hasattr(sub.X, 'toarray') else sub.X,
        index=sub.obs_names,
        columns=sub.var_names,
    ).T

    metadata = sub.obs[['group']].copy()
    metadata.index = sub.obs_names

    dds = DeseqDataSet(counts=counts, metadata=metadata, design_factors='group')
    dds.deseq2()
    stats = DeseqStats(dds, contrast=['group', 'perturbed', 'control'])
    stats.summary()
    return stats.results_df
