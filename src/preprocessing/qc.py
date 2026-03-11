"""Cell and gene QC filters for Perturb-seq data."""
import scanpy as sc
import numpy as np

# Default QC thresholds (override per dataset as needed)
MIN_CELLS_PER_PERT = 20
MIN_UMIS_PER_CELL = 1000
MIN_GENES_PER_CELL = 500
MAX_PCT_MITO = 20
CONTROL_FRAC_MIN = 0.10
CONTROL_FRAC_MAX = 0.30

CONTROL_LABELS = ['control', 'non-targeting', 'safe-harbor', 'CTRL', 'NT']


def run_qc(adata, min_umis=MIN_UMIS_PER_CELL, min_genes=MIN_GENES_PER_CELL,
           max_pct_mito=MAX_PCT_MITO):
    """Apply standard QC filters to an AnnData object.

    Parameters
    ----------
    adata : AnnData
    min_umis : int
    min_genes : int
    max_pct_mito : float

    Returns
    -------
    AnnData filtered copy
    """
    sc.pp.calculate_qc_metrics(adata, percent_top=None, log1p=False, inplace=True)
    adata.var['mt'] = adata.var_names.str.startswith('MT-')
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

    mask = (
        (adata.obs['total_counts'] >= min_umis) &
        (adata.obs['n_genes_by_counts'] >= min_genes) &
        (adata.obs['pct_counts_mt'] <= max_pct_mito)
    )
    return adata[mask].copy()


def check_control_fraction(adata, control_labels=None, obs_key='perturbation'):
    """Warn if control cell fraction is outside expected range."""
    if control_labels is None:
        control_labels = CONTROL_LABELS
    ctrl_mask = adata.obs[obs_key].isin(control_labels)
    frac = ctrl_mask.sum() / len(adata)
    if frac < CONTROL_FRAC_MIN or frac > CONTROL_FRAC_MAX:
        import warnings
        warnings.warn(f"Control fraction {frac:.2%} outside expected range "
                      f"[{CONTROL_FRAC_MIN:.0%}, {CONTROL_FRAC_MAX:.0%}]")
    return frac
