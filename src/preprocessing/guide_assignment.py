"""Guide RNA assignment utilities for Perturb-seq data."""


def assign_guides_from_obs(adata, guide_col='guide_id', sep='_'):
    """Extract perturbation gene name from guide ID column.

    Parameters
    ----------
    adata : AnnData
    guide_col : str
        Column in adata.obs containing the guide RNA identifier.
    sep : str
        Separator between gene name and guide number in guide_id.

    Returns
    -------
    AnnData with adata.obs['perturbation'] populated.
    """
    adata.obs['perturbation'] = adata.obs[guide_col].str.split(sep).str[0]
    return adata


def filter_by_cells_per_perturbation(adata, min_cells=20, obs_key='perturbation'):
    """Remove perturbations with fewer than min_cells cells."""
    counts = adata.obs[obs_key].value_counts()
    keep = counts[counts >= min_cells].index
    return adata[adata.obs[obs_key].isin(keep)].copy()
