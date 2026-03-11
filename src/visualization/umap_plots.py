"""UMAP visualization utilities."""
import scanpy as sc
import matplotlib.pyplot as plt


def plot_perturbation_umap(adata, obs_key='perturbation', highlight=None,
                           figsize=(8, 6), save=None):
    """Plot UMAP colored by perturbation label.

    Parameters
    ----------
    adata : AnnData with adata.obsm['X_umap']
    obs_key : str
    highlight : list of str, optional  — perturbations to highlight
    figsize : tuple
    save : str, optional  — path to save figure
    """
    fig, ax = plt.subplots(figsize=figsize)
    sc.pl.umap(adata, color=obs_key, ax=ax, show=False)
    if save:
        fig.savefig(save, dpi=150, bbox_inches='tight')
    return fig
