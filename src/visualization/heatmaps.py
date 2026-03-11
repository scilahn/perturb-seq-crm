"""Heatmap visualization utilities."""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


def perturbation_similarity_heatmap(similarity_matrix, figsize=(12, 10),
                                     cmap='RdBu_r', title='Perturbation Similarity',
                                     save=None):
    """Plot a heatmap of perturbation similarity scores.

    Parameters
    ----------
    similarity_matrix : pd.DataFrame  (square, perturbation x perturbation)
    figsize : tuple
    cmap : str
    title : str
    save : str, optional
    """
    fig, ax = plt.subplots(figsize=figsize)
    sns.heatmap(similarity_matrix, cmap=cmap, center=0, ax=ax,
                xticklabels=True, yticklabels=True)
    ax.set_title(title)
    plt.tight_layout()
    if save:
        fig.savefig(save, dpi=150, bbox_inches='tight')
    return fig
