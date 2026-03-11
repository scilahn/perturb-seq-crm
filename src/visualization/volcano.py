"""Volcano plot utilities."""
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd


def volcano_plot(de_results, lfc_col='log2FoldChange', pval_col='padj',
                 lfc_thresh=1.0, pval_thresh=0.05, top_n=10,
                 title='', figsize=(8, 6), save=None):
    """Generate a volcano plot from DE results.

    Parameters
    ----------
    de_results : pd.DataFrame  (genes as index)
    lfc_col, pval_col : str
    lfc_thresh : float
    pval_thresh : float
    top_n : int  — label top N significant genes
    title : str
    figsize : tuple
    save : str, optional
    """
    df = de_results[[lfc_col, pval_col]].dropna().copy()
    df['-log10p'] = -np.log10(df[pval_col].clip(lower=1e-300))

    fig, ax = plt.subplots(figsize=figsize)
    ax.scatter(df[lfc_col], df['-log10p'], alpha=0.4, s=10, color='grey')

    sig = df[(df[pval_col] < pval_thresh) & (df[lfc_col].abs() > lfc_thresh)]
    ax.scatter(sig[lfc_col], sig['-log10p'], alpha=0.7, s=15, color='red')

    for gene in sig.nsmallest(top_n, pval_col).index:
        ax.annotate(gene, (sig.loc[gene, lfc_col], sig.loc[gene, '-log10p']),
                    fontsize=7)

    ax.axhline(-np.log10(pval_thresh), color='blue', lw=0.8, ls='--')
    ax.axvline(lfc_thresh, color='black', lw=0.8, ls='--')
    ax.axvline(-lfc_thresh, color='black', lw=0.8, ls='--')
    ax.set_xlabel('log2 Fold Change')
    ax.set_ylabel('-log10(padj)')
    ax.set_title(title)

    if save:
        fig.savefig(save, dpi=150, bbox_inches='tight')
    return fig
