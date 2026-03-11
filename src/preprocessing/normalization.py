"""Normalization utilities."""
import scanpy as sc

N_HVG = 2000
N_PCA_COMPONENTS = 50


def normalize_log1p(adata, target_sum=1e4):
    """Library-size normalize and log1p-transform counts."""
    sc.pp.normalize_total(adata, target_sum=target_sum)
    sc.pp.log1p(adata)
    return adata


def select_hvg_pca(adata, n_hvg=N_HVG, n_components=N_PCA_COMPONENTS):
    """Select highly variable genes and compute PCA."""
    sc.pp.highly_variable_genes(adata, n_top_genes=n_hvg, subset=True)
    sc.pp.scale(adata, max_value=10)
    sc.tl.pca(adata, n_comps=n_components)
    return adata
