"""Signature scoring and GSEA wrappers."""
import scanpy as sc
import decoupler as dc


def score_gene_sets(adata, gene_sets: dict, method='gsva'):
    """Score gene sets using decoupler.

    Parameters
    ----------
    adata : AnnData (log-normalized)
    gene_sets : dict  {set_name: [gene1, gene2, ...]}
    method : str  'gsva' | 'aucell' | 'ulm'

    Returns
    -------
    AnnData with scores in obsm
    """
    net = []
    for name, genes in gene_sets.items():
        for g in genes:
            net.append({'source': name, 'target': g, 'weight': 1.0})
    import pandas as pd
    net_df = pd.DataFrame(net)

    dc.run_ulm(mat=adata, net=net_df, source='source', target='target',
               weight='weight', verbose=True, use_raw=False)
    return adata


def run_gsea(de_results, gene_set_db='h.all', organism='human'):
    """Run GSEA on DE results using GSEApy.

    Parameters
    ----------
    de_results : pd.DataFrame  with 'log2FoldChange' and index as gene names
    gene_set_db : str  MSigDB collection name
    organism : str

    Returns
    -------
    gseapy results object
    """
    import gseapy as gp
    ranked = de_results['log2FoldChange'].dropna().sort_values(ascending=False)
    res = gp.prerank(rnk=ranked, gene_sets=gene_set_db, organism=organism,
                     outdir=None, verbose=True)
    return res
