"""GWAS Catalog locus mapping."""
import requests
import pandas as pd

GWAS_REST_API = "https://www.ebi.ac.uk/gwas/rest/api"

CARDIAC_TRAITS = [
    'coronary artery disease',
    'heart failure',
    'atrial fibrillation',
    'myocardial infarction',
    'cardiac hypertrophy',
]

RENAL_TRAITS = [
    'chronic kidney disease',
    'estimated glomerular filtration rate',
    'end-stage renal disease',
]

METABOLIC_TRAITS = [
    'type 2 diabetes',
    'LDL cholesterol',
    'triglycerides',
    'HDL cholesterol',
    'body mass index',
]


def search_associations(trait_keyword, p_threshold=5e-8, max_results=1000):
    """Fetch GWAS associations for a trait keyword.

    Parameters
    ----------
    trait_keyword : str
    p_threshold : float
    max_results : int

    Returns
    -------
    pd.DataFrame with SNP, mapped_gene, trait, p_value columns
    """
    url = f"{GWAS_REST_API}/associations/search"
    params = {
        'q': trait_keyword,
        'size': max_results,
    }
    resp = requests.get(url, params=params, timeout=30)
    resp.raise_for_status()
    data = resp.json()

    rows = []
    for assoc in data.get('_embedded', {}).get('associations', []):
        p_val = assoc.get('pvalue', 1.0)
        if p_val > p_threshold:
            continue
        snp = assoc.get('strongestRiskAlleles', [{}])[0].get('riskAlleleName', '')
        genes = [g.get('geneName', '') for g in assoc.get('loci', [{}])[0].get('authorReportedGenes', [])]
        rows.append({
            'snp': snp,
            'mapped_genes': ','.join(genes),
            'trait': trait_keyword,
            'p_value': p_val,
        })
    return pd.DataFrame(rows)


def intersect_with_perturbed_genes(gwas_df, perturbed_genes, gene_col='mapped_genes'):
    """Find GWAS genes that overlap with perturbation screen targets.

    Parameters
    ----------
    gwas_df : pd.DataFrame  (output of search_associations)
    perturbed_genes : list of str
    gene_col : str

    Returns
    -------
    tuple (overlapping_genes, missing_genes)
    """
    gwas_genes = set()
    for entry in gwas_df[gene_col].dropna():
        gwas_genes.update(g.strip() for g in entry.split(',') if g.strip())

    perturbed = set(perturbed_genes)
    overlapping = gwas_genes & perturbed
    missing = gwas_genes - perturbed
    return sorted(overlapping), sorted(missing)
