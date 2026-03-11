"""Open Targets Genetics API queries."""
import requests

OT_GENETICS_API = "https://api.genetics.opentargets.org/graphql"


def get_genes_for_locus(variant_id, radius_kb=500):
    """Query Open Targets Genetics for genes near a variant.

    Parameters
    ----------
    variant_id : str  e.g. '1_109817590_G_T'
    radius_kb : int

    Returns
    -------
    list of gene symbols
    """
    query = """
    query genesForVariant($variantId: String!) {
      genesForVariantSchema(variantId: $variantId) {
        gene { symbol }
        overallScore
      }
    }
    """
    resp = requests.post(
        OT_GENETICS_API,
        json={'query': query, 'variables': {'variantId': variant_id}},
        timeout=30,
    )
    resp.raise_for_status()
    data = resp.json()
    entries = data.get('data', {}).get('genesForVariantSchema', [])
    return [e['gene']['symbol'] for e in entries]


def prioritize_gwas_genes(gene_list, disease_efo):
    """Score a list of genes against a disease EFO term using Open Targets.

    Parameters
    ----------
    gene_list : list of str  (Ensembl IDs or gene symbols)
    disease_efo : str  e.g. 'EFO_0000319' (coronary artery disease)

    Returns
    -------
    pd.DataFrame with gene, overallAssociationScore
    """
    import pandas as pd
    query = """
    query targets($ensemblIds: [String!]!, $efoId: String!) {
      targets(ensemblIds: $ensemblIds) {
        id
        approvedSymbol
        associatedDiseases(filter: {id: [$efoId]}) {
          rows { score }
        }
      }
    }
    """
    resp = requests.post(
        "https://api.platform.opentargets.org/api/v4/graphql",
        json={'query': query, 'variables': {'ensemblIds': gene_list, 'efoId': disease_efo}},
        timeout=30,
    )
    resp.raise_for_status()
    targets = resp.json().get('data', {}).get('targets', [])
    rows = []
    for t in targets:
        score = 0.0
        assoc = t.get('associatedDiseases', {}).get('rows', [])
        if assoc:
            score = assoc[0]['score']
        rows.append({'gene': t['approvedSymbol'], 'ensembl_id': t['id'], 'ot_score': score})
    return pd.DataFrame(rows).sort_values('ot_score', ascending=False)
