"""SCVI / CPA model training utilities."""
import scvi


def train_scvi(adata, batch_key=None, n_latent=30, n_epochs=400):
    """Train SCVI model on raw counts.

    Parameters
    ----------
    adata : AnnData  (raw counts in .X)
    batch_key : str, optional
    n_latent : int
    n_epochs : int

    Returns
    -------
    Trained SCVI model
    """
    scvi.model.SCVI.setup_anndata(adata, batch_key=batch_key)
    model = scvi.model.SCVI(adata, n_latent=n_latent)
    model.train(max_epochs=n_epochs)
    adata.obsm['X_scVI'] = model.get_latent_representation()
    return model


def train_cpa(adata, perturbation_key='perturbation', control_label='control',
              n_latent=64, n_epochs=500):
    """Train CPA model for perturbation effect decomposition.

    Parameters
    ----------
    adata : AnnData  (log-normalized)
    perturbation_key : str
    control_label : str
    n_latent : int
    n_epochs : int

    Returns
    -------
    Trained CPA model
    """
    scvi.external.CPA.setup_anndata(
        adata,
        perturbation_key=perturbation_key,
        control_group=control_label,
    )
    model = scvi.external.CPA(adata, n_latent=n_latent)
    model.train(max_epochs=n_epochs)
    return model
