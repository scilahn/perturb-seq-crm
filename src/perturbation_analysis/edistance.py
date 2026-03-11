"""E-distance / E-test utilities."""
import numpy as np
import pandas as pd


def compute_edistances(adata, obs_key='perturbation', rep_key='X_pca',
                       control_label='control'):
    """Compute pairwise E-distances between perturbations and control.

    Wraps scperturb.edist if available, otherwise falls back to manual computation.

    Parameters
    ----------
    adata : AnnData with adata.obsm[rep_key]
    obs_key : str
    rep_key : str
    control_label : str

    Returns
    -------
    pd.Series  index=perturbation names, values=E-distance to control
    """
    try:
        from scperturb import edist
        estats = edist(adata, obs_key=obs_key)
        return estats[control_label].drop(control_label)
    except ImportError:
        pass

    X = adata.obsm[rep_key]
    labels = adata.obs[obs_key].values
    ctrl_idx = np.where(labels == control_label)[0]
    X_ctrl = X[ctrl_idx]

    results = {}
    for pert in np.unique(labels):
        if pert == control_label:
            continue
        pert_idx = np.where(labels == pert)[0]
        X_pert = X[pert_idx]
        e = _edist(X_pert, X_ctrl)
        results[pert] = e
    return pd.Series(results)


def _edist(X, Y):
    """Compute E-distance between two sample matrices."""
    from scipy.spatial.distance import cdist
    n, m = len(X), len(Y)
    XY = cdist(X, Y).mean()
    XX = cdist(X, X).mean()
    YY = cdist(Y, Y).mean()
    return 2 * XY - XX - YY
