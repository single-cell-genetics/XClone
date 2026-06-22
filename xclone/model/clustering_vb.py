"""
Variational-Bayes subclonal clustering for XClone combine module.

Methods
-------
knn_em : hierarchical Ward + Bayesian EM refine (original; handled in combine wrap)
rdr    : GaussianMixtureVB on RDR log-ratio layer
baf    : BinomMixtureVB on phased BAF AD/DP bins
joint  : JointBinomGaussianVB (shared clone assignment)
"""

from __future__ import annotations

import os
from typing import Optional

import numpy as np
import pandas as pd
from anndata import AnnData

from ._vb_binomial import get_binom_mixture_vb_class
from ._vb_models import GaussianMixtureVB, JointBinomGaussianVB, impute_gene_medians, to_dense

CLUSTERING_METHODS = ("knn_em", "rdr", "baf", "joint")


def _resolve_n_clones(n_clones: Optional[int], default: int = 3) -> int:
    k = int(n_clones) if n_clones is not None else default
    if k < 2:
        raise ValueError(f"n_clones must be >= 2, got {k}")
    return k


def _assign_clone_id(ID_prob: np.ndarray, prefix: str = "clone") -> np.ndarray:
    labels = ID_prob.argmax(axis=1)
    return np.array([f"{prefix}_{i:02d}" for i in labels])


def _compute_clone_gene_prob(
    prob_layer: np.ndarray,
    labels: np.ndarray,
    n_clones: int,
) -> np.ndarray:
    """Compute clone-level average CNA probability: (K, n_genes, n_states)."""
    n_genes, n_states = prob_layer.shape[1], prob_layer.shape[2]
    clone_gene_prob = np.zeros((n_clones, n_genes, n_states), dtype=np.float32)
    for cid in range(n_clones):
        mask = labels == cid
        if mask.any():
            clone_gene_prob[cid] = np.mean(prob_layer[mask], axis=0)
    return clone_gene_prob


def _write_vb_posteriors(adata: AnnData, ID_prob: np.ndarray, out_dir: str, sample_name: str):
    os.makedirs(out_dir, exist_ok=True)
    k = ID_prob.shape[1]
    cols = [f"C{i+1:02d}" for i in range(k)]
    df = pd.DataFrame(ID_prob, index=adata.obs_names, columns=cols)
    df.index.name = "cell_barcode"
    df["clone_id_vb"] = _assign_clone_id(ID_prob)
    df["max_posterior"] = ID_prob.max(axis=1)
    path = os.path.join(out_dir, f"{sample_name}_vb_cell_clone_posteriors.tsv")
    df.to_csv(path, sep="\t")
    print(f"[XClone VB] Saved posteriors -> {path}")


def xclone_vb_clustering(
    combined_adata: AnnData,
    baf_adata: AnnData,
    *,
    method: str,
    n_clones: Optional[int] = None,
    rdr_layer: str = "WMA_smoothed_log_ratio_ab_dynamic",
    baf_ad_layer: str = "ad_bin_phased_rev",
    baf_dp_layer: str = "dp_bin",
    n_init: int = 50,
    min_iter: int = 30,
    random_seed: int = 0,
    out_dir: str = ".",
    sample_name: str = "sample",
    run_refinement: bool = True,
) -> AnnData:
    """
    Run VB-based subclonal clustering and optionally Bayesian EM refinement.

    Parameters
    ----------
    method : {'rdr', 'baf', 'joint'}
    n_clones : int
        Number of clones K (required for VB methods).
    """
    if method not in CLUSTERING_METHODS:
        raise ValueError(f"method must be one of {CLUSTERING_METHODS}, got {method!r}")
    if method == "knn_em":
        raise ValueError("knn_em is handled by xclone_subclonal_analysis + refine_clones_bayesian")

    k = _resolve_n_clones(n_clones)
    max_iter_pre = max(15, min_iter // 2)
    verbose = True

    print(f"[XClone VB clustering] method={method} K={k} n_init={n_init}")

    elbo_last = np.nan

    if method == "rdr":
        if rdr_layer not in combined_adata.layers:
            raise KeyError(f"RDR layer {rdr_layer!r} missing from combined AnnData")
        Y = impute_gene_medians(to_dense(combined_adata.layers[rdr_layer]))
        model = GaussianMixtureVB(n_cell=combined_adata.n_obs, n_var=Y.shape[1], n_donor=k)
        model.fit(
            Y,
            n_init=n_init,
            min_iter=min_iter,
            max_iter=min_iter,
            max_iter_pre=max_iter_pre,
            random_seed=random_seed,
            verbose=verbose,
        )
        ID_prob = model.ID_prob
        elbo_last = float(model.ELBO_iters[-1])

    elif method == "baf":
        if baf_adata is None:
            raise ValueError("baf_adata required for method='baf'")
        for layer in (baf_ad_layer, baf_dp_layer):
            if layer not in baf_adata.layers:
                raise KeyError(f"BAF layer {layer!r} missing")
        ad = to_dense(baf_adata.layers[baf_ad_layer])
        dp = to_dense(baf_adata.layers[baf_dp_layer])
        BinomMixtureVB = get_binom_mixture_vb_class()
        model = BinomMixtureVB(n_var=baf_adata.n_vars, n_cell=baf_adata.n_obs, n_donor=k)
        model.fit(ad.T, dp.T, n_init=n_init, min_iter=min_iter, max_iter=min_iter, max_iter_pre=max_iter_pre, random_seed=random_seed, verbose=verbose)
        ID_prob = model.ID_prob
        elbo_last = float(model.ELBO_iters[-1])

    elif method == "joint":
        if baf_adata is None:
            raise ValueError("baf_adata required for method='joint'")
        if rdr_layer not in combined_adata.layers:
            raise KeyError(f"RDR layer {rdr_layer!r} missing from combined AnnData")
        for layer in (baf_ad_layer, baf_dp_layer):
            if layer not in baf_adata.layers:
                raise KeyError(f"BAF layer {layer!r} missing")
        Y = to_dense(combined_adata.layers[rdr_layer])
        joint = JointBinomGaussianVB(
            n_cell=combined_adata.n_obs,
            n_var_baf=baf_adata.n_vars,
            n_var_rdr=Y.shape[1],
            n_donor=k,
        )
        joint.fit(
            baf_adata.layers[baf_ad_layer],
            baf_adata.layers[baf_dp_layer],
            Y,
            n_init=n_init,
            max_iter=min_iter,
            max_iter_pre=max_iter_pre,
            random_seed=random_seed,
            verbose=verbose,
        )
        ID_prob = joint.ID_prob
        elbo_last = float(joint.ELBO_iters[-1])
    else:
        raise ValueError(method)

    adata_out = combined_adata.copy()
    labels = ID_prob.argmax(axis=1)
    adata_out.obs["clone_id"] = _assign_clone_id(ID_prob)
    adata_out.obs["clone_posterior_max"] = ID_prob.max(axis=1)
    adata_out.obsm["X_clone_posteriors_vb"] = ID_prob.astype("float32")
    adata_out.uns["xclone_clustering_method"] = method
    adata_out.uns["xclone_clustering_n_clones"] = k
    clone_gene_prob = None
    if "prob1_merge" in adata_out.layers:
        prob_layer = to_dense(adata_out.layers["prob1_merge"])
        if prob_layer.ndim == 3 and prob_layer.shape[2] == 4:
            clone_gene_prob = _compute_clone_gene_prob(prob_layer, labels, k)
            adata_out.uns["clone_gene_prob"] = clone_gene_prob
            adata_out.uns["clone_gene_state"] = np.argmax(clone_gene_prob, axis=2).astype("int8")
            adata_out.uns["clone_ids"] = [f"clone_{i:02d}" for i in range(k)]

    _write_vb_posteriors(adata_out, ID_prob, out_dir, sample_name)
    print(f"[XClone VB] ELBO last = {elbo_last:.2f}")

    if run_refinement:
        from .clustering import refine_clones_bayesian

        print("[XClone VB] Running Bayesian EM refinement on VB initial labels...")
        adata_out = refine_clones_bayesian(
            adata_out,
            initial_col="clone_id",
            prob_layer="prob1_merge",
            n_iter=15,
            alpha=20.0,
            min_cells=50,
            n_clones=k,
            out_dir=out_dir,
            sample_name=sample_name,
        )
    elif clone_gene_prob is not None:
        # Keep combine plotting behavior consistent: always provide a clone-level
        # per-cell layer even when EM refinement is disabled.
        adata_out.obs["clone_id_refined"] = adata_out.obs["clone_id"].astype("category")
        adata_out.layers["prob1_merge_refined"] = clone_gene_prob[labels].astype("float32")
        adata_out.layers["state_merge_refined"] = np.argmax(
            adata_out.layers["prob1_merge_refined"], axis=2
        ).astype("int8")
        print("[XClone VB] EM refinement disabled; exported clone-level refined heatmap layers.")

    return adata_out


def run_post_combine_clustering(
    combine_Xdata: AnnData,
    baf_adata: Optional[AnnData],
    config,
    out_dir: str,
    sample_name: str,
) -> AnnData:
    """
    Dispatch clustering after combine based on ``config.clustering_method``.

    knn_em -> Ward hierarchical + Bayesian EM (original).
    rdr / baf / joint -> VB models in this module.
    """
    method = getattr(config, "clustering_method", "knn_em")
    if method not in CLUSTERING_METHODS:
        raise ValueError(f"Unknown clustering_method={method!r}; choose from {CLUSTERING_METHODS}")

    if method == "knn_em":
        from .clustering import refine_clones_bayesian, xclone_subclonal_analysis

        print("[XClone clustering] method=knn_em (Ward + Bayesian EM)")
        combine_Xdata = xclone_subclonal_analysis(
            combined_adata=combine_Xdata,
            baf_adata=baf_adata,
            method="combined",
            n_clones=config.n_clones,
            out_dir=out_dir,
            sample_name=sample_name,
        )
        combine_Xdata = refine_clones_bayesian(
            adata=combine_Xdata,
            initial_col="clone_id",
            prob_layer="prob1_merge",
            n_iter=15,
            alpha=20.0,
            min_cells=50,
            n_clones=config.n_clones,
            out_dir=out_dir,
            sample_name=sample_name,
        )
        combine_Xdata.uns["xclone_clustering_method"] = "knn_em"
        return combine_Xdata

    return xclone_vb_clustering(
        combine_Xdata,
        baf_adata,
        method=method,
        n_clones=config.n_clones,
        rdr_layer=getattr(config, "vb_rdr_layer", "WMA_smoothed_log_ratio_ab_dynamic"),
        baf_ad_layer=getattr(config, "vb_baf_ad_layer", "ad_bin_phased_rev"),
        baf_dp_layer=getattr(config, "vb_baf_dp_layer", "dp_bin"),
        n_init=getattr(config, "vb_n_init", 50),
        min_iter=getattr(config, "vb_min_iter", 30),
        random_seed=getattr(config, "vb_random_seed", 0),
        out_dir=out_dir,
        sample_name=sample_name,
        run_refinement=getattr(config, "vb_run_refinement", True),
    )
