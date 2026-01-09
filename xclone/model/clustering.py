# tumor non-tumor classification and clustering post CNA calling
# Author: Jiamu James Qiao
# Date: 12/09/2025

import anndata as ad
import pandas as pd
import os
import numpy as np
from scipy.cluster.hierarchy import linkage, fcluster

from typing import Optional, Union
from anndata import AnnData
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import silhouette_score

# tumor non-tumor classification based on modified expression matrix
def tumor_classify(adata, layer_name, out_dir):
    """
    Cluster cells based on expression matrix using hierarchical clustering with Ward linkage
    and Euclidean distance, forcing 2 clusters. Label the cluster with the lowest aberration score
    (mean absolute deviation from median expression) as 'normal' and the other as 'tumor'.
    Modified from bcd
    Parameters:
    -----------
    adata : anndata.AnnData
    AnnData object containing the expression matrix in adata.layers[layer_name].
    layer_name : str
    Name of the layer containing the expression matrix (shape: n_cells, n_genes).
    Returns:
    --------
    anndata.AnnData
    Updated AnnData object with 'tumor_pred' column in obs ('normal', 'tumor').
    """
    # Step 1: Extract the expression matrix
    expr_matrix = adata.layers[layer_name] # Shape: (n_cells, n_genes)
    
    # Verify matrix shape
    n_cells, n_genes = expr_matrix.shape
    if expr_matrix.shape[0] != adata.n_obs:
        raise ValueError("Cell count mismatch between matrix and obs")

    # Step 2: Apply hierarchical clustering with Ward linkage and Euclidean distance
    Z = linkage(expr_matrix, method='ward', metric='euclidean')
    cluster_labels = fcluster(Z, t=2, criterion='maxclust') - 1 # Labels: 0 or 1

    # Step 3: Assign cluster labels based on aberration score
    # Compute aberration score as mean absolute deviation from median expression per gene
    aberration_scores = np.mean(np.abs(expr_matrix - 0), axis=1) # xclone: processed log ratio
    
    # Compute mean aberration scores per cluster
    mean_scores = [np.mean(aberration_scores[cluster_labels == i]) for i in [0, 1]]
    normal_cluster = np.argmin(mean_scores)
    tumor_pred = np.where(cluster_labels == normal_cluster, 'normal', 'tumor')

    # Step 4: Add predictions to adata.obs
    adata.obs['tumor_pred'] = tumor_pred
    # Print summary
    print(f"Processed {n_cells} cells and {n_genes} genes.")
    print(f"Mean aberration scores per cluster: {mean_scores}")
    print(f"Cluster labels: {['normal', 'tumor']}")
    for label in np.unique(tumor_pred):
        count = np.sum(tumor_pred == label)
    print(f"Number of cells in {label}: {count}")

    # Save predictions TSV
    predictions_df = pd.DataFrame({
        'barcode': adata.obs_names,
        'prediction': adata.obs['tumor_pred']
    })
    predictions_path = os.path.join(out_dir, 'xclone_tumor_predictions.tsv')
    predictions_df.to_csv(predictions_path, sep='\t', index=False)
    print(f"Tumor predictions saved to {predictions_path}")

    return adata

# clustering
def _to_dense(arr):
    """Convert sparse matrix or array to dense numpy array safely."""
    if hasattr(arr, "toarray"):
        return arr.toarray()
    elif hasattr(arr, "A"):  # scipy sparse
        return arr.A
    else:
        return np.asarray(arr)


def xclone_subclonal_analysis(
    combined_adata: AnnData,
    baf_adata: Optional[AnnData] = None,
    method: str = "combined",                    # 'prob', 'rdr_baf', 'combined'
    n_clones: Optional[int] = None,
    auto_n_clones: bool = True,
    max_clones: int = 12,
    pca_components: Union[int, float] = 50,
    random_state: int = 42,
    out_dir: str = "xclone_subclonal",
    sample_name: str = "sample"
) -> AnnData:
    """
    XClone clustering
    """
    np.random.seed(random_state)
    os.makedirs(out_dir, exist_ok=True)

    # 1. Extract prob1_merge (cell × gene × 4)
    prob_layer = _to_dense(combined_adata.layers["prob1_merge"])
    if prob_layer.shape[2] != 4:
        raise ValueError(f"prob1_merge must have 4 states, got {prob_layer.shape[2]}")

    n_cells, n_genes, _ = prob_layer.shape
    state_names = ["copy_loss", "loh", "copy_neutral", "copy_gain"]
    gene_names = combined_adata.var["GeneName"].tolist()

    # 2. Build clustering matrix
    prob_mean = np.nanmean(prob_layer, axis=1)  # (n_cells, 4)

    if method == "prob":
        X = prob_mean.copy()

    elif method == "rdr_baf":
        if baf_adata is None:
            raise ValueError("baf_adata required for method='rdr_baf'")
        rdr = _to_dense(combined_adata.layers["WMA_smoothed_log_ratio_ab_dynamic"])
        baf = _to_dense(baf_adata.layers["BAF_phased_WMA"])
        X_raw = StandardScaler().fit_transform(np.hstack([rdr, baf]))
        X = PCA(n_components=pca_components, random_state=random_state).fit_transform(X_raw)

    elif method == "combined":
        X = prob_mean.copy()
        if baf_adata is not None:
            rdr = _to_dense(combined_adata.layers["WMA_smoothed_log_ratio_ab_dynamic"])
            baf = _to_dense(baf_adata.layers["BAF_phased_WMA"])
            X_raw = StandardScaler().fit_transform(np.hstack([rdr, baf]))
            X_pcs = PCA(n_components=pca_components, random_state=random_state).fit_transform(X_raw)
            X = np.hstack([X, X_pcs])
    else:
        raise ValueError("method must be 'prob', 'rdr_baf', or 'combined'")

    # 3. Determine number of clones
    if n_clones is None and auto_n_clones:
        scores = []
        ks = range(2, max_clones + 1)
        for k in ks:
            lbl = fcluster(linkage(X, "ward"), k, criterion="maxclust") - 1
            scores.append(silhouette_score(X, lbl) if len(np.unique(lbl)) >= 2 else -1)
        final_k = ks[np.argmax(scores)]
        print(f"Auto-selected {final_k} clones")
    else:
        final_k = n_clones or 3

    # 4. Final clustering
    labels = fcluster(linkage(X, "ward"), final_k, criterion="maxclust") - 1
    clone_ids = np.array([f"clone_{i:02d}" for i in labels])

    # 5. Clone-level consensus CNA
    print("Computing clone-level consensus CNA states...")
    clone_consensus_prob = np.zeros((final_k, n_genes, 4), dtype=np.float32)
    for cid in range(final_k):
        mask = labels == cid
        if mask.any():
            clone_consensus_prob[cid] = np.mean(prob_layer[mask], axis=0)

    dominant_idx = np.argmax(clone_consensus_prob, axis=2)  # (final_k, n_genes)
    dominant_states = np.array(state_names)[dominant_idx]   # string array

    # 6. Cell-level layers
    cell_cna_codes = dominant_idx[labels]                    # (n_cells, n_genes)
    cell_cna_prob = np.zeros((n_cells, n_genes, 4), dtype=np.float32)
    for s in range(4):
        cell_cna_prob[:, :, s] = (cell_cna_codes == s)

    # 7. Save results
    adata_out = combined_adata.copy()
    adata_out.obs["clone_id"] = clone_ids

    # Layers
    adata_out.layers["clone_level_cna"] = cell_cna_codes.astype("int8")
    adata_out.layers["clone_level_cna_prob_for_plot"] = cell_cna_prob
    adata_out.uns["clone_level_cna_states"] = state_names
    '''
    # Save AnnData
    adata_out.write_h5ad(os.path.join(out_dir, f"{sample_name}.h5ad"))

    # TSVs
    pd.DataFrame({"cell_barcode": adata_out.obs_names, "clone_id": clone_ids}) \
        .to_csv(os.path.join(out_dir, f"{sample_name}_cell_to_clone.tsv"), sep="\t", index=False)

    pd.DataFrame({"clone_id": [f"clone_{i:02d}" for i in range(final_k)],
                  "n_cells": np.bincount(labels)}) \
        .to_csv(os.path.join(out_dir, f"{sample_name}_clone_summary.tsv"), sep="\t", index=False)

    # clone × gene: first column = clone_id, rest = gene names
    clone_gene_df = pd.DataFrame(
        dominant_states,
        index=[f"clone_{i:02d}" for i in range(final_k)],
        columns=gene_names
    )
    clone_gene_df = clone_gene_df.reset_index().rename(columns={"index": "clone_id"})
    clone_gene_df.to_csv(os.path.join(out_dir, f"{sample_name}_clone_gene_cna.tsv"), sep="\t", index=False)

    print(f"\nAll results saved to: {out_dir}")
    print(f"   • {sample_name}.h5ad")
    print(f"   • {sample_name}_clone_gene_cna.tsv  ← clone_id | gene1 | gene2 | ...")
    print(f"   • Layers: clone_level_cna + clone_level_cna_prob_for_plot")
    '''

    return adata_out

# refine initial clones using Bayesian EM with Dirichlet prior     
def refine_clones_bayesian(
    adata: AnnData,
    initial_col: str = "clone_id",
    prob_layer: str = "prob1_merge",           # (n_cells, n_genes, 4)
    n_iter: int = 15,
    alpha: float = 20.0,                       # Dirichlet prior strength
    min_cells: int = 50,
    n_clones: Optional[int] = None,            # fixed number of clones
    out_dir: str = "refined_clones",
    sample_name: str = "sample"
) -> AnnData:
    """
    Bayesian refinement of XClone subclones using EM + Dirichlet prior.
    - Uses real gene names from adata.var["GeneName"]
    - Saves posterior probability per cell to TSV
    - Maximizes: log P(data, z | θ) + log P(θ)  →  ELBO-like objective
    """

    # 1. Extract real gene names
    if "GeneName" in adata.var.columns:
        gene_names = adata.var["GeneName"].tolist()
    else:
        gene_names = adata.var_names.tolist()  # fallback

    # Extract probabilities
    prob = adata.layers[prob_layer]
    # if hasattr(prob, "toarray"):
    #     prob = prob.toarray()  # (n_cells, n_genes, 4)
    prob = _to_dense(prob)

    initial_labels = adata.obs[initial_col].astype(str).values
    unique_clones = np.unique(initial_labels[initial_labels != "nan"])

    print(f"Starting Bayesian refinement with {len(unique_clones)} initial clones...")

    # Initial clone profiles: mean posterior probability per clone
    temp_profiles = []
    temp_sizes = []
    temp_clone_names = []

    for clone in unique_clones:
        mask = initial_labels == clone
        size = mask.sum()
        if size < min_cells:
            continue
        profile = np.mean(prob[mask], axis=0) + 1e-6
        profile = profile / profile.sum(axis=1, keepdims=True)
        temp_profiles.append(profile)
        temp_sizes.append(size)
        temp_clone_names.append(clone)

    if len(temp_profiles) < 2:
        print("Not enough clones for refinement.")
        return adata

    # NEW: handle fixed n_clones
    if n_clones is not None:
        if n_clones > len(temp_profiles):
            print(f"Warning: requested n_clones={n_clones} > available valid clones ({len(temp_profiles)}). Using all available.")
            selected_idx = np.arange(len(temp_profiles))
        else:
            # select the n_clones largest initial clones by cell count
            selected_idx = np.argsort(temp_sizes)[-n_clones:]
        print(f"Using fixed n_clones = {n_clones} (selected largest initial clones)")
    else:
        selected_idx = np.arange(len(temp_profiles))
        print(f"Using all {len(selected_idx)} valid initial clones (n_clones=None)")

    clone_profiles = np.stack([temp_profiles[i] for i in selected_idx])
    K = len(clone_profiles)
    '''
    clone_profiles = []
    valid_clones = []
    for clone in unique_clones:
        mask = initial_labels == clone
        if mask.sum() < min_cells:
            continue
        profile = np.mean(prob[mask], axis=0) + 1e-6
        profile = profile / profile.sum(axis=1, keepdims=True)
        clone_profiles.append(profile)
        valid_clones.append(clone)

    if len(clone_profiles) < 2:
        print("Not enough clones for refinement.")
        return adata

    clone_profiles = np.stack(clone_profiles)  # (K, n_genes, 4), clone K's CNA profile
    K = len(clone_profiles)
    '''
    n_cells = prob.shape[0]

    # Prior over clone proportions (Dirichlet)
    prior_counts = np.ones(K) * alpha

    # EM loop
    responsibilities = np.zeros((n_cells, K))  # γ_ik = P(z_i=k | x_i)
    log_likelihoods = []

    for it in range(n_iter):
        # E-step: posterior responsibility
        log_resp = np.zeros((n_cells, K))
        for k in range(K):
            # log P(x_i | z_i=k) = sum_g log P(state_g | gene_g, clone_k)
            log_p = np.log(clone_profiles[k] + 1e-15)
            log_resp[:, k] = np.sum(prob * log_p, axis=(1, 2))

        # Add log prior P(z_i=k)
        clone_sizes = responsibilities.sum(axis=0) + prior_counts
        log_resp += np.log(clone_sizes)

        # Normalize
        log_resp -= log_resp.max(axis=1, keepdims=True)
        resp = np.exp(log_resp)
        responsibilities = resp / resp.sum(axis=1, keepdims=True)

        # M-step: update clone profiles
        for k in range(K):
            weight = responsibilities[:, k][:, None, None]
            clone_profiles[k] = (weight * prob).sum(axis=0) + 1e-6
            clone_profiles[k] /= clone_profiles[k].sum(axis=1, keepdims=True)

        # Track likelihood
        ll = np.sum(responsibilities * log_resp)
        log_likelihoods.append(ll)
        print(f"  Iter {it+1:2d} | Log-likelihood = {ll:,.0f}")

    # Final assignment
    final_assign = np.argmax(responsibilities, axis=1)
    final_clone_ids = [f"C{i+1:02d}" for i in final_assign]  # C01, C02, ...

    # Add results to AnnData
    adata.obs["clone_id_refined"] = pd.Categorical(final_clone_ids)
    adata.obs["clone_posterior_max"] = responsibilities.max(axis=1)
    adata.obsm["X_clone_posteriors"] = responsibilities.astype("float32")

    # XClone-compatible CNA probability reconstruction
    cell_cna_prob = responsibilities[:, :, None, None] * clone_profiles[None, :, :, :]
    cell_cna_prob = cell_cna_prob.sum(axis=1)
    cell_cna_prob = np.nan_to_num(cell_cna_prob)
    cell_cna_prob /= cell_cna_prob.sum(axis=2, keepdims=True)
    cell_cna_prob = np.nan_to_num(cell_cna_prob, nan=0.25)

    adata.layers["prob1_merge_refined"] = cell_cna_prob.astype("float32")
    adata.layers["state_merge_refined"] = np.argmax(cell_cna_prob, axis=2).astype("int8")

    print(f"\nRefinement complete! {K} high-quality clones.")
    print(f"XClone-compatible layers added:")
    print(f"   • layers['prob1_merge_refined']  ← (n_cells, n_genes, 4)")
    print(f"   • layers['state_merge_refined']  ← (n_cells, n_genes)")
    print(f"All results saved to: {out_dir}")

    # Save posterior probabilities per cell
    post_df = pd.DataFrame(
        responsibilities,
        index=adata.obs_names,
        columns=[f"C{i+1:02d}" for i in range(K)]
    )
    post_df.index.name = "cell_barcode"
    post_df["clone_id_refined"] = final_clone_ids
    post_df["max_posterior"] = responsibilities.max(axis=1)
    post_df = post_df[["clone_id_refined", "max_posterior"] + [f"C{i+1:02d}" for i in range(K)]]
    post_df.to_csv(os.path.join(out_dir, f"{sample_name}_cell_clone_posteriors.tsv"), sep="\t")
    print(f"Saved cell × clone posteriors → {sample_name}_cell_clone_posteriors.tsv")

    # Save refined clone × gene CNA states
    dominant_states = np.argmax(clone_profiles, axis=2)
    state_map = {0: "copy_loss", 1: "loh", 2: "copy_neutral", 3: "copy_gain"}
    state_names = np.vectorize(state_map.get)(dominant_states)

    clone_gene_df = pd.DataFrame(
        state_names,
        index=[f"C{i+1:02d}" for i in range(K)],
        columns=gene_names
    )
    clone_gene_df = clone_gene_df.reset_index().rename(columns={"index": "clone_id"})
    clone_gene_df.to_csv(os.path.join(out_dir, f"{sample_name}_refined_clone_gene_cna.tsv"), sep="\t", index=False)

    print(f"\nRefinement complete! Final clones: {K}")
    print(f"Top confident cell per clone:")
    for k in range(K):
        best_cell = adata.obs_names[responsibilities[:, k].argmax()]
        conf = responsibilities[:, k].max()
        print(f"   {f'C{k+1:02d}'} ← {best_cell} (p = {conf:.3f})")

    return adata