"""Base functions for XClone RDR Gaussian model.
"""

# Author: Jiamu James Qiao
# Date: 2025-06-10
# update: 

import numpy as np
import anndata as an
from sklearn.neighbors import KDTree
from joblib import Parallel, delayed
from scipy.ndimage import gaussian_filter1d, median_filter
from sklearn.decomposition import PCA
from sklearn.neighbors import NearestNeighbors

import xclone
from ..plot.smooth_plot import raw_ratio_visualization
from .._logging import get_logger
from datetime import datetime, timezone

# adaptive baseline for RDR Gaussian model
def preprocess_adaptive_baseline(
    RDR_Xdata,
    cell_anno_key='spot_anno',
    ref_celltype='normal',
    k=5,
    pseudo_count=1e-8,
    verbose=False,
    plot=False
):
    """
    Compute adaptive baseline log ratio for each cell and gene using k-nearest normal cells as reference.

    Steps:
    1. Normalize raw expression by library size.
    2. Add a pseudocount to avoid log(0).
    3. Log-transform the normalized expression.
    4. For each cell, find its k nearest normal cells (using KDTree in log space).
    5. For each cell and gene, compute the mean expression of its k nearest normal cells as the adaptive baseline.
    6. Subtract the adaptive baseline from the log-normalized expression to obtain the log ratio.
    7. Store the result in a new layer 'log_ratio_ab'.
    8. Optionally, visualize the log ratio.

    Parameters:
    - RDR_Xdata: AnnData object with 'raw_expr' layer (sparse matrix)
    - k: int, number of nearest neighbors (default: 50)
    - pseudo_count: float, value to add before log (default: 1e-8)
    - verbose: bool, print intermediate results (default: False)
    - plot: bool, plot the log ratio (default: True)

    Returns:
    - tmp_adata: AnnData object with new layer 'log_ratio_ab'

    Notes:
    - use this version to avoid shift after WMA
    - for this simulated dataset, pseudocount 1e-6 works better that 1e-8
    - k=5 work better than k=50
    """
    # from sklearn.neighbors import KDTree

    print(f"Generating adaptive baseline with k = {k}, pseudo_count = {pseudo_count}")

    tmp_adata = RDR_Xdata.copy()
    Y = tmp_adata.layers['raw_expr'].A
    Y_32 = Y.astype(np.float32)

    # Step 1: Normalize by library size
    library_sizes = np.sum(Y_32, axis=1, keepdims=True)
    normalized = Y_32 / library_sizes

    # Step 2: Add pseudocount
    normalized = normalized + pseudo_count

    # Step 3: Log transformation
    log_normalized = np.log(normalized)

    if verbose:
        print('log_normalized:')
        print(log_normalized)

    # Identify normal cells
    if isinstance(ref_celltype, list):
        # Multiple reference cell types
        is_not_normal = ~tmp_adata.obs[cell_anno_key].isin(ref_celltype)
    else:
        # Single reference cell type
        is_not_normal = tmp_adata.obs[cell_anno_key].values != ref_celltype

    # Prepare data for KDTree
    normal_cells = log_normalized[~is_not_normal, :]
    all_cells = log_normalized

    # test if number of normal cells is less than k
    if normal_cells.shape[0] < k:
        # If not enough normal cells, use all available normal cells
        k = normal_cells.shape[0]
        print(f"Warning: Only {normal_cells.shape[0]} normal cells available, using k={k} for adaptive baseline.")


    # Build KDTree and query k nearest normal cells for each cell
    kdtree = KDTree(normal_cells)
    distances, indices = kdtree.query(all_cells, k=k)

    # Calculate adaptive baseline for each cell
    adaptive_baseline = np.zeros_like(all_cells)
    for i, neighbor_indices in enumerate(indices):
        nearest_neighbors = normal_cells[neighbor_indices, :]
        # Scale mean by library size ratio
        adaptive_baseline[i, :] = nearest_neighbors.mean(axis=0) * (
            np.sum(all_cells[i, :]) / (np.sum(nearest_neighbors, axis=1, keepdims=True)).mean()
        )

    log_ratio_ab = all_cells - adaptive_baseline

    if verbose:
        print('adaptive_baseline:')
        print(adaptive_baseline)
        print('log_ratio_ab:')
        print(log_ratio_ab)

    tmp_adata.layers['log_ratio_ab'] = log_ratio_ab

    if plot:
        # plot log ratio
        raw_ratio_visualization(
            tmp_adata,
            cell_anno_key=cell_anno_key,
            Xlayer="log_ratio_ab",
            vmin=-0.7,
            vmax=0.7,
            colorbar_name="RDR raw ratio (log)"
        )

    return tmp_adata


# smoothing, WMA is used
def WMA_smoothing(
    RDR_Xdata,
    start_layer='log_ratio_ab',
    wma_varp='WMA_connect',
    knn_obsp='connectivities',
    method='wma',  # 'wma', 'knn', or 'both'
    ):
    """
    Smooth an AnnData layer using WMA and/or KNN smoothing.

    Parameters:
    - RDR_Xdata: AnnData object
    - start_layer: str, layer to smooth
    - wma_varp: str, key in .varp for WMA matrix
    - knn_obsp: str, key in .obsp for KNN connectivities
    - method: str, one of 'wma', 'knn', or 'both'

    Returns:
    - tmp_adata: AnnData object with new smoothed layers
    """
    tmp_adata = RDR_Xdata.copy()
    suffix = start_layer

    print(f"Smoothing method: {method}")
    print(f"Using start layer: {start_layer}")

    if method in ['wma', 'both']:
        # WMA smoothing
        WMA_mat = tmp_adata.varp[wma_varp]
        start_mtx = tmp_adata.layers[start_layer]
        WMA_smoothed = start_mtx @ WMA_mat

        # Save to adata
        wma_layer_name = f'WMA_smoothed_{suffix}'
        tmp_adata.layers[wma_layer_name] = WMA_smoothed

        print(f"WMA smoothed layer saved: '{wma_layer_name}'")

    if method in ['knn', 'both']:
        # Use WMA smoothed as input if available, else start_layer
        if method == 'knn' and f'WMA_smoothed_{suffix}' in tmp_adata.layers:
            start_mtx = tmp_adata.layers[f'WMA_smoothed_{suffix}']
        elif method == 'both':
            start_mtx = tmp_adata.layers[f'WMA_smoothed_{suffix}']
        else:
            start_mtx = tmp_adata.layers[start_layer]

        connectivities = tmp_adata.obsp[knn_obsp]
        KNN_smoothed = connectivities @ start_mtx

        # Save
        knn_layer_name = f'KNN_smoothed_{suffix}'
        tmp_adata.layers[knn_layer_name] = KNN_smoothed

        print(f"KNN smoothed layer saved: '{knn_layer_name}'")

    return tmp_adata

# denoise function
def denoise(
    RDR_Xdata,
    layer='WMA_smoothed_log_ratio_ab',
    method='dynamic',
    cell_anno_key='spot_anno',
    ref_celltype='normal',
    sd_amplifier=1.5
):
    """
    Denoise an AnnData layer using dynamic thresholding or logistic method.

    Parameters:
    - RDR_Xdata: AnnData object
    - layer: str, input layer name
    - method: 'dynamic' or 'logistic'
    - cell_anno_key: str, obs key for cell annotation
    - ref_celltype: str, value in obs to select reference cells
    - sd_amplifier: float, multiplier for standard deviation (dynamic method)

    Returns:
    - tmp_adata: AnnData object with denoised layer saved as layer+method

    Note:
    - recommended: method='dynamic'
    """
    # Copy AnnData
    tmp_adata = RDR_Xdata.copy()
    data = tmp_adata.layers[layer]
    
    # Select reference cells, can be multiple cell types
    if isinstance(ref_celltype, list):
        ref_mask = tmp_adata.obs[cell_anno_key].isin(ref_celltype)
    else:
        # Single reference cell type
        ref_mask = tmp_adata.obs[cell_anno_key] == ref_celltype
    ref_data = data[ref_mask, :]
    
    # Compute reference mean and std
    mean_ref = np.mean(ref_data)
    std_ref = np.mean(np.std(ref_data, axis=0))
    threshold = sd_amplifier * std_ref

    upper_bound = mean_ref + threshold
    lower_bound = mean_ref - threshold

    # Denoising
    denoised = data.copy()
    if method == 'dynamic':
        mask = (denoised > lower_bound) & (denoised < upper_bound)
        denoised[mask] = mean_ref
    elif method == 'logistic':
        # Logistic squashing around mean_ref
        # You can adjust the steepness (k) as needed
        k = 10 / threshold if threshold != 0 else 1
        denoised = mean_ref + (denoised - mean_ref) / (1 + np.exp(-k * (denoised - mean_ref)))
    else:
        raise ValueError("method must be 'dynamic' or 'logistic'")

    # Save to new layer
    layer_name = layer + '_' + method
    tmp_adata.layers[layer_name] = denoised
    print(f"Denoised layer saved to '{layer_name}', using method '{method}'")

    return tmp_adata, layer_name

# Compute Gaussian probabilities
def compute_gaussian_probabilities(
    RDR_Xdata,
    layer='WMA_smoothed_log_ratio_ab',
    cell_anno_key='spot_anno',
    ref_celltype='normal',
    c_k=np.array([0.5, 1, 1.5])
):
    """
    Computes Gaussian log-likelihoods and probabilities for each cell and gene across 3 states,
    using reference cells to estimate mean and variance.

    Parameters:
    - RDR_Xdata: AnnData object
    - layer: str, layer to use for input data
    - cell_anno_key: str, obs key for cell annotation
    - ref_celltype: str, value in obs to select reference cells
    - c_k: array-like, state multipliers

    Returns:
    - tmp_adata: AnnData object with new layers 'gaussian_loglike', 'gaussian_likelihood', 'gaussian_probabilities'
    """
    # Copy AnnData
    tmp_adata = RDR_Xdata.copy()

    # Step 1: Identify reference cells
    if isinstance(ref_celltype, list):
        # Multiple reference cell types
        ref_cells_mask = tmp_adata.obs[cell_anno_key].isin(ref_celltype)
    else:
        # Single reference cell type
        ref_cells_mask = tmp_adata.obs[cell_anno_key] == ref_celltype

    # Step 2: Extract the relevant data from the specified layer
    if layer not in tmp_adata.layers:
        raise ValueError(f"Layer '{layer}' not found in the AnnData object.")

    ref_data = tmp_adata.layers[layer][ref_cells_mask, :]

    # Convert to dense array if the data is sparse
    # if not isinstance(ref_data, np.ndarray):
    #     ref_data_dense = ref_data.toarray()
    # else:
    #     ref_data_dense = ref_data

    # Calculate the mean and variance for each gene
    mean = np.mean(ref_data, axis=0)
    gene_variance = np.var(ref_data, axis=0)

    # debug 
    print("variance of genes:", gene_variance)

    # add a small value to avoid division by zero when necessary
    gene_variance = np.where(gene_variance < 1e-6, 1e-6, gene_variance)
    print("variance of genes after clip:", gene_variance)

    # Prepare data for all cells
    all_data = tmp_adata.layers[layer]
    # if not isinstance(all_data, np.ndarray):
    #     all_data = all_data.toarray()
    n_cells, n_genes = all_data.shape

    # Function to calculate log-likelihood for a single cell
    def calculate_log_likelihood(cell_data, c_k, mean, gene_variance):
        log_likelihoods = np.zeros((len(mean), len(c_k)))  # (n_genes, 3 states)
        for k, c in enumerate(c_k):
            state_mean = np.log(c) + mean
            log_likelihoods[:, k] = -0.5 * np.log(2 * np.pi * gene_variance) - ((cell_data - state_mean) ** 2 / (2 * gene_variance))
        return log_likelihoods

    # Parallel computation for all cells
    log_likelihoods_all_cells = Parallel(n_jobs=-1)(
        delayed(calculate_log_likelihood)(all_data[cell_idx, :], c_k, mean, gene_variance)
        for cell_idx in range(n_cells)
    )

    # Convert results to a 3D array (n_cells, n_genes, 3 states)
    log_likelihoods_all_cells = np.array(log_likelihoods_all_cells)
    tmp_adata.layers['gaussian_loglike'] = log_likelihoods_all_cells

    # Calculate likelihoods
    gaussian_likelihood = np.exp(log_likelihoods_all_cells)
    tmp_adata.layers['gaussian_likelihood'] = gaussian_likelihood

    # Normalize to get probabilities
    state_probabilities = gaussian_likelihood / gaussian_likelihood.sum(axis=2, keepdims=True)
    tmp_adata.layers['gaussian_probabilities'] = state_probabilities

    layer_name = 'gaussian_probabilities'

    print(f"Computed Gaussian probabilities saved to layer '{layer_name}'")
    print(f"Computed Gaussian log-likelihoods saved to layer 'gaussian_loglike'")
    print(f"Computed Gaussian likelihoods saved to layer 'gaussian_likelihood'")

    return tmp_adata, layer_name

# WMA smoothing and normalization on probabilities
# probably unnecessary, but keep it for now
def smooth_and_normalize(RDR_Xdata, WMA_mat_key='WMA_connect', layer='post_HMM_prob_2'):
    """
    Smooth each state in RDR_Xdata.layers['post_HMM_prob_2'] using WMA_mat from tmp_adata.varp['WMA_connect'],
    and normalize so that the sum of the 3 states equals 1 for each cell and gene.

    Parameters:
    - RDR_Xdata: AnnData object containing the posterior probabilities and smoothing matrix.
    - WMA_mat_key: Key in `varp` for the weighted moving average matrix (default: 'WMA_connect').
    - post_HMM_prob_key: Key in `layers` for the posterior probabilities (default: 'post_HMM_prob_2').

    Returns:
    - RDR_Xdata: Updated AnnData object with smoothed and normalized probabilities.
    """
    # Step 1: Extract the posterior probabilities and WMA matrix
    if layer not in RDR_Xdata.layers:
        raise ValueError(f"Layer '{layer}' not found in RDR_Xdata.layers.")
    if WMA_mat_key not in RDR_Xdata.varp:
        raise ValueError(f"Matrix '{WMA_mat_key}' not found in RDR_Xdata.varp.")

    post_HMM_prob = RDR_Xdata.layers[layer]  # Shape: (n_cells, n_genes, 3 states)
    WMA_mat = RDR_Xdata.varp[WMA_mat_key]  # Shape: (n_genes, n_genes)


    # Step 2: Smooth each state using WMA_mat
    n_states = post_HMM_prob.shape[2]
    smoothed_prob = np.zeros_like(post_HMM_prob)  # Shape: (n_cells, n_genes, 3 states)

    for state_idx in range(n_states):
        # Extract the matrix for the current state
        start_mtx = post_HMM_prob[:, :, state_idx]  # Shape: (n_cells, n_genes)

        # Smooth using WMA_mat
        smoothed_prob[:, :, state_idx] = start_mtx @ WMA_mat  # Shape: (n_cells, n_genes)

        print("finished for state", state_idx)

    # Step 3: Normalize so that the sum of the 3 states equals 1
    smoothed_prob_sum = smoothed_prob.sum(axis=2, keepdims=True)  # Sum across states
    normalized_prob = smoothed_prob / smoothed_prob_sum  # Normalize each state

    # Step 4: Store the normalized probabilities back in the AnnData object
    tmp_adata = RDR_Xdata.copy()
    layer_name = layer + '_WMA'
    tmp_adata.layers[layer_name] = normalized_prob

    return tmp_adata, layer_name


# Smooth an AnnData layer using Gaussian, median, or minimum segment smoothing.
def smooth_anndata_layer(
    RDR_Xdata,
    layer='post_HMM_prob_2_WMA',
    method='gauss',  # 'gauss', 'median', or 'minseg'
    gauss_window=3,
    median_window=3,
    min_segment_size=5
):
    """
    Perform a selected smoothing method on an AnnData layer.

    Parameters:
    - RDR_Xdata: AnnData object
    - layer: str, layer to smooth
    - method: str, one of 'gauss', 'median', 'minseg'
    - gauss_window: int, window size for Gaussian smoothing
    - median_window: int, window size for median filtering
    - min_segment_size: int, minimum segment size to enforce

    Returns:
    - tmp_adata: AnnData object with the new smoothed layer
    """
    tmp_adata = RDR_Xdata.copy()
    data = tmp_adata.layers[layer]
    if not isinstance(data, np.ndarray):
        data = data.toarray()
    n_cells, n_genes, n_states = data.shape

    print(f"Smoothing layer '{layer}' using method '{method}'")

    if method == 'gauss':
        gauss_smoothed = np.empty_like(data)
        for cell in range(n_cells):
            for state in range(n_states):
                gauss_smoothed[cell, :, state] = gaussian_filter1d(
                    data[cell, :, state], sigma=gauss_window/2, mode='nearest'
                )
        layer_name = layer + '_gauss'
        tmp_adata.layers[layer_name] = gauss_smoothed
        print(f"Smoothed layer saved: '{layer_name}'")

    elif method == 'median':
        median_smoothed = np.empty_like(data)
        for cell in range(n_cells):
            for state in range(n_states):
                median_smoothed[cell, :, state] = median_filter(
                    data[cell, :, state], size=median_window, mode='nearest'
                )
        layer_name = layer + '_median'
        tmp_adata.layers[layer_name] = median_smoothed
        print(f"Smoothed layer saved: '{layer_name}'")

    elif method == 'minseg':
        assigned_states = np.argmax(data, axis=2)  # shape: (n_cells, n_genes)

        def enforce_min_segment(arr, min_size):
            arr = arr.copy()
            start = 0
            while start < len(arr):
                curr_state = arr[start]
                end = start
                while end + 1 < len(arr) and arr[end + 1] == curr_state:
                    end += 1
                segment_len = end - start + 1
                if segment_len < min_size:
                    # Merge with previous if possible, else with next
                    if start > 0:
                        arr[start:end+1] = arr[start-1]
                    elif end + 1 < len(arr):
                        arr[start:end+1] = arr[end+1]
                start = end + 1
            return arr

        minseg_states = np.empty_like(assigned_states)
        for cell in range(n_cells):
            minseg_states[cell, :] = enforce_min_segment(assigned_states[cell, :], min_segment_size)
        # Convert to one-hot encoding (probabilities)
        minseg_smoothed = np.zeros((n_cells, n_genes, n_states), dtype=data.dtype)
        for state in range(n_states):
            minseg_smoothed[:, :, state] = (minseg_states == state)
        layer_name = layer + '_minseg'
        tmp_adata.layers[layer_name] = minseg_smoothed
        print(f"Smoothed layer saved: '{layer_name}'")

    else:
        raise ValueError("method must be one of 'gauss', 'median', or 'minseg'")

    return tmp_adata, layer_name

# Low-rank approximation using PCA
# for denoising
def low_rank_approximation(
    RDR_Xdata,
    layer='post_HMM_prob_2_WMA',
    n_components=10
):
    """
    Perform low-rank approximation (PCA) for each state in a probability tensor.

    Parameters:
    - RDR_Xdata: AnnData object
    - layer: str, layer to smooth (shape: n_cells x n_genes x n_states)
    - n_components: int, number of principal components to keep

    Returns:
    - tmp_adata: AnnData object with new layer '<layer>_lowrank'
    """
    tmp_adata = RDR_Xdata.copy()
    data = tmp_adata.layers[layer]
    if not isinstance(data, np.ndarray):
        data = data.toarray()
    n_cells, n_genes, n_states = data.shape

    print(f"Performing low-rank approximation on layer '{layer}' with {n_components} components")

    lowrank = np.empty_like(data)
    for state in range(n_states):
        X = data[:, :, state]
        pca = PCA(n_components=min(n_components, min(X.shape)), svd_solver='full')
        X_reduced = pca.fit_transform(X)
        X_approx = pca.inverse_transform(X_reduced)
        lowrank[:, :, state] = X_approx

    layer_name = layer + '_lowrank'
    tmp_adata.layers[layer_name] = lowrank
    print(f"Low-rank approximation layer saved: '{layer_name}'")
    return tmp_adata, layer_name


# KNN smoothing of a probability layer
def knn_smooth_prob_layer(
    RDR_Xdata,
    prob_layer='post_HMM_prob_2_lowrank',
    k=10,
    knn_obsp_name='knn_connectivities',
    smoothed_layer_name='KNN_smoothed_post_HMM_prob_2_lowrank'
):
    """
    Construct a KNN connectivity matrix using a probability layer and smooth the matrix.
    Normalizes so the sum of probabilities across 3 states is 1 for each cell/gene.

    Parameters:
    - RDR_Xdata: AnnData object
    - prob_layer: str, probability layer to use (shape: [cell, gene, 3_states])
    - k: int, number of neighbors for KNN
    - knn_obsp_name: str, key to save the KNN connectivities in .obsp
    - smoothed_layer_name: str, key to save the smoothed layer in .layers

    Returns:
    - tmp_adata: AnnData object with KNN connectivities and smoothed probability layer
    """
    tmp_adata = RDR_Xdata.copy()
    prob = tmp_adata.layers[prob_layer]
    if not isinstance(prob, np.ndarray):
        prob = prob.toarray()
    n_cells, n_genes, n_states = prob.shape

    # Flatten last two axes for KNN (cell x (gene*3))
    X = prob.reshape(n_cells, -1)

    # Compute KNN connectivities
    nbrs = NearestNeighbors(n_neighbors=k, metric='euclidean').fit(X)
    connectivities = nbrs.kneighbors_graph(X, mode='connectivity')
    tmp_adata.obsp[knn_obsp_name] = connectivities

    # Smooth each state separately, then stack
    smoothed = np.empty_like(prob)
    for state in range(n_states):
        smoothed[:, :, state] = connectivities @ prob[:, :, state]

    # Normalize so sum of 3 states = 1 for each cell/gene
    smoothed_sum = smoothed.sum(axis=2, keepdims=True)
    smoothed_sum[smoothed_sum == 0] = 1  # Avoid division by zero
    smoothed = smoothed / smoothed_sum

    tmp_adata.layers[smoothed_layer_name] = smoothed

    return tmp_adata, smoothed_layer_name

def run_RDR_gaussian_plot(RDR_adata,
            dataset_name,
            plot_cell_anno_key,
            log_ratio_layer = "log_ratio_ab",
            prob_layer = "post_HMM_prob_2_gauss",
            set_figtitle = True,
            rdr_plot_vmin = -0.7,
            rdr_plot_vmax = 0.7,
            out_dir = None,
            **kwargs):
    """
    modified from Rongting's run_RDR_plot function.
    support user-defined log_ratio_layer and prob_layer.
    """
    ## Result output prepare
    if out_dir is None:
        cwd = os.getcwd()
        out_dir = cwd + "/XCLONE_OUT/"
    
    out_plot_dir = str(out_dir) + "/plot/"
    xclone.al.dir_make(out_plot_dir)

    # default:XClone
    fig_title = ""
    rdr_smooth_fig = out_plot_dir + dataset_name + "_RDR_smooth.png"
    rdr_final_fig = out_plot_dir + dataset_name + "_RDR_CNV.png"

    sub_logger = get_logger("RDR plot module")
    sub_logger.info("RDR plot module started")
    start_time = datetime.now(timezone.utc)
    
    if set_figtitle:
        fig_title = dataset_name + " RDR_smooth (log scale ratio)"

    xclone.pl.RDR_smooth_visualization(RDR_adata, 
                                       Xlayer = log_ratio_layer, 
                                       cell_anno_key = plot_cell_anno_key,
                                       change_colorbar = False,
                                       vmin = rdr_plot_vmin, vmax = rdr_plot_vmax,
                                       title = fig_title,
                                       save_file = True, 
                                       out_file = rdr_smooth_fig,
                                       **kwargs)
    if set_figtitle:
        fig_title = dataset_name + " RDR module"
    
    xclone.pl.CNV_visualization(RDR_adata, 
                                Xlayer = prob_layer,
                                states_weight = np.array([1,2,3]), 
                                weights = True, 
                                cell_anno_key = plot_cell_anno_key, 
                                title = fig_title,
                                save_file = True, 
                                out_file = rdr_final_fig,
                                **kwargs)
    
    end_time = datetime.now(timezone.utc)
    time_passed = end_time - start_time
    sub_logger.info("RDR plot module finished (%d seconds)" % (time_passed.total_seconds()))
    return None