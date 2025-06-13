"""HMM for XClone RDR Gaussian model.
"""

# Author: Jiamu James Qiao, logdotexp and HMM by Yuanhua Huang
# Date: 2025-06-10
# update: 

import numpy as np
from joblib import Parallel, delayed


# HMM
# HMM by yh
import numpy as np
from scipy.special import logsumexp

def logdotexp(A_log, B_log):
    """
    Extension of logsumexp to logdotexp
    Note, this is not memory efficient and will have peak size: (n, k, m) 
    for (n, k) * (k, m)
    
    Examples
    --------
    >>> A = np.arange(1, 7).reshape(3, 2)
    >>> B = np.arange(1, 9).reshape(2, 4)
    >>> print(np.log(A @ B))
    >>> print(np.log(np.dot(A, B)))
    >>> print(logdotexp(np.log(A), np.log(B)))
    """
    if len(A_log.shape) == 1:
        A_log = A_log.reshape(1, -1)
    if len(B_log.shape) == 1:
        B_log = B_log.reshape(-1, 1)
        
    AB_log = logsumexp(
        np.expand_dims(A_log, 2) + np.expand_dims(B_log, 0),
        axis=1
    )
    
    return AB_log

def HMM(emm_p_log, pi=None, A=None, diag_prob=0.99):
    """
    Hidden Markov Model
    -------------------
    Forward-Backward algorithm for calculating the 
    posterior of the state assignment via HMM
    """
    n_index = emm_p_log.shape[0]
    n_state = emm_p_log.shape[1]

    # prior distribution
    if pi is None:
        pi = np.ones(n_state) / n_state

    # transition matrix
    if A is None:
        diag_prob_adj = (diag_prob - 1 / n_state) * n_state / (n_state - 1)
        A = np.eye(n_state) * diag_prob_adj + (1 - diag_prob_adj) / n_state
    A_log = np.log(A)
    
    # Forward-Backward algorithm for updating posterior
    fw_p_log = emm_p_log.copy()
    bw_p_log = np.zeros_like(emm_p_log)

    # Forward part of the algorithm
    fw_p_log[0, :] += np.log(pi)
    for i in range(1, n_index):
        fw_p_log[i, :] += logdotexp(fw_p_log[i - 1, :], A_log).reshape(-1)

    # Backward part of the algorithm
    for i in range(n_index - 1, 0, -1):
        bw_p_log[i - 1, :] = logdotexp(
            bw_p_log[i, :] + emm_p_log[i, :], A_log).reshape(-1)

    # Update posterior of state assignment
    z_post_log = fw_p_log + bw_p_log
    z_post_log -= logsumexp(z_post_log, axis=1, keepdims=True)
    z_post = np.exp(z_post_log)
    z_post = z_post / np.sum(z_post, axis=1, keepdims=True)
    
    return z_post

from joblib import Parallel, delayed
def calculate_z_post_parallel(RDR_Xdata, start_prob, trans_prob, layer='gaussian_loglike', n_jobs=-1):
    """
    Calculate z_post for each cell along each chromosome using HMM with parallel computing.
    parallel computing: require joblib
    """
    # Extract emission probabilities
    if layer not in RDR_Xdata.layers:
        raise ValueError("Layer not found in the AnnData object.")
    emm_p_log = RDR_Xdata.layers[layer]  # Shape: (n_cells, n_genes, n_states)
    print("using layer", layer, "as emm_p_log")

    # Ensure dense matrix if sparse
    if not isinstance(emm_p_log, np.ndarray):
        emm_p_log = emm_p_log.toarray()
    # 3min 10s if dense

    # Order genes along chromosomes
    chr_order = RDR_Xdata.var.sort_values(by=['chr', 'start', 'stop']).index
    RDR_Xdata = RDR_Xdata[:, chr_order]

    # Initialize results
    n_cells = emm_p_log.shape[0]
    n_genes = emm_p_log.shape[1]
    n_states = emm_p_log.shape[2]
    z_post_all = np.zeros((n_cells, n_genes, n_states))

    # Function to process a single chromosome
    def process_chromosome(chr_name):
        chr_mask = RDR_Xdata.var['chr'] == chr_name
        chr_emm_p_log = emm_p_log[:, chr_mask, :]  # Subset emission probabilities for this chromosome

        # Calculate z_post for each cell along this chromosome
        z_post_chr = np.zeros((n_cells, chr_mask.sum(), n_states))
        for cell_idx in range(n_cells):
            z_post_chr[cell_idx, :, :] = HMM(chr_emm_p_log[cell_idx, :, :], start_prob, trans_prob)
        return chr_name, z_post_chr

    # Parallel processing for each chromosome
    results = Parallel(n_jobs=n_jobs)(
        delayed(process_chromosome)(chr_name) for chr_name in RDR_Xdata.var['chr'].unique()
    )

    # Combine results
    for chr_name, z_post_chr in results:
        chr_mask = RDR_Xdata.var['chr'] == chr_name
        z_post_all[:, chr_mask, :] = z_post_chr

    return z_post_all