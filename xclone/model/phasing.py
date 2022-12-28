# This file contains methods for phasing alleles on SNPs or genomic regions by 
# spatial correlations on allelic ratio
# Author: Yuanhua Huang
# Date: 12/11/2020

import numpy as np
from scipy.stats import entropy
from scipy.optimize import linear_sum_assignment
from scipy.special import logsumexp, digamma, betaln, binom

from .base_utils import normalize, loglik_amplify


def Local_Phasing(AD, DP, min_iter=10, max_iter=1000, epsilon_conv=1e-2,
                  init_mode='warm', verbose=False):
    """
    Phase the small blocks into a medium sized bin by assuming the allelic
    ratio is the same for all blocks. This is equavilent to a binary clustering,
    whose likelihood can be maximised by an EM alogrithm

    # Note that AD DP should be in in Compressed Sparse Column format.
    # or other sparse format.
    
    Parameters
    ----------
    AD : sparse matrix of integers
        Read counts for ALT allele in N blocks and M cells
    DP : sparse matrix of integers
        Read counts for REF allele in N blocks and M cells
    
    Returns
    -------
    ad_sum, dp_sum, Z, thetas, logLik
    """
    N, M = AD.shape
    BD = DP - AD
    
    ## Initialization matters!!!
    if init_mode == 'warm':
        # Option 1: warm initialization (usually good?)
        Z = np.zeros((N, 2))
        Z[:, 0] = (AD.sum(1) / DP.sum(1)).reshape(-1)
        Z[:, 1] = 1 - Z[:, 0]
    elif init_mode == 'current':
        # Option 2: initializing with no flipping (can return poor local optimal)
        Z = np.zeros((N, 2))
        Z[:, 0] = 1 ## high prob means no flipping
    else:
        # Option 3: random initialization (may need large number of trials)
        Z = np.random.rand(N, 2)
        Z[:, 1] = 1 - Z[:, 0]
    
    # allele ratio parameters
    # thetas = np.array((AD.T * Z + BD.T * (1 - Z)) / (DP.T.sum(1)))
    thetas = np.array((AD.T @ Z + BD.T @ (1 - Z)) / (DP.T.sum(1))) 
    # works for both sparse matrix and nunmpy array
    
    # thetas = np.zeros((M, 2))
    # thetas[:, 0:1] = (AD.T * Z[:, 0:1] + BD.T * Z[:, 1:2]) / (DP.sum(0).T)
    # thetas[:, 1:2] = 1 - thetas[:, 0:1]
    
    # likelihood
    _logLik_mat = (AD * np.log(thetas) + BD * np.log(1 - thetas))
    _logLik_new = np.sum(logsumexp(_logLik_mat, axis=1))
    
    for it in range(max_iter):
        _logLik_old = _logLik_new + 0.0
        
        # E step: calculate the expecation
        Z = normalize(np.exp(loglik_amplify(np.array(_logLik_mat))))
        
        # M step: maximise the likihood over thetas
        thetas = np.array((AD.T * Z + BD.T * (1 - Z)) / (DP.T.sum(1)))
        
        # Likelihood
        _logLik_mat = (AD * np.log(thetas) + BD * np.log(1 - thetas))
        _logLik_new = np.sum(logsumexp(_logLik_mat, axis=1))
        
        # print(it, _logLik_old, _logLik_new)
        # convergence
        if it >= min_iter and _logLik_new - _logLik_old < epsilon_conv:
            if verbose:
                print("EM finished in %d iterations with logLik %.4e" 
                      %(it, _logLik_new))
            break
        elif _logLik_new < _logLik_old:
            if verbose:
                print("Warning: likelihood decreases in EM algorithm!")
                print(it, _logLik_old, print(it, _logLik_new))
    # if soft_phasing:
    ## soft phasing bins counts: `ad_sum`
    ad_sum = Z.T * AD + (1 - Z.T) * BD 
    # else:
    Z_argmax = np.argmax(Z, axis = 1)
    Z_hard_assign = np.vstack((1- Z_argmax, Z_argmax)).T
    ## hard phasing bins counts:`ad_sum1``
    ad_sum1 = Z_hard_assign.T * AD + (1 - Z_hard_assign.T) * BD 
    dp_sum = DP.sum(axis=0)
    return ad_sum, ad_sum1, dp_sum, Z, thetas, _logLik_new
        
    
    
def Global_Phasing(thetas, step_size=1, n_prefix=1):
    """
    Phase the medium-sized bins into large scale regions, e.g., whole 
    chromosome. The assumption here is that the allelic ratio is similar within
    a neibougdhood. A dynamic programming for recercive optimization is used for
    this.
    """
    
    def _Dyna_Programming(thetas):
        """recursive optimization
        """
        is_flips = np.zeros(thetas.shape[0], bool)
        distances = np.zeros(thetas.shape[0])
        thetas_new = thetas + 0
        
        if thetas.shape[0] == 1:
            # print("Bottom:", thetas.shape[0])
            return is_flips, distances, thetas_new
        else:
            is_filps_pre, distances_pre, thetas_pre = _Dyna_Programming(thetas[:-1, :])
            
            theta_tmp0 = 0 + thetas[-1, :]
            theta_tmp1 = 1 - thetas[-1, :]
            
            _dist_tmp0 = (theta_tmp0 - thetas_pre[-1, :])**2
            _dist_tmp1 = (theta_tmp1 - thetas_pre[-1, :])**2
            _dist_tmp0 = np.sum(_dist_tmp0[_dist_tmp0 >=0 ])
            _dist_tmp1 = np.sum(_dist_tmp1[_dist_tmp1 >=0 ])

            # print(thetas.shape[0], _dist_tmp0, _dist_tmp1)
            
            is_flips[:-1]  = is_filps_pre
            distances[:-1] = distances_pre
            if _dist_tmp0 < _dist_tmp1:
                is_flips[-1] = False
                distances[-1]  = _dist_tmp0 + 0
                thetas_new[-1, :] = theta_tmp0 + 0
            else:
                is_flips[-1] = True
                distances[-1]  = _dist_tmp1 + 0
                thetas_new[-1, :] = theta_tmp1 + 0
            
            return is_flips, distances, thetas_new
    
    ## Return Dynamical programming results
    return _Dyna_Programming(thetas)
