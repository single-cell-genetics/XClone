# Utility functions as heritages from vireoSNP

import numpy as np
import scipy as sp
from scipy.stats import entropy
from scipy.optimize import linear_sum_assignment
from scipy.special import logsumexp, digamma, betaln, binom


def normalize(X, axis=-1):
    """
    Normalization of tensor with sum to 1.
    X should be numpy.array, otherwise will be transformed to it.
    
    Example
    -------
    X = np.random.rand(3, 5, 8)
    tensor_normalize(X, axis=1)
    """
    shape2 = list(X.shape)
    shape2[axis] = 1
    
    if type(X) == np.matrix:
        X = X.A
    if sp.sparse.issparse(X):
        X = X.A
    
    X_sum = np.sum(X, axis=axis, keepdims=True)
    return X / X_sum


def loglik_amplify(X, axis=-1):
    """
    Amplify the log likelihood matrix by subtract the maximum.
    X should be numpy.array, otherwise will be transformed to it.
    
    Example
    -------
    X = np.random.rand(3, 5, 8)
    loglik_amplify(X, axis=1)
    """
    shape2 = list(X.shape)
    shape2[axis] = 1
    
    if type(X) == np.matrix:
        X = X.A
    
    X_max = np.max(X, axis=axis, keepdims=True)
    return X - X_max

def cal_log_lik(emm_prob_log, posterior_mtx_log):
    """
    Function:

    INPUT:
    -------
    c * g * states
    c * g * states
    c can be celltype based | can be cell based
    """
    emm_dim = emm_prob_log.ndim
    pos_dim = posterior_mtx_log.ndim
    ## check the dims of the two numpy array
    if emm_dim != 3 | pos_dim != 3:
        print("[XClone]: cal_log_lik dim mismatch.")
        return None

    log_lik_ = logsumexp(emm_prob_log + posterior_mtx_log, axis = 2).sum(axis=1).sum()
    return log_lik_