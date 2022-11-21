"""base functions for XClone BAF and RDR module (Variational inference)
"""
# Author: Rongting Huang
# Date: 2021/01/26
# update: 2021/10/06

import numpy as np
from scipy.stats import entropy
from scipy.optimize import linear_sum_assignment
from scipy.special import logsumexp, digamma, betaln, binom

from scipy.special import gammaln
from numpy import log as ln
from scipy.special import factorial


def logfact(n):
    """
    for np.array
    calculate ln(n!)
    """
    return ln(factorial(n))

def gamma_entropy(X, X_prior=None):
    """
    Get the entropy for gamma distributions. If X_prior is not None, return the
    Kullback-Leibler divergence
    See: https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.entropy.html
    Example
    -------
    theta_shapes1 = np.array([[0.3, 29.7], [3, 3], [29.7, 0.3]])
    theta_shapes2 = np.array([[364, 24197], [5886, 7475], [6075, 397]])
    gamma_entropy(theta_shapes2)
    gamma_entropy(theta_shapes2, theta_shapes1)
    """
    RV1 = 0
    if X_prior is None:
        X_prior = X.copy()
    else:
        RV1 = (X[:,0] * ln(X[:,1]) - gammaln(X[:,0]) +
        (X[:, 0] - 1) * (digamma(X[:,0]) - ln(X[:,1])) +
        X[:,0]
        )

    RV2 = (X_prior[:,0] * ln(X_prior[:,1]) - gammaln(X_prior[:,0]) +
    (X_prior[:, 0] - 1) * (digamma(X[:,0]) - ln(X[:,1])) +
    X_prior[:,1] * (X[:,0] / X[:,1])
    )

    return np.sum(RV1) - np.sum(RV2)


def gamma_entropy1(alpha, beta, alpha_prior=None, beta_prior=None):
    """
    updated version for new element-wise calculation.
    Get the entropy for gamma distributions. If X_prior is not None, return the
    Kullback-Leibler divergence
    See: https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.entropy.html
    Example
    -------

    """
    RV1 = 0
    if alpha_prior is None:
        alpha_prior = alpha.copy()
        beta_prior = beta.copy()
    else:
        RV1 = alpha * ln(beta) - gammaln(alpha) + (alpha - 1) * (digamma(alpha) - ln(beta)) + alpha

    RV2 = alpha_prior * ln(beta_prior) - gammaln(alpha_prior) + (alpha_prior - 1) * (digamma(alpha) - ln(beta)) + beta_prior * (alpha/beta)

    return np.sum(RV1) - np.sum(RV2)


# from vireo_base
def get_binom_coeff(AD, DP, max_val=700, is_log=True):
    """Get the binomial coefficients
    """
    # Since binom can't give log value, the maximum value in 64bit is
    # around e**700, close to binom(1000, 500)
    idx = DP > 0
    _AD = AD[idx].astype(np.int64)
    _DP = DP[idx].astype(np.int64)

    binom_coeff = np.log(binom(_DP, _AD))
    binom_coeff[binom_coeff > max_val] = max_val
    binom_coeff = binom_coeff.astype(np.float32)

    return binom_coeff

# from vireo_base
def normalize(X, axis=-1):
    """
    Normalization of tensor with sum to 1.

    Example
    -------
    X = np.random.rand(3, 5, 8)
    tensor_normalize(X, axis=1)
    """
    # shape2 = list(X.shape)
    # shape2[axis] = 1
    X_sum = np.sum(X, axis=axis, keepdims=True)
    return X / X_sum

# from vireo_base
def tensor_normalize(X, axis=1):
    return normalize(X, axis)

# from vireo_base
def loglik_amplify(X, axis=-1):
    """
    Amplify the log likelihood matrix by subtract the maximum.

    Example
    -------
    X = np.random.rand(3, 5, 8)
    loglik_amplify(X, axis=1)
    """
    shape2 = list(X.shape)
    shape2[axis] = 1
    X_max = np.max(X, axis=axis, keepdims=True)
    return X - X_max

# from vireo_base
def beta_entropy(X, X_prior=None):
    """
    Get the entropy for beta distributions. If X_prior is not None, return the
    Kullback-Leibler divergence
    See: https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.entropy.html
    Example
    -------
    theta_shapes1 = np.array([[0.3, 29.7], [3, 3], [29.7, 0.3]])
    theta_shapes2 = np.array([[364, 24197], [5886, 7475], [6075, 397]])
    beta_entropy(theta_shapes2)
    beta_entropy(theta_shapes2, theta_shapes1)
    """
    RV1 = 0
    if X_prior is None:
        X_prior = X.copy()
    else:
        RV1 = (- betaln(X[:, 0], X[:, 1]) +
               (X[:, 0] - 1) * digamma(X[:, 0]) +
               (X[:, 1] - 1) * digamma(X[:, 1]) -
               (X.sum(axis=1) - 2) * digamma(X.sum(axis=1)))

    RV2 = (- betaln(X_prior[:, 0], X_prior[:, 1]) +
           (X_prior[:, 0] - 1) * digamma(X[:, 0]) +
           (X_prior[:, 1] - 1) * digamma(X[:, 1]) -
           (X_prior.sum(axis=1) - 2) * digamma(X.sum(axis=1)))

    return np.sum(RV1) - np.sum(RV2)
