"""
Variational mixture models for XClone subclonal clustering (RDR + joint).

BAF binomial VB lives in ``_vb_binomial.py``.
"""

from __future__ import annotations

import numpy as np
from scipy.sparse import csc_matrix

from .base_utils import loglik_amplify, normalize
from ._vb_binomial import get_binom_mixture_vb_class


def to_dense(x) -> np.ndarray:
    if hasattr(x, "toarray"):
        return np.asarray(x.toarray(), dtype=float)
    return np.asarray(x, dtype=float)


def impute_gene_medians(Y: np.ndarray) -> np.ndarray:
    Y = np.asarray(Y, dtype=float)
    med = np.nanmedian(Y, axis=0)
    med = np.where(np.isfinite(med), med, 0.0)
    out = Y.copy()
    bad = ~np.isfinite(out)
    if bad.any():
        out[bad] = med[np.where(bad)[1]]
    return out


class GaussianMixtureVB:
    """Diagonal Gaussian mixture VB on RDR log-ratio matrix (cells x genes)."""

    def __init__(
        self,
        n_cell: int,
        n_var: int,
        n_donor: int,
        *,
        min_inv_sigma2: float = 1e-2,
        max_inv_sigma2: float = 1e4,
    ):
        self.n_cell = n_cell
        self.n_var = n_var
        self.n_donor = n_donor
        self.min_inv_sigma2 = min_inv_sigma2
        self.max_inv_sigma2 = max_inv_sigma2
        self.set_prior()
        self.set_initial()
        self.ELBO_iters = np.array([])

    def set_prior(self, ID_prior=None):
        self.ID_prior = (
            normalize(np.asarray(ID_prior, dtype=float), axis=-1)
            if ID_prior is not None
            else normalize(np.ones((self.n_cell, self.n_donor)))
        )

    def set_initial(self, mu_init=None, inv_sigma2_init=None, ID_prob_init=None):
        self.mu = (
            np.asarray(mu_init, dtype=float)
            if mu_init is not None
            else np.zeros((self.n_var, self.n_donor))
        )
        self.inv_sigma2 = (
            np.clip(np.asarray(inv_sigma2_init, dtype=float), self.min_inv_sigma2, self.max_inv_sigma2)
            if inv_sigma2_init is not None
            else np.ones((self.n_var, self.n_donor))
        )
        self.ID_prob = (
            normalize(ID_prob_init, axis=1)
            if ID_prob_init is not None
            else normalize(np.random.rand(self.n_cell, self.n_donor))
        )

    def get_E_logLik(self, Y: np.ndarray) -> np.ndarray:
        is2, mu = self.inv_sigma2, self.mu
        t3 = np.sum(mu * mu * is2, axis=0)
        quad = (
            np.einsum("cv,vk->ck", Y * Y, is2)
            - 2.0 * np.einsum("cv,vk->ck", Y, mu * is2)
            + t3
        )
        log_det = np.sum(np.log(2 * np.pi / is2), axis=0)
        return -0.5 * quad - 0.5 * log_det

    def update_profiles(self, Y: np.ndarray):
        w = self.ID_prob
        sw = w.sum(axis=0) + 1e-12
        self.mu = (Y.T @ w) / sw
        var = np.empty_like(self.mu)
        for k in range(self.n_donor):
            diff = Y - self.mu[:, k]
            var[:, k] = (w[:, k] @ (diff * diff)) / sw[k] + 1e-4
        self.inv_sigma2 = np.clip(1.0 / var, self.min_inv_sigma2, self.max_inv_sigma2)

    def update_ID_prob(self, Y=None, logLik_ID=None):
        if logLik_ID is None:
            logLik_ID = self.get_E_logLik(Y)
        self.ID_prob = normalize(
            np.exp(loglik_amplify(logLik_ID + np.log(self.ID_prior + 1e-300)))
        )

    def get_ELBO(self, logLik_ID: np.ndarray) -> float:
        from scipy.stats import entropy

        return float(np.sum(logLik_ID * self.ID_prob) - np.sum(entropy(self.ID_prob, self.ID_prior, axis=-1)))

    def _fit_vb(self, Y, max_iter=200, min_iter=20, epsilon_conv=1e-2, verbose=False):
        elbo = np.zeros(max_iter)
        for it in range(max_iter):
            self.update_profiles(Y)
            ll = self.get_E_logLik(Y)
            self.update_ID_prob(logLik_ID=ll)
            elbo[it] = self.get_ELBO(ll)
            if it > min_iter and elbo[it] - elbo[it - 1] < epsilon_conv:
                break
        if verbose and it >= max_iter - 1:
            print("Warning: RDR VB did not converge!")
        self.ELBO_iters = np.append(self.ELBO_iters, elbo[: it + 1])

    def fit(self, Y, *, n_init=10, max_iter=200, max_iter_pre=100, random_seed=None, verbose=False, **kwargs):
        Y = to_dense(Y)
        if Y.shape[0] == self.n_var and Y.shape[1] == self.n_cell:
            Y = Y.T
        if Y.shape != (self.n_cell, self.n_var):
            raise ValueError(f"Y must be ({self.n_cell}, {self.n_var}), got {Y.shape}")
        Y = impute_gene_medians(Y)
        if random_seed is not None:
            np.random.seed(random_seed)
        global_mu = np.nanmedian(Y, axis=0)
        self.ELBO_inits = []
        best = None
        for _ in range(n_init):
            self.set_initial(
                mu_init=global_mu[:, np.newaxis] + np.random.randn(self.n_var, self.n_donor) * 0.05
            )
            self._fit_vb(Y, max_iter=max_iter_pre, verbose=verbose, **kwargs)
            self.ELBO_inits.append(float(self.ELBO_iters[-1]))
            if best is None or self.ELBO_iters[-1] >= max(self.ELBO_inits[:-1]):
                best = (self.ID_prob.copy(), self.mu.copy(), self.inv_sigma2.copy(), self.ELBO_iters.copy())
        self.ID_prob, self.mu, self.inv_sigma2 = best[0], best[1], best[2]
        self.ELBO_iters = best[3]
        self._fit_vb(Y, max_iter=max_iter, verbose=verbose, **kwargs)
        return self


class JointBinomGaussianVB:
    """Shared clone assignment: BAF Binomial + RDR Gaussian emissions."""

    def __init__(self, n_cell, n_var_baf, n_var_rdr, n_donor, *, rdr_weight=1.0, fix_beta_sum=False):
        BinomMixtureVB = get_binom_mixture_vb_class()
        self.n_cell = n_cell
        self.n_var_baf = n_var_baf
        self.n_var_rdr = n_var_rdr
        self.n_donor = n_donor
        self.rdr_weight = rdr_weight
        self._baf = BinomMixtureVB(n_var=n_var_baf, n_cell=n_cell, n_donor=n_donor, fix_beta_sum=fix_beta_sum)
        self._rdr = GaussianMixtureVB(n_cell=n_cell, n_var=n_var_rdr, n_donor=n_donor)
        self.ID_prob = self._baf.ID_prob
        self.ELBO_iters = np.array([])

    def _sync_id(self):
        self._baf.ID_prob = self.ID_prob
        self._rdr.ID_prob = self.ID_prob

    def get_loglik_ID(self, AD, DP, Y):
        self._sync_id()
        return self._baf.get_E_logLik(AD, DP) + self.rdr_weight * self._rdr.get_E_logLik(Y)

    def _fit_joint(self, AD, DP, Y, max_iter=200, min_iter=20, epsilon_conv=1e-2, verbose=False):
        elbo = np.zeros(max_iter)
        for it in range(max_iter):
            self._sync_id()
            self._baf.update_theta_size(AD, DP)
            self._rdr.update_profiles(Y)
            ll = self.get_loglik_ID(AD, DP, Y)
            self.ID_prob = normalize(np.exp(loglik_amplify(ll + np.log(self._baf.ID_prior + 1e-300))))
            self._sync_id()
            elbo[it] = float(np.sum(ll * self.ID_prob))
            if it > min_iter and elbo[it] - elbo[it - 1] < epsilon_conv:
                break
        if verbose and it >= max_iter - 1:
            print("Warning: joint VB did not converge!")
        self.ELBO_iters = np.append(self.ELBO_iters, elbo[: it + 1])

    def fit(self, AD, DP, Y, *, n_init=10, max_iter=200, max_iter_pre=100, random_seed=None, verbose=False):
        AD, DP = to_dense(AD), to_dense(DP)
        if AD.shape[0] == self.n_var_baf and AD.shape[1] == self.n_cell:
            AD, DP = AD.T, DP.T
        Y = impute_gene_medians(to_dense(Y))
        if Y.shape[0] == self.n_var_rdr and Y.shape[1] == self.n_cell:
            Y = Y.T
        if isinstance(DP, np.ndarray) and np.mean(DP > 0) < 0.3:
            AD, DP = csc_matrix(AD).T, csc_matrix(DP).T
        else:
            AD, DP = AD.T, DP.T
        if random_seed is not None:
            np.random.seed(random_seed)
        self.ELBO_inits = []
        best_id, best_elbo = None, -np.inf
        for _ in range(n_init):
            self._baf.set_initial()
            self._rdr.set_initial()
            self.ID_prob = normalize(np.random.rand(self.n_cell, self.n_donor))
            self._sync_id()
            self._fit_joint(AD, DP, Y, max_iter=max_iter_pre, verbose=verbose)
            el = float(self.ELBO_iters[-1])
            self.ELBO_inits.append(el)
            if el > best_elbo:
                best_elbo, best_id = el, self.ID_prob.copy()
        self.ID_prob = best_id
        self._sync_id()
        self.ELBO_iters = np.array([])
        self._fit_joint(AD, DP, Y, max_iter=max_iter, verbose=verbose)
        return self
