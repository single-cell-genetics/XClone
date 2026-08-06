"""
Lightweight Binomial mixture VB for BAF bin-level clustering.

API-compatible with ``vireoSNP.BinomMixtureVB`` (variants x cells layout for AD/DP).
Uses ``xclone.model.base_utils``; falls back to vireoSNP when installed.
"""

from __future__ import annotations

import numpy as np
from scipy.sparse import csc_matrix
from scipy.special import betaln, binom, digamma
from scipy.stats import entropy

from .base_utils import loglik_amplify, normalize

try:
    from vireoSNP.utils.vireo_base import beta_entropy, get_binom_coeff
except ImportError:  # pragma: no cover
    def beta_entropy(alpha, alpha_prior):
        """
        KL divergence between Beta posteriors and priors.

        Matches vireoSNP.utils.vireo_base.beta_entropy behavior for (N, 2)
        shapes used by BinomMixtureVB.
        """
        if len(alpha.shape) == 1 and alpha.shape[0] == 2:
            alpha = alpha.reshape(-1, 2)
        if len(alpha_prior.shape) == 1 and alpha_prior.shape[0] == 2:
            alpha_prior = alpha_prior.reshape(-1, 2)

        def _beta_cross_entropy(xp, xq):
            return (
                betaln(xq[:, 0], xq[:, 1])
                - (xq[:, 0] - 1) * digamma(xp[:, 0])
                - (xq[:, 1] - 1) * digamma(xp[:, 1])
                + (xq.sum(axis=1) - 2) * digamma(xp.sum(axis=1))
            )

        rv = _beta_cross_entropy(alpha, alpha_prior) - _beta_cross_entropy(alpha, alpha)
        return np.sum(rv)

    def get_binom_coeff(AD, DP, max_val=700, is_log=True):
        """Get log binomial coefficients matching vireo fallback behavior."""
        idx = DP > 0
        _AD = AD[idx].astype(np.int64)
        _DP = DP[idx].astype(np.int64)
        coeff = np.log(binom(_DP, _AD))
        coeff[coeff > max_val] = max_val
        return coeff.astype(np.float32)


def get_binom_mixture_vb_class():
    """Return BinomMixtureVB class (vireoSNP if available, else built-in)."""
    try:
        from vireoSNP import BinomMixtureVB

        return BinomMixtureVB
    except ImportError:
        return XCloneBinomMixtureVB


class XCloneBinomMixtureVB:
    """
    Binomial mixture model with variational inference.

    Designed as a close in-repo replication of vireoSNP BinomMixtureVB for
    environments where vireoSNP is unavailable.
    """

    def __init__(
        self,
        n_cell,
        n_var,
        n_donor,
        fix_beta_sum=False,
        beta_mu_init=None,
        beta_sum_init=None,
        ID_prob_init=None,
    ):
        self.n_cell = n_cell
        self.n_var = n_var
        self.n_donor = n_donor
        self.fix_beta_sum = fix_beta_sum
        self.beta_mu_init = beta_mu_init
        self.beta_sum_init = beta_sum_init
        self.ID_prob_init = ID_prob_init
        self.set_prior()
        self.set_initial(beta_mu_init, beta_sum_init, ID_prob_init)
        self.ELBO_iters = np.array([])

    def set_prior(self, ID_prior=None, beta_mu_prior=None, beta_sum_prior=None):
        if beta_mu_prior is None:
            beta_mu_prior = np.ones((self.n_var, self.n_donor)) * 0.5
        if beta_sum_prior is None:
            beta_sum_prior = np.ones(beta_mu_prior.shape) * 2.0
        self.theta_s1_prior = beta_mu_prior * beta_sum_prior
        self.theta_s2_prior = (1 - beta_mu_prior) * beta_sum_prior
        if ID_prior is not None:
            if len(ID_prior.shape) == 1:
                ID_prior = np.expand_dims(ID_prior, axis=0)
            self.ID_prior = ID_prior
        else:
            self.ID_prior = normalize(np.ones((self.n_cell, self.n_donor)))

    def set_initial(self, beta_mu_init=None, beta_sum_init=None, ID_prob_init=None):
        self.beta_mu = (
            beta_mu_init
            if beta_mu_init is not None
            else np.ones((self.n_var, self.n_donor)) * 0.5
        )
        self.beta_sum = (
            beta_sum_init
            if beta_sum_init is not None
            else np.ones(self.beta_mu.shape) * 30
        )
        if ID_prob_init is not None:
            self.ID_prob = normalize(ID_prob_init, axis=1)
        else:
            self.ID_prob = normalize(np.random.rand(self.n_cell, self.n_donor))

    @property
    def theta_s1(self):
        return self.beta_mu * self.beta_sum

    @property
    def theta_s2(self):
        return (1 - self.beta_mu) * self.beta_sum

    def get_E_logLik(self, AD, DP):
        bd = DP - AD
        return (
            AD.T @ digamma(self.theta_s1)
            + bd.T @ digamma(self.theta_s2)
            - DP.T @ digamma(self.theta_s1 + self.theta_s2)
        )

    def update_theta_size(self, AD, DP):
        bd = DP - AD
        s1 = AD @ self.ID_prob + self.theta_s1_prior
        s2 = bd @ self.ID_prob + self.theta_s2_prior
        self.beta_mu = s1 / (s1 + s2)
        if not self.fix_beta_sum:
            self.beta_sum = s1 + s2

    def update_ID_prob(self, AD=None, DP=None, logLik_ID=None):
        if logLik_ID is None:
            logLik_ID = self.get_E_logLik(AD, DP)
        self.ID_prob = normalize(
            np.exp(loglik_amplify(logLik_ID + np.log(self.ID_prior + 1e-300)))
        )

    def get_ELBO(self, logLik_ID):
        lb = np.sum(logLik_ID * self.ID_prob)
        kl_id = np.sum(entropy(self.ID_prob, self.ID_prior, axis=-1))
        kl_theta = beta_entropy(
            np.append(
                np.expand_dims(self.theta_s1, 1),
                np.expand_dims(self.theta_s2, 1),
                axis=1,
            ),
            np.append(
                np.expand_dims(self.theta_s1_prior, 1),
                np.expand_dims(self.theta_s2_prior, 1),
                axis=1,
            ),
        )
        return lb - kl_id - kl_theta

    def _fit_BV(self, AD, DP, max_iter=200, min_iter=20, epsilon_conv=1e-2, verbose=True):
        ELBO = np.zeros(max_iter)
        for it in range(max_iter):
            self.update_theta_size(AD, DP)
            ll = self.get_E_logLik(AD, DP)
            self.update_ID_prob(logLik_ID=ll)
            ELBO[it] = self.get_ELBO(ll)
            if it > min_iter:
                delta = ELBO[it] - ELBO[it - 1]
                if delta < -1e-6:
                    if verbose:
                        print("Warning: ELBO decreases %.8f to %.8f!" % (ELBO[it - 1], ELBO[it]))
                elif it == max_iter - 1:
                    if verbose:
                        print("Warning: VB did not converge!")
                elif delta < epsilon_conv:
                    break

        # Keep behavior aligned with vireo code path.
        self.ELBO_iters = np.append(self.ELBO_iters, ELBO[:it])

    def fit(self, AD, DP, n_init=10, max_iter=200, max_iter_pre=100, random_seed=None, **kwargs):
        if random_seed is not None:
            np.random.seed(random_seed)
        if type(DP) is np.ndarray and np.mean(DP > 0) < 0.3:
            AD, DP = csc_matrix(AD), csc_matrix(DP)
        binom_coeff = np.sum(get_binom_coeff(AD, DP, is_log=True))
        self.ELBO_inits = []
        best = None
        for i in range(n_init):
            self.set_initial(self.beta_mu_init, self.beta_sum_init, self.ID_prob_init)
            self._fit_BV(AD, DP, max_iter=max_iter_pre, **kwargs)
            self.ELBO_inits.append(float(self.ELBO_iters[-1]))
            if i == 0 or self.ELBO_iters[-1] >= max(self.ELBO_inits[:-1]):
                best = (self.ID_prob.copy(), self.beta_mu.copy(), self.beta_sum.copy(), self.ELBO_iters.copy())
        self.set_initial(best[1], best[2], best[0])
        self.ELBO_iters = best[3]
        self._fit_BV(AD, DP, max_iter=max_iter, **kwargs)
        self.ELBO_iters = self.ELBO_iters + binom_coeff
        self.ELBO_inits = np.array(self.ELBO_inits) + binom_coeff
        return self
