"""XClone BAF and RDR module (Variational inference)
BAF: binomial mixture model (Updated from Vireo)
RDR: possion mixture model (simple one)
"""
# Author: Rongting Huang
# Date: 2021/01/26
# update: 2021/10/06

import itertools
import numpy as np
from scipy.stats import entropy
from scipy.sparse import issparse
# from scipy.sparse import isspmatrix
from scipy.optimize import minimize
from scipy.special import logsumexp, digamma, betaln, binom, gammaln
from .XClone_VB_base import normalize, loglik_amplify, beta_entropy, get_binom_coeff
# from .base_utils import normalize, loglik_amplify

from .XClone_VB_base import logfact, gamma_entropy

from numpy import log as ln

__docformat__ = "restructuredtext en"


def _check_X(X,):
    """Check the input data X.
    ensure the data format of X.
    """
    pass


class XCloneBase():
    """XClone model: Variational Inference model for calling CNVs and
    clustering cells.
    XCloneVBase provide the base functions.
    """
    pass


class BAFXClone(XCloneBase):
    """XClone BAF module.
    base model: binomial distribution

    Key properties
    --------------
    beta_mu: numpy array (1,n_CNVs) or (n_features, n_CNVs)
        Beta mean parameter of theta's posterior
    beta_sum: numpy array (1, n_CNVs) or (n_features, n_CNVs), same as beta_mu
        Beta concetration parameter of theta's posterior
        Note that params of theta (beta distribution) includes alpha, beta
            numpy array (1, n_CNVs) or (n_features, n_CNVs)
            alpha_: theta_s1 = beta_mu * beta_sum
            beta_: theta_s2 = (1-beta_mu) * beta_sum
    ID_prob: numpy array (n_samples, n_components)
        Posterior cell assignment probability to each clone/cell cluster
    CNV_prob: numpy array (n_features, n_components, n_CNVs)
        Posterior CNV_state probability per gene/block per clone/cell cluster
    """

    def __init__(self, n_features=1, n_samples=1, n_components=1, n_CNVs=16,
        learn_CNV=True, learn_params=True, fix_pi = False, fix_beta_sum=False, FSP_mode=True):
        """Initializing the XClone model
        Note, multiple initializations are highly recomended to avoid local
        optima.

        Parameters
        ----------
        n_samples : int.
            Number of cells
        n_features : int.
            Number of blocks/genes
        n_components : int.
            Number of clone/cell cluster
        n_CNVs : int.
            Number of CNV states (categories)
        learn_CNV: bool.
            Whether updating `CNV_prob`; otherwise using the initial
        FSP_mode: bool.h
            Whether setting params `theta`/`m` to be feature(block/gene) specific
        fix_beta_sum: bool
            Whether fix the concetration parameter of theta's posterior (beta distribution)
        fix_pi: bool
            Whether fix the pi to the initial value
        """
        self.n_features = n_features
        self.n_samples = n_samples
        self.K = n_components  # for short cut
        self.n_CNVs = n_CNVs
        self.learn_CNV = learn_CNV
        self.learn_params = learn_params

        self.fix_pi = fix_pi
        self.fix_beta_sum = fix_beta_sum

        self.FSP_mode = FSP_mode
        self.theta_len = n_features if FSP_mode else 1
        # initial key params by random initial
        self.set_init()
        self.set_prior()
        self.ELBO_ = np.zeros((0))

    def set_init(self, beta_mu_init=None, beta_sum_init=None, ID_prob_init=None,
        CNV_prob_init=None,random_state=None):
        '''
        params_init: numpy array (1, n_CNVs) or (n_features, n_CNVs)
            Initial value of beta_mu, the mean parameter of theta
        ID_prob_init: numpy array (n_samples, n_components)
            Initial value of ID_prob, cell assignment probability to each donor
        CNV_prob_init: numpy array (n_features, n_components, n_CNVs)
            Initial value of CNV_prob, CNV state probability per gene/block in a clone/cell cluster
        '''
        # inital key parameters
        if random_state is not None:
            np.random.seed(seed =  random_state)

        if beta_mu_init is not None:
            self.beta_mu = beta_mu_init
        else:
            self.beta_mu = (np.ones((self.theta_len, self.n_CNVs)) *
                np.linspace(0.01, 0.99, self.n_CNVs).reshape(1, -1))

        if beta_sum_init is not None:
            self.beta_sum = beta_sum_init
        else:
            self.beta_sum = np.ones((self.theta_len, self.n_CNVs)) * 50

        if ID_prob_init is not None:
            self.ID_prob = normalize(ID_prob_init, axis=1)
        else:
            self.ID_prob = normalize(np.random.rand(self.n_samples, self.K))

        if CNV_prob_init is not None:
            self.CNV_prob = normalize(CNV_prob_init)## check
        else:
            self.CNV_prob = normalize(np.random.rand(self.n_features, self.K, self.n_CNVs))

    def set_prior(self, CNV_prior=None, ID_prior=None, beta_mu_prior=None,
        beta_sum_prior=None, min_CNVS=0.00001):
        """Set prior for key variables: theta, CNV_prob and ID_prob.
        The priors are in the same shape as its according variables.
        min_CNV: float. Minimun genotype probability in CNV_prior.
        """
        if beta_mu_prior is None:
            beta_mu_prior = np.expand_dims(
                np.linspace(0.01, 0.99, self.beta_mu.shape[1]), axis=0)
        if beta_sum_prior is None:
            beta_sum_prior = np.ones(beta_mu_prior.shape) * 50.0
        self.theta_s1_prior = beta_mu_prior * beta_sum_prior
        self.theta_s2_prior = (1 - beta_mu_prior) * beta_sum_prior

        if ID_prior is not None:
            if len(ID_prior.shape) == 1:
                ID_prior = np.expand_dims(ID_prior, axis=0)
            self.ID_prior = ID_prior
        else:
            self.ID_prior = normalize(np.ones(self.ID_prob.shape))

        if CNV_prior is not None:
            if len(CNV_prior.shape) == 2:
                CNV_prior = np.expand_dims(CNV_prior, axis=0)
            CNV_prior[CNV_prior < min_GP] = min_GP
            CNV_prior[CNV_prior > 1 - min_GP] = 1 - min_GP
            CNV_prior = normalize(CNV_prior)
            self.CNV_prior = CNV_prior
        else:
            self.CNV_prior = normalize(np.ones(self.CNV_prob.shape))

    @property
    def theta_s1(self):
        """Beta concetration1 parameter for theta posterior"""
        return self.beta_mu * self.beta_sum

    @property
    def theta_s2(self):
        """Beta concetration2 parameter for theta posterior"""
        return (1 - self.beta_mu) * self.beta_sum

    @property
    def digamma1_(self):
        """Digamma of Beta concetration1 parameter"""
        return np.expand_dims(digamma(self.theta_s1), 1)

    @property
    def digamma2_(self):
        """Digamma of Beta concetration2 parameter"""
        return np.expand_dims(digamma(self.theta_s2), 1)

    @property
    def digammas_(self):
        """Digamma of Beta concetration summary parameter"""
        return np.expand_dims(digamma(self.theta_s1 + self.theta_s2), 1)

    def update_ID_prob(self, data):
        """Coordinate ascent for updating assignment probability
        """
        AD, DP = data[0], data[1]
        BD = DP - AD
        logLik_ID = np.zeros((AD.shape[1], self.K))
        for t in range(self.n_CNVs):
            S1 = AD.T @ (self.CNV_prob[:, :, t] * self.digamma1_[:, :, t])
            S2 = BD.T @ (self.CNV_prob[:, :, t] * self.digamma2_[:, :, t])
            SS = DP.T @ (self.CNV_prob[:, :, t] * self.digammas_[:, :, t])
            logLik_ID += (S1 + S2 - SS)

        self.ID_prob = normalize(np.exp(loglik_amplify(
            logLik_ID + np.log(self.ID_prior))))

        return logLik_ID

    def update_CNV_prob(self, data):
        """Coordinate ascent for updating CNV probability
        """
        AD, DP = data[0], data[1]
        # BD = DP - AD
        S1_cnv = AD @ self.ID_prob
        SS_cnv = DP @ self.ID_prob
        S2_cnv = SS_cnv - S1_cnv

        logLik_CNV = np.zeros(self.CNV_prior.shape)
        for t in range(self.n_CNVs):
            logLik_CNV[:, :, t] = (
                S1_cnv * self.digamma1_[:, :, t] +
                S2_cnv * self.digamma2_[:, :, t] -
                SS_cnv * self.digammas_[:, :, t])

        self.CNV_prob = normalize(np.exp(loglik_amplify(
            logLik_CNV + np.log(self.CNV_prior))))

    def update_params(self, data):
        """Coordinate ascent for updating theta posterior parameters
        """
        AD, DP = data[0], data[1]
        BD = DP - AD
        S1_ik = AD @ self.ID_prob  #(n_features, n_components)
        S2_ik = BD @ self.ID_prob  #(n_features, n_components)

        _theta_s1 = np.zeros(self.beta_mu.shape)
        _theta_s2 = np.zeros(self.beta_mu.shape)
        _theta_s1 += self.theta_s1_prior.copy()
        _theta_s2 += self.theta_s2_prior.copy()
        for t in range(self.n_CNV):
            _axis = 1 if self.FSP_mode else None
            _theta_s1[:, t:(t+1)] += np.sum(
                S1_ik * self.CNV_prob[:, :, t], axis=_axis, keepdims=True)
            _theta_s2[:, t:(t+1)] += np.sum(
                S2_ik * self.CNV_prob[:, :, t], axis=_axis, keepdims=True)

        self.beta_mu = _theta_s1 / (_theta_s1 + _theta_s2)
        if self.fix_beta_sum == False:
            self.beta_sum = _theta_s1 + _theta_s2

    def get_ELBO(self,logLik_ID, data=None):
        """calculate variational evidence lower bound with current parameters.

        The lower bound is used to detect the convergence and has to increase at
        each iteration.

        For BAF module-

        parameters
        ----------

        Returns
        -------
        lower_bound : float
        """
        AD, DP = data[0], data[1]
        if logLik_ID is None:
            BD = DP - AD
            logLik_ID = np.zeros((AD.shape[1],self.n_components))
            for t in range(self.n_CNV):
                S0 = get_binom_coeff(AD, DP).T @ self.CNV_prob[:, :, t] # vireo has an error here, maybe
                S1 = AD.T @ (self.CNV_prob[:, :, t] * self.digamma1_[:, :, t])
                S2 = BD.T @ (self.CNV_prob[:, :, t] * self.digamma2_[:, :, t])
                SS = DP.T @ (self.CNV_prob[:, :, t] * self.digammas_[:, :, t])
                logLik_ID += (S0 + S1 + S2 - SS)

        LB_p = np.sum(logLik_ID * self.ID_prob)
        KL_ID = -np.sum(entropy(self.ID_prob, self.ID_prior, axis=-1))
        KL_CNV = -np.sum(entropy(self.CNV_prob, self.CNV_prior, axis=-1))
        KL_paras = -beta_entropy(
            np.append(
                np.expand_dims(self.theta_s1, 1),
                np.expand_dims(self.theta_s2, 1), axis = 1),
            np.append(
                np.expand_dims(self.theta_s1_prior, 1),
                np.expand_dims(self.theta_s2_prior, 1), axis = 1))
        return LB_p + KL_ID + KL_CNV + KL_paras # different from vireo

    def fit(self,data,max_iter=200, min_iter=5,epsilon_conv=1e-2,
        delay_fit_theta=0, verbose=False):
        """Fit XClone model (BAF module) with coordinate ascent
        Parameters
        ----------
        AD : scipy.sparse.csc_matrix (n_genes, n_cells)
            Sparse count matrix for alternative allele
        DP : scipy.sparse.csc_matrix (n_genes, n_cells)
            Sparse count matrix for depths, alternative + refeerence alleles
        max_iter : int
            Maximum number of iterations
        min_iter :
            Minimum number of iterations
        epsilon_conv : float
            Threshold for detecting convergence
        delay_fit_theta : int
            Number of steps to delay updating theta. This can be very useful
            for common genetics when there is good prior on allelic ratio.
        verbose : bool
            Whether print out log info
        """
        AD, DP = data[0], data[1]
         # check AD and DP's type: scipy.sparse or numpy.matrix
        if type(DP) is np.ndarray and np.mean(DP > 0) < 0.3:
            print("Warning: input matrices is %.1f%% sparse, "
                  %(100 - np.mean(DP > 0) * 100) +
                  "change to scipy.sparse.csc_matrix" )
            AD = csc_matrix(AD)
            DP = csc_matrix(DP)

        # _binom_coeff = np.sum(get_binom_coeff([AD, DP]))
        ELBO = np.zeros(max_iter)
        for it in range(max_iter):
            if self.learn_theta and it >= delay_fit_theta:
                self.update_params([AD, DP])
            if self.learn_CNV:
                self.update_CNV_prob([AD, DP])

            _logLik_ID = self.update_ID_prob([AD, DP])
            # ELBO[it] = self.get_ELBO(_logLik_ID) + _binom_coeff
            ELBO[it] = self.get_ELBO(_logLik_ID)

            if it > min_iter:
                if ELBO[it] < ELBO[it - 1]:
                    if verbose:
                        print("Warning: Lower bound decreases!\n")
                elif it == max_iter - 1:
                    if verbose:
                        print("Warning: VB did not converge!\n")
                elif ELBO[it] - ELBO[it - 1] < epsilon_conv:
                    break
        self.ELBO_ = np.append(self.ELBO_, ELBO[:it])


class RDRXClone_test1(XCloneBase):
    """XClone RDR module.  first version--not tested
    base model: Poisson mixture model
    Base model of X[i, j] in CNV state t:
        Poisson(X[i,j] | self.params[i, t] * self.mean_factors[j])

    Key properties
    --------------
    beta_mu: numpy array (1,n_CNVs) or (n_features, n_CNVs)
        Beta mean parameter of theta's posterior
    beta_sum: numpy array (1, n_CNVs) or (n_features, n_CNVs), same as beta_mu
        Beta concetration parameter of theta's posterior
        Note that params of theta (beta distribution) includes alpha, beta
            numpy array (1, n_CNVs) or (n_features, n_CNVs)
            alpha_: theta_s1 = beta_mu * beta_sum
            beta_: theta_s2 = (1-beta_mu) * beta_sum
    ID_prob: numpy array (n_samples, n_components)
        Posterior cell assignment probability to each clone/cell cluster
    CNV_prob: numpy array (n_features, n_components, n_CNVs)
        Posterior CNV_state probability per gene/block per clone/cell cluster
    """

    def __init__(self, n_features=1, n_samples=1, n_components=1, n_CNVs=16,
        learn_CNV=True, learn_params=True, FSP_mode=True,err1=0.01):
        """Initializing the XClone model
        Note, multiple initializations are highly recomended to avoid local
        optima.

        Parameters
        ----------
        n_samples : int.
            Number of cells
        n_features : int.
            Number of blocks/genes
        n_components : int.
            Number of clone/cell cluster
        n_CNVs : int.
            Number of CNV states (categories)
        learn_CNV: bool.
            Whether updating `CNV_prob`; otherwise using the initial
        FSP_mode: bool.h
            Whether setting params `theta`/`m` to be feature(block/gene) specific
        fix_beta_sum: bool
            Whether fix the concetration parameter of theta's posterior (beta distribution)
        fix_pi: bool
            Whether fix the pi to the initial value
        """
        self.n_features = n_features
        self.n_samples = n_samples
        self.K = n_components  # for short cut
        self.n_CNVs = n_CNVs
        self.learn_CNV = learn_CNV
        self.learn_params = learn_params

        self.FSP_mode = FSP_mode
        self.params_len = n_features if FSP_mode else 1

        self.err1_ = err1
        # initial key params by random initial
        self.set_init()
        self.set_prior()
        self.ELBO_ = np.zeros((0))

    def set_init(self, alpha_s1_init=None, beta_s2_init=None, ID_prob_init=None,
        CNV_prob_init=None,random_state=None):
        '''
        alpha_s1_init: numpy array (1, n_CNVs) or (n_features, n_CNVs)
            Initial value of alpha_s1, the alpha parameter of gamma distribution
        beta_s2_init: numpy array (1, n_CNVs) or (n_features, n_CNVs)
            Initial value of beta_s2, the beta parameter of gamma distribution
        ID_prob_init: numpy array (n_samples, n_components)
            Initial value of ID_prob, cell assignment probability to each donor
        CNV_prob_init: numpy array (n_features, n_components, n_CNVs)
            Initial value of CNV_prob, CNV state probability per gene/block in a clone/cell cluster
        '''
        # inital key parameters
        if random_state is not None:
            np.random.seed(seed =  random_state)

        if alpha_s1_init is not None:
            self.alpha_s1 = alpha_s1_init
        else:
            self.alpha_s1 = (np.ones((self.params_len, self.n_CNVs)) *
                np.linspace(0.01, 0.99, self.n_CNVs).reshape(1, -1))

        if beta_s2_init is not None:
            self.beta_s2 = beta_s2_init
        else:
            self.beta_s2 = np.ones((self.params_len, self.n_CNVs)) * 50

        if ID_prob_init is not None:
            self.ID_prob = normalize(ID_prob_init, axis=1)
        else:
            self.ID_prob = normalize(np.random.rand(self.n_samples, self.K))

        if CNV_prob_init is not None:
            self.CNV_prob = normalize(CNV_prob_init)## check
        else:
            self.CNV_prob = normalize(np.random.rand(self.n_features, self.K, self.n_CNVs))

    def set_prior(self, CNV_prior=None, ID_prior=None, alpha_s1_prior=None,
        beta_s2_prior=None, min_CNVS=0.00001): ## to do
        """Set prior for key variables: theta, CNV_prob and ID_prob.
        The priors are in the same shape as its according variables.
        min_CNV: float. Minimun genotype probability in CNV_prior.
        """
        if alpha_s1_prior is None:
            alpha_s1_prior = np.expand_dims(
                np.linspace(0.01, 0.99, self.alpha_s1.shape[1]), axis=0)
        if beta_s2_prior is None:
            beta_s2_prior = np.ones(alpha_s1_prior.shape) * 50.0
        self.alpha_s1_prior = alpha_s1_prior
        self.beta_s2_prior = beta_s2_prior

        if ID_prior is not None:
            if len(ID_prior.shape) == 1:
                ID_prior = np.expand_dims(ID_prior, axis=0)
            self.ID_prior = ID_prior
        else:
            self.ID_prior = normalize(np.ones(self.ID_prob.shape))

        if CNV_prior is not None:
            if len(CNV_prior.shape) == 2:
                CNV_prior = np.expand_dims(CNV_prior, axis=0)
            CNV_prior[CNV_prior < min_GP] = min_GP
            CNV_prior[CNV_prior > 1 - min_GP] = 1 - min_GP
            CNV_prior = normalize(CNV_prior)
            self.CNV_prior = CNV_prior
        else:
            self.CNV_prior = normalize(np.ones(self.CNV_prob.shape))

    # @property
    # def alpha_s1(self):
    #     """Beta concetration1 parameter for theta posterior"""
    #     return self.beta_mu * self.beta_sum ## TO DO
    #
    # @property
    # def beta_s2(self):
    #     """Beta concetration2 parameter for theta posterior"""
    #     return (1 - self.beta_mu) * self.beta_sum ## TO DO

    @property
    def gamma_exp1_(self):
        """Digamma of Beta concetration1 parameter"""
        return np.expand_dims(digamma(self.alpha_s1)-ln(self.beta_s2), 1)

    @property
    def gamma_exp2_(self):
        """Digamma of Beta concetration2 parameter"""
        return np.expand_dims(self.alpha_s1/self.beta_s2, 0)

    @property
    def tau_(self):
        ## define paras for different CNV states
        ct1_ = np.array([0,1,0,1,2,0,2,1,2,3,0,3,1,3,2,3])
        ct2_ = np.array([0,0,1,1,0,2,1,2,2,0,3,1,3,2,3,3])
        tau_value = (ct1_ + ct2_) / 2
        tau_value[0] = self.err1_ #TO DO- GIVE default error value
        return tau_value

    # @property
    # def omega_(self):
    #     ## define paras for different CNV states
    #     omega_value = [0,]
    #     for i in range(1,16):
    #         omega_value.append(1)
    #     omega_value = np.array(omega_value)
    #     return omega_value

    def update_ID_prob(self, data):## to do
        """Coordinate ascent for updating assignment probability
        """
        RDR = data
        EX_per_cell = np.expand_dims(RDR.mean(axis=0).T,1)
        logLik_ID = np.zeros((RDR.shape[1], self.K))
        for t in range(self.n_CNVs):
            S1 = (RDR.T * ln(EX_per_cell)) @ self.CNV_prob[:, :, t]
            S3 = logfact(RDR).T @ self.CNV_prob[:, :, t]
            if t == 0:
                S2 = RDR.T @ (self.CNV_prob[:, :, t] * (ln(self.tau_[t])))
                S4 = self.tau_[t] * (EX_per_cell @ np.ones(self.gamma_exp2_.shape)) @ self.CNV_prob[:, :, t]
            S2 = RDR.T @ (self.CNV_prob[:, :, t] * (ln(self.tau_[t]) + self.gamma_exp1_))
            S4 = self.tau_[t] * (EX_per_cell @ self.gamma_exp2_) @ self.CNV_prob[:, :, t]
            logLik_ID += (S1 + S2 - S3 - S4)
        self.ID_prob = normalize(np.exp(loglik_amplify(
            logLik_ID + np.log(self.ID_prior))))

        return logLik_ID

    def update_CNV_prob(self, data):## to do
        """Coordinate ascent for updating CNV probability
        """
        RDR = data
        EX_per_cell = np.expand_dims(RDR.mean(axis=0).T,1)

        S1_cnv = RDR * ln(EX_per_cell) @ self.ID_prob
        S3_cnv = logfact(RDR) @ self.ID_prob

        logLik_CNV = np.zeros(self.CNV_prior.shape)
        for t in range(self.n_CNVs):
            if t == 0:
                S2_cnv = RDR @ self.ID_prob * (ln(self.tau_[t]))
                S4_cnv =  self.tau_[t] * (EX_per_cell @ np.ones(self.gamma_exp2_.shape)) @ self.ID_prob
            else:
                S2_cnv = RDR @ self.ID_prob * (ln(self.tau_[t]) + self.gamma_exp1_)
                S4_cnv = self.tau_[t] * (EX_per_cell @ self.gamma_exp2_) @ self.ID_prob
            logLik_CNV[:,:,t] = S1_cnv + S2_cnv - S3_cnv - S4_cnv

        self.CNV_prob = normalize(np.exp(loglik_amplify(
            logLik_CNV + np.log(self.CNV_prior))))

    def update_params(self, data):## to do
        """Coordinate ascent for updating theta posterior parameters
        """
        RDR = data
        EX_per_cell = np.expand_dims(RDR.mean(axis=0).T,1)
        S1_ik = RDR @ self.ID_prob  #(n_features, n_components)
        S2_ik = EX_per_cell @ self.ID_prob  #(n_features, n_components)

        _alpha_s1 = np.zeros(self.alpha_s1.shape)
        _beta_s2 = np.zeros(self.beta_s2.shape)
        _alpha_s1 += self.alpha_s1_prior.copy()
        _beta_s2 += self.beta_s2_prior.copy()
        for t in range(1,self.n_CNV):
            _axis = 1 if self.FSP_mode else None
            _alpha_s1[:, t:(t+1)] += np.sum(
                S1_ik * self.CNV_prob[:, :, t], axis=_axis, keepdims=True)
            _beta_s2[:, t:(t+1)] += np.sum(
                S2_ik * self.tau_[t] * self.CNV_prob[:, :, t], axis=_axis, keepdims=True)

        self.alpha_s1 = _alpha_s1
        self.beta_s2 = _beta_s2


    def get_ELBO(self,logLik_ID, data):
        """calculate variational evidence lower bound with current parameters.

        The lower bound is used to detect the convergence and has to increase at
        each iteration.

        For RDR module-

        parameters
        ----------

        Returns
        -------
        lower_bound : float
        """
        RDR = data
        EX_per_cell = np.expand_dims(RDR.mean(axis=0).T,1)
        if logLik_ID is None:
            logLik_ID = np.zeros((RDR.shape[1],self.n_components))
            for t in range(self.n_CNV):
                S1 = (RDR.T * ln(EX_per_cell)) @ self.CNV_prob[:, :, t]
                S3 = logfact(RDR).T @ self.CNV_prob[:, :, t]
                if t == 0:
                    S2 = RDR.T @ (self.CNV_prob[:, :, t] * (ln(self.tau_[t])))
                    S4 = self.tau_[t] * (EX_per_cell @ np.ones(self.gamma_exp2_.shape)) @ self.CNV_prob[:, :, t]
                #S2 = RDR.T @ (self.CNV_prob[:, :, t] * (ln(self.tau_[t]) + self.omega_[t] * self.gamma_exp1_))
                S2 = RDR.T @ (self.CNV_prob[:, :, t] * (ln(self.tau_[t]) + self.gamma_exp1_))
                S4 = self.tau_[t] * (EX_per_cell @ self.gamma_exp2_) @ self.CNV_prob[:, :, t]
                logLik_ID += (S1 + S2 - S3 - S4)

        LB_p = np.sum(logLik_ID * self.ID_prob)
        KL_ID = -np.sum(entropy(self.ID_prob, self.ID_prior, axis=-1))
        KL_CNV = -np.sum(entropy(self.CNV_prob, self.CNV_prior, axis=-1))
        KL_paras = -gamma_entropy(
        np.append(
            np.expand_dims(self.alpha_s1, 1),
            np.expand_dims(self.beta_s2, 1), axis = 1),
        np.append(
            np.expand_dims(self.alpha_s1_prior, 1),
            np.expand_dims(self.beta_s2_prior, 1), axis = 1))
        return LB_p + KL_ID + KL_CNV + KL_paras # FRAMEWORK

    def fit(self,data,max_iter=200, min_iter=5,epsilon_conv=1e-2,
        delay_fit_theta=0, verbose=False):## to do
        """Fit XClone model (RDR module) with coordinate ascent
        Parameters
        ----------
        RDR : scipy.sparse.csc_matrix (n_genes, n_cells)
            Sparse count matrix for expression data
        max_iter : int
            Maximum number of iterations
        min_iter :
            Minimum number of iterations
        epsilon_conv : float
            Threshold for detecting convergence
        delay_fit_theta : int
            Number of steps to delay updating theta. This can be very useful
            for common genetics when there is good prior on allelic ratio.
        verbose : bool
            Whether print out log info
        """
        RDR = data
         # check RDR's type: scipy.sparse or numpy.matrix
        if type(RDR) is np.ndarray and np.mean(RDR > 0) < 0.3:
            print("Warning: input matrices is %.1f%% sparse, "
                  %(100 - np.mean(DP > 0) * 100) +
                  "change to scipy.sparse.csc_matrix" )
            RDR = csc_matrix(RDR)

        ELBO = np.zeros(max_iter)
        for it in range(max_iter):
            if self.learn_theta and it >= delay_fit_theta:## to do
                self.update_params(RDR)
            if self.learn_CNV:
                self.update_CNV_prob(RDR)

            _logLik_ID = self.update_ID_prob(RDR)
            ELBO[it] = self.get_ELBO(_logLik_ID)# to do

            if it > min_iter:
                if ELBO[it] < ELBO[it - 1]:
                    if verbose:
                        print("Warning: Lower bound decreases!\n")
                elif it == max_iter - 1:
                    if verbose:
                        print("Warning: VB did not converge!\n")
                elif ELBO[it] - ELBO[it - 1] < epsilon_conv:
                    break
        self.ELBO_ = np.append(self.ELBO_, ELBO[:it])

class RDRXClone(XCloneBase):
    """XClone RDR module.  
    new version 2021-10-06 updated 2021-10-27
    TO DO: efficiency part-- cellbased not so fast
    Variational Inference for reconstruction of CNV state.
    
    The prior can be set via set_prior() before fitting the model.
    
    base model: Poisson mixture model
    Base model of X[i, j] in CNV state t:
        Poisson(X[i,j] | self.params[i, j, t] * self.mean_factors[j])
    
    The prior can be set via set_prior() before fitting the model.
    
    Key properties
    --------------
    alpha_s1: numpy array
    beta_s2: numpy array

    CNV_prob: numpy array (n_features, n_components, n_CNVs)
        Posterior CNV_state probability per gene/block per clone/cell cluster
    
    TEST PART:
    ----------
        i genes j cells  t CNV states
        # def get_libratio(self, ):

        libratio = np.random.rand(1,8)
        ref_counts = np.random.randint(10, size = (10,1))
        r_counts = np.random.randint(10, size = (10,8))

        alpha_ =  np.random.rand(3,10,8)
        beta_ = np.random.rand(3,10,8)

        update 
        alpha_ =  np.random.rand(10,8,3)
        beta_ = np.random.rand(10,8,3)
        ## input format, or maybe need reassign 
        libratio[:,:,np.newaxis]
        ref_counts[:,:,np.newaxis]
        r_counts[:,:,np.newaxis]

    Example:
    --------
        ## previous record
        # libratio = np.random.rand(1,self.n_samples)
        # ref_counts = np.random.randint(10, size = (self.n_features,1))

        _model_rdr_vb = xclone.XClone_VB_model.RDRXClone1(n_features=1, n_samples=1, n_CNVs=3, libratio=libratio, ref_counts=ref_counts) # , libratio=libratio, ref_counts=ref_counts
        min_iter = 20
        _model_rdr_vb.fit(r_counts, min_iter=min_iter)
    """

    def __init__(self, n_features=1, n_samples=1, n_CNVs=3,
        ref_counts=None, libratio=None, err1=0.01,
        learn_CNV=True, learn_params=True, FSP_mode=True,
        mu_init = np.array([0.5, 1, 1.5]), var_init = 0.01,
        alpha_s1_init=None, beta_s2_init=None, CNV_prob_init=None, 
        random_state=None, verbose = False):
        """Initializing the XClone model
        Note, multiple initializations are highly recomended to avoid local
        optima.

        Parameters
        ----------
        n_features: int.
            Number of genes/blocks
        n_samples: int.
            Number of cells/celltypes
        n_CNVs: int.
            Number of CNV states (categories)
        ref_counts: numpy array.
            reference counts(each feature)
        libratio: numpy array.
            libsize ratio(each cell/celltype) for rescaling the r_counts
        learn_CNV: bool.
            Whether updating `CNV_prob`; otherwise using the initial
        learn_params: bool.
            Whether updating params `alpha_s1`, `beta_s2`; otherwise using the initial
        FSP_mode: bool.
            Whether setting params `alpha_s1`/`beta_s2` to be feature(gene/block) specific
        alpha_s1_init: numpy array (n_features, n_components, n_CNVs)
            Initial value of alpha_s1, the gamma parameter
        beta_s2_init: numpy array (n_features, n_components, n_CNVs), same as alpha_s1_init
            Initial value of beta_s2, the gamma parameter
        CNV_prob_init: numpy array (n_features, n_components, n_CNVs)
            Initial value of CNV_prob, assignment probability to each CNV state per feature and sample
        random_state: int.
            random seed for the model init.
        """
        self.verbose = verbose
        self.n_features = n_features
        self.n_samples = n_samples
        self.n_CNVs = n_CNVs

        ## libratio check
        try:
            if libratio.ndim == 3:
                print("[XClone] libratio is provided")
        except AttributeError:
            # reason-"AttributeError-'NoneType' object has no attribute 'ndim'"
            raise ValueError("[XClone] ValueError: pls set libratio for XClone RDR module!")
        else:
            if libratio.ndim == 3:
                self.libratio = libratio
            elif libratio.ndim == 2:
                self.libratio = libratio[:,:,np.newaxis]
            else:
                print("[XClone] pls check libratio format for XClone RDR module!")
        
        ## ref_counts check
        try:
            if ref_counts.ndim == 3:
                print("[XClone] ref_counts is provided")
        except AttributeError:
            raise ValueError("[XClone] ValueError: pls set ref_counts for XClone RDR module!")
        else:
            if ref_counts.ndim == 3:
                self.ref_counts = ref_counts
            elif ref_counts.ndim == 2:
                self.ref_counts = ref_counts[:,:,np.newaxis]
            else:
                print("[XClone] pls check ref_counts format for XClone RDR module!")

        self.learn_CNV = learn_CNV
        self.learn_params = learn_params

        # self.FSP_mode = FSP_mode
        # self.params_len = n_features if FSP_mode else 1

        # self.err1_ = err1
        self.mu_init = mu_init
        self.var_init = var_init

        self.ELBO_ = np.zeros((0))

        # initial key params of RDRXClone by random initial
        self.set_init(alpha_s1_init, beta_s2_init, CNV_prob_init, random_state)
        # set hyper parameters for the prior
        self.set_prior()

    def set_init(self, alpha_s1_init=None, beta_s2_init=None,
        CNV_prob_init=None, random_state=None):
        """
        alpha_s1_init: numpy array (n_features, n_samples, n_CNVs)
            Initial value of alpha_s1, the alpha parameter of gamma distribution
        beta_s2_init: numpy array (n_features, n_samples, n_CNVs)
            Initial value of beta_s2, the beta parameter of gamma distribution
        CNV_prob_init: numpy array (n_features, n_samples, n_CNVs)
            Initial value of CNV_prob, CNV state probability per gene/block per sample/bulk sample
        """
        # inital key parameters
        if random_state is not None:
            np.random.seed(seed =  random_state)
            print("random seed:", random_state)

        if alpha_s1_init is not None:
            self.alpha_s1 = alpha_s1_init
        else:
            self.alpha_s1 = (np.ones((self.n_features, self.n_samples, self.n_CNVs)) *
                (np.power(self.mu_init,2)/self.var_init))

        if beta_s2_init is not None:
            self.beta_s2 = beta_s2_init
        else:
            self.beta_s2 = np.ones((self.n_features, self.n_samples, self.n_CNVs)) * (self.mu_init/self.var_init)

        if CNV_prob_init is not None:
            self.CNV_prob = normalize(CNV_prob_init, axis = 2)
        else:
            self.CNV_prob = normalize(np.random.rand(self.n_features, self.n_samples, self.n_CNVs), axis=2)
        print(self.CNV_prob)

    def set_prior(self, CNV_prior=None, 
        alpha_s1_prior=None, beta_s2_prior=None, min_CNVs=0.00001):
        """Set prior for key variables: CNV_prob.
        The priors are in the same shape as its according variables.
        
        min_CNVs: float. 
            Minimun genotype probability in CNV_prior.
        """

        if CNV_prior is not None:
            CNV_prior[CNV_prior < min_CNVs] = min_CNVs
            CNV_prior[CNV_prior > 1 - min_CNVs] = 1 - min_CNVs
            CNV_prior = normalize(CNV_prior)
            self.CNV_prior = CNV_prior
        else:
            self.CNV_prior = normalize(np.ones(self.CNV_prob.shape))

        if alpha_s1_prior is None:
            alpha_s1_prior =(np.ones((self.n_features, self.n_samples, self.n_CNVs)) *
                (np.power(self.mu_init,2)/self.var_init))
        if beta_s2_prior is None:
            beta_s2_prior = np.ones(alpha_s1_prior.shape) * (self.mu_init/self.var_init)
        
        self.alpha_s1_prior = alpha_s1_prior
        self.beta_s2_prior = beta_s2_prior

    @property
    def gamma_exp1_(self):
        """Digamma of Beta concetration1 parameter"""
        if self.verbose:
            print("gamma_exp1_", digamma(self.alpha_s1)-ln(self.beta_s2))
        return digamma(self.alpha_s1)-ln(self.beta_s2)

    @property
    def gamma_exp2_(self):
        """Digamma of Beta concetration2 parameter"""
        if self.verbose:
            print("gamma_exp2_", self.alpha_s1/self.beta_s2)
        return self.alpha_s1/self.beta_s2
    
    def update_CNV_prob(self, data):
        """Coordinate ascent for updating CNV probability
        """
        r_counts = data

        # logLik_CNV = r_counts * ln(self.libratio) + r_counts * ln(self.ref_counts) + self.gamma_exp1_ * r_counts  - self.gamma_exp2_ * self.libratio * self.ref_counts - logfact(r_counts)
        logLik_CNV = r_counts * ln(self.libratio) + r_counts * ln(self.ref_counts) + self.gamma_exp1_ * r_counts  - self.gamma_exp2_ * self.libratio * self.ref_counts - gammaln(r_counts+1)
        if self.verbose:
            print("logLik_CNV 1:", logLik_CNV)
        # logLik_CNV = logLik_CNV.sum(axis=0,keepdims=True).sum(axis=1,keepdims=True) ## check
        
        if self.verbose:
            print("logLik_CNV 2:", logLik_CNV)
        self.CNV_prob = normalize(np.exp(loglik_amplify(
            logLik_CNV + np.log(self.CNV_prior))))
        self.CNV_prob_log = logLik_CNV + np.log(self.CNV_prior)
        
        return logLik_CNV

    def update_params(self, data):
        """Coordinate ascent for updating theta posterior parameters
        """
        r_counts = data

        _alpha_s1 = self.CNV_prob * r_counts
        _beta_s2 = self.CNV_prob * self.libratio * self.ref_counts
        
        _alpha_s1 += self.alpha_s1_prior.copy()
        _beta_s2 += self.beta_s2_prior.copy()
 
        self.alpha_s1 = _alpha_s1
        self.beta_s2 = _beta_s2


    def get_ELBO(self, logLik_CNV, data):## todo
        """calculate variational evidence lower bound with current parameters.

        The lower bound is used to detect the convergence and has to increase at
        each iteration.

        parameters
        ----------

        Returns
        -------
        lower_bound : float
        """
        r_counts = data
        if logLik_CNV is None:
            logLik_CNV = r_counts * ln(self.libratio) + r_counts * ln(self.ref_counts) + self.gamma_exp1_ * r_counts  - self.gamma_exp2_ * self.libratio * self.ref_counts  - gammaln(r_counts+1)

        LB_p = np.sum(logLik_CNV * self.CNV_prob)
        KL_CNV = -np.sum(entropy(self.CNV_prob, self.CNV_prior, axis=-1)) ## check
        
        KL_paras = -gamma_entropy(
        np.append(
            np.expand_dims(self.alpha_s1, 1),
            np.expand_dims(self.beta_s2, 1), axis = 1),
        np.append(
            np.expand_dims(self.alpha_s1_prior, 1),
            np.expand_dims(self.beta_s2_prior, 1), axis = 1))
        return LB_p + KL_CNV + KL_paras # FRAMEWORK

    def fit(self, data, max_iter=200, min_iter=5, epsilon_conv=1e-2,
        delay_fit_theta=0, verbose=False):## to do
        """Fit XClone model (RDR module) with coordinate ascent
        Parameters
        ----------
        r_counts : scipy.sparse.csc_matrix (n_genes, n_cells)
            Sparse count matrix for expression data
        max_iter : int
            Maximum number of iterations
        min_iter :
            Minimum number of iterations
        epsilon_conv : float
            Threshold for detecting convergence
        delay_fit_theta : int
            Number of steps to delay updating theta. This can be very useful
            for common genetics when there is good prior on allelic ratio.
        verbose : bool
            Whether print out log info
        """
        r_counts = data[:,:,np.newaxis]
         # check RDR's type: scipy.sparse or numpy.matrix
        # if type(r_counts) is np.ndarray and np.mean(r_counts > 0) < 0.3:
        #     print("Warning: input matrices is %.1f%% sparse, "
        #           %(100 - np.mean(r_counts > 0) * 100) +
        #           "change to scipy.sparse.csc_matrix" )
        #     r_counts = csc_matrix(r_counts)

        # TypeError: expected dimension <= 2 array or matrix
        # so scipy sparse matrix can not support here

        ELBO = np.zeros(max_iter)
        for it in range(max_iter):
            if self.learn_params and it >= delay_fit_theta:## to do
                self.update_params(r_counts)
            # if self.learn_CNV:
            #     self.update_CNV_prob(r_counts)
            
            _logLik_CNV = self.update_CNV_prob(r_counts)
            ELBO[it] = self.get_ELBO(_logLik_CNV, r_counts)# to do

            if it > min_iter:
                if ELBO[it] < ELBO[it - 1]:
                    if verbose:
                        print("Warning: Lower bound decreases!\n")
                elif it == max_iter - 1:
                    if verbose:
                        print("Warning: VB did not converge!\n")
                elif ELBO[it] - ELBO[it - 1] < epsilon_conv:
                    break
        self.ELBO_ = np.append(self.ELBO_, ELBO[:it])


class LinkedXClone(XCloneBase):
    """docstring for LinkedXClone."""
    def xxx():
        """func
        """
        pass
