"""Base class for mixture models and EM algorithm."""

# This file contains a few mixture models and its combinations
# Author: Yuanhua Huang
# Date: 02/12/2020

import numpy as np
from scipy.stats import entropy
from scipy.optimize import linear_sum_assignment
from scipy.special import logsumexp, digamma, betaln, binom
from scipy.sparse import issparse

from .base_utils import normalize, loglik_amplify

class MixtureBase():
    """Base object for mixture model
    """
    def __init__(self, n_components=1, n_samples=1, learn_pi=False):
        """Initializing the model
        """
        self.n_samples = n_samples
        self.n_components = n_components
        self.params = np.array([]) # empty parameters
        self.pi = np.ones(self.n_components) / self.n_components
        self.learn_pi = learn_pi
        
        self.random_initial()
        
    def random_initial(self, random_state=None):
        """Random initialization of the assignment probability and 
        key parameters
        """
        if random_state is not None:
            np.random.seed(seed = random_state)
            
        self.proba_ = np.random.dirichlet(np.ones(self.n_components), 
                                          size=self.n_samples)
        
    def update_logLik(self, data):
        """Caculate log likelihood with current parameters
        
        Require re-implementation for specific model
        """
        print("Warning: Base model's update_logLik() hasn't been defined!")
        self.logLik_mat_ = np.zeros((self.n_samples, self.n_components))
        self.logLik = 0
    
    def E_step(self, data):
        """E step: get the expectation of the assignment
        """
        self.update_logLik(data)
        self.proba_ = normalize(np.exp(loglik_amplify(self.logLik_mat_)))
        
    def M_step(self, data):
        """M ste: update the parameters to maximizing the logLikelihood
        
        Require re-implementation for specific model
        """
        print("Warning: Base model's M_step() hasn't been defined!")
        # self.params = None
        
    def _EM_fit(self, data, min_iter=10, max_iter=1000, epsilon_conv=1e-2):
        """General framework to fit the model with EM algorithm
        """
        ## EM algorithm
        logLik_all = np.zeros(max_iter)
        for it in range(max_iter):        
            self.M_step(data)
            self.E_step(data)
            
            logLik_all[it] = self.logLik + 0
            
            # convergence
            if it >= min_iter:
                if logLik_all[it] - logLik_all[it - 1] < epsilon_conv:
                    break
                elif logLik_all[it] < logLik_all[it - 1]:
                    print("Warning: loglik decrease from %.2e to %.2e" 
                          %(logLik_all[it-1] < logLik_all[it]))
                    
        self.logLik_all = logLik_all[:it]
        
        
    def fit(self, data, n_init=4, **kwargs):
        """Fit the model with multiple initialization to avoid local optioma
        """
        ## Fit the first initialization
        self.random_initial()
        self._EM_fit(data, **kwargs)
        
        _proba_best = self.proba_ + 0
        _logLik_best = self.logLik + 0
        _logLik_all_best = self.logLik_all + 0
        self.logLik_inits = [_logLik_best]
        
        for i in range(n_init - 1):
            self.random_initial()
            self._EM_fit(data, **kwargs)
            self.logLik_inits.append(self.logLik + 0)
            
            if self.logLik > _logLik_best:
                _proba_best = self.proba_ + 0
                _logLik_best = self.logLik + 0
                _logLik_all_best = self.logLik_all + 0
                
        ## Re-fit with best parameters
        self.random_initial()
        self.proba_ = _proba_best + 0
        self._EM_fit(data)
        # self.logLik_all = np.append(_logLik_all_best, self.logLik_all)
    
    @property
    def n_parameters(self):
        """Calculate number of parameters
        """
        _n_parameters = self.params.size 
        if self.learn_pi:
            _n_parameters += self.n_components
        return _n_parameters

    def get_BIC(self):
        """Get BIC = n_parameters * log(n_samples) - 2 * logLik
        """
        if hasattr(self, 'logLik') == False:
            print("Error: no logLik. Please fit the model and update_logLik.")
            return None
        elif self.logLik == None:
            print("Error: logLik is None. Please check update_logLik.")
            return None

        return self.n_parameters * np.log(self.n_samples) - 2 * self.logLik
    


class BinomialMixture(MixtureBase):
    """Binomial mixture model:
    
    Base model of AD[i, j] and DP[i,j] in component k: 
        Binomial(AD[i,j] | DP[i,j], self.params[i, k])
    """
    def update_logLik(self, data):
        """Get loglikelihood matrix
        
        Args
        ----
        data: a list of two elements, [AD, DP]. Both are numpy.array or 
            scipy.sparse matrix, (n_dimensions, n_samples)
        """
        AD, DP = data[0], data[1]
        BD = DP - AD
        self.logLik_mat_ = (
            AD.T @ np.log(self.params) + 
            BD.T @ np.log(1 - self.params)
        )
        self.logLik_mat_ += np.log(self.pi).reshape(1, -1)

        self.logLik = np.sum(logsumexp(self.logLik_mat_, axis=-1))
        
        ### binomial coefficient as a constant
        # from .XClone_VB_base import get_binom_coeff
        # self.logLik += get_binom_coeff(AD, DP, is_log=True).sum()
        
        
    def M_step(self, data):
        """M ste: update the parameters to maximizing the logLikelihood
        
        Args
        ----
        data: a list of two elements, [AD, DP]. Both are numpy.array or 
            scipy.sparse matrix, (n_dimensions, n_samples)
        
        Updates
        -------
        self.params: numpy.array, (n_dimension, n_components)
        """
        AD, DP = data[0], data[1]
        self.params = (AD @ self.proba_ + 0.01) / (DP @ self.proba_ + 0.02)
        if self.learn_pi:
            self.pi = self.proba_.mean(axis=0)
        

        
class PoissonMixture(MixtureBase):
    """Poisson mixture model:
    
    Base model of X[i, j] in component k: 
        Poisson(X[i,j] | self.params[i, k] * self.size_factors[j])
    """        
    def random_initial(self, random_state=None):
        """Random initialization of the assignment probability and 
        key parameters
        """
        if random_state is not None:
            np.random.seed(seed = random_state)
            
        self.proba_ = np.random.dirichlet(np.ones(self.n_components), 
                                          size=self.n_samples)
        self.size_factors = np.ones((self.n_samples, 1))
    
    def update_logLik(self, X):
        """Get loglikelihood matrix
        
        Args
        ----
        X: numpy.array or scipy.sparse matrix, (n_dimensions, n_samples) 
        """
        if type(X) == np.matrix:
            X = X.A
            
        if issparse(X):
            X_dot_size_factors = X.T.multiply(np.log(self.size_factors))
            X_dot_size_factors_sum = X_dot_size_factors.sum(1).A
        else:
            X_dot_size_factors = np.multiply(X.T, np.log(self.size_factors))
            X_dot_size_factors_sum = X_dot_size_factors.sum(1, keepdims=True)
            
        # logLik_mat_: (n_samples, n_components)
        self.logLik_mat_ = (
            X.T @ np.log(self.params) -
            self.size_factors @ self.params.sum(0, keepdims=True) +
            X_dot_size_factors_sum)
        self.logLik_mat_ += np.log(self.pi).reshape(1, -1)

        self.logLik = np.sum(logsumexp(self.logLik_mat_, axis=-1))
        
    def M_step(self, X):
        """M ste: update the parameters to maximizing the logLikelihood
        
        Args
        ----
        X: numpy.array or scipy.sparse matrix, (n_dimensions, n_samples)
            
        Updates
        -------
        self.params: numpy.array, (n_dimension, n_components)
        self.size_factors: numpy.array, (n_samples, 1)
        """
        if type(X) == np.matrix:
            X = X.A
            
        if issparse(X):
            X_colsum = X.sum(0).A
        else:
            X_colsum = X.sum(0, keepdims=True)
        
        self.params = (X @ self.proba_) / (self.size_factors.T @ self.proba_)
        
        self.size_factors = (
            X_colsum.T / (self.proba_ @ self.params.sum(0, keepdims=True).T)
        )

        if self.learn_pi:
            self.pi = self.proba_.mean(axis=0)


        
class MultinomialMixture(MixtureBase):
    """Multinomial mixture model:
    
    Base model of X[:, j] in component k:
        Multinomial(X[:, j] | self.params[:, k])
    """
    def update_logLik(self, X):
        """Get loglikelihood matrix
        
        Args
        ----
        X: numpy.array or scipy.sparse matrix, (n_dimensions, n_samples)
        """
        self.logLik_mat_ = X.T @ np.log(self.params)
        self.logLik_mat_ += np.log(self.pi).reshape(1, -1)
        self.logLik = np.sum(logsumexp(self.logLik_mat_, axis=-1))
        
    def M_step(self, X):
        """M ste: update the parameters to maximizing the log Likelihood
        
        Args
        ----
        X: numpy.array or scipy.sparse matrix, (n_dimensions, n_samples)
            
        Updates
        -------
        self.params: numpy.array, (n_dimension, n_components)
        """
        self.params = normalize(X @ self.proba_, axis=0)
        if self.learn_pi:
            self.pi = self.proba_.mean(axis=0)
        

        
class LinkedMixture(MixtureBase):
    """A general framework for mixture models with multiple modules.
    """ 
    def __init__(self, n_components=1, n_samples=1, learn_pi=False, models=[]):
        """Initializing the model
        """
        self.n_samples = n_samples
        self.n_components = n_components
        self.pi = np.ones(self.n_components) / self.n_components
        self.learn_pi = learn_pi

        self.models = models
        self.check_models()

    def check_models(self):
        """Check models are consistent
        """
        for im, _model in enumerate(self.models):
            if _model.n_samples != self.n_samples:
                print("Warning: model %d has different n_samples." %(im))
            if _model.n_components != self.n_components:
                print("Warning: model %d has different n_components." %(im))
            if _model.learn_pi != self.learn_pi:
                print("Warning: model %d has different learn_pi." %(im))
    
    def random_initial(self, random_state=None):
        """Random initialization of the assignment probability and 
        key parameters
        """
        if random_state is not None:
            np.random.seed(seed = random_state)
            
        for i in range(len(self.models)):
            self.models[i].random_initial()
            
        self.proba_ = np.random.dirichlet(np.ones(self.n_components), 
                                          size=self.n_samples)
        
    def update_logLik(self, datasets):
        """Get loglikelihood matrix
        
        Args
        ----
        datasets: a list of data, each element is for a model in same order
        """
        self.logLik_mat_ = np.zeros((self.n_samples, self.n_components))
        for i in range(len(self.models)):
            self.models[i].update_logLik(datasets[i])
            self.logLik_mat_ += self.models[i].logLik_mat_

        # deduct over counted pi; should only count once
        self.logLik_mat_ -= (
            np.log(self.pi).reshape(1, -1) * (len(self.models) - 1)
        )
            
        self.logLik = np.sum(logsumexp(self.logLik_mat_, axis=-1))
        
    def M_step(self, datasets):
        """M ste: update the parameters to maximizing the log Likelihood
        
        Args
        ----
        datasets: a list of data, each element is for a model in same order
            
        Updates
        -------
        self.models[i].params
        """
        for i in range(len(self.models)):
            self.models[i].M_step(datasets[i])        
        
        if self.learn_pi:
            self.pi = self.proba_.mean(axis=0)
        
    def E_step(self, datasets):
        """E step: get the expectation of the assignment
        
        Args
        ----
        datasets: a list of data, each element is for a model in same order
            
        Updates
        -------
        self.models[i].proba_
        """
        self.update_logLik(datasets)
        self.proba_ = normalize(np.exp(loglik_amplify(self.logLik_mat_)))
        for i in range(len(self.models)):
            self.models[i].proba_ = self.proba_ + 0

    @property
    def n_parameters(self):
        """Calculate number of parameters
        """
        self.n_parameters = 0
        for _model in self.models:
            self.n_parameters += _model.params.size

        if np.mean(self.learn_pi) > 0:
            self.n_parameters += self.n_components


class HMMixture(MixtureBase):
    """A general framework for mixture models with multiple modules.
    Rongting
    """
    pass
