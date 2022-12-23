"""Base functions for XClone HMM smoothing.
develop version.
"""

import numpy as np
import statsmodels.api as sm
from scipy.stats import nbinom
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
    
    AB_log = logsumexp(
        np.expand_dims(A_log, 2) + np.expand_dims(B_log, 0),
        axis=1
    )
    
    return AB_log


class Model_Base():
    pass
    
class Model_NB():
    def __init__(self, n_obs, n_var, n_sta, 
                 gene_specific = True, dispersion=None,
                 ref_value_init=None, lib_ratio_init=None, 
                 cnv_ratio_init=None):        
        """Initialise Model
        """
        self.n_obs = n_obs
        self.n_var = n_var
        self.n_sta = n_sta
        
        self.ref_value = np.ones((1, n_var, 1))
        self.lib_ratio = np.ones((n_obs, 1, 1))
        self.cnv_ratio = np.ones((1, 1, n_sta))
        
        if ref_value_init is not None:
            self.ref_value[0, :, 0] = ref_value_init
            
        if lib_ratio_init is not None:
            self.lib_ratio[:, 0, 0] = lib_ratio_init
            
        if cnv_ratio_init is not None:
            self.cnv_ratio[0, 0, :] = cnv_ratio_init
        
        self.dispersion = np.ones((1, n_var, 1)) * 0.1
        if dispersion is not None:
            if gene_specific:
                self.dispersion[0, :, 0] = dispersion
            else:
                self.dispersion[:, :, :] = dispersion
                # NOTE: only use one dispersion for all genes here
                    
    @property
    def mean(self):
        return self.ref_value * self.lib_ratio * self.cnv_ratio
        
    
    def get_prob_log(self, data):
        _var = self.mean + self.dispersion * self.mean**2
    
        _nb_prob  = self.mean / _var
        _nb_total = self.mean * _nb_prob / (1 - _nb_prob)
        
        emm_prob_log = nbinom.logpmf(np.expand_dims(data, 2), 
                                     _nb_total, _nb_prob)
        
        return emm_prob_log
    
    
    def fit(self, data, weights=None, options={'mode': 'cnv_ratio'}):
        """Fitting function with supporting multiple purposes
        """
        NB_kwargs = {'disp': True, 'skip_hessian': True}
        
        if options['mode'] == 'dispersion':
            print('Not supported yet.')
            
        elif options['mode'] == 'lib_ratio':
            print('Not supported yet.')
            
        elif options['mode'] == 'cnv_ratio':
            sm_results = []
            new_cnv_ratio = []
            for k in range(self.cnv_ratio.shape[2]):
                obs_y = data.reshape(-1)
                freq_w = weights[:, :, k].reshape(-1)
                feature_x = np.ones(len(obs_y))
                exposure_x = (self.ref_value * self.lib_ratio).reshape(-1)
                
                # print(obs_y.shape, feature_x.shape, exposure_x.shape, freq_w.shape)
                
                _dispersion = self.dispersion[0, 0, 0]
                # if k == 0:
                #     _dispersion = 5.0
                # elif k == 1:
                #     _dispersion = 1.0
                
                W_NB_glm = sm.GLM(
                    obs_y, feature_x,
                    exposure = exposure_x,
                    freq_weights = freq_w,
                    family=sm.families.NegativeBinomial(alpha=_dispersion)
                )

                NB_results = W_NB_glm.fit(
                    start_params = self.cnv_ratio[0, 0, k : (k+1)],
                    maxiter=100, **NB_kwargs
                )
                
                sm_results.append(NB_results)
                new_cnv_ratio.append(np.exp(NB_results.params[0]))
                
            self.cnv_ratio[0, 0, :] = new_cnv_ratio
            print(new_cnv_ratio)
            
        return sm_results
        

class HMM_Frame:
    def __init__(self, data, pi, A, model=None, var_groups=None):
        """
        data: a tuple of matrix or matrices e.g., ( X ) or (AD, DP)
        """
        self.A = A
        self.pi = pi
        self.data = data
        self.A_log = np.log(self.A)
        self.pi_log = np.log(self.pi)
        
        self.n_sta = len(pi)
        self.n_obs = data.shape[0]
        self.n_var = data.shape[1]
        
        if var_groups is None:
            self.var_groups = np.ones(self.n_var, dtype=int)
        else:
            self.var_groups = var_groups
            
        self.model = model        
        self.emm_p_log = np.zeros((self.n_obs, self.n_var, self.n_sta))
        
        
    def update_emm_prob(self):
        if self.model is None:
            print("No model is available, please specify first.")
        self.emm_p_log = self.model.get_prob_log(self.data)
        
    
    def M_step(self):
        if self.model is None:
            print("No model is available, please specify first.")
        self.model.fit(self.data, weights=np.exp(self.z_post_log))
        
    
    def update_posterior(self):
        """
        Forward-Backward algorithm will be used to calculate the 
        posterior of the state assignment via HMM
        """
        self.fw_p_log = self.emm_p_log.copy()
        self.bw_p_log = self.emm_p_log.copy()
        
        # Forward part of the algorithm
        for g in range(self.n_var):
            if g == 0 or self.var_groups[g] != self.var_groups[g - 1]:
                self.fw_p_log[:, g, :] += self.pi_log.reshape(1, -1)
            else:
                self.fw_p_log[:, g, :] += logdotexp(
                    self.fw_p_log[:, g-1, :], self.A_log
                )
        
        # Backward part of the algorithm
        for g in range(self.n_var):
            rg = self.n_var - 1 - g
            if rg == self.n_var - 1 or self.var_groups[rg] != self.var_groups[rg+1]:
                self.bw_p_log[:, rg, :] = 0
            else:
                self.bw_p_log[:, rg, :] = logdotexp(
                    self.bw_p_log[:, rg + 1, :] + self.emm_p_log[:, rg + 1, :], 
                    self.A_log
                )
                
        # Update posterior of state assignment
        self.z_post_log  = self.fw_p_log + self.bw_p_log
        self.z_post_log -= logsumexp(self.z_post_log, axis=2, keepdims=True)
        
        # LogLikehood to monitor
        self.logLik_mat = logsumexp(self.fw_p_log, axis=2)
        self.logLik = np.sum(self.logLik_mat)
        
    
    def EM_fit(self, KNN_connect=None, min_iter=1, max_iter=100, epsilon_conv=1e-2):
        """EM algorithm to fit HMM
        """
        logLik_all = np.zeros(max_iter)
        
        for it in range(max_iter):
            # E step:
            self.update_emm_prob()
            self.update_posterior()
            logLik_all[it] = self.logLik + 0.0
            
            # KNN smoothing across obs
            if KNN_connect is not None:
                for k in range(self.z_post_log.shape[2]):
                    # NOTE: may consider logdotexp
                    self.z_post_log[:, :, k] = KNN_connect @ self.z_post_log[:, :, k]
            
            # M step:
            self.M_step()
            
            # Convergence check
            if it >= min_iter:
                if logLik_all[it] < logLik_all[it - 1]:
                    print("Warning: step %d, loglik decrease from %.2e to %.2e" 
                          %(it, logLik_all[it-1], logLik_all[it]))         
                elif logLik_all[it] - logLik_all[it - 1] < epsilon_conv:
                    break
                    
        self.logLik_all = logLik_all[:it]