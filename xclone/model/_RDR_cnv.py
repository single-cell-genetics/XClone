import datetime
import numpy as np
import scipy as sp

from ._RDR_CNVratio  import hard_assign_prob

## developing
def fit_CNV_ratio_update(Xdata, 
                  init_prob, hard_assign = True, 
                  n_sample_cells = 1000,
                  n_top_genes = 15,
                  feature_X = None, 
                  depth_key = "library_ratio", 
                  dispersion_set = None, random_seed = None, verbose=True,
                  NB_kwargs={'disp': True, 'skip_hessian': True}):
    """
    Function:
    version 0.0.2-try to fit from previous identified CNV regions.

    Based on selected group genes and GLM model to get the estimate CNV ratio
    for different state, default state is copy loss, neutral, copy gain.
    GLM loop is for different gene groups. 

    We estimate CNV state-specific ratio for each gene group, for each gene we
    sample the cells from the cells pool.
 
    Here use weighted GLM.
    https://www.statsmodels.org/devel/examples/notebooks/generated/glm_weights.html#
    https://www.statsmodels.org/dev/examples/notebooks/generated/glm.html
    https://www.statsmodels.org/dev/_modules/statsmodels/discrete/discrete_model.html#NegativeBinomial

    Parameters:
    -----------
    Xdata: anndata.
    init_prob: np.array.
               probability for each gene across cells for all CNV states.
    hard_assign: bool. 
                 preprocessing for init_prob to get the hard assign prob or not.
    n_sample_cells: int.
    dispersion_set: float.
    random_seed: int.      
    verbose : bool
            Whether print out log info
    test_func: bool.
    
    Return:
    ------
    Xdata: anndata.


    Example:
    -------

    """
    import statsmodels.api as sm
    import scipy as sp
    print("[XClone] fit CNV ratio")

    if sp.sparse.issparse(Xdata.X):
        Xmtx = Xdata.X.A
    else:
        Xmtx = Xdata.X
    
    if random_seed is not None:
        rvgen = np.random.RandomState(random_seed)
    else:
        rvgen = np.random
    
    # time stamp
    start_time_ = datetime.datetime.now()
    
    ## for the first round and no region selection
    ## hard_assign should be true
    if hard_assign == True:
        if verbose:
            print("[XClone] fit_CNV_ratio: hard_assign prob")
        init_prob = hard_assign_prob(init_prob)

    # check input
    if Xdata.shape[0] == init_prob.shape[0]:
        if Xdata.shape[1] == init_prob.shape[1]:
            if verbose:
                print("input data is in matched dimension")
        else:
            raise ValueError("[XClone] gene dimension is not matched! Pls check!")
    else:
        raise ValueError("[XClone] cell dimension is not matched! Pls check!")
    
    ## data by different gene_expression groups
    group_genes_lst = Xdata.uns["group_genes"]
    n_gene_group = len(group_genes_lst)
    n_states = init_prob.shape[2]
    n_cells = Xdata.shape[0]

    model_results = []
    CNV_ratio = []
    ref_bio_ratio = []

    ## GLM model
    for i in range(n_gene_group):
        flag_ = Xdata.var["GeneName"].isin(group_genes_lst[i])

        ref_vec = np.array(Xdata.var["ref_avg"][flag_])
        tmp_mtx = Xmtx[:,flag_]
        tmp_prob = init_prob[:,flag_, :]
        
        libratio = Xdata.obs[depth_key]
        
        if dispersion_set is None:
            # dispersion_ = Xdata.var["dispersion"][flag_].median()
            dispersion_ = Xdata.var["dispersion"][flag_]
            print(dispersion_)
        else:
            dispersion_ = dispersion_set

        for state_idx in range(n_states):
            # summary of the high prob states across all cells
            prob_sum_ = (np.argmax(tmp_prob, axis = -1) == state_idx).sum(axis=0)
            rank_idx = prob_sum_.argsort()[::-1][1:n_top_genes+1] 
            # the 1st can be false positive, removed
            gene_idx = np.repeat(rank_idx, n_sample_cells)

            ref_exp_ = ref_vec[gene_idx]

            sample_cell_idx = rvgen.randint(0, n_cells, n_sample_cells*n_top_genes)

            gene_exp_ = tmp_mtx[sample_cell_idx, gene_idx]
            freq_prob_ = tmp_prob[sample_cell_idx, gene_idx, :]
        
            if dispersion_set is None:
                dispersion_used = dispersion_[gene_idx]
            else:
                dispersion_used = dispersion_set
            
            libratio_used = libratio[sample_cell_idx]

            obs_y = gene_exp_.reshape(-1)

            freq_w = freq_prob_[:, state_idx]

            if feature_X is None:
                feature_x = np.ones(len(obs_y))
            else:
                feature_x1 = np.ones(len(obs_y))
                feature_x2 = feature_X[sample_cell_idx,gene_idx].reshape(-1)
                feature_x = np.vstack((feature_x1, feature_x2)).T
            #should take libratio into account
            exposure_x = np.array(ref_exp_ * libratio_used).reshape(-1) 

            try:
            # if True:
                W_NB_glm = sm.GLM(obs_y, 
                                  feature_x,
                                  exposure = exposure_x,
                                  family=sm.families.NegativeBinomial(alpha=dispersion_used),
                                  freq_weights = freq_w)
            
                W_NB_res = W_NB_glm.fit(**NB_kwargs)
            except:
                print("[XClone] GLM error:")
                print(obs_y, feature_x, exposure_x, freq_w)
                model_results.append("Error in this gene group")
                tmp_ratio = np.NaN
                CNV_ratio.append(tmp_ratio)
                if feature_X is None:
                    pass
                else:
                    tmp_ratio1 = np.NaN
                    ref_bio_ratio.append(tmp_ratio1)
            else:
                print("[XClone] GLM success:")
                print(obs_y, feature_x,exposure_x, freq_w)
                model_results.append(W_NB_res)
                tmp_ratio = np.exp(W_NB_res.params[0])
                CNV_ratio.append(tmp_ratio)
                if feature_X is None:
                    pass
                else:
                    tmp_ratio1 = np.exp(W_NB_res.params[1])
                    ref_bio_ratio.append(tmp_ratio1)
                

    CNV_ratio = np.array(CNV_ratio).reshape(-1, n_states)

    # time stamp
    end_time_ = datetime.datetime.now()
    time_used = end_time_ - start_time_
    print("Time used", time_used.seconds, "seconds")
    
    if feature_X is None:
        return CNV_ratio, model_results
    else:
        ref_bio_ratio = np.array(ref_bio_ratio).reshape(-1, n_states)

        return CNV_ratio, ref_bio_ratio, model_results