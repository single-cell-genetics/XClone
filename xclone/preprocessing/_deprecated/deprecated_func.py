## from _preprocessing


## Functions for phasing
## to do list-reconstructure the anndata
import multiprocessing
import datetime

import xclone
xclone.__version__

def local_phasing_wrap(AD, DP, regions, regions_mode, chr_lst, phasing_len = 100, nproc=1):
    start_t = datetime.datetime.now()
    if nproc > 1:
        result = []
        pool = multiprocessing.Pool(processes = nproc)
        for chr_ in chr_lst:
            ad_idx, AD_region = select_chr_region(AD, regions, "AD", regions_mode, [chr_], output_format="sp.sparse_mtx")
            dp_idx, DP_region = select_chr_region(DP, regions, "DP", regions_mode, [chr_], output_format="sp.sparse_mtx")
            result.append(pool.apply_async(process_region,(chr_, AD_region, DP_region),callback = None))
        pool.close()
        pool.join()
        result = [res.get() for res in result] 
    else:
        result = []
        for chr_ in chr_lst:
            ad_idx, AD_region = select_chr_region(AD, regions, "AD", regions_mode, [chr_], output_format="sp.sparse_mtx")
            dp_idx, DP_region = select_chr_region(DP, regions, "DP", regions_mode, [chr_], output_format="sp.sparse_mtx")
            RV = process_region(chr_, AD_region, DP_region, phasing_len = phasing_len, nproc=nproc)
            result.append(RV)
    # process the data from the list to a dataset
    ## need add the region annotation for the bins data
    for i , RV in zip(range(len(result)), result):
        if i == 0:
            ad_phase = RV['ad']
            dp_phase = RV['dp']
            phase_region = RV['id_record_bin']
        else:
            ad_phase = np.vstack((ad_phase, RV['ad']))
            dp_phase = np.vstack((dp_phase, RV['dp']))
            phase_region = np.append(phase_region, RV['id_record_bin'])
    
    end_t = datetime.datetime.now()
    elapsed_sec = (end_t - start_t).total_seconds()
    print("[Local_phasing] time_used: " + "{:.2f}".format(elapsed_sec) + "seconds")
    return ad_phase, dp_phase, phase_region

def process_region(ir, AD_region, DP_region, phasing_len = 100, nproc=1):
    n_bins = int(AD_region.shape[0] / phasing_len)
    ad_bin = np.zeros((n_bins, AD_region.shape[1]))
    dp_bin = np.zeros((n_bins, AD_region.shape[1]))
    theta_bin = np.zeros((n_bins, AD_region.shape[1]))
    id_record_bin = np.zeros((n_bins))
        
    if nproc > 1:
        result = []
        pool = multiprocessing.Pool(processes = nproc)
        for ib in range(n_bins):
            idx = range(ib*phasing_len, (ib+1)*phasing_len)
            result.append(pool.apply_async(process_bin,
                (ib, AD_region[idx, :], DP_region[idx, :]), 
                callback = None))
        pool.close()
        pool.join()
        result = [res.get() for res in result]  
    else:
        result = []
        for ib in range(n_bins):
            idx = range(ib*phasing_len, (ib+1)*phasing_len)
            RV_bin = process_bin(ib, AD_region[idx, :], DP_region[idx, :])
            result.append(RV_bin)
            
    for i, RV_bin in zip(range(len(result)), result):
        ad_bin[i, :], dp_bin[i, :] = RV_bin['ad_bin'], RV_bin['dp_bin']
        theta_bin[i, :] = RV_bin['theta_bin']
        id_record_bin[i] = RV_bin['idx']
    RV = {}
    RV['ir'] = ir
    RV['ad'] = ad_bin
    RV['dp'] = dp_bin
    RV['theta'] = dp_bin
    RV['id_record_bin'] = id_record_bin
    return RV

def process_bin(idx, AD_use, DP_use):
    ad_sum, dp_sum, Z, thetas, _logLik_new = xclone.model.Local_Phasing(AD_use, DP_use)
    RV = {}
    RV['idx'] = idx
    RV['ad_bin'] = ad_sum[0, :]
    RV['dp_bin'] = dp_sum[0, :]
    RV['theta_bin'] = np.array(thetas)[:, 0]
    RV['Z'] = Z
    RV['logLik'] = _logLik_new
    return RV
