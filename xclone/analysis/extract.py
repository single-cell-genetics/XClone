"""Base functions for XClone result extraction.
"""

def dir_make(out_dir):
    import os
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

import pandas as pd
import numpy as np
def extract_xclone_matrix(Xdata, Xlayer="prob1_merge", states = None, region_lst = ["chr", "start", "stop"], 
                          index = False, header = False, out_dir = None):
    """
    extract cell by features matrix from xclone
    state = 0: copy loss
    state = 1: loh
    state = 2: copy neutral
    state = 3: copy gain
    state = None: argmax CNV state.
    """
    if out_dir is not None:
        dir_make(out_dir)
        ## mtx file
        mtx_file = out_dir + "matrix.csv"
        if states is not None:
            pd.DataFrame(Xdata.layers[Xlayer][:,:,states]).to_csv(mtx_file, index = index, header = header)     
        else:
            cnv_state = np.argmax(Xdata.layers[Xlayer], axis =2)
            pd.DataFrame(cnv_state).to_csv(mtx_file, index = index, header = header)
        ## cells file
        cells_file = out_dir + "cells.csv"
        pd.Series(Xdata.obs.index).to_csv(cells_file, index = index, header = header)
        ## feature file
        feature_file = out_dir + "Features.csv"
        if region_lst is not None:
            Xdata.var[region_lst].to_csv(feature_file, index = index)
        else:
            Xdata.var.to_csv(feature_file, index = index)
    else:
        print("pls provide out_dir")
    return None