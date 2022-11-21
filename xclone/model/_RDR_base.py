"""Base functions for XClone RDR processing
"""

# Author: Rongting Huang
# Date: 2022/07/12
# update: 2022/07/12

import numpy as np


def list_update(lst, bool_flag, reverse = False):
    """
    Function:
    update lst based on the bool_flag.

    Parameters:
    ----------
    lst: list.
    bool_flag: np.array.
    reverse: bool.

    Return:
    ------
    lst_update: list

    Example:
    -------
    update_model_lst = list_update(Xdata.uns["fit_lib_ratio_model"].copy(), ~FLAG_)
    in FUNC `remove_cells`.

    """
    if reverse:
        bool_flag = ~bool_flag
    lst_update = []
    lst_index = np.array(range(len(bool_flag)))[bool_flag]
    for idx_ in lst_index:
        lst_update.append(lst[idx_])
    return lst_update