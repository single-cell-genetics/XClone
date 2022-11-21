"""Base functions for XClone efficiency improvement
"""
# Author: Rongting Huang
# Date: 2022-03-17
# update: 2022-03-17

## Part I: 

import multiprocessing

def efficiency_preview():
    """
    Example:
    xclone.pp.efficiency_preview()
    """
    print("[XClone efficiency] multiprocessing cpu total count in your device", multiprocessing.cpu_count())
    return None

def packages_info():
    """
    print packages version used in XClone.
    """
    pass
