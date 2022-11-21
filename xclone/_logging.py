"""Logging and Profiling in XClone

Author: Rongting Huang
"""

import logging
# import anndata.logging

def set_main_logger(settings, module_name, level = logging.INFO):
    root_logger = logging.getLogger(module_name)
    root_logger.setLevel(level)
    
    def _set_log_file(settings, level):
        file = settings.logfile
        name = settings.logpath

        h = logging.StreamHandler(file) if name is None else logging.FileHandler(name)
        formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s - (%(filename)s).%(funcName)s(%(lineno)d)')
        h.setFormatter(formatter)
        h.setLevel(level)
        return h
    
    # create the logging file handler
    fh = _set_log_file(settings, level)
    
    # fh = logging.FileHandler("xclone.log")
    # formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    # fh.setFormatter(formatter)
    # add handler to logger object
    root_logger.addHandler(fh)
    return root_logger

def set_sub_logger(module_name):
    sub_logger = logging.getLogger(module_name)
    return sub_logger


import logging

_log_format = f"%(asctime)s - [%(levelname)s] - %(name)s  - %(message)s - (%(filename)s).%(funcName)s(%(lineno)d)"

def get_file_handler():
    file_handler = logging.FileHandler("XClone.log")
    file_handler.setLevel(logging.INFO)
    file_handler.setFormatter(logging.Formatter(_log_format))
    return file_handler

def get_stream_handler():
    stream_handler = logging.StreamHandler()
    stream_handler.setLevel(logging.INFO)
    stream_handler.setFormatter(logging.Formatter(_log_format))
    return stream_handler

def get_logger(name):
    logger = logging.getLogger(name)
    logger.setLevel(logging.INFO)
    logger.addHandler(get_file_handler())
    # logger.addHandler(get_stream_handler())
    return logger