"""XClone: Inference of clonal Copy Number Alterations in single cells.
"""

from .version import __version__
from time import gmtime, strftime
print (f'(Running XClone {__version__})')
print (strftime("%Y-%m-%d %H:%M:%S", gmtime()))

# set simplified alias
from . import preprocessing as pp
from . import plot as pl
from . import analysis as al
from . import model
from . import data

from .model import phasing as phasing
from .model import mixture as mixture

from ._logging import get_logger
from ._config import PreprocessingConfig, XCloneConfig, settings
from ._config import XCloneGeneral_config, RDR_General_config, BAF_General_config, Combine_General_config
from ._config import HMM_Configs, Smartseq_Config, Spatial_Config