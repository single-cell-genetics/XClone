from .version import __version__
print (f'(Running XClone {__version__})')

# set simplified alias
from . import preprocessing as pp
from . import plot as pl
from . import analysis as al
from . import model
from . import data

from .model import phasing as phasing
from .model import mixture as mixture

from ._logging import get_logger
from ._config import XCloneConfig, settings