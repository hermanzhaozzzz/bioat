"""Doc of bioat.

"""
from __future__ import annotations, absolute_import
from .about import about
from .fastx import Fastx
from .hic import HiC
from .logger import set_logging_level
from .mgi import Mgi
from .system import System
from .table import Table
from .version import (__version__, __author__, __doc_format__)

__name__ = "bioat"
