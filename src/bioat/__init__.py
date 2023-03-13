"""Doc of bioat.

"""
from __future__ import annotations, absolute_import
from .version import (__version__, __author__, __doc_format__)
from .logger import set_logging_level
from .fastx import Fastx
from .table import Table
from .hic import HiC
from .mgi import Mgi
from .system import System

__name__ = "bioat"
