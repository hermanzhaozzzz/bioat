"""BioAT package.
BioAT can be a package to import.
BioAT also can be a command-line tool.
It is a bioinformatic tool/pkg bundle for python.
"""

from bioat.bamtools import BamTools
from bioat.bedtools import BedTools
from bioat.crisprtools import CrisprTools
from bioat.exceptions import *  # !must be the first
from bioat.fastxtools import FastxTools
from bioat.hictools import HiCTools
from bioat.logger import __PKG_NAME__  # !must be the second
from bioat.metatools import MetaTools
from bioat.searchtools import SearchTools
from bioat.systemtools import SystemTools
from bioat.tabletools import TableTools
from bioat.target_seq import TargetSeq
from bioat.version import __author__, __doc_format__, __version__
