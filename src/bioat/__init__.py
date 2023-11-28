"""Doc of bioat.

"""
from bioat.exceptions import BioatFileFormatError, BioatFileNameError, BioatParameterFormatError # must be the first
from bioat.logger import get_logger # must be the first
from bioat.about import about
from bioat.bamtools import BamTools
from bioat.bedtools import BedTools
from bioat.crisprtools import CrisprTools
from bioat.fastxtools import FastxTools
from bioat.hictools import HiCTools
from bioat.metatools import MetaTools
from bioat.searchtools import SearchTools
from bioat.systemtools import SystemTools
from bioat.tabletools import TableTools
from bioat.target_seq import TargetSeq
from bioat.version import (__version__, __author__, __doc_format__)

__name__ = "bioat"
