"""BioAT package.
BioAT can be a package to import.
BioAT also can be a command-line tool.
It is a bioinformatic tool/pkg bundle for python.
"""
# !use "# ruff: isort: skip_file" annotation to skip this file when sorting imports
# ruff: isort: skip_file
from bioat._meta import (
    __PKG_NAME__,
    __AUTHOR__,
    __AUTHOR_EMAIL__,
    __DESCRIPTION__,
    __DOC_FORMAT__,
    __DOC_PAGE__,
    __HOME_PAGE__,
    __ISSUE_PAGE__,
    __LICENSE__,
    __VERSION__,
)

from bioat.bamtools import BamTools
from bioat.bedtools import BedTools
from bioat.crisprtools import CrisprTools

from bioat.fastxtools import FastxTools
from bioat.foldtools import FoldTools
from bioat.hictools import HiCTools

from bioat.metatools import MetaTools
from bioat.searchtools import SearchTools
from bioat.systemtools import SystemTools
from bioat.tabletools import TableTools
from bioat.target_seq import TargetSeq


__all__ = [
    "__PKG_NAME__",
    "__AUTHOR__",
    "__AUTHOR_EMAIL__",
    "__DESCRIPTION__",
    "__DOC_FORMAT__",
    "__DOC_PAGE__",
    "__HOME_PAGE__",
    "__ISSUE_PAGE__",
    "__LICENSE__",
    "__VERSION__",
    "BamTools",
    "BedTools",
    "CrisprTools",
    "FastxTools",
    "FoldTools",
    "HiCTools",
    "MetaTools",
    "SearchTools",
    "SystemTools",
    "TableTools",
    "TargetSeq",
]
