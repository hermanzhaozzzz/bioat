"""BioAT package.
BioAT can be a package to import.
BioAT also can be a command-line tool.
It is a bioinformatic tool/pkg bundle for python.
"""

from importlib import import_module

from bioat._meta import (
    __AUTHOR__,
    __AUTHOR_EMAIL__,
    __DESCRIPTION__,
    __DOC_FORMAT__,
    __DOC_PAGE__,
    __HOME_PAGE__,
    __ISSUE_PAGE__,
    __LICENSE__,
    __PKG_NAME__,
    __VERSION__,
)

_LAZY_ATTRS = {
    "BamTools": ("bioat.bamtools", "BamTools"),
    "BedTools": ("bioat.bedtools", "BedTools"),
    "CrisprTools": ("bioat.crisprtools", "CrisprTools"),
    "FastxTools": ("bioat.fastxtools", "FastxTools"),
    "FoldTools": ("bioat.foldtools", "FoldTools"),
    "HiCTools": ("bioat.hictools", "HiCTools"),
    "MetaTools": ("bioat.metatools", "MetaTools"),
    "SearchTools": ("bioat.searchtools", "SearchTools"),
    "SystemTools": ("bioat.systemtools", "SystemTools"),
    "TableTools": ("bioat.tabletools", "TableTools"),
    "TargetSeq": ("bioat.target_seq", "TargetSeq"),
}


def __getattr__(name: str):
    if name not in _LAZY_ATTRS:
        msg = f"module {__name__!r} has no attribute {name!r}"
        raise AttributeError(msg)

    module_name, attr_name = _LAZY_ATTRS[name]
    value = getattr(import_module(module_name), attr_name)
    globals()[name] = value
    return value


def __dir__():
    return sorted(set(globals()) | set(__all__))


__all__ = [
    "__AUTHOR_EMAIL__",
    "__AUTHOR__",
    "__DESCRIPTION__",
    "__DOC_FORMAT__",
    "__DOC_PAGE__",
    "__HOME_PAGE__",
    "__ISSUE_PAGE__",
    "__LICENSE__",
    "__PKG_NAME__",
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
