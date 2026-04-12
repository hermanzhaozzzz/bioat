"""bioat.cli.

This module provides the command line interface (CLI) for the BioAT (Bioinformatic Analysis Tools) toolkit.
The BioAT toolkit is designed for bioinformatic analyses, allowing users to handle various biological data formats
including BED, BAM, FASTA, FASTQ, and VCF.

Usage Examples:
    Example 1:
        In shell:
        $ bioat list

        In Python console:
        >>> from bioat.cli import Cli
        >>> bioat = Cli()
        >>> bioat.list()
        >>> print(bioat.list())

    Example 2:
        In shell:
        $ bioat about
        In Python console:
        >>> from bioat.cli import Cli
        >>> bioat = Cli()
        >>> bioat.about()
        >>> print(bioat.about())



Description:
    The CLI interface allows users to access the functionalities of the BioAT package through
    terminal commands or import it directly in Python scripts. The module provides various
    tools for tasks such as mining CRISPRs, downloading metagenomic data, and retrieving
    search results from Google Scholar.

Commands:
    - about: Displays information about the BioAT toolkit.
    - list: Returns available groups and commands in a structured format.
    - version: Returns the version information of the BioAT toolkit.

Copyright:
    For researchers: Freely applicable for academic research; citation is appreciated but not mandatory.
    For commercial use: Not permitted without prior permission from the author.
"""

import ast
from functools import lru_cache
from importlib import import_module
from pathlib import Path

from bioat.about import __ABOUT__
from bioat._meta import __VERSION__

_TOOL_SPECS = {
    "bam": ("bioat.bamtools", "BamTools"),
    "bed": ("bioat.bedtools", "BedTools"),
    "crispr": ("bioat.crisprtools", "CrisprTools"),
    "fastx": ("bioat.fastxtools", "FastxTools"),
    "fold": ("bioat.foldtools", "FoldTools"),
    "hic": ("bioat.hictools", "HiCTools"),
    "meta": ("bioat.metatools", "MetaTools"),
    "search": ("bioat.searchtools", "SearchTools"),
    "table": ("bioat.tabletools", "TableTools"),
    "target_seq": ("bioat.target_seq", "TargetSeq"),
}


@lru_cache(maxsize=None)
def _get_public_tool_methods(module_name: str, class_name: str) -> tuple[str, ...]:
    module_path = Path(__file__).resolve().parent / f"{module_name.rsplit('.', 1)[1]}.py"
    tree = ast.parse(module_path.read_text(encoding="utf-8"))

    for node in tree.body:
        if isinstance(node, ast.ClassDef) and node.name == class_name:
            return tuple(
                child.name
                for child in node.body
                if isinstance(child, ast.FunctionDef) and not child.name.startswith("_")
            )

    return ()


class Cli:
    """Cli interface of BioAT.

    Brief:
        - "bioat" is short for "Bioinformatic Analysis Tools." It is a command-line toolkit and a Python package that can be used through this CLI interface in a terminal or via the `import` method in Python code.
        - "bioat" has many subcommands to handle different bio-formats: BED, BAM, FASTA, FASTQ, VCF, etc.
        - "bioat" can be used for mining CRISPRs, downloading metagenomes, and even reporting Google Scholar search results!
        - For more information, run `$ bioat about`.

    Copyright:
        For researchers: freely applied to academic research. Please cite my work:
            1. bibtex
            2. software copyright
        For commercial use:
            NOT PERMITTED unless permission is obtained from the author.

    ---

    COMMAND:
        Such as `bioat version`. This represents a direct command.

    GROUPS:
        Such as `bioat crispr`. This represents a group/bundle name for commands related to a specific topic.

    ---

    Usage:
        A demo to understand `bioat`:

            bioat about
                <in shell>:
                    $ bioat about
                <in python console>:
                    >>> from bioat.cli import Cli
                    >>> bioat = Cli()
                    >>> bioat.list()
                    >>> print(bioat.list())
    """

    # - Citation:
    #     - Citation is not forcible and I would be appreciated if you wann to
    #  do this.
    #     - citation reference format:
    #         [1] Huanan ZHAO. bioat, a Bioinformatic Analysis Tool kit (2023),
    # https://github.com/hermanzhaozzzz/bioat

    # """

    def __init__(self):
        self._tool_cache = {}

    def _load_tool(self, name: str):
        if name not in self._tool_cache:
            module_name, class_name = _TOOL_SPECS[name]
            tool_class = getattr(import_module(module_name), class_name)
            self._tool_cache[name] = tool_class()
        return self._tool_cache[name]

    @property
    def bam(self):
        return self._load_tool("bam")

    @property
    def bed(self):
        return self._load_tool("bed")

    @property
    def crispr(self):
        return self._load_tool("crispr")

    @property
    def fastx(self):
        return self._load_tool("fastx")

    @property
    def fold(self):
        return self._load_tool("fold")

    @property
    def hic(self):
        return self._load_tool("hic")

    @property
    def meta(self):
        return self._load_tool("meta")

    @property
    def search(self):
        return self._load_tool("search")

    @property
    def table(self):
        return self._load_tool("table")

    @property
    def target_seq(self):
        return self._load_tool("target_seq")

    @classmethod
    def about(cls):
        """Returns information about the `bioat` tool.

        This class method provides a description or metadata regarding
        the `bioat` application, which may include version information,
        author details, or usage instructions.

        Returns:
            str: Information about the `bioat` tool.
        """
        return __ABOUT__

    def list(
        self,
    ):
        """Returns a formatted string of GROUPS and COMMANDS.

        This method retrieves and formats the attributes and sub-attributes
        of the instance, excluding private members (those starting with "_").

        Returns:
            str: A tree formatted string representing the available
                 subcommands and their attributes.
        """
        out = ""
        for att in dir(self):
            if not att.startswith("_"):
                out += f"{att}\n"
                if att in _TOOL_SPECS:
                    module_name, class_name = _TOOL_SPECS[att]
                    for sub_att in _get_public_tool_methods(module_name, class_name):
                        out += f"  ├── {sub_att}\n"
        return out

    @classmethod
    def version(cls):
        """Return version information.

        Returns:
            str: The version information.
        """
        return __VERSION__
