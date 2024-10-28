"""
bioat.cli

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

from bioat import (
    __VERSION__,
    BamTools,
    BedTools,
    CrisprTools,
    FastxTools,
    FoldTools,
    HiCTools,
    MetaTools,
    SearchTools,
    # SystemTools,
    TableTools,
    TargetSeq,
)
from bioat.about import __ABOUT__


class Cli(object):
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
        self.bam = BamTools()
        self.bed = BedTools()
        self.crispr = CrisprTools()
        self.fastx = FastxTools()
        self.fold = FoldTools()
        self.hic = HiCTools()
        self.meta = MetaTools()
        self.search = SearchTools()
        # self.system = SystemTools()
        self.table = TableTools()
        self.target_seq = TargetSeq()

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
        # ! do not change OBJ method to CLASS method
        # ! because some attr and methods do not exit before obj instantiation
        out = ""
        for att in dir(self):
            if not att.startswith("_"):
                out += f"{att}\n"

                for sub_att in dir(getattr(self, att)):
                    if not sub_att.startswith("_"):
                        out += f"  ├── {sub_att}\n"
        return out

    @classmethod
    def version(cls):
        """Return version information.

        Returns:
            str: The version information.
        """
        return __VERSION__
