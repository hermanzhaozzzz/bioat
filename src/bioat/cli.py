"""The Cli bioat command entry module

author: Herman Huanan Zhao
email: hermanzhaozzzz@gmail.com
homepage: https://github.com/hermanzhaozzzz

>>> example 1:
    in shell:
    $ bioat list

    in python consolo
    >>> from bioat.cli import Cli
    >>> bioat = Cli()
    >>> bioat.list()
    >>> print(bioat.list())
>>> example 2:
    
"""

from bioat import (  # SystemTools,
    BamTools,
    BedTools,
    CrisprTools,
    FastxTools,
    HiCTools,
    MetaTools,
    SearchTools,
    TableTools,
    TargetSeq,
    __version__,
)
from bioat.about import __ABOUT__


class Cli(object):
    """Cli interface of BioAT

    <Brief>
    - "bioat" is short for "Bioinformatic Analysis Tools" and is a commandline toolkit and a python package. It can be used through this cli interface in terminal or the `import` way in python codes.
    - "bioat" has many subcommand to deal with different bio-format: BED, BAM, FASTA, FASTQ, VCF, et al.
    - "bioat" can be used for mining CRISPRs, downloading metagenome and even report a Google Scholar search result!
    - for more information, run `$ bioat about`

    <Copyright>
    For researchers: freely applied to academic research, please cite my work:
            1. #TODO bibtex
            2. #TODO software copyright
    For commercial use:
        NOT PERMITTED unless permission is obtained from the author

    ---

    <COMMAND>
    Such as `bioat version`, it means a directly command

    <GROUPS>
    Such as `bioat crispr`, it is a group/bundle name for commands with a specific topic

    ---

    <Usage>
    A Demo to understand `bioat`:

        bioat about
            <in shell>:
                $ bioat about
            <in python consolo>:
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
        self.hic = HiCTools()
        self.meta = MetaTools()
        self.search = SearchTools()
        # self.system = SystemTools()
        self.table = TableTools()
        self.target_seq = TargetSeq()

    @classmethod
    def about(cls):
        """Return information about `bioat`."""
        return __ABOUT__

    def list(
        self,
    ):
        """Return GROUPS and COMMANDS.

        :return: tree formatted subcommands that can be used
        :rtype: string
        """
        # ! do not change obj method to class method
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

        :return: version info
        :rtype: string
        """
        return __version__
