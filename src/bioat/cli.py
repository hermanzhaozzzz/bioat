from __future__ import absolute_import
from bioat import (
    __version__,
    about,
    BamTools,
    BedTools,
    # FastxTools,
    HiCTools,
    MgiTools,
    # SystemTools,
    TableTools,
    TargetSeq
)


class Cli(object):
    """Cli interface of bioat

    - `bioat` is short for "Bioinformatic Analysis Tools" and is a commandline toolkit and a python package.
        It can be used through this cli interface in terminal or the `import` way in python codes.
    - `bioat` has many subcommand to deal with different bio-format:
        BED, BAM, FASTA, FASTQ, VCF, et al.

    """

    # - Citation:
    #     - Citation is not forcible and I would be appreciated if you wann to do this.
    #     - citation reference format:
    #         [1] Huanan ZHAO. bioat, a Bioinformatic Analysis Tool kit (2023), https://github.com/hermanzhaozzzz/bioat

    # """

    def __init__(self):
        self.bam = BamTools()
        self.bed = BedTools()
        # self.fastx = FastxTools()
        self.hic = HiCTools()
        self.mgi = MgiTools()
        # self.system = SystemTools()
        self.table = TableTools()
        self.target_seq = TargetSeq()

    @classmethod
    def about(self):
        """Print information about `bioat`."""
        return about

    def list(self):  # not class method because some attr and methods do exit when obj instantiation
        """Print commands and subcommands."""
        out = ""
        for att in dir(self):
            if not att.startswith('_'):
                out += f'{att}\n'

                for sub_att in dir(getattr(self, att)):
                    if not sub_att.startswith('_'):
                        out += f'  ├── {sub_att}\n'
        return out

    @classmethod
    def version(self):
        """Print version information."""
        return __version__
