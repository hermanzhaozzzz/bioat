from __future__ import absolute_import
import fire
from bioat import (
    __version__,
    about,
    BamTools,
    BedTools,
    FastxTools,
    HiCTools,
    MgiTools,
    SystemTools,
    TableTools,
    TargetedDeepSequencing
)


class Cli(object):
    """Cli interface of python package <bioat>

    - <bioat> is a commandline tool and a python package, which is short for
        Bioinformatic Analysis Tools.
    - <bioat> has many subcommand to deal with different bio-format:
        BED, BAM, FASTA, FASTQ, VCF, et al.
    """

    def __init__(self):
        self.bam = BamTools()
        self.bed = BedTools()
        # self.fastx = FastxTools()
        self.hic = HiCTools()
        self.mgi = MgiTools()
        # self.system = SystemTools()
        self.table = TableTools()
        self.targeted_deep_sequencing = TargetedDeepSequencing()

    @classmethod
    def about(self):
        """Print information about <bioat>.
        """
        return about

    def list(self):
        """Print commands and subcommands.
        """
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
        """Print version information.
        """
        return __version__


def main() -> int:
    calculator = Cli()
    fire.Fire(calculator, name='bioat')


if __name__ == '__main__':
    main()
