from __future__ import absolute_import
import fire
from .about import about
from .fastx import Fastx
from .table import Table
from .hic import HiC
from .mgi import Mgi
from .system import System


class Cli(object):
    """Cli interface of python package <bioat>

    - <bioat> is a commandline tool and a python package, which is short for
        Bioinformatic Analysis Tools.
    - <bioat> has many subcommand to deal with different bio-format:
        BED, BAM, FASTA, FASTQ, VCF, et al.
    """

    def __init__(self):
        self.fastx = Fastx()
        self.hic = HiC()
        self.mgi = Mgi()
        self.system = System()
        self.table = Table()

    def about(self):
        """Print information about <bioat>
        """
        return about

    def list(self):
        return "all cmds"


def main() -> int:
    calculator = Cli()
    fire.Fire(calculator, name='bioat')


if __name__ == '__main__':
    main()
