import dataclasses
import gzip
import os
import sys
from bioat import get_logger
from bioat.lib.libdataclasses import Bed


class BedTools:
    """Bed toolbox."""

    def __init__(self):
        pass

    def random(self, length=200, number=100, fasta=None, seed=0, rm_N=False) -> Bed:
        """Get random region

        :param length: random region length.
        :param number: random region number.
        :param fasta: get random region from which genome, must defined when rm_n=True
        :param seed: random seed
        :param rm_N: remove random region contains N base (It's often the telomere region)
        :return: bed file
        """
        pass
