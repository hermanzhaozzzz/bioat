import bioat
import os
import sys
import random
import string
from io import StringIO
from unittest.mock import patch
from bioat.cli import Cli

f_mpileup = f'{bioat.__path__[0]}/../../data/bam/test_sorted.mpileup'
f_bam = f'{bioat.__path__[0]}/../../data/bam/test_sorted.bam'


def test_bam_mpileup_to_table():
    cli = Cli()
    cli.bam.mpileup_to_table(
        mpileup=f_mpileup,
        output=sys.stdout,
        threads=os.cpu_count() - 1,
        mutation_number_threshold=0,
        temp_dir=f"/tmp/bioat_{''.join(random.sample(string.ascii_letters + string.digits, 16))}",
        remove_temp=True,
        log_level='INFO'
    )
    cli.bam.mpileup_to_table(
        mpileup=f_mpileup,
        output='/dev/null',
        threads=os.cpu_count() - 1,
        mutation_number_threshold=3,
        temp_dir=f"/tmp/bioat_{''.join(random.sample(string.ascii_letters + string.digits, 16))}",
        remove_temp=True,
        log_level='INFO'
    )

def test_bam_mpileup_to_table():
    cli = Cli()
    cli.bam.remove_clip(
        input=f_bam
    )