import bioat
import os
import sys
import random
import string
import subprocess
from bioat.cli import Cli
from settings import DATA_PATH

cli = Cli()

f_mpileup = os.path.join(DATA_PATH, 'bam/test_sorted.mpileup.gz')
f_bam_sortp = os.path.join(DATA_PATH, 'bam/test_sorted.bam')
f_bam_sortn = os.path.join(DATA_PATH, 'bam/test_sorted_n.bam')
temp_dir = f"/tmp/bioat_{''.join(random.sample(string.ascii_letters + string.digits, 16))}"


def test_mpileup_to_table():
    cli.bam.mpileup_to_table(
        mpileup=f_mpileup,
        output='/dev/null',
        threads=os.cpu_count() - 1,
        mutation_number_threshold=3,
        temp_dir=temp_dir,
        remove_temp=True,
        log_level='WARNING'
    )


def test_cli_mpileup_to_table():
    args = [
        'bioat', 'bam', 'mpileup_to_table',
        '--mpileup', f_mpileup,
        '--output', '/dev/null',
        '--threads', str(os.cpu_count() - 1),
        '--mutation_number_threshold', '0',
        '--temp_dir', temp_dir,
        '--remove_temp', 'True',
        '--log_level', 'WARNING'
    ]
    assert subprocess.check_call(args) == 0


def test_remove_clip():
    cli.bam.remove_clip(
        input=f_bam_sortn,
        output='/dev/null',
        threads=os.cpu_count() - 1,
        output_fmt='SAM',
        remove_as_paired=True,
        max_clip=0,
        log_level='WARNING'
    )


def test_cli_remove_clip():
    args = [
        'bioat', 'bam', 'remove_clip',
        '--input', f_bam_sortn,
        '--output', '/dev/null',
        '--threads', str(os.cpu_count() - 1),
        '--output_fmt', 'SAM',
        '--remove_as_paired', 'True',
        '--max_clip', '0',
        '--log_level', 'WARNING'
    ]
    assert subprocess.check_call(args) == 0
