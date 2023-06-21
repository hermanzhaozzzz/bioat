import bioat
import os
import sys
import random
import string
import subprocess
from bioat.cli import Cli
from settings import DATA_PATH, TESTOUT_PATH

cli = Cli()

f_bmat = os.path.join(DATA_PATH, 'target_seq/test_sorted.mpileup.info.tsv')
f_refseq = os.path.join(DATA_PATH, 'target_seq/HK4-AOut-1.ref.upper.fa')
f_output_fig = os.path.join(TESTOUT_PATH, 'test_target_seq.pdf')


def test_region_heatmap_1():
    # test no reference_seq
    cli.target_seq.region_heatmap(
        input_table=f_bmat,
        output_fig=f_output_fig,
        target_seq="GGCACTGCGGCTGGAGGTGG^NGG^",  # HEK4
        log_level='WARNING'
    )


def test_region_heatmap_2():
    # test reference_seq as fasta file
    cli.target_seq.region_heatmap(
        input_table=f_bmat,
        output_fig=f_output_fig,
        target_seq="GGCACTGCGGCTGGAGGTGG^NGG^",  # HEK4
        reference_seq=f_refseq,  # target locus
        log_level='WARNING'
    )


def test_region_heatmap_3():
    # test reference_seq as str
    cli.target_seq.region_heatmap(
        input_table=f_bmat,
        output_fig=f_output_fig,
        target_seq="GGCACTGCGGCTGGAGGTGG^NGG^",  # HEK4
        reference_seq="GTCACGACCCTTCCACAAGAGGAGATACTGACAGTGGGAACACTGTCCACCTTTCTGCCTCTGGAGAGGGAGGAGGGGCCTCTTCTTGCCAGCCA"
                      "TAAGGGCTGGCATATGCAGGGGGGTAAAAATAGAACCTCAGCGGCGCTGTCTCCCCTTCCAGTGAGGGCGGGCACAGGCTGGGGATGGAGGAGCT"
                      "GGCTGTGAGGGCAGTCTGCAGCAAACGTGATCTTCCCCGCTCTGAACCTCTGTGCCTGATACAACCCACCAGGGTC",
        # target locus
        log_level='WARNING'
    )


def test_cli_region_heatmap():
    args = [
        'bioat', 'target_seq', 'region_heatmap',
        '--input_table', f_bmat,
        '--output_fig', f_output_fig,
        '--target_seq', "GGCACTGCGGCTGGAGGTGG^NGG^",
        '--reference_seq', f_refseq,
        '--log_level', 'WARNING'
    ]
    assert subprocess.check_call(args) == 0

# f"""{PYTHON} ./program/multiplot-targetseq-bmat-heatmap-V3.py \
#     -i \
#     target_seq/test_sorted.mpileup.info.tsv,\
#     target_seq/test_sorted.mpileup.info.tsv,\
#     target_seq/test_sorted.mpileup.info.tsv \
#     -l \
#     test1,\
#     test2,\
#     test3 \
#     --to_base A,G,C,T,Ins,Del \
#     --count_ratio all \
#     -o 'test_countratio.pdf' \
#     --mut_direction "CT,GA" \
#     --plot_heatmap 'test_heatmap.pdf' \
#     --region_extend_length 5 \
#     --sgRNA GGCACTGCGGCTGGAGGTGG """
