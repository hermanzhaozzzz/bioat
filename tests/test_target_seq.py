"""
Tests for `target_seq` package.
"""

import os
import subprocess
import time

import pytest

from bioat.cli import Cli

from .settings import DATA_PATH, TESTOUT_PATH

cli = Cli()

f_bmat = os.path.join(DATA_PATH, "target_seq/test_sorted.mpileup.info.tsv")
f_refseq = os.path.join(DATA_PATH, "target_seq/HK4-AOut-1.ref.upper.fa")
f_output_fig = os.path.join(TESTOUT_PATH, "test_target_seq.pdf")
f_output_fig_count_ratio = os.path.join(
    TESTOUT_PATH, "test_target_seq_compare_count-ratio.pdf"
)
f_output_table_count_ratio = os.path.join(
    TESTOUT_PATH, "test_target_seq_compare_count-ratio.csv"
)
f_output_fig_heatmap = os.path.join(
    TESTOUT_PATH, "test_target_seq_compare_heatmap.pdf"
)
f_output_table_heatmap = os.path.join(
    TESTOUT_PATH, "test_target_seq_compare_heatmap.csv"
)
output_files = [
    f_output_fig,
    f_output_fig_count_ratio,
    f_output_table_count_ratio,
    f_output_fig_heatmap,
    f_output_table_heatmap,
]


def test_region_heatmap_1():
    # test no reference_seq
    cli.target_seq.region_heatmap(
        input_table=f_bmat,
        output_fig=f_output_fig,
        target_seq="GGCACTGCGGCTGGAGGTGG^NGG^",  # HEK4
        log_level="WARNING",
    )


def test_region_heatmap_2():
    # test reference_seq as fasta file
    cli.target_seq.region_heatmap(
        input_table=f_bmat,
        output_fig=f_output_fig,
        target_seq="GGCACTGCGGCTGGAGGTGG^NGG^",  # HEK4
        reference_seq=f_refseq,  # target locus
        log_level="WARNING",
    )


def test_region_heatmap_3():
    # test reference_seq as str
    cli.target_seq.region_heatmap(
        input_table=f_bmat,
        output_fig=f_output_fig,
        target_seq="GGCACTGCGGCTGGAGGTGG^NGG^",  # HEK4
        reference_seq="GTCACGACCCTTCCACAAGAGGAGATACTGACAGTGGGAACACTGTCCACCTTTCTGCCTCTGGAGAGGGAGGAGGGGCCTCTTCTTGCCAGCCA"
        "TAAGGGCTGGCATATGCAGGGGGGTAAAAATAGAACCTCAGCGGCGCTGTCTCCCCTTCCAGTGAGGGCGGGCACAGGCTGGGGATGGAGGAGCT"
        "GGCTGTGAGGGCAGTCTGCAGCAAACGTGATCTTCCCCGCTCTGAACCTCTGTGCCTGATACAACCCACCAGGGTC",  # target locus
        log_level="WARNING",
    )


def test_cli_region_heatmap():
    args = [
        "bioat",
        "target_seq",
        "region_heatmap",
        "--input_table",
        f_bmat,
        "--output_fig",
        f_output_fig,
        "--target_seq",
        "GGCACTGCGGCTGGAGGTGG^NGG^",
        "--reference_seq",
        f_refseq,
        "--log_level",
        "WARNING",
    ]
    assert subprocess.run(args).returncode == 0


def test_region_heatmap_compare_1():
    cli.target_seq.region_heatmap_compare(
        input_tables=f"{f_bmat},{f_bmat},{f_bmat}",
        labels=("test1", "test2", "test3"),
        to_base=("A", "G", "C", "T", "Ins", "Del"),
        heatmap_mut_direction=("AT", "GC"),
        count_ratio="all",
        output_fig_count_ratio=f_output_fig_count_ratio,
        output_fig_heatmap=f_output_fig_heatmap,
        output_table_count_ratio=f_output_table_count_ratio,
        output_table_heatmap=f_output_table_heatmap,
        region_extend_length=5,
        target_seq="GGCACTGCGGCTGGAGGTGG^NGG^",  # HEK4
        # reference_seq=f_refseq,  # target locus
        log_level="WARNING",
    )


def test_cli_region_heatmap_compare_1():
    args = [
        "bioat",
        "target_seq",
        "region_heatmap_compare",
        "--input_tables",
        f"{f_bmat},{f_bmat},{f_bmat}",
        "--labels",
        ",".join(("test1", "test2", "test3")),
        "--to_base",
        ",".join(("A", "G", "C", "T", "Ins", "Del")),
        "--heatmap_mut_direction",
        ",".join(("AT", "GC")),
        "--count_ratio",
        "all",
        "--output_fig_count_ratio",
        f_output_fig_count_ratio,
        "--output_fig_heatmap",
        f_output_fig_heatmap,
        "--output_table_count_ratio",
        f_output_table_count_ratio,
        "--output_table_heatmap",
        f_output_table_heatmap,
        "--region_extend_length",
        "5",
        "--target_seq",
        "GGCACTGCGGCTGGAGGTGG^NGG^",
        "--log_level",
        "WARNING",
    ]
    assert subprocess.run(args, check=True)


def test_region_heatmap_compare_2():
    cli.target_seq.region_heatmap_compare(
        input_tables=f"{f_bmat},{f_bmat},{f_bmat}",
        labels=("test1", "test2", "test3"),
        to_base=("A", "G", "C", "T", "Ins", "Del"),
        heatmap_mut_direction=("AT", "GC"),
        count_ratio="all",
        output_fig_count_ratio=f_output_fig_count_ratio,
        output_fig_heatmap=f_output_fig_heatmap,
        output_table_count_ratio=f_output_table_count_ratio,
        output_table_heatmap=f_output_table_heatmap,
        region_extend_length=5,
        target_seq="GGCACTGCGGCTGGAGGTGG^NGG^",  # HEK4
        reference_seq=f_refseq,  # target locus
        log_level="WARNING",
    )


def test_cli_region_heatmap_compare_2():
    args = [
        "bioat",
        "target_seq",
        "region_heatmap_compare",
        "--input_tables",
        f"{f_bmat},{f_bmat},{f_bmat}",
        "--labels",
        ",".join(("test1", "test2", "test3")),
        "--to_base",
        ",".join(("A", "G", "C", "T", "Ins", "Del")),
        "--heatmap_mut_direction",
        ",".join(("AT", "GC")),
        "--count_ratio",
        "all",
        "--output_fig_count_ratio",
        f_output_fig_count_ratio,
        "--output_fig_heatmap",
        f_output_fig_heatmap,
        "--output_table_count_ratio",
        f_output_table_count_ratio,
        "--output_table_heatmap",
        f_output_table_heatmap,
        "--region_extend_length",
        "5",
        "--target_seq",
        "GGCACTGCGGCTGGAGGTGG^NGG^",
        "--reference_seq",
        f_refseq,
        "--log_level",
        "WARNING",
    ]
    assert subprocess.run(args, check=True)


def test_region_heatmap_compare_3():
    cli.target_seq.region_heatmap_compare(
        input_tables=f"{f_bmat},{f_bmat},{f_bmat}",
        labels=("test1", "test2", "test3"),
        to_base=("A", "G", "C", "T", "Ins", "Del"),
        heatmap_mut_direction=("AT", "GC"),
        count_ratio="all",
        output_fig_count_ratio=f_output_fig_count_ratio,
        output_fig_heatmap=f_output_fig_heatmap,
        output_table_count_ratio=f_output_table_count_ratio,
        output_table_heatmap=f_output_table_heatmap,
        region_extend_length=5,
        target_seq="GGCACTGCGGCTGGAGGTGG^NGG^",  # HEK4
        reference_seq="GTCACGACCCTTCCACAAGAGGAGATACTGACAGTGGGAACACTGTCCACCTTTCTGCCTCTGGAGAGGGAGGAGGGGCCTCTTCTTGCCAGCCA"
        "TAAGGGCTGGCATATGCAGGGGGGTAAAAATAGAACCTCAGCGGCGCTGTCTCCCCTTCCAGTGAGGGCGGGCACAGGCTGGGGATGGAGGAGCT"
        "GGCTGTGAGGGCAGTCTGCAGCAAACGTGATCTTCCCCGCTCTGAACCTCTGTGCCTGATACAACCCACCAGGGTC",  # target locus
        log_level="WARNING",
    )


def test_cli_region_heatmap_compare_3():
    str_refseq = "GTCACGACCCTTCCACAAGAGGAGATACTGACAGTGGGAACACTGTCCACCTTTCTGCCTCTGGAGAGGGAGGAGGGGCCTCTTCTTGCCAGCCA"
    "TAAGGGCTGGCATATGCAGGGGGGTAAAAATAGAACCTCAGCGGCGCTGTCTCCCCTTCCAGTGAGGGCGGGCACAGGCTGGGGATGGAGGAGCT"
    "GGCTGTGAGGGCAGTCTGCAGCAAACGTGATCTTCCCCGCTCTGAACCTCTGTGCCTGATACAACCCACCAGGGTC"

    args = [
        "bioat",
        "target_seq",
        "region_heatmap_compare",
        "--input_tables",
        f"{f_bmat},{f_bmat},{f_bmat}",
        "--labels",
        ",".join(("test1", "test2", "test3")),
        "--to_base",
        ",".join(("A", "G", "C", "T", "Ins", "Del")),
        "--heatmap_mut_direction",
        ",".join(("AT", "GC")),
        "--count_ratio",
        "all",
        "--output_fig_count_ratio",
        f_output_fig_count_ratio,
        "--output_fig_heatmap",
        f_output_fig_heatmap,
        "--output_table_count_ratio",
        f_output_table_count_ratio,
        "--output_table_heatmap",
        f_output_table_heatmap,
        "--region_extend_length",
        "5",
        "--target_seq",
        "GGCACTGCGGCTGGAGGTGG^NGG^",
        "--reference_seq",
        str_refseq,
        "--log_level",
        "WARNING",
    ]
    assert subprocess.run(args, check=True)


def test_finish_test_and_clean_files():
    for i in output_files:
        try:
            os.remove(i)
        except:
            pass
