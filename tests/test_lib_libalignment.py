"""_summary_.

author: Herman Huanan Zhao
email: hermanzhaozzzz@gmail.com
homepage: https://github.com/hermanzhaozzzz

_description_

example 1:
    bioat list
        <in shell>:
            $ bioat list
        <in python consolo>:
            >>> from bioat.cli import Cli
            >>> bioat = Cli()
            >>> bioat.list()
            >>> print(bioat.list())

example 2:
    _example_
"""

import pytest
from Bio.Seq import Seq

from bioat.lib.libalignment import (
    get_aligned_seq,
    get_alignment_info,
    get_best_alignment,
    instantiate_pairwise_aligner,
)


@pytest.fixture
def test_seq_pair():
    return (
        Seq("GGCACTGCGGCTGGAAAAAAAAAAAAAAAGT"),
        Seq("GGCAGNNGCTGGAAAAAAAAANNNAAAAGGT"),
    )


def test_instantiate_pairwise_aligner(test_seq_pair):
    aligner = instantiate_pairwise_aligner(
        scoring_match=1,
        penalty_mismatch=-0.8,
        penalty_gap_open=-5,
        penalty_gap_extension=-2,
        penalty_query_left_gap_score=0,
        penalty_query_right_gap_score=0,
        mode="global",
        log_level="DEBUG",
    )
    alignment = get_best_alignment(
        *test_seq_pair, aligner=aligner, consider_strand=True
    )
    assert alignment.score is not None
    assert alignment.aligned is not None


def test_instantiate_pairwise_aligner_letn_match(test_seq_pair):
    aligner = instantiate_pairwise_aligner(
        scoring_match=1,
        penalty_mismatch=-0.8,
        penalty_gap_open=-5,
        penalty_gap_extension=-2,
        penalty_query_left_gap_score=0,
        penalty_query_right_gap_score=0,
        mode="global",
        letn_match=True,
        log_level="DEBUG",
    )
    alignment = get_best_alignment(
        *test_seq_pair, aligner=aligner, consider_strand=True
    )
    assert alignment.score is not None


def test_get_aligned_seq(test_seq_pair):
    aligner = instantiate_pairwise_aligner()
    alignment = get_best_alignment(
        *test_seq_pair, aligner=aligner, consider_strand=True
    )
    result_normal = get_aligned_seq(alignment, reverse=False)
    result_reverse = get_aligned_seq(alignment, reverse=True)
    assert "reference_seq" in result_normal
    assert "target_seq" in result_reverse


def test_get_aligned_seq_letn_match(test_seq_pair):
    aligner = instantiate_pairwise_aligner()
    alignment = get_best_alignment(
        *test_seq_pair, aligner=aligner, consider_strand=True
    )
    result = get_aligned_seq(alignment, reverse=False, letn_match=True)
    assert set(result.keys()) >= {"reference_seq", "target_seq", "aln_info"}


def test_get_alignment_info(test_seq_pair):
    aligner = instantiate_pairwise_aligner()
    alignment = get_best_alignment(
        *test_seq_pair, aligner=aligner, consider_strand=True
    )
    result = get_alignment_info(alignment, reverse=False)
    assert "alignment" in result
    assert set(result["alignment"].keys()) >= {
        "reference_seq",
        "target_seq",
        "aln_info",
    }


def test_get_alignment_info_letn_match():
    aligner = instantiate_pairwise_aligner(
        scoring_match=1,
        penalty_mismatch=-0.8,
        penalty_gap_open=-5,
        penalty_gap_extension=-2,
        penalty_query_left_gap_score=0,
        penalty_query_right_gap_score=0,
        mode="global",
        letn_match=True,
        log_level="DEBUG",
    )
    alignment = get_best_alignment(
        Seq("AGGCCCCCCCT"),
        Seq("AGGNNGGGT"),
        aligner=aligner,
        consider_strand=True,
    )
    result = get_alignment_info(
        alignment,
        reverse=False,
        letn_match=True,
        log_level="DEBUG",
    )
    assert set(result["alignment"].keys()) >= {
        "reference_seq",
        "target_seq",
        "aln_info",
    }
