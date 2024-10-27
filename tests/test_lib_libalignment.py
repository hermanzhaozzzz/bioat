"""_summary_

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

from Bio.Seq import Seq

from bioat.devtools import profile
from bioat.lib.libalignment import (
    get_aligned_seq,
    get_alignment_info,
    get_best_alignment,
    instantiate_pairwise_aligner,
)


def test_instantiate_pairwise_aligner():
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
        Seq("GGCACTGCGGCTGGAAAAAAAAAAAAAAAGT"),
        Seq("GGCAGNNGCTGGAAAAAAAAANNNAAAAGGT"),
        aligner=aligner,
        consider_strand=True,
    )
    print(alignment.score)
    print(alignment)


def test_instantiate_pairwise_aligner_letn_match():
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
    print(aligner.substitution_matrix)
    alignment = get_best_alignment(
        Seq("GGCACTGCGGCTGGAAAAAAAAAAAAAAAGT"),
        Seq("GGCAGNNGCTGGAAAAAAAAANNNAAAAGGT"),
        aligner=aligner,
        consider_strand=True,
    )
    print(alignment.score)
    print(alignment)


@profile(num_iterations=10000)
def test_instantiate_pairwise_aligner_10000x():
    aligner = instantiate_pairwise_aligner(
        scoring_match=1,
        penalty_mismatch=-0.8,
        penalty_gap_open=-5,
        penalty_gap_extension=-2,
        penalty_query_left_gap_score=0,
        penalty_query_right_gap_score=0,
        mode="global",
        log_level="WARNING",
    )
    alignments = aligner.align(
        Seq("GGCACTGCGGCTGGAAAAAAAAAAAAAAAGT"),
        Seq("GGCAGCGGCTGGAAAAAAAAAAAAAAAAGGT"),
    )
    alignment = get_best_alignment(
        Seq("GGCACTGCGGCTGGAAAAAAAAAAAAAAAGT"),
        Seq("GGCAGNNGCTGGAAAAAAAAANNNAAAAGGT"),
        aligner=aligner,
        consider_strand=True,
    )
    # print(alignment.aligned)


def test_get_aligned_seq():
    aligner = instantiate_pairwise_aligner(log_level="DEBUG")
    alignment = get_best_alignment(
        Seq("GGCACTGCGGCTGGAAAAAAAAAAAAAAAGT"),
        Seq("GGCAGNNGCTGGAAAAAAAAANNNAAAAGGT"),
        aligner=aligner,
        consider_strand=True,
    )
    res = get_aligned_seq(alignment, reverse=False)
    print(res)
    res = get_aligned_seq(alignment, reverse=True)
    print(res)


@profile(num_iterations=10000)
def test_get_aligned_seq10000x():
    aligner = instantiate_pairwise_aligner(log_level="WARNING")
    alignment = get_best_alignment(
        Seq("GGCACTGCGGCTGGAAAAAAAAAAAAAAAGT"),
        Seq("GGCAGNNGCTGGAAAAAAAAANNNAAAAGGT"),
        aligner=aligner,
        consider_strand=True,
    )
    res = get_aligned_seq(alignment, reverse=False)


def test_get_aligned_seq_letn_match():
    aligner = instantiate_pairwise_aligner(log_level="DEBUG")
    alignment = get_best_alignment(
        Seq("GGCACTGCGGCTGGAAAAAAAAAAAAAAAGT"),
        Seq("GGCAGNNGCTGGAAAAAAAAANNNAAAAGGT"),
        aligner=aligner,
        consider_strand=True,
    )
    res = get_aligned_seq(alignment, reverse=False, letn_match=True)
    print(res)
    print(res["reference_seq"])
    print(res["aln_info"])
    print(res["target_seq"])


@profile(num_iterations=10000)
def test_get_aligned_seq_letn_match10000x():
    aligner = instantiate_pairwise_aligner(log_level="WARNING")
    alignment = get_best_alignment(
        Seq("GGCACTGCGGCTGGAAAAAAAAAAAAAAAGT"),
        Seq("GGCAGNNGCTGGAAAAAAAAANNNAAAAGGT"),
        aligner=aligner,
        consider_strand=True,
    )
    res = get_aligned_seq(alignment, reverse=False, letn_match=True)


def test_get_alignment_info():
    aligner = instantiate_pairwise_aligner(log_level="DEBUG")
    alignment = get_best_alignment(
        Seq("GGCACTGCGGCTGGAAAAAAAAAAAAAAAGT"),
        Seq("GGCAGNNGCTGGAAAAAAAAANNNAAAAGGT"),
        aligner=aligner,
        consider_strand=True,
    )
    res = get_alignment_info(alignment, reverse=False)
    print(res)
    res = get_alignment_info(alignment, reverse=True)
    print(res)


@profile(num_iterations=10000)
def test_get_alignment_info10000x():
    aligner = instantiate_pairwise_aligner(log_level="WARNING")
    alignment = get_best_alignment(
        Seq("GGCACTGCGGCTGGAAAAAAAAAAAAAAAGT"),
        Seq("GGCAGNNGCTGGAAAAAAAAANNNAAAAGGT"),
        aligner=aligner,
        consider_strand=True,
    )
    res = get_alignment_info(alignment, reverse=False)


def test_get_alignment_info_letn_match(letn_match=False, log_level="DEBUG"):
    aligner = instantiate_pairwise_aligner(
        scoring_match=1,
        penalty_mismatch=-0.8,
        penalty_gap_open=-5,
        penalty_gap_extension=-2,
        penalty_query_left_gap_score=0,
        penalty_query_right_gap_score=0,
        mode="global",
        letn_match=letn_match,
        log_level=log_level,
    )
    alignment = get_best_alignment(
        Seq("AGGCCCCCCCT"),
        Seq("AGGNNGGGT"),
        aligner=aligner,
        consider_strand=True,
    )
    print(alignment.__dir__())
    print("reverse=False, letn_match=False")
    res = get_alignment_info(
        alignment, reverse=False, letn_match=letn_match, log_level=log_level
    )
    print(res)
    print(res["alignment"]["reference_seq"])
    print(res["alignment"]["aln_info"])

    print(res["alignment"]["target_seq"])
    print(hasattr(alignment, "is_a_reverse_complement_alignment"))


if __name__ == "__main__":
    # test_instantiate_pairwise_aligner()
    # test_instantiate_pairwise_aligner_letn_match()
    # test_instantiate_pairwise_aligner_10000x()
    # test_get_aligned_seq()
    # test_get_aligned_seq10000x()
    # test_get_aligned_seq_letn_match()
    # test_get_aligned_seq_letn_match10000x()
    # test_get_alignment_info()
    # test_get_alignment_info10000x()
    test_get_alignment_info_letn_match()
    pass
