"""Library for alignment, depends on Biopython.

author: Herman Huanan Zhao
email: hermanzhaozzzz@gmail.com
homepage: https://github.com/hermanzhaozzzz

Library for alignment, depends on Biopython

example 1:
    bioat.lib.libalignment
        <in python consolo>:
            >>> from bioat.lib.libalignment import instantiate_pairwise_aligner
            >>> aligner = instantiate_pairwise_aligner()
            >>> alignment = get_best_alignment(
            >>>     Seq("GGCACTGCGGCTGGAAAAAAAAAAAAAAAGT"),
            >>>     Seq("GGCAGCGGCTGGAAAAAAAAAAAAAAAAGGT"),
            >>>     aligner=aligner,
            >>>     consider_strand=True
            >>> )
            >>> res = get_aligned_seq(alignment, reverse=False)
            >>> print(res)
            >>> res = get_aligned_seq(alignment, reverse=True)
            >>> print(res)
            >>> res = get_alignment_info(alignment, reverse=False)
            >>> print(res)
            >>> res = get_alignment_info(alignment, reverse=True)
            >>> print(res)
"""

import sys

from Bio.Align import Alignment, PairwiseAligner, substitution_matrices
from Bio.Seq import Seq

from bioat.logger import LoggerManager

lm = LoggerManager(mod_name="bioat.lib.libalignment")


def instantiate_pairwise_aligner(
    scoring_match: float = 1,
    penalty_mismatch: float = 0.8,
    penalty_gap_open: float = -5,
    penalty_gap_extension: float = -2,
    penalty_query_left_gap_score: float = 0,
    penalty_query_right_gap_score: float = 0,
    mode="global",
    letn_match=False,
    score_matrix_dict: None | dict = None,
    log_level="DEBUG",
) -> PairwiseAligner:
    """Returns a PairwiseAligner object.

    Args:
        scoring_match (float, optional): Scoring for matches. Defaults to 1.
        penalty_mismatch (float, optional): Penalty for mismatches. Defaults to 0.8.
        penalty_gap_open (float, optional): Penalty for opening a gap. Defaults to -5.
        penalty_gap_extension (float, optional): Penalty for extending a gap. Defaults to -2.
        penalty_query_left_gap_score (float, optional): Penalty for left query gap. Defaults to 0.
        penalty_query_right_gap_score (float, optional): Penalty for right query gap. Defaults to 0.
        mode (str, optional): Alignment mode, either "global" or "local". Defaults to "global".
        letn_match (bool, optional): Whether to treat base 'N' as a match. Defaults to False.
        log_level (str, optional): Logging level. Defaults to "DEBUG".
        score_matrix_dict (dict, optional): Custom scoring matrix for bases or amino acids. Defaults to None.

    Returns:
        PairwiseAligner: An object instantiated from PairwiseAligner.

    Details for `score_matrix_dict`:
        You can use any base and any score, such as AGCTN, or 20AA + 2 abnormal AAs,
        or any characters you want. For example:

            score_matrix_dict = {

                ("A", "A"): scoring_match,
                ("A", "G"): penalty_mismatch,
                ("A", "C"): penalty_mismatch,
                ("A", "T"): penalty_mismatch,
                ("A", "N"): scoring_match,
                ("G", "A"): penalty_mismatch,
                ("G", "G"): scoring_match,
                ("G", "C"): penalty_mismatch,
                ("G", "T"): penalty_mismatch,
                ("G", "N"): scoring_match,
                ("C", "A"): penalty_mismatch,
                ("C", "G"): penalty_mismatch,
                ("C", "C"): scoring_match,
                ("C", "T"): penalty_mismatch,
                ("C", "N"): scoring_match,
                ("T", "A"): penalty_mismatch,
                ("T", "G"): penalty_mismatch,
                ("T", "C"): penalty_mismatch,
                ("T", "T"): scoring_match,
                ("T", "N"): scoring_match,
                ("N", "A"): scoring_match,
                ("N", "G"): scoring_match,
                ("N", "C"): scoring_match,
                ("N", "T"): scoring_match,
                ("N", "N"): scoring_match

            }
    """
    lm.set_names(func_name="instantiate_pairwise_aligner")
    lm.set_level(log_level)

    aligner = PairwiseAligner()

    if not letn_match and score_matrix_dict is None:
        # !如果 letn_match为False，
        # !并且未指定 score_matrix_dict
        # !则考虑本流程
        # !否则进入else中的自定义流程
        aligner.mode = mode
        aligner.match = scoring_match
        aligner.mismatch = penalty_mismatch
        aligner.open_gap_score = penalty_gap_open
        aligner.extend_gap_score = penalty_gap_extension
        aligner.query_left_gap_score = penalty_query_left_gap_score
        aligner.query_right_gap_score = penalty_query_right_gap_score
        aligner.wildcard = None
    else:
        if not score_matrix_dict:
            score_matrix_dict = {
                ("A", "A"): scoring_match,
                ("A", "G"): penalty_mismatch,
                ("A", "C"): penalty_mismatch,
                ("A", "T"): penalty_mismatch,
                ("A", "N"): scoring_match,
                ("G", "A"): penalty_mismatch,
                ("G", "G"): scoring_match,
                ("G", "C"): penalty_mismatch,
                ("G", "T"): penalty_mismatch,
                ("G", "N"): scoring_match,
                ("C", "A"): penalty_mismatch,
                ("C", "G"): penalty_mismatch,
                ("C", "C"): scoring_match,
                ("C", "T"): penalty_mismatch,
                ("C", "N"): scoring_match,
                ("T", "A"): penalty_mismatch,
                ("T", "G"): penalty_mismatch,
                ("T", "C"): penalty_mismatch,
                ("T", "T"): scoring_match,
                ("T", "N"): scoring_match,
                ("N", "A"): scoring_match,
                ("N", "G"): scoring_match,
                ("N", "C"): scoring_match,
                ("N", "T"): scoring_match,
                ("N", "N"): scoring_match,
            }

        substitution_matrix = substitution_matrices.Array(
            data=score_matrix_dict,
        )
        aligner.mode = mode
        aligner.substitution_matrix = substitution_matrix
        # !自定义矩阵时，不可定义match和mismatch的打分，不然会出问题，在矩阵中给分即可
        # !aligner.match = scoring_match
        # !aligner.mismatch = penalty_mismatch
        aligner.open_gap_score = penalty_gap_open
        aligner.extend_gap_score = penalty_gap_extension
        aligner.query_left_gap_score = penalty_query_left_gap_score
        aligner.query_right_gap_score = penalty_query_right_gap_score
        aligner.wildcard = None
        # !自定义的打分矩阵中如果有"N"，则wildcard最好注释掉
        # !wildcard是一个特殊字符，用于在序列比对过程中表示一个可以与任何字符匹配的通配符
        # !设置wildcard的目的通常是为了在处理包含未知或可变字符的序列时提高比对的灵活性
        # !使用wildcard时需要谨慎，因为它可能会影响比对的准确性和生物学意义
        # !可以将wildcard设置为'N'，这样在比对过程中，字符'N'将被视为与任何字符都匹配

    lm.logger.debug(aligner)
    return aligner


def get_best_alignment(
    seq_a: Seq,
    seq_b: Seq,
    aligner: PairwiseAligner,
    consider_strand=True,
):
    """Get best alignment by heightest alignment score.

    Get best alignment by heightest alignment score,
    default to consider strand, that means, it will attempt to use
    seq_b.reverse_complement() to align

    :param seq_a: seq_a
    :type seq_a: Seq
    :param seq_b: seq_b
    :type seq_b: Seq
    :param aligner: aligner returns from instantiate_pairwise_aligner
    :type aligner: PairwiseAligner
    :param consider_strand: consider to use seq_b.reverse_complement() or not, defaults to True
    :type consider_strand: bool, optional
    :return: Alignment with a new attribute: alignment.is_a_reverse_complement_alignment (True / False)
    :rtype: Alignment
    """
    alignments = aligner.align(seq_a, seq_b)
    # ! // sort alignments by score
    # ! old
    # alignment = sorted(alignments, key=lambda x: x.score, reverse=True)[0]
    # ! new
    alns = []
    scores = []
    for aln in alignments:
        alns.append(aln)
        scores.append(aln.score)
        if len(alns) >= 100:
            # ! break the loop if more than 100 alignments
            # ! this is a bug fix, in case of too many alignments
            # ! with too many memory occupied
            break
    # find max score index
    max_index = scores.index(max(scores))
    # get max score alignment
    alignment = alns[max_index]
    # ! sort alignments by score //

    is_a_reverse_complement_alignment = None
    final_alignment = None

    if consider_strand:
        alignments_rc = aligner.align(seq_a, seq_b.reverse_complement())
        alignment_rc = sorted(
            alignments_rc,
            key=lambda x: x.score,
            reverse=True,
        )[0]

        if alignment.score >= alignment_rc.score:
            final_alignment = alignment
            is_a_reverse_complement_alignment = False
        else:
            final_alignment = alignment_rc
            is_a_reverse_complement_alignment = True
    else:
        final_alignment = alignment
        is_a_reverse_complement_alignment = False

    # obj add attr
    final_alignment.is_a_reverse_complement_alignment = (
        is_a_reverse_complement_alignment
    )
    return final_alignment


def get_aligned_seq(
    alignment: Alignment,
    reverse: bool = False,
    letn_match=False,
) -> dict:
    """Retrieve aligned sequence information from a Bio.Align.Alignment object.

    Args:
        alignment (Bio.Align.Alignment):
            A `Bio.Align.Alignment` object containing the alignment information.
        reverse (bool, optional):
            Whether to return reversed sequence information. Defaults to False.

    Returns:
        dict: A dictionary containing the alignment details with the following keys:
            - 'reference_seq' (str): The aligned reference sequence.
            - 'aln_info' (str): The alignment information string, indicating matches ('|'),
              mismatches ('.'), and gaps ('-').
            - 'target_seq' (str): The aligned target sequence.

    Example:
        >>> result = get_aligned_seq(alignment)
        >>> print(result)
        {
            'reference_seq': 'GGCACTGCGGCTGGAAAAAAAAAAAAAAA--GT',
            'aln_info':      '--..|.|||||||||||||||||||||||--||',
            'target_seq':    '--GGCAGCGGCTGGAAAAAAAAAAAAAAAAGGT'
        }

    """
    ls = []

    for x in str(alignment).strip().split("\n"):
        if x:
            ls.append(x[20:].split(" ")[0])
    if not reverse:
        res = {
            "reference_seq": "".join(ls[::3]),
            "aln_info": "".join(ls[1::3]),
            "target_seq": "".join(ls[2::3]),
        }
    else:
        res = {
            "reference_seq": "".join(ls[::3])[::-1],
            "aln_info": "".join(ls[1::3])[::-1],
            "target_seq": "".join(ls[2::3])[::-1],
        }
    if letn_match:
        length = len(res["reference_seq"])
        # change base N to match
        i_n_a = [
            i
            for i in range(length)
            if (res["reference_seq"][i] == "N") and (res["target_seq"][i] != "-")
        ]
        i_n_b = [
            i
            for i in range(length)
            if (res["target_seq"][i] == "N") and (res["reference_seq"][i] != "-")
        ]
        i_n = list(set(i_n_a) | set(i_n_b))
        aln_list = list(res["aln_info"])

        for i in sorted(i_n):
            aln_list[i] = "★"
        res["aln_info"] = "".join(aln_list)

    return res


def get_alignment_info(
    alignment: Alignment,
    reverse: bool = False,
    letn_match=False,
    log_level="WARNING",
):
    """Analyze a Bio.Align.Alignment object to extract alignment metrics.

    This function calculates the match count, mismatch count, gap count, alignment score,
    and start and end indices for a given alignment.

    Args:
        alignment (Bio.Align.Alignment): A `Bio.Align.Alignment` object containing the alignment information.
        reverse (bool, optional): Whether to return reversed sequence information. Defaults to False.

    Returns:
        list: A list containing the following metrics:

            - match_count (int): The number of matches in the alignment.
            - mismatch_count (int): The number of mismatches in the alignment.
            - gap_count (int): The number of gaps within the sgRNA alignment region. Note: This counts only gaps within the alignment region, not total gaps.
            - alignment_score (float): The alignment score.
            - start_index (int): The start index of the sgRNA alignment region.
            - end_index (int): The end index of the sgRNA alignment region.

    Notes:
        The gap count specifically counts gaps within the sgRNA alignment region, not the total gaps in the entire sequence.
        For example:

            >>> Reference: AGTGGTAAGAAGAAGACGAGACATAATGAG
            >>>             ------||||||||||||||.|----||||
            >>> Target: ------AAGAAGAAGACGAGCC - ---TGAG

        In this example: Gap count = 4 (within sgRNA alignment region), not 10 (total gaps).


    Example:
        The returned list will be as follows: return_list = [match_count, mismatch_count, gap_count, alignment_score, start_index, end_index]

    """
    lm.set_names(func_name="get_alignment_info")
    lm.set_level(log_level)

    if reverse:
        aln_res = get_aligned_seq(
            alignment,
            reverse=True,
            letn_match=letn_match,
        )
    else:
        aln_res = get_aligned_seq(alignment, letn_match=letn_match)

    reference_seq = aln_res["reference_seq"]
    aln_info = aln_res["aln_info"]
    target_seq = aln_res["target_seq"]
    aln_score = alignment.score
    target_seq_length = len(target_seq)
    target_seq_length_rm_gap = len(target_seq.replace("-", ""))
    aligned_coordinates = alignment.aligned[::-1]
    if hasattr(alignment, "is_a_reverse_complement_alignment"):
        is_a_reverse_complement_alignment = alignment.is_a_reverse_complement_alignment
    else:
        is_a_reverse_complement_alignment = None

    lm.logger.debug(f"reference_seq = {reference_seq}")
    lm.logger.debug(f"aln_info      = {aln_info}")
    lm.logger.debug(f"target_seq    = {target_seq}")
    lm.logger.debug(f"aln_score     = {aln_score}")
    lm.logger.debug(f"target_seq_length (remove gap) = {target_seq_length_rm_gap}")
    lm.logger.debug(
        f"is_a_reverse_complement_alignment = {is_a_reverse_complement_alignment}",
    )
    lm.logger.debug(f"aligned_coordinates =\n{aligned_coordinates}")

    # define params
    ref_aln_start = target_seq_length - len(target_seq.lstrip("-"))
    ref_aln_end = ref_aln_start + len(target_seq.strip("-"))
    match_count = 0
    mismatch_count = 0
    gap_count = 0
    parsed_target_bases = 0
    encounter_target_seq = False

    for idx, base_info in enumerate(aln_info):
        """
        ---TTCTGCCTCTGGAGAGGG-GGAGGGGCCT---
        -------|.|.|||..|-.||-||.||--------
        -------GGCACTGCGG-TGGAGGTGG--------
        """
        if not encounter_target_seq:
            if target_seq[idx] == "-":
                # if start with -
                # not encounter with target_seq yet!
                # print(idx, base_info)
                continue
            # find first aligned base
            # encounter with target_seq!
            encounter_target_seq = True

        if encounter_target_seq and (parsed_target_bases < target_seq_length_rm_gap):
            # from the same loop of [else: catch_start = True]
            # parse start base and later base
            # print(parsed_target_bases, target_seq_length_rm_gap, base_info)
            if base_info in {"|", "★"}:
                # match
                match_count += 1
                parsed_target_bases += 1
            elif base_info == ".":  # . instead of X after Biopython 1.8.0
                # mismatch
                mismatch_count += 1
                parsed_target_bases += 1
            elif base_info == "-":
                # gap on reference_seq or target_seq
                gap_count += 1

                if reference_seq[idx] == "-":
                    # gap on reference_seq
                    # value base on target_seq
                    parsed_target_bases += 1
                    pass
                else:
                    # gap on target_seq
                    # no value base on target_seq
                    parsed_target_bases += 0
                    pass
            else:
                lm.logger.critical(
                    f"find not surpported alignment symbol: {base_info} ({idx} on the reference), check"
                    f"whether Biopython >= 1.8.0 or not",
                )
                sys.exit(1)
        else:
            # parsed_target_bases = target_seq_length_rm_gap
            # full target_seq bases were processed
            break

    return {
        "match_count": match_count,
        "mismatch_count": mismatch_count,
        "gap_count": gap_count,
        "aln_score": aln_score,
        "ref_aln_start": ref_aln_start,
        "ref_aln_end": ref_aln_end,
        "is_a_reverse_complement_alignment": is_a_reverse_complement_alignment,
        "alignment": aln_res,
    }
