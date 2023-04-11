from __future__ import absolute_import
import logging
import sys
from Bio.Seq import Seq
from Bio.Align import PairwiseAligner
from bioat.lib.libalignment import get_alignment_info

"""TARGET_REGIONS_LIB.

bioat targeted_deep_sequencing region_heatmap test_sorted.mpileup.info.tsv test.pdf --get_built_in_target_region

INFO:root:You can use <key> in built-in <target_regions> to represent your target_region:
        <key>   <target_region>
        EMX1    GAGTCCGAGCAGAAGAAGAA^GGG^
        HEK3    GGCCCAGACTGAGCACGTGA^TGG^
        HEK4    GGCACTGCGGCTGGAGGTGG^GGG^
        ...
"""
TARGET_SEQ_LIB = {
    "EMX1": "GAGTCCGAGCAGAAGAAGAA^GGG^",
    "HEK3": "GGCCCAGACTGAGCACGTGA^TGG^",
    "HEK4": "GGCACTGCGGCTGGAGGTGG^GGG^",
    "RNF2": "GTCATCTTAGTCATTACCTG^AGG^",
    "VEGFA": "GACCCCCTCCACCCCGCCTC^CGG^",
    "CDKN2A": "^TTTA^GCCCCAATAATCCCCACATGTCA",
    "DYRK1A": "^TTTA^GAAGCACATCAAGGACATTCTAA",
    "ND4": "TGCTAGTAACCACGTTCTCCTGATCAAATATCACTCTCCTACTTACAGGA",
    "ND5.1": "TAGCATTAGCAGGAATACCTTTCCTCACAGGTTTCTACTCCAAAGA",
    "ND6": "TGACCCCCATGCCTCAGGATACTCCTCAATAGCCATCG",
}


def find_seq_PAM_index(query_seq):
    """Find seq PAM index.

    :param query_seq: e.g. query_seq = "AAAGAGAG"
    :return: index list = [1,3,5], which means search NAG, NGG at the same time
    """
    pam_index_list = []

    for index, base in enumerate(query_seq[:-2]):
        if query_seq[index + 1:index + 3] == "AG":
            pam_index_list.append((index, "NAG"))

        elif query_seq[index + 1:index + 3] == "GG":
            pam_index_list.append((index, "NGG"))

    return (pam_index_list)


def sign_value(x):
    """
    HELP
        sign function
    """
    if x < 0:
        return (-1)
    elif x == 0:
        return (0)
    else:
        return (1)


def cmp_align_list(align_a, align_b):
    """
    INPUT
        like [17, 3, 0, 73.0, 1, 22 'AAGAAGAAGACGAGTCTGCA', '||||||||||||||X|||XX', 'AAGAAGAAGACGAGCCTGAG']

    HELP
        compare function for align_list
    """
    align_alphabet = {"-": 0, "X": 1, "|": 2}
    sort_index_list = [3, 2, 1, 0]
    sort_rev_state_list = [True, False, False, True]

    for order_index, align_index in enumerate(sort_index_list):
        if (align_a[align_index] - align_b[align_index]) == 0:
            continue
        else:
            if sort_rev_state_list[order_index]:
                return (-1 *
                        sign_value(align_a[align_index] - align_b[align_index]))
            else:
                return (
                    sign_value(
                        align_a[align_index] -
                        align_b[align_index]))

    for index, char_a in enumerate(align_a[7]):
        if index <= (len(align_b[7]) - 1):
            value_a = align_alphabet[char_a]
            value_b = align_alphabet[align_b[7][index]]

            if value_a > value_b:
                return 1
            elif value_a < value_b:
                return -1

    return 0


def run_sgRNA_alignment(align_ref_seq, align_sgRNA, aligner, extend_len=3):
    """
    INPUT
        <ref_seq>

        <align_sgRNA>
            sgRNA seq without PAM

        <possible_sgRNA_region>

    RETURN
        <final_align_res_list>
    """
    align_sgRNA_rev = align_sgRNA[::-1]

    # find all PAM
    PAM_info_list = find_seq_PAM_index(align_ref_seq)
    if len(PAM_info_list) == 0:
        return ([])

    # forward alignment
    final_align_res_list = []

    for PAM_start_idx, PAM_type in PAM_info_list:

        # filter 5' end PAM
        if (PAM_start_idx - extend_len) < len(align_sgRNA):
            continue

        # select PAM and sgRNA in the possible region
        region_seq_start = PAM_start_idx - len(align_sgRNA) - extend_len
        region_seq_end = PAM_start_idx + 3

        # alignment part
        region_seq = align_ref_seq[region_seq_start: PAM_start_idx]
        alignments = aligner.align(region_seq[::-1], align_sgRNA_rev)

        # parse alignment
        # if contain multiple alignment result with the same score, keep the best one;
        # sort reason score -> gap -> mismatch -> match
        align_res_list = []
        for alignment in alignments:
            align_analysis_res = get_alignment_info(alignment, reverse=True)
            align_info_list = [x[::-1] for x in str(alignment).strip().split("\n")]
            align_analysis_res += align_info_list
            align_analysis_res += [PAM_start_idx, PAM_type]
            align_res_list.append(align_analysis_res)

        if len(align_res_list) == 1:
            final_align_res_list.append(align_res_list[0])

        else:
            align_res_list.sort(cmp=cmp_align_list)
            final_align_res_list.append(align_res_list[0])

    # sort final alignment
    final_align_res_list.sort(cmp=cmp_align_list)

    return (final_align_res_list)


def run_no_PAM_sgRNA_alignment_no_chop(
        ref_seq: Seq,
        target_seq: Seq,
        aligner: PairwiseAligner
) -> list:
    """global alignment for target_seq

    :param ref_seq: a Seq object from BioPython, reference sequence for the targeted deep sequencing
    :param target_seq: a Seq object from BioPython, target region sequence for the editing window, usually the sgRNA
        sequencing without PAM
    :param aligner: an obj from BioPython pairwise alignment
    :return: a list: list contains some dict, each one is an alignment result.
    [
        {
            'match_count': 26,
            'mismatch_count': 3,
            'gap_count': 4,
            'alignment_score': 96.0,
            'ref_align_start_index': 0,
            'ref_align_end_index': 32,
            'alignment': {
                'aligned_reference_seq': 'GGCACTGCGGCTGGAAAAAAAAAAAAAAA--GT',
                'aligned_seq_info':      '--..|.|||||||||||||||||||||||--||',
                'aligned_target_seq':    '--GGCAGCGGCTGGAAAAAAAAAAAAAAAAGGT'
            },
            'PAM_start_index': None,
            'PAM_type': None
        },
        ...
    ]

    """
    # set logger
    lib_name = __name__
    function_name = sys._getframe().f_code.co_name
    logger = logging.getLogger(f'{lib_name}.{function_name} ==> ')

    # alignment part
    alignments = aligner.align(ref_seq, target_seq)
    # parse alignment
    # if contain multiple alignment result with the same score, keep the best one;
    # sort reason score -> gap -> mismatch -> match
    alignment_results = []

    for alignment in alignments:
        # logging.debug(type(align))  # Bio.Align.Alignment object
        alignment_analysis_result = get_alignment_info(alignment, reverse=False)
        logger.debug(f'alignment_analysis_result:\n{alignment_analysis_result}')
        alignment_results.append(alignment_analysis_result)

        if len(alignment_results) == 1:
            logger.debug(f'find 1 alignment result:\n{alignment_results}')
        return alignment_results

    else:
        alignment_results.sort(cmp=cmp_align_list)
        logger.debug(f'find more than 1 alignment results:\n{alignment_results}')
        return alignment_results
