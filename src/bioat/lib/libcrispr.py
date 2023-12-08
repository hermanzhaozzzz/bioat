import sys

import pandas as pd
from Bio.Seq import Seq
from Bio.Align import PairwiseAligner
from bioat.lib.libalignment import get_alignment_info
from bioat import get_logger

__module_name__ = 'bioat.lib.libcrispr'

# TARGET_REGIONS_LIB
#     for target_seq alignment
#
# bioat target_seq region_heatmap test_sorted.mpileup.info.tsv test.pdf --get_built_in_target_region
#
# INFO:root:You can use <key> in built-in <target_regions> to represent your target_region:
#         <key>   <target_region>
#         EMX1    GAGTCCGAGCAGAAGAAGAA^NGG^
#         HEK3    GGCCCAGACTGAGCACGTGA^NGG^
#         HEK4    GGCACTGCGGCTGGAGGTGG^NGG^
#         ...

TARGET_SEQ_LIB = {
    "EMX1": "GAGTCCGAGCAGAAGAAGAA^NGG^",
    "HEK3": "GGCCCAGACTGAGCACGTGA^NGG^",
    "HEK4": "GGCACTGCGGCTGGAGGTGG^NGG^",
    "RNF2": "GTCATCTTAGTCATTACCTG^NGG^",
    "VEGFA": "GACCCCCTCCACCCCGCCTC^NGG^",
    "CDKN2A": "^TTTA^GCCCCAATAATCCCCACATGTCA",
    "DYRK1A": "^TTTA^GAAGCACATCAAGGACATTCTAA",
    "ND4": "TGCTAGTAACCACGTTCTCCTGATCAAATATCACTCTCCTACTTACAGGA",
    "ND5.1": "TAGCATTAGCAGGAATACCTTTCCTCACAGGTTTCTACTCCAAAGA",
    "ND6": "TGACCCCCATGCCTCAGGATACTCCTCAATAGCCATCG",
}


def cmp_align_list(aln_a: dict, aln_b: dict):
    """
    INPUT
        like [17, 3, 0, 73.0, 1, 22 'AAGAAGAAGACGAGTCTGCA', '||||||||||||||X|||XX', 'AAGAAGAAGACGAGCCTGAG']
        {
            'match_count': 26,
            'mismatch_count': 3,
            'gap_count': 4,
            'aln_score': 96.0,
            'ref_aln_start': 0,
            'ref_aln_end': 32,
            'alignment': {
                'reference_seq': 'GGCACTGCGGCTGGAAAAAAAAAAAAAAA--GT',
                'aln_info':      '--..|.|||||||||||||||||||||||--||',
                'target_seq':    '--GGCAGCGGCTGGAAAAAAAAAAAAAAAAGGT'
            }
        }
    HELP
        compare function for align_list
    """

    # sort reason score -> gap -> mismatch -> match
    def sign_value(x):
        """
        HELP
            sign function
        """
        if x > 0:
            return 1
        elif x < 0:
            return -1
        else:
            return 0

    align_alphabet_score = {"-": 0, ".": 1, "|": 2}
    sort_keys = ['alignment_score', 'gap_count', 'mismatch_count', 'match_count']
    sort_rev_state_list = [True, False, False, True]

    for idx, key in enumerate(sort_keys):
        if aln_a[key] - aln_b[key] == 0:
            continue
        else:
            if sort_rev_state_list[idx]:
                return -1 * sign_value(aln_a[key] - aln_b[key])
            else:
                return sign_value(aln_a[key] - aln_b[key])

    aln_info_a = aln_a['alignment']['aln_info']
    aln_info_b = aln_b['alignment']['aln_info']

    for idx, char_a in enumerate(aln_info_a):
        if idx <= (len(aln_info_b) - 1):
            value_a = align_alphabet_score[char_a]
            value_b = align_alphabet_score[aln_info_a[idx]]
            return sign_value(value_a - value_b)
        else:
            pass


def run_target_seq_align(
        ref_seq: Seq,
        target_seq: Seq,
        aligner: PairwiseAligner,
        PAM: dict = None,
        log_level='WARNING'
) -> list:
    """global alignment for target_seq

    :param ref_seq: a Seq object from BioPython, reference sequence for the targeted deep sequencing
    :param target_seq: a Seq object from BioPython, target region sequence for the editing window, usually the sgRNA
        sequencing without PAM
    :param aligner: an obj from BioPython pairwise alignment
    :param PAM: a dict like {'PAM': 'AGG', 'position': 20, 'weight': 1.0},
        position is the insertion site of target_seq,
        weight is the priority, PAM alignment score will multiply by this weight.
    :return: a list: list contains some dict, each one is an alignment result.
    [
        {
            'match_count': 26,
            'mismatch_count': 3,
            'gap_count': 4,
            'aln_score': 96.0,
            'ref_aln_start': 0,
            'ref_aln_end': 32,
            'alignment': {
                'reference_seq': 'GGCACTGCGGCTGGAAAAAAAAAAAAAAA--GT',
                'aln_info':      '--..|.|||||||||||||||||||||||--||',
                'target_seq':    '--GGCAGCGGCTGGAAAAAAAAAAAAAAAAGGT'
            }
        },
        ...
    ]

    """
    # set logger
    logger = get_logger(level=log_level, module_name=__module_name__, func_name=sys._getframe().f_code.co_name)

    if PAM:
        # considering PAM
        # DNA info https://www.bioinformatics.org/sms/iupac.html
        alignments_pam = aligner.align(ref_seq, Seq(PAM['PAM']))
        # use pam to find aim sequence on ref and then use target_seq to align them and give it a higher priority
        dt_PAMs = {}

        for idx, alignment in enumerate(alignments_pam):
            alignment_analysis_result = get_alignment_info(alignment, reverse=False)
            logger.debug(f'alignment_analysis_result for PAM:\n{alignment_analysis_result}')
            ref_aln_start = alignment_analysis_result['ref_aln_start']
            ref_aln_end = alignment_analysis_result['ref_aln_end']
            aln_score = alignment_analysis_result['aln_score'] * PAM['weight']
            dt_PAMs[f'{ref_aln_start}_{ref_aln_end}'] = dict(aln_score=aln_score * PAM['weight'])

        logger.debug(f'find PAMs on reference: {dt_PAMs}')

        # add PAM left or right or any position
        target_seq = target_seq[:PAM['position']] + Seq(PAM['PAM']) + target_seq[PAM['position']:]
        logger.info(f'find PAM on target_seq, target_seq is updated to: {str(target_seq)}')

    # alignment part
    alignments = aligner.align(ref_seq, target_seq)
    # parse alignment
    alignment_results = []

    for alignment in alignments:
        alignment_analysis_result = get_alignment_info(alignment, reverse=False)
        logger.debug(f'alignment_analysis_result:\n{alignment_analysis_result}')
        alignment_results.append(alignment_analysis_result)

    if len(alignment_results) == 1:
        logger.debug(f'find 1 alignment result:\n{alignment_results}')
        return alignment_results[0]
    else:
        logger.debug(f'find {len(alignment_results)} alignment results:\n{alignment_results}')
        df = pd.DataFrame(alignment_results)
        # if contain multiple alignment result with the same score, keep the best one;
        # sort reason score -> gap -> mismatch -> match
        df.sort_values(
            by=['aln_score', 'gap_count', 'mismatch_count', 'match_count'],
            ascending=[False, True, True, False],
            inplace=True
        )
        return df.iloc[0, :].to_dict()

