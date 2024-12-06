
import pandas as pd
from Bio.Align import PairwiseAligner
from Bio.Seq import Seq

from bioat.lib.libalignment import get_alignment_info
from bioat.logger import LoggerManager

lm = LoggerManager(mod_name="bioat.lib.libcrispr")

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
    """Compare function for alignment lists.

    This function processes alignment data and performs comparisons based on the input structure.

    Args:
        input_data (list | dict): 
            The alignment data to be processed. Can be provided as:
            - A list, e.g.: [17, 3, 0, 73.0, 1, 22, 'AAGAAGAAGACGAGTCTGCA', '||||||||||||||X|||XX', 'AAGAAGAAGACGAGCCTGAG']
            - A dictionary, e.g.: {'match_count': 26, 'mismatch_count': 3, 'gap_count': 4, 'aln_score': 96.0, 'ref_aln_start': 0, 'ref_aln_end': 32, 'alignment': {'reference_seq': 'GGCACTGCGGCTGGAAAAAAAAAAAAAAA--GT','aln_info': '--..|.|||||||||||||||||||||||--||', 'target_seq': '--GGCAGCGGCTGGAAAAAAAAAAAAAAAAGGT'}}

    Returns:
        None: The function does not return a value but performs comparisons.

    Example:
        >>> compare_alignments([
        ...     17, 3, 0, 73.0, 1, 22, 
        ...     'AAGAAGAAGACGAGTCTGCA', 
        ...     '||||||||||||||X|||XX', 
        ...     'AAGAAGAAGACGAGCCTGAG'
        ... ])
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
    """Perform global alignment for the target sequence.

    This function performs global alignment between a reference sequence and a target sequence (usually the sgRNA sequence without PAM) using the specified aligner.

    Args:
        ref_seq (Seq): 
            A Seq object from BioPython representing the reference sequence 
            for the targeted deep sequencing.
        target_seq (Seq): 
            A Seq object from BioPython representing the target region sequence 
            for the editing window (usually the sgRNA sequence without PAM).
        aligner (object): 
            A PairwiseAligner object from BioPython used for performing the alignment.
        PAM (dict, optional): 
            A dictionary containing PAM information with the following keys:
            - 'PAM' (str): The PAM sequence (e.g., "AGG").
            - 'position' (int): The insertion site of `target_seq`.
            - 'weight' (float): The priority weight for PAM alignment. The alignment score will be multiplied by this weight.
            - Example: {'PAM': 'AGG', 'position': 20, 'weight': 1.0}.

    Returns:
        list[dict]: 
            A list of dictionaries, each containing alignment results. Each dictionary has the following keys:
            - 'match_count' (int): The number of matching bases in the alignment.
            - 'mismatch_count' (int): The number of mismatched bases.
            - 'gap_count' (int): The number of gaps.
            - 'aln_score' (float): The alignment score.
            - 'ref_aln_start' (int): The starting position of the alignment on the reference sequence.
            - 'ref_aln_end' (int): The ending position of the alignment on the reference sequence.
            - 'alignment' (dict): The alignment details containing:'reference_seq' (str): The aligned reference sequence.'aln_info' (str): The alignment information (e.g., matches, mismatches, and gaps).'target_seq' (str): The aligned target sequence.

    Example:
        An example alignment result might look like this:[{'match_count': 26,'mismatch_count': 3,'gap_count': 4,'aln_score': 96.0,'ref_aln_start': 0,'ref_aln_end': 32,'alignment': {'reference_seq': 'GGCACTGCGGCTGGAAAAAAAAAAAAAAA--GT','aln_info': '--..|.|||||||||||||||||||||||--||','target_seq': '--GGCAGCGGCTGGAAAAAAAAAAAAAAAAGGT'}}]

    """

    lm.set_names(func_name="run_target_seq_align")
    lm.set_level(log_level)


    if PAM:
        # considering PAM
        # DNA info https://www.bioinformatics.org/sms/iupac.html
        alignments_pam = aligner.align(ref_seq, Seq(PAM['PAM']))
        # use pam to find aim sequence on ref and then use target_seq to align them and give it a higher priority
        dt_PAMs = {}

        for idx, alignment in enumerate(alignments_pam):
            alignment_analysis_result = get_alignment_info(alignment, reverse=False)
            lm.logger.debug(f'alignment_analysis_result for PAM:\n{alignment_analysis_result}')
            ref_aln_start = alignment_analysis_result['ref_aln_start']
            ref_aln_end = alignment_analysis_result['ref_aln_end']
            aln_score = alignment_analysis_result['aln_score'] * PAM['weight']
            dt_PAMs[f'{ref_aln_start}_{ref_aln_end}'] = dict(aln_score=aln_score * PAM['weight'])

        lm.logger.debug(f'find PAMs on reference: {dt_PAMs}')

        # add PAM left or right or any position
        target_seq = target_seq[:PAM['position']] + Seq(PAM['PAM']) + target_seq[PAM['position']:]
        lm.logger.info(f'find PAM on target_seq, target_seq is updated to: {str(target_seq)}')

    # alignment part
    alignments = aligner.align(ref_seq, target_seq)
    # parse alignment
    alignment_results = []

    for alignment in alignments:
        alignment_analysis_result = get_alignment_info(alignment, reverse=False)
        lm.logger.debug(f'alignment_analysis_result:\n{alignment_analysis_result}')
        alignment_results.append(alignment_analysis_result)

    if len(alignment_results) == 1:
        lm.logger.debug(f'find 1 alignment result:\n{alignment_results}')
        return alignment_results[0]
    else:
        lm.logger.debug(f'find {len(alignment_results)} alignment results:\n{alignment_results}')
        df = pd.DataFrame(alignment_results)
        # if contain multiple alignment result with the same score, keep the best one;
        # sort reason score -> gap -> mismatch -> match
        df.sort_values(
            by=['aln_score', 'gap_count', 'mismatch_count', 'match_count'],
            ascending=[False, True, True, False],
            inplace=True
        )
        return df.iloc[0, :].to_dict()
