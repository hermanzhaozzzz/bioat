from Bio.Align import Alignment, PairwiseAligner
import sys
import logging


def instantiate_pairwise_aligner(
        scoring_match,
        penalty_mismatch,
        penalty_gap_open,
        penalty_gap_extension,
        # penalty_query_left_gap
        mode='global'
):
    # set logger
    lib_name = __name__
    function_name = sys._getframe().f_code.co_name
    logger = logging.getLogger(f'{lib_name}.{function_name} ==> ')

    aligner = PairwiseAligner()
    aligner.match = scoring_match
    aligner.mismatch = penalty_mismatch
    aligner.open_gap_score = penalty_gap_open
    aligner.extend_gap_score = penalty_gap_extension
    aligner.query_left_gap_score = 0
    aligner.query_right_gap_score = 0
    aligner.mode = mode
    logger.info(aligner)
    return aligner

def get_aligned_seq(alignment: Alignment, reverse: bool = False) -> dict:
    """Get_aligned_seq.

    :param alignment: Bio.Align.Alignment object
    :param reverse: return reversed seq info or not
    :return: dict, like this:
         {
            'reference_seq': 'GGCACTGCGGCTGGAAAAAAAAAAAAAAA--GT',
            'aln_info':      '--..|.|||||||||||||||||||||||--||',
            'target_seq':    '--GGCAGCGGCTGGAAAAAAAAAAAAAAAAGGT'
        }
    """
    ls = []

    for x in str(alignment).strip().split("\n"):
        if x:
            ls.append(x[20:].split(' ')[0])
    if not reverse:
        res = dict(
            reference_seq=''.join(ls[::3]),
            aln_info=''.join(ls[1::3]),
            target_seq=''.join(ls[2::3])
        )
    else:
        res = dict(
            reference_seq=''.join(ls[::3])[::-1],
            aln_info=''.join(ls[1::3])[::-1],
            target_seq=''.join(ls[2::3])[::-1]
        )
    return res


def get_alignment_info(alignment: Alignment, reverse: bool = False):
    """

    :param alignment: Bio.Align.Alignment object
    :param reverse: return reversed seq info or not
    :return:
            1. match count
            2. mismatch count
            3. gap count
            4. alignment.score
    :doc:
            The gap count should be the num of gap contain in sgRNA alignment region.

            e.g.

            AGTGGTAAGAAGAAGACGAGACATAATGAG
            ------||||||||||||||X|----||||
            ------AAGAAGAAGACGAGCC----TGAG

            gap count should be 4, rather than 10.

            add return info, start_index, end_index, now the retrun list will be

            return_list = [
                match_count,
                mismatch_count,
                gap_count,
                alignment.score,
                start_index,
                end_index
            ]

            The <start_index> and <end_index> are index related to sgRNA alignment string
    """
    # set logger
    lib_name = __name__
    function_name = sys._getframe().f_code.co_name
    logger = logging.getLogger(f'{lib_name}.{function_name} ==> ')

    if reverse:
        aln_res = get_aligned_seq(alignment, reverse=True)
    else:
        aln_res = get_aligned_seq(alignment)

    reference_seq = aln_res['reference_seq']
    aln_info = aln_res['aln_info']
    target_seq = aln_res['target_seq']
    aln_score = alignment.score
    target_seq_length = len(target_seq)
    target_seq_length_rm_gap = len(target_seq.replace('-', ''))
    aligned_coordinates = alignment.aligned[::-1]

    logger.debug(f'reference_seq = {reference_seq}')
    logger.debug(f'aln_info      = {aln_info}')
    logger.debug(f'target_seq    = {target_seq}')
    logger.debug(f'aln_score     = {aln_score}')
    logger.debug(f'target_seq_length (remove gap) = {target_seq_length_rm_gap}')
    logger.debug(f'aligned_coordinates =\n{aligned_coordinates}')

    # define params
    ref_aln_start = target_seq_length - len(target_seq.lstrip('-'))
    ref_aln_end = ref_aln_start + len((target_seq.strip('-'))) - 1
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
            else:
                # find first aligned base
                # encounter with target_seq!
                encounter_target_seq = True

        if encounter_target_seq and (parsed_target_bases < target_seq_length_rm_gap):
            # from the same loop of [else: catch_start = True]
            # parse start base and later base
            # print(parsed_target_bases, target_seq_length_rm_gap, base_info)
            if base_info == "|":
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
                logger.critical(f'find not surpported alignment symbol: {base_info} ({idx} on the reference), check'
                                f'whether Biopython >= 1.8.0 or not')
                exit(1)
        else:
            # parsed_target_bases = target_seq_length_rm_gap
            # full target_seq bases were processed
            break

    res = dict(
        match_count=match_count,
        mismatch_count=mismatch_count,
        gap_count=gap_count,
        aln_score=alignment.score,
        ref_aln_start=ref_aln_start,
        ref_aln_end=ref_aln_end,
        alignment=aln_res
    )
    return res
