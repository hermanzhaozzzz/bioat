from Bio.Align import Alignment


def parse_alignment_str(alignment: Alignment) -> dict:
    """parse_alignment_str.

    :param alignment: Bio.Align.Alignment object
    :return: dict, like this:
         {
            'aligned_reference_seq': 'GTCACCTGCCTCTGGAGAGGGAGGAGGGGCCTCT',
            'aligned_seq_info':      '-------|.|.|||..|..|||||.||-------',
            'aligned_target_seq':    '-------GGCACTGCGGCTGGAGGTGG-------'
        }
    """
    ls = []

    for x in str(alignment).strip().split("\n"):
        if x:
            ls.append(x[20:].split(' ')[0])
    res = dict(
        aligned_reference_seq=''.join(ls[::3]),
        aligned_seq_info=''.join(ls[1::3]),
        aligned_target_seq=''.join(ls[2::3])
    )
    return res
