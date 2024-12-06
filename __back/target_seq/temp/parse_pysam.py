# import logging
# import sys

# import pysam

# u"""说明文档
# 本文件中构造了用于解析pysam.AlignedSegment对象中
# """


# def back_indel_shift(info_index_list, cur_index) -> int:
#     """return acc_shift(back_indel_shift)

#     Args:
#         info_index_list (list/tuple): list or tuples generated from align.cigar tuples
#         cur_index (int): index related to MD tag in BAM file

#     Returns:
#         int: acc_shift index
#     """
#     # parse soft clip and insertion
#     if len(info_index_list) == 0:
#         return 0

#     acc_shift = 0
#     for info_start_index, info_len in info_index_list:
#         if info_start_index >= cur_index:
#             return acc_shift
#         else:
#             acc_shift += info_len
#     return acc_shift


# def get_align_mismatch_pairs(align, ref_genome_dict=None) -> list:
#     """input a pysam AlignedSegment object

#     Args:
#         align (pysam.AlignedSeqment object): pysam.AlignedSeqment object
#         ref_genome_dict (dict, optional): returned dict from load_reference_fasta_as_dict(). Defaults to None.

#     Returns:
#         list/None:
#             it returns mismatch_pair_list, just like [ref_index, align_index, ref_base, align_base];
#             and the "ref_index" is the same coordinate with UCSC genome browser;
#             When NM == 0, it returns None.
#     """
#     # No mismatch
#     try:
#         if align.get_tag("NM") == 0:
#             return None
#     except:
#         return None

#     MD_tag_state = align.has_tag("MD")

#     if MD_tag_state:
#         # parse softclip, insertion and deletion
#         info_index_list = []
#         accu_index = 0

#         for cigar_type, cigar_len in align.cigartuples:
#             if cigar_type == 1 or cigar_type == 4:
#                 info_index_list.append((accu_index + 1, cigar_len))

#             elif cigar_type == 2:
#                 info_index_list.append((accu_index + 1, -cigar_len))

#             accu_index += cigar_len

#         # parse MD tag
#         mismatch_pair_list = []
#         cur_base = ""
#         cur_index = 0
#         bases = align.get_tag("MD")

#         i = 0
#         while i < len(bases):
#             base = bases[i]

#             if base.isdigit():
#                 cur_base += base
#                 i += 1

#             else:
#                 cur_index += int(cur_base)
#                 cur_base = ""

#                 if base == "^":
#                     i += 1
#                     del_str = ""

#                     while (bases[i].isalpha()) and (i < len(bases)):
#                         del_str += bases[i]
#                         i += 1

#                     cur_index += len(del_str)
#                     del_str = ""

#                 elif base.isalpha():
#                     cur_index += 1
#                     ref_base = base
#                     i += 1

#                     # add into list
#                     fix_index = cur_index + back_indel_shift(info_index_list, cur_index)

#                     if fix_index < len(align.query_sequence):
#                         mismatch_pair_list.append(
#                             [
#                                 cur_index + align.reference_start,
#                                 cur_index - 1,
#                                 ref_base,
#                                 align.query_sequence[fix_index - 1],
#                             ]
#                         )
#                     else:
#                         return None

#         return mismatch_pair_list
#     else:
#         mismatch_pair_list = []
#         for align_idx, ref_idx in align.get_aligned_pairs():
#             if (align_idx is not None) and (ref_idx is not None):
#                 align_base = align.query_sequence[align_idx]
#                 ref_base = ref_genome_dict[align.reference_name][ref_idx]

#                 if align_base != ref_base:
#                     mismatch_pair_list.append(
#                         [ref_idx + 1, align_idx, ref_base, align_base]
#                     )

#         return mismatch_pair_list


# def get_query_site_mpileup_info(
#     site_idx_list,
#     bam_obj,
#     genome_obj,
#     ignore_overlaps=False,
#     min_base_quality=20,
#     min_mapping_quality=50,
#     mpileup_extend_length=10,
#     dist_cutoff=200,
#     max_depth=4000,
#     log_verbose=logging.DEBUG,
# ) -> dict:
#     """get each base mutation count info from the input list (site_idx_list)

#     Args:
#         site_idx_list (list):
#             It is a list format like:
#                 ["chr1_125178576_CT", "chr1_125178578_CT", "chr1_125178580_CT", "chr1_125178588_CA"]
#                 or
#                 ["chr1_125178576", "chr1_125178578", "chr1_125178580", "chr1_125178588"]
#             The site coordinate is related to the UCSC genoem browser index.
#             Site index have to be sorted as ascending order.
#         bam_obj (pysam.AlignmentFile object):
#             pysam.AlignmentFile() object
#         genome_obj (pysam.FastaFile object):
#             genome object, which is created by pysam.FastaFile()
#         ignore_overlaps (bool, optional):
#             Whether to ignore overlapping regions.
#             Defaults to False.
#         min_base_quality (int, optional):
#             Minimum base quality.
#             Defaults to 20.
#         min_mapping_quality (int, optional):
#             Minimum mapping quality.
#             Defaults to 50.
#         mpileup_extend_length (int, optional):
#             To prevent reads from being missed, set an appropriate number to extend fetching region.
#             Defaut number is fun.
#             Defaults to 10.
#         dist_cutoff (int, optional):
#             The number of bases to be processed at one time.
#             An appropriate number, can improve performance.
#             Defaults to 200.
#         max_depth (int, optional):
#             Pileup file max depth.
#             Defaults to 4000.
#         log_verbose (logging object, optional):
#             A logging states object, which can be logging.DEBUG, logging.INFO, logging.WARNING, logging.ERROR, logging.CRITICAL.
#             Defaults to logging.DEBUG.
#     Raises:
#         IOError: [description]
#         IOError: [description]

#     Returns:
#         dict/None:
#             key:
#                 site_index only key chr_name and chr_pos like "chr1_125178576" rather than "chr1_125178576_CT"
#             value:
#                 A, T, C, G, N, other count
#             If an empty list is given, it will return None.
#     """
#     # log setting
#     logging.basicConfig(
#         level=log_verbose,
#         format="%(levelname)-5s @ %(asctime)s: %(message)s ",
#         datefmt="%Y-%m-%d %H:%M:%S",
#         stream=sys.stderr,
#         filemode="w",
#         force=True,
#     )
#     # check site index list
#     if len(site_idx_list) == 0:
#         return None

#     chr_name = site_idx_list[0].split("_")[0]
#     region_start = region_end = int(site_idx_list[0].split("_")[1])
#     site_dist = 0

#     for site_idx in site_idx_list:
#         if site_idx.split("_")[0] != chr_name:
#             logging.error("<site_idx_list> error: %s!" % site_idx)
#             raise IOError(
#                 "<site_idx_list> all site index have to on the same chromosome."
#             )

#     if len(site_idx_list) >= 2:
#         region_start = int(site_idx_list[0].split("_")[1])
#         region_end = int(site_idx_list[-1].split("_")[1])
#         site_dist = region_end - region_start

#     if site_dist >= dist_cutoff:
#         raise IOError(
#             "<site_idx_list> mpileup region is too large, which could take a lot of time!"
#         )

#     # make raw dict
#     site_dict = {}
#     for site_idx in site_idx_list:
#         key = "_".join(site_idx.split("_")[:2])
#         site_dict[key] = {"A": 0, "T": 0, "C": 0, "G": 0, "N": 0, "total": 0}

#     # run mpileup
#     mpileup_iter = bam_obj.pileup(
#         chr_name,
#         region_start - mpileup_extend_length,
#         region_end + mpileup_extend_length,
#         fastafile=genome_obj,
#         max_depth=max_depth,
#         ignore_overlaps=ignore_overlaps,
#         min_base_quality=min_base_quality,
#         min_mapping_quality=min_mapping_quality,
#     )

#     for pileup in mpileup_iter:
#         run_index = "%s_%s" % (chr_name, pileup.reference_pos + 1)

#         if (pileup.reference_pos + 1) > region_end:
#             break

#         if site_dict.get(run_index) is not None:
#             base_list = [i.upper() for i in pileup.get_query_sequences()]
#             site_dict[run_index]["A"] = base_list.count("A")
#             site_dict[run_index]["T"] = base_list.count("T")
#             site_dict[run_index]["C"] = base_list.count("C")
#             site_dict[run_index]["G"] = base_list.count("G")
#             site_dict[run_index]["N"] = base_list.count("N")
#             site_dict[run_index]["total"] = len(base_list)

#     return site_dict


# if __name__ == "__main__":
#     pt_bam = "/Users/zhaohuanan/zhaohn_HD/3.project/2021_DdCBE_topic/20210224_DetectSeq_all_bams/bam/293T-DdCBE-ND6-All-PD_rep2_hg38.MAPQ20.bam"
#     pt_ref = "/Users/zhaohuanan/zhaohn_HD/2.database/db_genomes/genome_fa/genome_ucsc_hg38/chr20.fa"

#     # load bam
#     bam_file = pysam.AlignmentFile(pt_bam, "rb")
#     # load genome
#     genome_ref = pysam.FastaFile(filename=pt_ref)
#     get_query_site_mpileup_info(
#         site_idx_list=["chr20_53025829", "chr20_53025830"],
#         bam_obj=bam_file,
#         genome_obj=genome_ref,
#     )
