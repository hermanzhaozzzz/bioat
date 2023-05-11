# from __future__ import absolute_import, division, print_function
#
# import logging
# import sys
# import os
# import gzip
# import random
# import string
# import json
# import pysam
# import time
# import math
# import multiprocessing
# from scipy import stats
# from Bio import SeqIO
#
# from statsmodels.stats import proportion
# from statsmodels.stats.multitest import multipletests
# import numpy as np
#
#
# class mpmatLine(object):
#     """
#     INPUT:
#         <mpmat_line_list>
#             list, info like
#             [
#                 'chr1',
#                 '629627',
#                 '629632',
#                 '4',
#                 '4',
#                 '0',
#                 'chr1_629627_CT,chr1_629628_CT,chr1_629631_CT,chr1_629632_CT',
#                 '3,5,3,4',
#                 '37,38,37,38',
#                 '0.08108,0.13158,0.08108,0.10526',
#                 'False,False,False,False',
#                 '0,0,0,0',
#                 'Pass,Pass,Pass,Pass'
#             ]
#
#         <filter_info_index>
#             int, col index describe filter info, default=None, start at 1
#
#         <block_info_index>
#             int, col index describe block info, default=None, start at 1
#
#     RETURN:
#         <mpmatLine obj>
#     """
#
#     def __init__(self, mpmat_line_list, filter_info_index=None, block_info_index=None):
#         # standard .mpmat file
#         self.chr_name = mpmat_line_list[0]
#         self.chr_start = int(mpmat_line_list[1])
#         self.chr_end = int(mpmat_line_list[2])
#
#         self.site_num = int(mpmat_line_list[3])
#         self.mut_site_num = int(mpmat_line_list[4])
#         self.SNP_site_num = int(mpmat_line_list[5])
#
#         self.site_index_list = mpmat_line_list[6].split(",")
#         self.mut_count_list = list(map(int, mpmat_line_list[7].split(",")))
#         self.cover_count_list = list(map(int, mpmat_line_list[8].split(",")))
#         self.mut_ratio_list = list(map(float, mpmat_line_list[9].split(",")))
#         self.SNP_ann_list = list(map(eval, mpmat_line_list[10].split(",")))
#         self.tandem_info_list = mpmat_line_list[11].split(",")
#
#         # region str
#         self.region = "%s:%s-%s" % (mpmat_line_list[0], mpmat_line_list[1], mpmat_line_list[2])
#
#         # get mut type
#         self.mut_type = mpmat_line_list[6].split(",")[0].split("_")[-1]
#
#         # mutation key
#         mut_key_list = ["N"] * self.site_num
#         for index, site_index in enumerate(self.site_index_list):
#             if self.SNP_ann_list[index]:
#                 mut_key_list[index] = "S"
#
#         self.mut_key = "-".join(mut_key_list)
#         self.mut_key_list = mut_key_list
#
#         # load filter
#         if filter_info_index is not None:
#             try:
#                 self.filter_state_list = mpmat_line_list[filter_info_index - 1].split(",")
#             except:
#                 raise IOError("Parsing error occur at <filter_info_index>")
#
#         if block_info_index is not None:
#             try:
#                 self.block_info_list = list(map(eval, mpmat_line_list[block_info_index - 1].split(",")))
#                 self.block_site_num = self.block_info_list.count(True)
#             except:
#                 raise IOError("Parsing error occur at <block_info_index>")
#
#
# class siteIndex(object):
#     """
#     INPUT
#         <site_index>
#             str, like chr1_10000_CT
#
#     RETURN
#         siteIndex obj
#     """
#
#     def __init__(self, _site_index):
#         site_index_split = _site_index.split("_")
#
#         self.chr_name = site_index_split[0]
#         self.site_pos = int(site_index_split[1])
#         self.site_index = site_index_split[0] + "_" + site_index_split[1]
#         self.site_index_raw = _site_index
#
#         if len(site_index_split) > 2:
#             self.mut_type = site_index_split[2]
#
#
# class bmatLine(object):
#     """
#     INPUT:
#         <bmat_line_list>
#             list, info like
#             [
#                'chr1',
#                '10013',
#                'T',
#                '0',
#                '0',
#                '0',
#                '6',
#                '0',
#                '0',
#                '0',
#                '.',
#                '.',
#                '.',
#                '0'
#             ]
#
#     RETURN:
#         A bmatLine obj
#     """
#
#     def __init__(self, bmat_line_list):
#         self.chr_name = bmat_line_list[0]
#         self.chr_index = int(bmat_line_list[1])
#         self.ref_base = bmat_line_list[2]
#
#         self.count_del = int(bmat_line_list[7])
#         self.count_insert = int(bmat_line_list[8])
#         self.count_ambis = int(bmat_line_list[9])
#
#         self.deletion = bmat_line_list[10]
#         self.ambiguous = bmat_line_list[11]
#         self.insertion = bmat_line_list[12]
#         self.mut_num = int(bmat_line_list[13])
#
#         # make dict
#         self.count_dict = {
#             "A": int(bmat_line_list[3]),
#             "G": int(bmat_line_list[4]),
#             "C": int(bmat_line_list[5]),
#             "T": int(bmat_line_list[6])
#         }
#
#
# #################################################################################
# # FUN
# #################################################################################
# def check_mpmat_test_state_by_block_info(mpmat_line_obj, block_site_min_num_cutoff=1, block_num_ratio_cutoff=0.8,
#                                          block_num_ratio_check_min_num=5, block_num_max_num_cutoff=15):
#     """
#     INPUT:
#         <mpmat_line_obj>
#             obj, mpmatLine obj contain BlockInfo annotation
#
#         <block_site_min_num_cutoff>
#             int, non-block site num less than this cutoff will return False
#
#         <block_num_ratio_check_min_num>
#             int, non-block site num less than this will not check ratio
#
#         <block_num_ratio_cutoff>
#             float, when non-block site num >= <block_num_ratio_check_min_num> check ratio.
#                 If block site num >= int(<block_num_ratio_check_min_num> * <block_num_ratio_cutoff>), return False.
#
#         <block_num_max_num_cutoff>
#             int, block site num >= this value, return False
#
#     RETURN:
#         <test_state>
#             bool, TRUE means region need to be processed Poisson test
#                   False means region will be omitted
#
#         <state_reason>
#             str, 'NoSignalSite' 'BlockSiteRatio' 'BlockSiteNum'
#
#     """
#     # load block number
#     block_site_num = mpmat_line_obj.block_state_list.count(True)
#     non_block_site_num = mpmat_line_obj.block_state_list.count(False)
#
#     if non_block_site_num < block_site_min_num_cutoff:
#         return False, "NoSignalSite"
#
#     if block_site_num >= block_num_ratio_check_min_num:
#         block_site_ratio = block_site_num / 1.0 / len(mpmat_line_obj.block_state_list)
#         if block_site_ratio >= block_num_ratio_cutoff:
#             return False, "BlockSiteRatio"
#
#     if block_site_num >= block_num_max_num_cutoff:
#         return False, "BlockSiteNum"
#
#     return True, None
#
#
# #################################################################################
# # FUN
# #################################################################################
# def back_indel_shift(info_index_list, cur_index):
#     """
#     INPUT:
#         <info_index_list>
#             generated from align.cigar tuples
#
#         <cur_index>
#             index related to MD tag in BAM file
#
#     RETURN
#         <acc_shift>
#     """
#
#     # parse soft clip and insertion
#     if len(info_index_list) == 0:
#         return 0
#
#     acc_shift = 0
#     for info_start_index, info_len in info_index_list:
#         if info_start_index >= cur_index:
#             return acc_shift
#
#         else:
#             acc_shift += info_len
#
#     return acc_shift
#
#
# #################################################################################
# # FUN
# #################################################################################
# def get_align_mismatch_pairs(align):
#     """
#     INPUT
#         <align>
#             pysam AlignedSegment object
#
#     RETURN
#         <mismatch_pair_list>
#             [ref_index, align_index, ref_base, align_base]
#
#             ref_index is the same coordinate with UCSC genome browser
#
#             When NM == 0, return None
#     """
#     # Hisat-3n aligner NM will be set as 0, but with several conversions.
#     # So I have to comment out this part
#     # 2022-07-29
#
#     # No mismatch
#     # try:
#     #     if align.get_tag("NM") == 0:
#     #         return None
#     # except:
#     #     return None
#
#     # parse soft clip, insertion and deletion
#     info_index_list = []
#     accu_index = 0
#
#     for cigar_type, cigar_len in align.cigartuples:
#         if cigar_type == 1 or cigar_type == 4:
#             info_index_list.append((accu_index + 1, cigar_len))
#
#         elif cigar_type == 2:
#             info_index_list.append((accu_index + 1, -cigar_len))
#
#         accu_index += cigar_len
#
#     # parse MD tag
#     mismatch_pair_list = []
#     cur_base = ""
#     cur_index = 0
#     bases = align.get_tag("MD")
#
#     i = 0
#     while i < len(bases):
#         base = bases[i]
#
#         if base.isdigit():
#             cur_base += base
#             i += 1
#
#         else:
#             if cur_base != "":
#                 cur_index += int(cur_base)
#
#             cur_base = ""
#
#             if base == "^":
#                 i += 1
#                 del_str = ""
#
#                 while (bases[i].isalpha()) and (i < len(bases)):
#                     del_str += bases[i]
#                     i += 1
#
#                 cur_index += len(del_str)
#
#             elif base.isalpha():
#                 cur_index += 1
#                 ref_base = base
#                 i += 1
#
#                 # add into list
#                 fix_index = cur_index + back_indel_shift(info_index_list, cur_index)
#
#                 if fix_index < len(align.query_sequence):
#                     mismatch_pair_list.append([cur_index + align.reference_start, cur_index - 1, ref_base,
#                                                align.query_sequence[fix_index - 1]])
#                 else:
#                     return None
#
#     if len(mismatch_pair_list) == 0:
#         return None
#     else:
#         return mismatch_pair_list
#
#
# #################################################################################
# # FUN
# #################################################################################
# def get_No_MD_align_mismatch_pairs(align, ref_genome_dict):
#     """
#     INPUT
#         <align>
#             pysam AlignedSegment object
#
#         <genome_dict>
#             key is like chr1, chr2, ...
#             value is chromosome sequence
#
#     RETURN
#         <mismatch_pair_list>
#             [ref_index, align_index, ref_base, align_base]
#
#             ref_index is the same coordinate with UCSC genome browser
#
#             When NM == 0, return None
#
#     """
#     # No mismatch
#
#     # Hisat-3n aligner NM will be set as 0, but with several conversions.
#     # So I have to comment out this part
#     # 2022-07-29
#
#     # try:
#     #     if align.get_tag("NM") == 0:
#     #         return None
#     # except:
#     #     return None
#
#     mismatch_pair_list = []
#     for align_idx, ref_idx in align.get_aligned_pairs():
#         if (align_idx is not None) and (ref_idx is not None):
#             align_base = align.query_sequence[align_idx]
#             ref_base = ref_genome_dict[align.reference_name][ref_idx]
#
#             if align_base != ref_base:
#                 mismatch_pair_list.append([
#                     ref_idx + 1,
#                     align_idx,
#                     ref_base,
#                     align_base
#                 ])
#
#     if len(mismatch_pair_list) == 0:
#         return None
#     else:
#         return mismatch_pair_list
#
#
# #################################################################################
# # FUN
# #################################################################################
# def analyse_align_mut_state(mpmat_info, align_mismatch_pairs, query_mut_type, site_index_dict, snp_index_dict,
#                             block_index_dict):
#     """
#     INPUT:
#         <mpmat_info>
#             obj, mpmat info
#
#         <align_mismatch_pairs>
#             list, format like
#                 [[2474644, 81, 'G', 'A'], [2474656, 93, 'G', 'A'], [2474679, 116, 'C', 'T']]
#
#         <query_mut_type>
#             str, like "GA"
#
#         <site_index_dict>
#             dict, key is pos as str like '2474644', value is site order in mpmat_info.site_index_list
#
#         <snp_index_dict>
#             dict, key is pos as str like '2474644', value is site order in mpmat_info.site_index_list
#
#         <block_index_dict>
#             dict, key is pos as str like '2474644', value is site order in mpmat_info.site_index_list
#
#     RETURN:
#         <align_mut_state> str
#             "M-M-M" means tandem mutation,
#             "M-B-M" means mut, block, mut
#             "M-S-M" mean mut, SNV, mut
#             "N-N-N" means no mut
#
#         <query_mut_count>
#             int, query type mutation count on whole alignment read
#
#         <other_mut_count>
#             int, other position or other type mutation in query position
#
#         <total_mut_count>
#             int, total number of mutation on whole alignment read
#
#         <region_query_mut_count>
#             int, total number of reads with query mutation in mpmat region
#
#     VERSION:
#         Final edition date: 2022-07-31
#     """
#
#     # var init
#     total_mut_count = 0
#     other_mut_count = 0
#     query_mut_count = 0
#     region_query_mut_count = 0
#
#     align_mut_state_list = mpmat_info.mut_key_list[:]
#
#     for site_mis_info in align_mismatch_pairs:
#         total_mut_count += 1
#
#         site_index, site_align_pos, from_base, to_base = site_mis_info
#
#         # check block
#         block_site_order = block_index_dict.get(str(site_index))
#
#         # check SNP
#         snp_site_order = snp_index_dict.get(str(site_index))
#
#         if block_site_order is not None:
#             continue
#
#         if snp_site_order is not None:
#             continue
#
#         # not in block info list
#         site_order = site_index_dict.get(str(site_index))
#
#         # count mutation num
#         if (from_base == query_mut_type[0]) and (to_base == query_mut_type[1]):
#             query_mut_count += 1
#
#             # check mutation site only in mpmat region
#             if site_order is not None:
#                 region_query_mut_count += 1
#                 align_mut_state_list[site_order] = "M"
#
#         else:
#             other_mut_count += 1
#
#     # align mut state
#     align_mut_state = "-".join(align_mut_state_list)
#
#     return align_mut_state, query_mut_count, other_mut_count, total_mut_count, region_query_mut_count
#
#
# #################################################################################
# # FUN
# #################################################################################
# # make a function
# def get_mpmat_region_count(in_bam_obj, mpmat_info, ref_genome_dict, query_mut_type, query_mut_min_cutoff=1,
#                            query_mut_max_cutoff=16, total_mut_max_cutoff=20, other_mut_max_cutoff=16):
#     """
#     INPUT:
#         <in_bam_obj>
#             obj, pysam.AlignmentFile
#
#         <mpmat_info>
#             obj, mpmatLine obj
#
#         <ref_genome_dict>
#             dict, key is chr_name, value is reference sequence
#
#         <query_mut_min_cutoff>
#             int, if mutation number >= query_mut_min_cutoff in the mpmat region,
#                 will be counted as 'region_mut_count'.
#
#         <query_mut_max_cutoff>
#             int, larger than this will be marked.
#
#         <total_mut_max_cutoff>
#             int, larger than this will be marked.
#
#         <other_mut_max_cutoff>
#             int, larger than this will be marked.
#
#     RETURN:
#         <align_count_dict>
#             dict, key and value like:
#                 {
#                     "all_align_count": 0,
#                     "region_mut_count": 0,
#                     "region_non_mut_count": 0,
#                     "total_high_mismatch_count": 0,
#                     "other_high_mismatch_count": 0,
#                     "query_high_mismatch_count": 0,
#                     "all_filter_count": 0
#                 }
#
#         <align_mut_tandem_dict>
#             dict, key is  <align_mut_state> (refer to FUN analyse_align_mut_state)
#                 "M-M-M" means tandem mutation,
#                 "M-B-M" means mut, block, mut
#                 "M-S-M" means mut, SNP, mut
#                 "N-N-N" means no mut
#
#         <align_mut_count_dict>
#             dict, key is mutation number like 0,1,2,3,4.... value is align reads count
#
#         <mpmat_info>
#             obj, mpmatLine obj, change <mut_key_list> and <mut_key> according to mpmat block info
#
#     VERSION:
#         Final edition date: 2022-07-31
#     """
#     # ---------------------------------------------------------->>>>>>>>>>
#     # init var
#     # ---------------------------------------------------------->>>>>>>>>>
#     # site index dict
#     site_index_dict = {}
#
#     # site SNP dict
#     snp_site_index_dict = {}
#
#     # site block dict
#     block_site_index_dict = {}
#
#     # fix mpmat_info
#     for index, site_index in enumerate(mpmat_info.site_index_list):
#         site_index_dict[site_index.split("_")[1]] = index
#
#         if mpmat_info.block_info_list[index]:
#             block_site_index_dict[site_index.split("_")[1]] = index
#
#         if mpmat_info.SNP_ann_list[index]:
#             snp_site_index_dict[site_index.split("_")[1]] = index
#
#     # make non mut key
#     non_mut_key = mpmat_info.mut_key
#     align_mut_tandem_dict = {non_mut_key: 0}
#
#     # define count dict
#     align_count_dict = {
#         "all_align_count": 0,
#         "region_mut_count": 0,
#         "region_non_mut_count": 0,
#         "total_high_mismatch_count": 0,
#         "other_high_mismatch_count": 0,
#         "query_high_mismatch_count": 0,
#         "all_filter_count": 0
#     }
#
#     # mutation count dict
#     align_mut_count_dict = {}
#     for mut_num in range(mpmat_info.site_num + 1):
#         align_mut_count_dict[mut_num] = 0
#
#     # ---------------------------------------------------------->>>>>>>>>>
#     # iter align info
#     # ---------------------------------------------------------->>>>>>>>>>
#     for align in in_bam_obj.fetch(reference=mpmat_info.chr_name,
#                                   start=mpmat_info.chr_start - 1,
#                                   end=mpmat_info.chr_end + 1):
#
#         # count total
#         align_count_dict["all_align_count"] += 1
#
#         # make sure MD state
#         MD_state = True
#         try:
#             MD_str_tag = align.get_tag("MD")
#         except:
#             MD_state = False
#
#         # get mismatch pairs
#         if MD_state:
#             align_mismatch_pairs = get_align_mismatch_pairs(align)
#         else:
#             align_mismatch_pairs = get_No_MD_align_mismatch_pairs(align, ref_genome_dict)
#
#         # analysis of mismatch pairs
#         if align_mismatch_pairs is None:
#             align_count_dict["region_non_mut_count"] += 1
#             align_mut_count_dict[0] += 1
#             align_mut_tandem_dict[non_mut_key] += 1
#
#         else:
#             align_mut_analyse_res = analyse_align_mut_state(mpmat_info=mpmat_info,
#                                                             align_mismatch_pairs=align_mismatch_pairs,
#                                                             query_mut_type=query_mut_type,
#                                                             site_index_dict=site_index_dict,
#                                                             snp_index_dict=snp_site_index_dict,
#                                                             block_index_dict=block_site_index_dict)
#
#             if align_mut_tandem_dict.get(align_mut_analyse_res[0]) is None:
#                 align_mut_tandem_dict[align_mut_analyse_res[0]] = 1
#             else:
#                 align_mut_tandem_dict[align_mut_analyse_res[0]] += 1
#
#             # load mut num
#             query_mut_count, other_mut_count, total_mut_count, region_query_mut_count = align_mut_analyse_res[1:]
#             align_mut_count_dict[region_query_mut_count] += 1
#
#             # filter high mismatch reads
#             if total_mut_count >= total_mut_max_cutoff:
#                 align_count_dict["total_high_mismatch_count"] += 1
#                 align_count_dict["all_filter_count"] += 1
#
#             elif other_mut_count >= other_mut_max_cutoff:
#                 align_count_dict["other_high_mismatch_count"] += 1
#                 align_count_dict["all_filter_count"] += 1
#
#             elif query_mut_count >= query_mut_max_cutoff:
#                 align_count_dict["query_high_mismatch_count"] += 1
#                 align_count_dict["all_filter_count"] += 1
#
#             else:
#                 if region_query_mut_count == 0:
#                     align_count_dict["region_non_mut_count"] += 1
#
#                 elif region_query_mut_count >= query_mut_min_cutoff:
#                     align_count_dict["region_mut_count"] += 1
#
#     return align_count_dict, align_mut_tandem_dict, align_mut_count_dict
#
#
# #################################################################################
# # Statistics Function part
# #################################################################################
# def _zstat_generic_parse(value, std_diff, alternative):
#     """
#     INPUT:
#         <value>
#             Stats
#
#         <std_diff>
#             Std
#
#         <alternative>
#             Alternative string
#
#     RETURN:
#         z-score, pvalue
#
#     HELP:
#         Copied and fixed from statsmodels.stats.weight stats
#     """
#
#     zstat = value / std_diff
#
#     if alternative in ['two-sided', '2-sided', '2s']:
#         pvalue = stats.norm.sf(np.abs(zstat)) * 2
#
#     elif alternative in ['larger', 'l']:
#         pvalue = stats.norm.sf(zstat)
#
#     elif alternative in ['smaller', 's']:
#         pvalue = stats.norm.cdf(zstat)
#
#     else:
#         raise ValueError('invalid alternative')
#
#     return zstat, pvalue
#
#
# def poisson_test(lambda_1, lambda_2, method='sqrt', alternative='two-sided'):
#     """
#     INPUT:
#
#         <lambda_1>, <lambda_2>
#             Poisson parameter in case1 and case2
#
#         <method>
#             Support method:
#                 'wald': method W1A, wald test, variance based on separate estimates
#                 'score': method W2A, score test, variance based on estimate under Null
#                 'sqrt': W5A, based on variance stabilizing square root transformation
#
#                 'exact-cond': exact conditional test based on binomial distribution
#                 'cond-midp': midpoint-pvalue of exact conditional test
#
#         <alternative>
#             'two-sided'
#                 means H0: lambda_1 = lambda_2
#
#             'larger'
#                 means H0: lambda_1 > lambda_2
#
#             'smaller'
#                 means H0: lambda_1 < lambda_2
#
#     HELP:
#
#         Reference paper:
#             Gu, Ng, Tang, Schucany 2008: Testing the Ratio of Two Poisson Rates,
#             Biometrical Journal 50 (2008) 2, 2008
#
#         Raw code is copied from html and then fixed realted to DETECT-Seq project requirement
#             https://stackoverflow.com/questions/33944914/implementation-of-e-test-for-poisson-in-python
#     """
#
#     # calculate stat
#     dist = ""
#     stat = 0
#
#     if method in ['score']:
#         stat = (lambda_1 - lambda_2) / np.sqrt((lambda_1 + lambda_2))
#         dist = 'normal'
#
#     elif method in ['wald']:
#         stat = (lambda_1 - lambda_2) / np.sqrt((lambda_1 + lambda_2))
#         dist = 'normal'
#
#     elif method in ['sqrt']:
#         stat = 2 * (np.sqrt(lambda_1 + 3 / 1.0 / 8) - np.sqrt((lambda_2 + 3 / 1.0 / 8)))
#         stat /= np.sqrt(2)
#         dist = 'normal'
#
#     elif method in ['exact-cond', 'cond-midp']:
#         lambda_total = lambda_1 + lambda_2
#         stat = None
#         pvalue = proportion.binom_test(lambda_1, lambda_total, prop=0.5, alternative=alternative)
#         return stat, pvalue
#
#     # return part
#     if dist == 'normal':
#         return _zstat_generic_parse(stat, 1, alternative)
#
#
# def run_mpmat_poisson_test(
#         mpmat_filename, out_poisson_filename, ctrl_bam_filename, treat_bam_filename, ref_genome_fa_filename,
#         select_chr_name_list, scale_factor_dict=None, normalize_scale_factor_dict=None, genome_bg_dict=None,
#         lambda_bg_method="ctrl_max", poisson_method="mutation", region_block_mut_num_cutoff=2,
#         reads_query_mut_min_cutoff=1, reads_query_mut_max_cutoff=16, reads_total_mut_max_cutoff=20,
#         reads_other_mut_max_cutoff=16, log_verbose=3, mpmat_filter_col_idx=-1, mpmat_block_col_idx=-1
# ):
#     """
#
#     RETURN:
#             0, means everything is okay.
#     """
#     # ---------------------------------------------------------->>>>>>>>>>>>>>>>>>>>
#     # log setting
#     # ---------------------------------------------------------->>>>>>>>>>>>>>>>>>>>
#     logging.basicConfig(level=(4 - log_verbose) * 10,
#                         format='%(levelname)-5s @ %(asctime)s: %(message)s ',
#                         datefmt='%Y-%m-%d %H:%M:%S',
#                         stream=sys.stderr,
#                         filemode="w")
#
#     process_cmd_str = """run_mpmat_poisson_test step with params:
#         mpmat_filename={mpmat_filename}
#         out_poisson_filename={out_poisson_filename}
#         ctrl_bam_filename={ctrl_bam_filename}
#         treat_bam_filename={treat_bam_filename}
#         ref_genome_fa_filename={ref_genome_fa_filename}
#         select_chr_name_list={select_chr_name_list}
#         scale_factor_dict={scale_factor_dict}
#         normalize_scale_factor_dict={normalize_scale_factor_dict}
#         genome_bg_dict={genome_bg_dict}
#         lambda_bg_method={lambda_bg_method}
#         poisson_method={poisson_method}
#         region_block_mut_num_cutoff={region_block_mut_num_cutoff}
#         reads_query_mut_min_cutoff={reads_query_mut_min_cutoff}
#         reads_query_mut_max_cutoff={reads_query_mut_max_cutoff}
#         reads_total_mut_max_cutoff={reads_total_mut_max_cutoff}
#         reads_other_mut_max_cutoff={reads_other_mut_max_cutoff}
#         log_verbose={log_verbose}""".format(
#         mpmat_filename=mpmat_filename,
#         out_poisson_filename=out_poisson_filename,
#         ctrl_bam_filename=ctrl_bam_filename,
#         treat_bam_filename=treat_bam_filename,
#         ref_genome_fa_filename=ref_genome_fa_filename,
#         select_chr_name_list=select_chr_name_list,
#         scale_factor_dict=scale_factor_dict,
#         normalize_scale_factor_dict=normalize_scale_factor_dict,
#         genome_bg_dict=genome_bg_dict,
#         lambda_bg_method=lambda_bg_method,
#         poisson_method=poisson_method,
#         region_block_mut_num_cutoff=region_block_mut_num_cutoff,
#         reads_query_mut_min_cutoff=reads_query_mut_min_cutoff,
#         reads_query_mut_max_cutoff=reads_query_mut_max_cutoff,
#         reads_total_mut_max_cutoff=reads_total_mut_max_cutoff,
#         reads_other_mut_max_cutoff=reads_other_mut_max_cutoff,
#         log_verbose=log_verbose)
#
#     # ---------------------------------------------------------->>>>>>>>>>>>>>>>>>>>
#     # open files
#     # ---------------------------------------------------------->>>>>>>>>>>>>>>>>>>>
#     try:
#         chr_mpmat_file = open(mpmat_filename, "r")
#         ctrl_bam = pysam.AlignmentFile(ctrl_bam_filename, "rb")
#         treat_bam = pysam.AlignmentFile(treat_bam_filename, "rb")
#
#         if out_poisson_filename == "stdout":
#             out_file = sys.stdout
#         else:
#             out_file = open(out_poisson_filename, "w")
#
#     except:
#         raise IOError("Open files error! Please check INPUT and OUTPUT filenames!")
#
#     logging.info("Starting to run Poisson test on \n\t%s" % mpmat_filename)
#
#     logging.debug(process_cmd_str)
#
#     # ---------------------------------------------------------->>>>>>>>>>>>>>>>>>>>
#     # load FASTA genome
#     # ---------------------------------------------------------->>>>>>>>>>>>>>>>>>>>
#     ref_genome_dict = load_reference_fasta_as_dict(ref_fasta_path=ref_genome_fa_filename,
#                                                    ref_name_list=select_chr_name_list)
#
#     # get pysam Fasta obj
#     ref_genome_fa_obj = pysam.FastaFile(ref_genome_fa_filename)
#
#     # ---------------------------------------------------------->>>>>>>>>>>>>>>>>>>>
#     # iter to run test
#     # ---------------------------------------------------------->>>>>>>>>>>>>>>>>>>>
#     if mpmat_filter_col_idx == -1:
#         parse_mpmat_filter_col_idx = None
#     else:
#         parse_mpmat_filter_col_idx = mpmat_filter_col_idx
#
#     if mpmat_block_col_idx == -1:
#         parse_mpmat_block_col_idx = None
#     else:
#         parse_mpmat_block_col_idx = mpmat_block_col_idx
#
#     for line_index, line in enumerate(chr_mpmat_file):
#         line_list = line.strip().split("\t")
#         mpmat_info = mpmatLine(line_list,
#                                block_info_index=parse_mpmat_block_col_idx,
#                                filter_info_index=parse_mpmat_filter_col_idx)
#
#         # make dict
#         mpmat_info_dict = {
#             "ctrl_all_count": None,
#             "treat_all_count": None,
#             "ctrl_mut_count": None,
#             "treat_mut_count": None,
#             "ctrl_all_count.norm": None,
#             "treat_all_count.norm": None,
#             "ctrl_mut_count.norm": None,
#             "treat_mut_count.norm": None,
#             "log2FC_all_count": None,
#             "log2FC_mut_count": None,
#             "region_count_info": None,
#             "test_state": None,
#             "pvalue": None,
#             "region_site_index": None,
#             "region_site_num": None,
#             "region_block_site_num": None,
#             "region_mut_site_num": None,
#             "region_block_state": None,
#             "region_highest_site_index": None,
#             "region_highest_site_mut_num": None,
#             "region_highest_site_cover_num": None,
#             "region_highest_site_mut_ratio": None,
#         }
#
#         # chr name
#         mpmat_chr_name = mpmat_info.chr_name
#
#         # log
#         if log_verbose == 0:
#             run_report_num = 10000
#         elif log_verbose == 1:
#             run_report_num = 1000
#         else:
#             run_report_num = 100
#
#         if line_index % run_report_num == 0:
#             logging.info("Running Poisson test on %s, processed line number %s" % (mpmat_chr_name, line_index + 1))
#
#         # check mpmat region state [old-block-methods]
#         # if parse_mpmat_block_col_idx is not None:
#         #     mpmat_check_res = check_mpmat_test_state_by_block_info(
#         #         mpmat_info,
#         #         block_site_min_num_cutoff=region_block_site_min_num_cutoff,
#         #         block_num_ratio_cutoff=region_block_num_ratio_cutoff,
#         #         block_num_ratio_check_min_num=region_block_num_ratio_check_min_num,
#         #         block_num_max_num_cutoff=region_block_num_max_num_cutoff
#         #     )
#         # else:
#         #     mpmat_check_res = [True]
#
#         # ---------------------------------------------------------->>>>>>>>>>
#         # new-add-step.a hard block and find highest signal
#         # ---------------------------------------------------------->>>>>>>>>>
#         mpmat_info_block = find_block_info_and_highest_signal(mpmat_info,
#                                                               ctrl_bam_obj=ctrl_bam,
#                                                               treat_bam_obj=treat_bam,
#                                                               ref_genome_obj=ref_genome_fa_obj,
#                                                               block_mut_num_cutoff=2)
#
#         # ---------------------------------------------------------->>>>>>>>>>
#         # step1. get align count info
#         # ---------------------------------------------------------->>>>>>>>>>
#         # pipeline
#         ctrl_count_res = get_mpmat_region_count(in_bam_obj=ctrl_bam,
#                                                 mpmat_info=mpmat_info_block,
#                                                 ref_genome_dict=ref_genome_dict,
#                                                 query_mut_type=mpmat_info.mut_type,
#                                                 query_mut_min_cutoff=reads_query_mut_min_cutoff,
#                                                 query_mut_max_cutoff=reads_query_mut_max_cutoff,
#                                                 total_mut_max_cutoff=reads_total_mut_max_cutoff,
#                                                 other_mut_max_cutoff=reads_other_mut_max_cutoff)
#
#         treat_count_res = get_mpmat_region_count(in_bam_obj=treat_bam,
#                                                  mpmat_info=mpmat_info_block,
#                                                  ref_genome_dict=ref_genome_dict,
#                                                  query_mut_type=mpmat_info.mut_type,
#                                                  query_mut_min_cutoff=reads_query_mut_min_cutoff,
#                                                  query_mut_max_cutoff=reads_query_mut_max_cutoff,
#                                                  total_mut_max_cutoff=reads_total_mut_max_cutoff,
#                                                  other_mut_max_cutoff=reads_other_mut_max_cutoff)
#
#         ctrl_mpmat_count_dict = ctrl_count_res[0]
#         treat_mpmat_count_dict = treat_count_res[0]
#
#         # ---------------------------------------------------------->>>>>>>>>>
#         # step2. align count normalization
#         # ---------------------------------------------------------->>>>>>>>>>
#         # region mut lambda
#         # try:
#         ctrl_mpmat_lambda = ctrl_mpmat_count_dict["region_mut_count"] / 1.0 \
#                             / scale_factor_dict["ctrl"][mpmat_chr_name]["all_align_count.scale_factor"]
#         treat_mpmat_lambda = treat_mpmat_count_dict["region_mut_count"] / 1.0 \
#                              / scale_factor_dict["treat"][mpmat_chr_name]["all_align_count.scale_factor"]
#
#         # all reads lambda
#         ctrl_mpmat_lambda_all = ctrl_mpmat_count_dict["all_align_count"] / 1.0 \
#                                 / scale_factor_dict["ctrl"][mpmat_chr_name]["all_align_count.scale_factor"]
#         treat_mpmat_lambda_all = treat_mpmat_count_dict["all_align_count"] / 1.0 \
#                                  / scale_factor_dict["treat"][mpmat_chr_name]["all_align_count.scale_factor"]
#
#         # ---------------------------------------------------------->>>>>>>>>>
#         # step3. Poisson test
#         # ---------------------------------------------------------->>>>>>>>>>
#         # global vars
#         treat_lambda = 0
#         bg_lambda = 0
#
#         # Poisson test
#         if lambda_bg_method == "raw":
#             if poisson_method == "mutation":
#                 bg_lambda = ctrl_mpmat_count_dict["region_mut_count"]
#                 treat_lambda = treat_mpmat_count_dict["region_mut_count"]
#
#             elif poisson_method == "all":
#                 bg_lambda = ctrl_mpmat_count_dict["all_align_count"]
#                 treat_lambda = treat_mpmat_count_dict["all_align_count"]
#
#             else:
#                 logging.error("Set wrong Poisson method!")
#                 raise ValueError("Set wrong Poisson method!")
#
#         else:
#             if poisson_method == "mutation":
#                 # make mutation lambda list
#                 ctrl_bg_lambda_list = [
#                     genome_bg_dict["ctrl"]["scale_mut_bg"]["genome_bg"],
#                     genome_bg_dict["ctrl"]["scale_mut_bg"][mpmat_chr_name],
#                     ctrl_mpmat_lambda
#                 ]
#
#                 treat_bg_lambda_list = [
#                     genome_bg_dict["treat"]["scale_mut_bg"]["genome_bg"],
#                     genome_bg_dict["treat"]["scale_mut_bg"][mpmat_chr_name]
#                 ]
#
#                 bg_lambda_list = [
#                     genome_bg_dict["ctrl"]["scale_mut_bg"]["genome_bg"],
#                     genome_bg_dict["ctrl"]["scale_mut_bg"][mpmat_chr_name],
#                     ctrl_mpmat_lambda,
#                     genome_bg_dict["treat"]["scale_mut_bg"]["genome_bg"],
#                     genome_bg_dict["treat"]["scale_mut_bg"][mpmat_chr_name]
#                 ]
#
#                 # set treat lambda value
#                 treat_lambda = treat_mpmat_lambda
#
#             elif poisson_method == "all":
#                 # make mutation lambda list
#                 ctrl_bg_lambda_list = [
#                     genome_bg_dict["ctrl"]["scale_all_bg"]["genome_bg"],
#                     genome_bg_dict["ctrl"]["scale_all_bg"][mpmat_chr_name],
#                     ctrl_mpmat_lambda_all
#                 ]
#
#                 treat_bg_lambda_list = [
#                     genome_bg_dict["treat"]["scale_all_bg"]["genome_bg"],
#                     genome_bg_dict["treat"]["scale_all_bg"][mpmat_chr_name]
#                 ]
#
#                 bg_lambda_list = [
#                     genome_bg_dict["ctrl"]["scale_all_bg"]["genome_bg"],
#                     genome_bg_dict["ctrl"]["scale_all_bg"][mpmat_chr_name],
#                     ctrl_mpmat_lambda_all,
#                     genome_bg_dict["treat"]["scale_all_bg"]["genome_bg"],
#                     genome_bg_dict["treat"]["scale_all_bg"][mpmat_chr_name]
#                 ]
#
#                 # set treat lambda value
#                 treat_lambda = treat_mpmat_lambda_all
#
#             else:
#                 logging.error("Set wrong Poisson method!")
#                 raise ValueError("Set wrong Poisson method!")
#
#             # select bg lambda
#             if lambda_bg_method == "ctrl_max":
#                 bg_lambda = max(ctrl_bg_lambda_list)
#
#             elif lambda_bg_method == "treat_max":
#                 bg_lambda = max(treat_bg_lambda_list)
#
#             elif lambda_bg_method == "max":
#                 bg_lambda = max(bg_lambda_list)
#
#             elif lambda_bg_method == "raw":
#                 pass
#
#             else:
#                 logging.error("Set wrong lambda method!")
#
#         # Poisson test
#         mut_pvalue = poisson_test(treat_lambda, bg_lambda, alternative="larger", method="sqrt")[1]
#
#         # record state
#         state_test = "TestOK"
#
#         # ---------------------------------------------------------->>>>>>>>>>
#         # step4. record signal and pvalue
#         # ---------------------------------------------------------->>>>>>>>>>
#         # make count string
#         mut_num_list = []
#         ctrl_count_list = []
#         treat_count_list = []
#
#         for mut_num in range(mpmat_info.site_num + 1):
#             mut_num_list.append(mut_num)
#             ctrl_count_list.append(ctrl_count_res[2][mut_num])
#             treat_count_list.append(treat_count_res[2][mut_num])
#
#         count_str_list = [mut_num_list, ctrl_count_list, treat_count_list]
#         count_str = " ".join([",".join(map(str, x)) for x in count_str_list])
#         mpmat_info_dict["region_count_info"] = count_str
#
#         # calculate normalized signal
#         mpmat_ctrl_mut_count_norm = ctrl_mpmat_count_dict["region_mut_count"] / 1.0 / \
#                                     normalize_scale_factor_dict["ctrl"][mpmat_chr_name]["all_align_count.scale_factor"]
#         mpmat_treat_mut_count_norm = treat_mpmat_count_dict["region_mut_count"] / 1.0 / \
#                                      normalize_scale_factor_dict["treat"][mpmat_chr_name][
#                                          "all_align_count.scale_factor"]
#
#         mpmat_ctrl_all_count_norm = ctrl_mpmat_count_dict["all_align_count"] / 1.0 / \
#                                     normalize_scale_factor_dict["ctrl"][mpmat_chr_name]["all_align_count.scale_factor"]
#         mpmat_treat_all_count_norm = treat_mpmat_count_dict["all_align_count"] / 1.0 / \
#                                      normalize_scale_factor_dict["treat"][mpmat_chr_name][
#                                          "all_align_count.scale_factor"]
#
#         # All -> calculate log2 fold change fix
#
#         if mpmat_ctrl_all_count_norm != 0:
#             all_count_FC = mpmat_treat_all_count_norm / 1.0 / mpmat_ctrl_all_count_norm
#         else:
#             all_count_FC = mpmat_treat_all_count_norm / 1.0 / genome_bg_dict["ctrl"]["norm_scale_all_bg"][
#                 mpmat_chr_name]
#
#         if all_count_FC > 0:
#             log2FC_all_count = math.log(all_count_FC, 2)
#         else:
#             log2FC_all_count = "NA"
#
#         # Mut -> calculate log2 fold change fix
#         if mpmat_ctrl_mut_count_norm != 0:
#             mut_count_FC = mpmat_treat_mut_count_norm / 1.0 / mpmat_ctrl_mut_count_norm
#         else:
#             mut_count_FC = mpmat_treat_mut_count_norm / 1.0 / genome_bg_dict["ctrl"]["norm_scale_mut_bg"][
#                 mpmat_chr_name]
#
#         if mut_count_FC > 0:
#             log2FC_mut_count = math.log(mut_count_FC, 2)
#         else:
#             log2FC_mut_count = "NA"
#
#         # add all info into dict
#         mpmat_info_dict["ctrl_all_count"] = ctrl_mpmat_count_dict["all_align_count"]
#         mpmat_info_dict["treat_all_count"] = treat_mpmat_count_dict["all_align_count"]
#
#         mpmat_info_dict["ctrl_mut_count"] = ctrl_mpmat_count_dict["region_mut_count"]
#         mpmat_info_dict["treat_mut_count"] = treat_mpmat_count_dict["region_mut_count"]
#
#         mpmat_info_dict["ctrl_all_count.norm"] = mpmat_ctrl_all_count_norm
#         mpmat_info_dict["treat_all_count.norm"] = mpmat_treat_all_count_norm
#
#         mpmat_info_dict["ctrl_mut_count.norm"] = mpmat_ctrl_mut_count_norm
#         mpmat_info_dict["treat_mut_count.norm"] = mpmat_treat_mut_count_norm
#
#         mpmat_info_dict["log2FC_all_count"] = log2FC_all_count
#         mpmat_info_dict["log2FC_mut_count"] = log2FC_mut_count
#
#         # add highest info
#         mpmat_info_dict["region_site_index"] = ",".join(mpmat_info_block.site_index_list)
#         mpmat_info_dict["region_site_num"] = mpmat_info_block.site_num
#         mpmat_info_dict["region_block_site_num"] = mpmat_info_block.block_site_num
#         mpmat_info_dict["region_mut_site_num"] = mpmat_info_block.site_num - mpmat_info_block.block_site_num
#         mpmat_info_dict["region_block_state"] = mpmat_info_block.mut_key
#         mpmat_info_dict["region_highest_site_index"] = mpmat_info_block.highest_site_dict["full_site_index"]
#         mpmat_info_dict["region_highest_site_mut_num"] = mpmat_info_block.highest_site_dict["mut_num"]
#         mpmat_info_dict["region_highest_site_cover_num"] = mpmat_info_block.highest_site_dict["total"]
#         mpmat_info_dict["region_highest_site_mut_ratio"] = mpmat_info_block.highest_site_dict["mut_ratio"]
#
#         # record pvalue
#         mpmat_info_dict["test_state"] = state_test
#         mpmat_info_dict["pvalue"] = mut_pvalue
#
#         # ---------------------------------------------------------->>>>>>>>>>
#         # step5. output part
#         # ---------------------------------------------------------->>>>>>>>>>
#         info_list = line_list[:3]
#
#         # region index
#         info_list.append("%s_%s_%s" % (line_list[0], line_list[1], line_list[2]))
#
#         info_list += [
#             mpmat_info_dict["region_site_num"],
#             mpmat_info_dict["region_block_site_num"],
#             mpmat_info_dict["region_mut_site_num"],
#             mpmat_info_dict["region_site_index"],
#             mpmat_info_dict["region_block_state"],
#             mpmat_info_dict["region_highest_site_index"],
#             mpmat_info_dict["region_highest_site_mut_num"],
#             mpmat_info_dict["region_highest_site_cover_num"],
#             mpmat_info_dict["region_highest_site_mut_ratio"],
#             mpmat_info_dict["ctrl_all_count"],
#             mpmat_info_dict["treat_all_count"],
#             mpmat_info_dict["ctrl_mut_count"],
#             mpmat_info_dict["treat_mut_count"],
#             mpmat_info_dict["ctrl_all_count.norm"],
#             mpmat_info_dict["treat_all_count.norm"],
#             mpmat_info_dict["ctrl_mut_count.norm"],
#             mpmat_info_dict["treat_mut_count.norm"],
#             mpmat_info_dict["region_count_info"],
#             mpmat_info_dict["log2FC_all_count"],
#             mpmat_info_dict["log2FC_mut_count"],
#             mpmat_info_dict["test_state"],
#             mpmat_info_dict["pvalue"]
#         ]
#
#         # write into file
#         info_list_str = "\t".join(map(str, info_list))
#         out_file.write(info_list_str + "\n")
#
#     # close files
#     chr_mpmat_file.close()
#     ctrl_bam.close()
#     treat_bam.close()
#
#     if out_poisson_filename != "stdout":
#         out_file.close()
#
#     # logging
#     logging.info("Run Poisson test successfully on .mpmat file \n\t%s" % mpmat_filename)
#
#     return 0
#
#
# def multi_run_mpmat_poisson_test(mpmat_block_split_dict, ctrl_bam_filename, treat_bam_filename, ref_genome_fa_filename,
#                                  scale_factor_dict=None, normalize_scale_factor_dict=None, genome_bg_dict=None,
#                                  lambda_bg_method="ctrl_max", poisson_method="mutation", log_verbose=3, thread=1,
#                                  **region_filter_args):
#     """
#     INPUT:
#
#         <**region_filter_args>
#             args list:
#                 <region_block_mut_num_cutoff>
#                 <reads_query_mut_min_cutoff>
#                 <reads_query_mut_max_cutoff>
#                 <reads_total_mut_max_cutoff>
#                 <reads_other_mut_max_cutoff>
#                 <mpmat_filter_info_col_index>
#                 <mpmat_block_info_col_index>
#
#     RETURN
#         <out_poisson_filename_dict>
#             dict, contains output mpmat files with block information, and format like
#
#             {
#                 'chr1': 'chr1.eqdA2zJkHCm0YW1o.AddBlockInfo',
#                 'chr19': 'chr19.ZIrU6xC2DFcBayE1.AddBlockInfo',
#                 'chr20': 'chr20.nK2goEB6xh9MzXpD.AddBlockInfo',
#                 'chr_name_order': ['chr1', 'chr19', 'chr20']
#             }
#
#     """
#     # ------------------------------------------------------------>>>>>>>>>>
#     # log setting
#     # ------------------------------------------------------------>>>>>>>>>>
#     logging.basicConfig(level=(4 - log_verbose) * 10,
#                         format='%(levelname)-5s @ %(asctime)s: %(message)s ',
#                         datefmt='%Y-%m-%d %H:%M:%S',
#                         stream=sys.stderr,
#                         filemode="w")
#
#     # ------------------------------------------------------------>>>>>>>>>>
#     # params setting
#     # ------------------------------------------------------------>>>>>>>>>>
#     run_params_dict = {
#         "region_block_mut_num_cutoff": 2,
#         "reads_query_mut_min_cutoff": 1,
#         "reads_query_mut_max_cutoff": 16,
#         "reads_total_mut_max_cutoff": 20,
#         "reads_other_mut_max_cutoff": 16,
#         "mpmat_filter_info_col_index": -1,
#         "mpmat_block_info_col_index": -1
#     }
#
#     for args_key in region_filter_args:
#         if region_filter_args.get(args_key) is not None:
#             run_params_dict[args_key] = region_filter_args[args_key]
#
#     # ------------------------------------------------------------>>>>>>>>>>
#     # check params
#     # ------------------------------------------------------------>>>>>>>>>>
#     if scale_factor_dict is None:
#         raise IOError("FUN <multi_run_mpmat_poisson_test> need set 'scale_factor_dict' !")
#
#     if normalize_scale_factor_dict is None:
#         raise IOError("FUN <multi_run_mpmat_poisson_test> need set 'normalize_scale_factor_dict' !")
#
#     if genome_bg_dict is None:
#         raise IOError("FUN <multi_run_mpmat_poisson_test> need set 'genome_bg_dict' !")
#
#     # ------------------------------------------------------------>>>>>>>>>>
#     # get chr name order
#     # ------------------------------------------------------------>>>>>>>>>>
#     chr_name_order_list = mpmat_block_split_dict["chr_name_order"]
#     out_poisson_filename_dict = {"chr_name_order": chr_name_order_list}
#
#     # ------------------------------------------------------------>>>>>>>>>>
#     # run part
#     # ------------------------------------------------------------>>>>>>>>>>
#     logging.info("-" * 80)
#     logging.info("Starting to run Poisson test...")
#
#     pool = multiprocessing.Pool(processes=thread)
#
#     run_return_info_list = []
#
#     for chr_name in chr_name_order_list:
#         # mpmat input and output
#         chr_mpmat_old_filename = mpmat_block_split_dict[chr_name]
#         chr_poisson_new_filename = chr_mpmat_old_filename + "." + "PoissonResult"
#         out_poisson_filename_dict[chr_name] = chr_poisson_new_filename
#
#         run_return_info_list.append(
#             pool.apply_async(
#                 func=run_mpmat_poisson_test,
#                 args=(
#                     chr_mpmat_old_filename,
#                     chr_poisson_new_filename,
#                     ctrl_bam_filename,
#                     treat_bam_filename,
#                     ref_genome_fa_filename,
#                     [chr_name],
#                     scale_factor_dict,
#                     normalize_scale_factor_dict,
#                     genome_bg_dict,
#                     lambda_bg_method,
#                     poisson_method,
#                     run_params_dict["region_block_mut_num_cutoff"],
#                     run_params_dict["reads_query_mut_min_cutoff"],
#                     run_params_dict["reads_query_mut_max_cutoff"],
#                     run_params_dict["reads_total_mut_max_cutoff"],
#                     run_params_dict["reads_other_mut_max_cutoff"],
#                     log_verbose,
#                     run_params_dict["mpmat_filter_info_col_index"],
#                     run_params_dict["mpmat_block_info_col_index"],
#                 )
#             )
#         )
#
#     pool.close()
#     pool.join()
#
#     # check run state
#     final_run_state = 0
#     for index, res in enumerate(run_return_info_list):
#         run_state = res.get()
#         if run_state != 0:
#             logging.error("Poisson test error occur with %s!" % chr_name_order_list[index])
#             if final_run_state == 0:
#                 final_run_state = 1
#
#     if final_run_state == 0:
#         logging.info("Calculation of Poisson test result. Done!")
#         return out_poisson_filename_dict
#
#     else:
#         logging.error("Something wrong with Poisson test step!")
#         raise RuntimeError()
#
#
# def make_qvalue_with_BH_method(pval_list):
#     """
#     INPUT:
#         <pval_list>
#             list, may contain 'NA' value
#
#     RETURN
#         <fdr_list>
#             list, return FDR list with BH method.
#
#     """
#
#     # init vars
#     raw_pval_index_dict = {}
#     rm_NA_pval_list = []
#     run_index = 0
#
#     for index, pval in enumerate(pval_list):
#         if pval != "NA":
#             rm_NA_pval_list.append(pval)
#             raw_pval_index_dict[run_index] = index
#             run_index += 1
#
#     FDR_qvalue = multipletests(np.array(rm_NA_pval_list), alpha=0.05, method="fdr_bh", is_sorted=False)
#     FDR_qvalue_vec = FDR_qvalue[1]
#
#     return_fdr_list = ["NA"] * len(pval_list)
#
#     for index, fdr in enumerate(FDR_qvalue_vec):
#         raw_index = raw_pval_index_dict[index]
#         return_fdr_list[raw_index] = fdr
#
#     return return_fdr_list
#
#
# def merge_split_files(split_file_dict, key_order_list=None, out_filename="stdout", header_list=None, in_sep="\t",
#                       out_sep="\t", log_verbose=3, return_col_index=None):
#     """
#     INPUT:
#         <split_file_dict>
#             dict, typical format like
#             {
#                 'chr1': 'chr1.eqdA2zJkHCm0YW1o.AddBlockInfo',
#                 'chr19': 'chr19.ZIrU6xC2DFcBayE1.AddBlockInfo',
#                 'chr20': 'chr20.nK2goEB6xh9MzXpD.AddBlockInfo',
#                 'chr_name_order': ['chr1', 'chr19', 'chr20']
#             }
#
#         <key_order_list>
#             list, if set 'None' as default, will merge all split files and make output
#
#
#         <out_filename>
#             str
#
#         <header_list>
#             list, length have to match with split files
#
#         <sep>
#             str, default '\t', can set by custom
#
#     RETURN:
#         0, everything is okay.
#
#     INFO:
#         Final date 2020-11-10
#     """
#     # log
#     logging.basicConfig(level=(4 - log_verbose) * 10,
#                         format='%(levelname)-5s @ %(asctime)s: %(message)s ',
#                         datefmt='%Y-%m-%d %H:%M:%S',
#                         stream=sys.stderr,
#                         filemode="w")
#
#     # return col
#     return_col_list = []
#
#     # open output file
#     if out_filename == "stdout":
#         out_file = sys.stdout
#     else:
#         out_file = open(out_filename, "w")
#
#     # init vars
#     ignore_key_list = ["chr_name_order"]
#
#     # check keys
#     if key_order_list is None:
#         key_order_list = split_file_dict.keys()
#
#     # header output
#     if header_list is not None:
#         out_file.write(out_sep.join(map(str, header_list)) + "\n")
#
#     # merge files
#     logging.info("Starting to merge files...")
#
#     for run_key in key_order_list:
#         if run_key not in ignore_key_list:
#             logging.debug("Merging files, processing on \n\t%s" % split_file_dict[run_key])
#
#             with open(split_file_dict[run_key], "r") as run_file:
#                 for line in run_file:
#                     if in_sep != out_sep:
#                         line_list = line.strip().split(in_sep)
#                         out_file.write(out_sep.join(line_list) + "\n")
#
#                         if return_col_index is not None:
#                             return_col_list.append(line_list[return_col_index])
#                     else:
#                         out_file.write(line)
#
#                         if return_col_index is not None:
#                             line_list = line.strip().split(in_sep)
#                             return_col_list.append(line_list[return_col_index])
#
#     # log
#     logging.info("Merging files, done!")
#
#     # close file
#     if out_filename != "stdout":
#         out_file.close()
#
#     return 0, return_col_list
#
#
# def clear_temp_files_by_dict(temp_file_dict, log_verbose=3):
#     """
#     INPUT:
#         <temp_file_dict>
#             dict, format like:
#                 {
#                     'chr1': 'chr1.eqdA2zJkHCm0YW1o.AddBlockInfo',
#                     'chr19': 'chr19.ZIrU6xC2DFcBayE1.AddBlockInfo',
#                     'chr20': 'chr20.nK2goEB6xh9MzXpD.AddBlockInfo',
#                     'chr_name_order': ['chr1', 'chr19', 'chr20']
#                 }
#
#     RETURN:
#         0, everything is okay.
#
#     INFO:
#         Final date 2020-11-10
#     """
#
#     # log
#     logging.basicConfig(level=(4 - log_verbose) * 10,
#                         format='%(levelname)-5s @ %(asctime)s: %(message)s ',
#                         datefmt='%Y-%m-%d %H:%M:%S',
#                         stream=sys.stderr,
#                         filemode="w")
#
#     run_state = 0
#
#     for run_key in temp_file_dict.keys():
#         if type(temp_file_dict[run_key]) is str:
#             if os.path.exists(temp_file_dict[run_key]):
#                 if os.path.isfile(temp_file_dict[run_key]):
#                     try:
#                         os.remove(temp_file_dict[run_key])
#                     except:
#                         logging.error("Removing error on \n\t%s" % temp_file_dict[run_key])
#                         run_state = 1
#
#     return run_state
#
#
# def cmp_site_index(site_index_a, site_index_b, ref_order_dict):
#     """
#     INPUT:
#         <site_index_a>, <site_index_b>
#             str, like chr1_629627_CT
#
#         <ref_order>
#             dict, format like:
#
#             ref_order_dict = {
#                 'chr1': 0,
#                 'chr19': 1,
#                 'chr20': 2
#             }
#
#     RETURN:
#         0, site_index_a == site_index_b at position level, ignore mutation info;
#         -1, site_a at upstream of site_b;
#         1, site_a at downstream of site_b
#     """
#
#     site_A = siteIndex(site_index_a)
#     site_B = siteIndex(site_index_b)
#
#     if site_A.chr_name == site_B.chr_name:
#         if site_A.site_pos == site_B.site_pos:
#             return 0
#
#         elif site_A.site_pos < site_B.site_pos:
#             return -1
#
#         elif site_A.site_pos > site_B.site_pos:
#             return 1
#
#     else:
#         site_A_chr_order = ref_order_dict.get(site_A.chr_name)
#         site_B_chr_order = ref_order_dict.get(site_B.chr_name)
#
#         if (site_A_chr_order is not None) and (site_B_chr_order is not None):
#             if site_A_chr_order < site_B_chr_order:
#                 return -1
#
#             elif site_A_chr_order > site_B_chr_order:
#                 return 1
#
#         else:
#             raise TypeError("Site index not in your reference!")
#
#
# def query_region_bmat_info(bmat_file, site_index_list, genome_order_dict):
#     """
#     INPUT:
#         <bmat_file>
#             file.obj, bmat file handle
#
#         <site_index_list>
#             list, like [chr1_20452_CT, chr1_20467_C., chr1_20474_CT]
#
#         <genome_order_dict>
#             dict, for FUN <cmp_site_index>
#
#     RETURN
#         <site_dict>
#             dict, key is site_index, value bmat line list
#
#     """
#
#     # init
#     site_index = siteIndex(site_index_list[0])
#     bmat_line = bmat_file.readline()
#
#     # define dict
#     query_total_num = len(site_index_list)
#     query_site_num = 0
#     query_res_dict = {
#         "site_index_list": []
#     }
#
#     # make init val
#     for raw_site_index in site_index_list:
#         raw_site_index_obj = siteIndex(raw_site_index)
#         query_res_dict["site_index_list"].append(raw_site_index_obj.site_index)
#         query_res_dict[raw_site_index_obj.site_index] = None
#
#     # query run
#     while bmat_line != "":
#         bmat_line_list = bmat_line.strip().split("\t")
#         bmat_site_index = "_".join(bmat_line_list[0:3])
#
#         cmp_res = cmp_site_index(site_index.site_index_raw, bmat_site_index, genome_order_dict)
#
#         if cmp_res == 0:
#             query_site_num += 1
#
#             # add info into dict
#             # query_res_dict["site_index_list"].append(site_index.site_index)
#             query_res_dict[site_index.site_index] = bmatLine(bmat_line_list)
#
#             if query_site_num >= query_total_num:
#                 break
#
#             else:
#                 # read new
#                 site_index = siteIndex(site_index_list[query_site_num])
#                 bmat_line = bmat_file.readline()
#
#         elif cmp_res == 1:
#             bmat_line = bmat_file.readline()
#
#         elif cmp_res == -1:
#             query_site_num += 1
#
#             # add info into dict
#             # query_res_dict["site_index_list"].append(site_index.site_index)
#             query_res_dict[site_index.site_index] = None
#
#             if query_site_num >= query_total_num:
#                 break
#
#             else:
#                 # read new
#                 site_index = siteIndex(site_index_list[query_site_num])
#                 bmat_line = bmat_file.readline()
#
#     return query_res_dict
#
#
# def site_binomial_test(bmat_line_obj, mut_type, background_pval, return_int=True):
#     """
#     INPUT:
#         <bmat_line_obj>
#             obj, bmatLine
#
#         <mut_type>
#             str, "CT" mean ref_base is "C" and mut_base is "T"
#
#         <background_pval>
#             float, binomial test background pvalue
#
#         <return_int>
#             bool, set True will return a close int after -10 * log(pvalue), OR return raw pvalue
#
#     RETURN
#         <pval>
#             int|float, binomial test pvalue
#     """
#
#     # count base info
#     ref_base_count = bmat_line_obj.count_dict.get(mut_type[0])
#     mut_base_count = bmat_line_obj.count_dict.get(mut_type[1])
#     total_base_count = ref_base_count + mut_base_count
#
#     try:
#         binom_pval = stats.binom_test(x=mut_base_count,
#                                       n=total_base_count,
#                                       p=background_pval,
#                                       alternative='greater')
#     except:
#         return None
#
#     if return_int:
#         if binom_pval > 1e-127:
#             binom_pval_int = int(round(-10 * math.log(binom_pval, 10), 0))
#         else:
#             binom_pval_int = 127
#
#         return binom_pval_int
#
#     else:
#         return binom_pval
#
#
# def site_hard_filter(site_bmat_obj, mut_type, mut_count_cutoff=5, other_mut_count_cutoff=5, other_mut_ratio_cutoff=0.25,
#                      cover_min_ratio_check_cutoff=6, cover_up_limit_cutoff=500):
#     """
#     INPUT:
#         <site_bmat_obj>
#             obj, bmatLine
#
#         <mut_type>
#             str, "CT" mean ref_base is "C" and mut_base is "T"
#
#         <mut_count_cutoff>
#             int, site larger than this will be blocked
#
#         <other_mut_count_cutoff>
#             int, site other type mutation larger than this will be blocked
#
#         <other_mut_ratio_cutoff>
#             float, site other type mutation ratio larger than this will be blocked
#
#         <cover_min_ratio_check_cutoff>
#             int, mutation count have to larger than this can process with <other_mut_ratio_cutoff>
#
#         <cover_up_limit_cutoff>
#             int, site cover larger than this cutoff will be blocked.
#
#     RETURN:
#         <block_state, block_reason>
#             tuple, contain two element
#
#             <block_state>
#                 True, site should be blocked
#
#                 False, site should not be blocked
#
#             <block_reason>
#                 MC, too much query mutation count;
#                 MR, too high query mutation ratio;
#                 OC, too much other mutation count;
#                 OR, too high other mutation ratio;
#                 TC, too much cover count
#                 None, when <block_state> is False, <block_reason> is None
#     """
#
#     # get number
#     ref_base_count = site_bmat_obj.count_dict[mut_type[0]]
#     mut_base_count = site_bmat_obj.count_dict[mut_type[1]]
#     other_mut_base_count = 0
#     total_cover_count = 0
#
#     # check ref base count
#     if ref_base_count < cover_min_ratio_check_cutoff:
#         return False, None
#
#     # check ref upper cutoff
#     if ref_base_count >= cover_up_limit_cutoff:
#         return True, "TC"
#
#     # check mut count
#     if mut_base_count >= mut_count_cutoff:
#         return True, "MC"
#
#     # count ratio
#     for base in site_bmat_obj.count_dict:
#         total_cover_count += site_bmat_obj.count_dict[base]
#         if base not in mut_type:
#             other_mut_base_count += site_bmat_obj.count_dict[base]
#
#     other_mut_ratio = other_mut_base_count / 1.0 / total_cover_count
#
#     # check other mut type
#     if other_mut_base_count >= other_mut_count_cutoff:
#         return True, "OC"
#
#     if other_mut_ratio >= other_mut_ratio_cutoff:
#         return True, "OR"
#
#     # final return
#     return False, None
#
#
# def add_block_info_to_mpmat(chr_mpmat_filename, chr_bmat_ctrl_filename, chr_bmat_treat_filename, query_mutation_type,
#                             chr_ctrl_mut_bg_pval, ref_order_dict, out_chr_mpmat_filename="stdout",
#                             hf_ctrl_mut_count_cutoff=3, hf_ctrl_other_mut_count_cutoff=5,
#                             hf_ctrl_other_mut_ratio_cutoff=0.25, hf_ctrl_cover_min_ratio_check_cutoff=6,
#                             hf_ctrl_cover_up_limit_cutoff=500, ctrl_binomial_cutoff_int=30,
#                             hf_treat_mut_count_cutoff=5000, hf_treat_other_mut_count_cutoff=50,
#                             hf_treat_other_mut_ratio_cutoff=0.6, hf_treat_cover_min_ratio_check_cutoff=10,
#                             hf_treat_cover_up_limit_cutoff=5000, treat_site_mut_min_cutoff=1, log_verbose=3):
#     """
#     HELP:
#
#     """
#
#     # ---------------------------------------------------------->>>>>>>>>>>>>>>>>>>>
#     # log setting
#     # ---------------------------------------------------------->>>>>>>>>>>>>>>>>>>>
#     logging.basicConfig(level=10,
#                         format='%(levelname)-5s @ %(asctime)s: %(message)s ',
#                         datefmt='%Y-%m-%d %H:%M:%S',
#                         stream=sys.stderr,
#                         filemode="w")
#
#     # ---------------------------------------------------------->>>>>>>>>>>>>>>>>>>>
#     # part I open file
#     # ---------------------------------------------------------->>>>>>>>>>>>>>>>>>>>
#     try:
#         chr_mpmat_file = open(chr_mpmat_filename, "r")
#
#         # open ctrl bmat file
#         if (chr_bmat_ctrl_filename[-3:] == ".gz") or (chr_bmat_ctrl_filename[-5:] == ".gzip"):
#             chr_bmat_ctrl_file = gzip.open(chr_bmat_ctrl_filename, "r")
#         else:
#             chr_bmat_ctrl_file = open(chr_bmat_ctrl_filename, "r")
#
#         # open treat bmat file
#         if (chr_bmat_treat_filename[-3:] == ".gz") or (chr_bmat_treat_filename[-5:] == ".gzip"):
#             chr_bmat_treat_file = gzip.open(chr_bmat_treat_filename, "r")
#         else:
#             chr_bmat_treat_file = open(chr_bmat_treat_filename, "r")
#
#         if out_chr_mpmat_filename == "stdout":
#             out_chr_mpmat_file = sys.stdout
#         else:
#             out_chr_mpmat_file = open(out_chr_mpmat_filename, "w")
#     except:
#         raise IOError("Load or Open files error!")
#
#     logging.info("Start to analysis sites block info. \n\t%s" % chr_mpmat_filename)
#     # ---------------------------------------------------------->>>>>>>>>>>>>>>>>>>>
#     # part II iteration and add block info
#     # ---------------------------------------------------------->>>>>>>>>>>>>>>>>>>>
#     for line_index, line in enumerate(chr_mpmat_file):
#
#         line_list = line.strip().split("\t")
#         mpmat_info = mpmatLine(line_list)
#
#         # log
#         if log_verbose == 0:
#             run_report_num = 1000
#         elif log_verbose == 1:
#             run_report_num = 100
#         elif log_verbose == 2:
#             run_report_num = 50
#         else:
#             run_report_num = 20
#
#         if line_index % run_report_num == 0:
#             logging.info(
#                 "Running mpmat block step on %s \n\tProcessed line number %s" % (chr_mpmat_filename, line_index + 1))
#
#         # check mut type
#         if mpmat_info.mut_type not in query_mutation_type:
#             raise ValueError("mpmatLine.mut_type doesn't match query mutation type!")
#
#         # ---------------------------------------------------------->>>>>>>>>>
#         # query site bmat info
#         # ---------------------------------------------------------->>>>>>>>>>
#         query_site_idx_list = mpmat_info.site_index_list
#
#         ctrl_bmat_res_dict = query_region_bmat_info(bmat_file=chr_bmat_ctrl_file,
#                                                     site_index_list=query_site_idx_list,
#                                                     genome_order_dict=ref_order_dict)
#
#         treat_bmat_res_dict = query_region_bmat_info(bmat_file=chr_bmat_treat_file,
#                                                      site_index_list=query_site_idx_list,
#                                                      genome_order_dict=ref_order_dict)
#
#         # ---------------------------------------------------------->>>>>>>>>>
#         # ctrl Hard filter
#         # ---------------------------------------------------------->>>>>>>>>>
#         ctrl_site_hard_filter_state = []
#         ctrl_site_hard_filter_reason = []
#
#         for site_idx in ctrl_bmat_res_dict["site_index_list"]:
#             if ctrl_bmat_res_dict[site_idx] is None:
#                 ctrl_site_hard_filter_state.append(False)
#
#                 # NS means 'Non sequencing data'
#                 ctrl_site_hard_filter_reason.append("NS")
#
#             else:
#                 hard_filter_res = site_hard_filter(site_bmat_obj=ctrl_bmat_res_dict[site_idx],
#                                                    mut_type=mpmat_info.mut_type,
#                                                    mut_count_cutoff=hf_ctrl_mut_count_cutoff,
#                                                    other_mut_count_cutoff=hf_ctrl_other_mut_count_cutoff,
#                                                    other_mut_ratio_cutoff=hf_ctrl_other_mut_ratio_cutoff,
#                                                    cover_min_ratio_check_cutoff=hf_ctrl_cover_min_ratio_check_cutoff,
#                                                    cover_up_limit_cutoff=hf_ctrl_cover_up_limit_cutoff)
#
#                 ctrl_site_hard_filter_state.append(hard_filter_res[0])
#                 ctrl_site_hard_filter_reason.append(hard_filter_res[1])
#
#         # ---------------------------------------------------------->>>>>>>>>>
#         # ctrl Binomial test
#         # ---------------------------------------------------------->>>>>>>>>>
#         binom_test_pval = []
#         binom_test_res = 0
#
#         for index, site_idx in enumerate(ctrl_bmat_res_dict["site_index_list"]):
#             if not ctrl_site_hard_filter_state[index]:
#                 if ctrl_site_hard_filter_reason[index] != "NS":
#                     binom_test_res = site_binomial_test(bmat_line_obj=ctrl_bmat_res_dict[site_idx],
#                                                         mut_type=mpmat_info.mut_type,
#                                                         background_pval=chr_ctrl_mut_bg_pval,
#                                                         return_int=True)
#
#             binom_test_pval.append(binom_test_res)
#
#         # ---------------------------------------------------------->>>>>>>>>>
#         # treat Hard filter
#         # ---------------------------------------------------------->>>>>>>>>>
#         treat_site_hard_filter_state = []
#         treat_site_hard_filter_reason = []
#
#         for site_idx in treat_bmat_res_dict["site_index_list"]:
#             if treat_bmat_res_dict[site_idx] is None:
#                 treat_site_hard_filter_state.append(True)
#
#                 # NS means 'Non sequencing data'
#                 treat_site_hard_filter_reason.append("NS")
#
#             else:
#                 hard_filter_res = site_hard_filter(site_bmat_obj=treat_bmat_res_dict[site_idx],
#                                                    mut_type=mpmat_info.mut_type,
#                                                    mut_count_cutoff=hf_treat_mut_count_cutoff,
#                                                    other_mut_count_cutoff=hf_treat_other_mut_count_cutoff,
#                                                    other_mut_ratio_cutoff=hf_treat_other_mut_ratio_cutoff,
#                                                    cover_min_ratio_check_cutoff=hf_treat_cover_min_ratio_check_cutoff,
#                                                    cover_up_limit_cutoff=hf_treat_cover_up_limit_cutoff)
#
#                 treat_site_hard_filter_state.append(hard_filter_res[0])
#                 treat_site_hard_filter_reason.append(hard_filter_res[1])
#
#         # ---------------------------------------------------------->>>>>>>>>>
#         # make block state
#         # ---------------------------------------------------------->>>>>>>>>>
#         block_state_list = [False] * mpmat_info.site_num
#         block_reason_list = ["None"] * mpmat_info.site_num
#
#         for index in range(mpmat_info.site_num):
#             try:
#                 if ctrl_site_hard_filter_state[index]:
#                     block_state_list[index] = True
#                     # block reason: ctrl hard filter
#                     block_reason_list[index] = "CHF"
#                     continue
#             except:
#                 print(ctrl_site_hard_filter_state)
#                 print(ctrl_bmat_res_dict)
#
#             if binom_test_pval[index] > ctrl_binomial_cutoff_int:
#                 block_state_list[index] = True
#                 # block reason: ctrl binomial test
#                 block_reason_list[index] = "CBT"
#                 continue
#
#             # block by treat info
#             site_idx = treat_bmat_res_dict["site_index_list"][index]
#
#             if treat_bmat_res_dict[site_idx] is None:
#                 block_state_list[index] = True
#                 # block reason: treat non cover
#                 block_reason_list[index] = "TNC"
#                 continue
#
#             else:
#                 site_mut_count = treat_bmat_res_dict[site_idx].count_dict[mpmat_info.mut_type[-1]]
#                 if site_mut_count < treat_site_mut_min_cutoff:
#                     block_state_list[index] = True
#                     # block reason: treat non mut
#                     block_reason_list[index] = "TNM"
#                     continue
#
#             if treat_site_hard_filter_state[index]:
#                 block_state_list[index] = True
#                 # block reason: treat hard filter
#                 block_reason_list[index] = "THF"
#                 continue
#
#         # ---------------------------------------------------------->>>>>>>>>>
#         # out block info
#         # ---------------------------------------------------------->>>>>>>>>>
#         out_block_info_list = [
#             "CHF,CHR,CBT,THF,THR,FBS,FBR"
#         ]
#
#         info_list = [
#             ctrl_site_hard_filter_state,
#             ctrl_site_hard_filter_reason,
#             binom_test_pval,
#             treat_site_hard_filter_state,
#             treat_site_hard_filter_reason,
#             block_state_list,
#             block_reason_list
#         ]
#
#         for info in info_list:
#             out_block_info_list.append(",".join(map(str, info)))
#
#         # out to file
#         line_list.append(" ".join(out_block_info_list))
#         out_chr_mpmat_file.write("\t".join(line_list) + "\n")
#
#     # ---------------------------------------------------------->>>>>>>>>>>>>>>>>>>>
#     # close file
#     # ---------------------------------------------------------->>>>>>>>>>>>>>>>>>>>
#     chr_mpmat_file.close()
#     chr_bmat_ctrl_file.close()
#     chr_bmat_treat_file.close()
#
#     if out_chr_mpmat_filename != "stdout":
#         out_chr_mpmat_file.close()
#
#     logging.info("Done! \n\t%s" % chr_mpmat_filename)
#
#     return 1
#
#
# def multi_mpmat_block_site(mpmat_split_dict, ctrl_bmat_split_dict, treat_bmat_split_dict, query_mutation_type,
#                            ctrl_mut_bg_pval_dict, ref_order_dict, thread=1, log_verbose=3, **filter_args):
#     """
#     INPUT:
#
#         <**filter_args>
#             args list:
#                 <hf_ctrl_mut_count_cutoff>
#                 <hf_ctrl_other_mut_count_cutoff>
#                 <hf_ctrl_other_mut_ratio_cutoff>
#                 <hf_ctrl_cover_min_ratio_check_cutoff>
#                 <hf_ctrl_cover_up_limit_cutoff>
#                 <ctrl_binomial_cutoff_int>
#                 <hf_treat_mut_count_cutoff>
#                 <hf_treat_other_mut_count_cutoff>
#                 <hf_treat_other_mut_ratio_cutoff>
#                 <hf_treat_cover_min_ratio_check_cutoff>
#                 <hf_treat_cover_up_limit_cutoff>
#                 <treat_site_mut_min_cutoff>
#
#     RETURN
#         <out_mpmat_filename_dict>
#             dict, contains output mpmat files with block information, and format like
#
#             {
#                 'chr1': 'chr1.eqdA2zJkHCm0YW1o.AddBlockInfo',
#                 'chr19': 'chr19.ZIrU6xC2DFcBayE1.AddBlockInfo',
#                 'chr20': 'chr20.nK2goEB6xh9MzXpD.AddBlockInfo',
#                 'chr_name_order': ['chr1', 'chr19', 'chr20']
#             }
#
#     """
#     # ------------------------------------------------------------>>>>>>>>>>
#     # log setting
#     # ------------------------------------------------------------>>>>>>>>>>
#     logging.basicConfig(level=(4 - log_verbose) * 10,
#                         format='%(levelname)-5s @ %(asctime)s: %(message)s ',
#                         datefmt='%Y-%m-%d %H:%M:%S',
#                         stream=sys.stderr,
#                         filemode="w")
#
#     # ------------------------------------------------------------>>>>>>>>>>
#     # params setting
#     # ------------------------------------------------------------>>>>>>>>>>
#     run_params_dict = {
#         "hf_ctrl_mut_count_cutoff": 3,
#         "hf_ctrl_other_mut_count_cutoff": 5,
#         "hf_ctrl_other_mut_ratio_cutoff": 0.25,
#         "hf_ctrl_cover_min_ratio_check_cutoff": 6,
#         "hf_ctrl_cover_up_limit_cutoff": 500,
#         "ctrl_binomial_cutoff_int": 30,
#         "hf_treat_mut_count_cutoff": 5000,
#         "hf_treat_other_mut_count_cutoff": 50,
#         "hf_treat_other_mut_ratio_cutoff": 0.6,
#         "hf_treat_cover_min_ratio_check_cutoff": 10,
#         "hf_treat_cover_up_limit_cutoff": 5000,
#         "treat_site_mut_min_cutoff": 1
#     }
#
#     for args_key in filter_args:
#         if filter_args.get(args_key) is not None:
#             run_params_dict[args_key] = filter_args[args_key]
#
#     # ------------------------------------------------------------>>>>>>>>>>
#     # var initial
#     # ------------------------------------------------------------>>>>>>>>>>
#     # check mpmat chr_name and bmat chr_name
#     chr_name_order_list = []
#     for mpmat_chr_name in mpmat_split_dict["chr_name_order"]:
#         if ref_order_dict.get(mpmat_chr_name) is not None:
#             chr_name_order_list.append(mpmat_chr_name)
#
#     out_mpmat_filename_dict = {"chr_name_order": chr_name_order_list}
#     # ------------------------------------------------------------>>>>>>>>>>
#     # run part
#     # ------------------------------------------------------------>>>>>>>>>>
#     sys.stderr.write("-" * 80 + "\n")
#     logging.info("Start to calculate genome block info...")
#
#     pool = multiprocessing.Pool(processes=thread)
#
#     run_return_info_list = []
#
#     for chr_name in chr_name_order_list:
#         # mpmat input and output
#         chr_mpmat_old_filename = mpmat_split_dict[chr_name]
#         chr_mpmat_new_filename = chr_mpmat_old_filename + "." + "AddBlockInfo"
#         out_mpmat_filename_dict[chr_name] = chr_mpmat_new_filename
#
#         # bmat
#         ctrl_bmat_filename = ctrl_bmat_split_dict[chr_name]
#         treat_bmat_filename = treat_bmat_split_dict[chr_name]
#
#         run_return_info_list.append(
#             pool.apply_async(
#                 func=add_block_info_to_mpmat,
#                 args=(
#                     chr_mpmat_old_filename,
#                     ctrl_bmat_filename,
#                     treat_bmat_filename,
#                     query_mutation_type,
#                     ctrl_mut_bg_pval_dict[chr_name]["query_mut_bg_pval"],
#                     ref_order_dict,
#                     chr_mpmat_new_filename,
#                     run_params_dict["hf_ctrl_mut_count_cutoff"],
#                     run_params_dict["hf_ctrl_other_mut_count_cutoff"],
#                     run_params_dict["hf_ctrl_other_mut_ratio_cutoff"],
#                     run_params_dict["hf_ctrl_cover_min_ratio_check_cutoff"],
#                     run_params_dict["hf_ctrl_cover_up_limit_cutoff"],
#                     run_params_dict["ctrl_binomial_cutoff_int"],
#                     run_params_dict["hf_treat_mut_count_cutoff"],
#                     run_params_dict["hf_treat_other_mut_count_cutoff"],
#                     run_params_dict["hf_treat_other_mut_ratio_cutoff"],
#                     run_params_dict["hf_treat_cover_min_ratio_check_cutoff"],
#                     run_params_dict["hf_treat_cover_up_limit_cutoff"],
#                     run_params_dict["treat_site_mut_min_cutoff"],
#                 )
#             )
#         )
#
#     pool.close()
#     pool.join()
#
#     # check run state
#     final_run_state = 0
#     for index, res in enumerate(run_return_info_list):
#         run_state = res.get()
#         if run_state != 1:
#             logging.error("Blocking error occur with %s!" % chr_name_order_list[index])
#             if final_run_state == 0:
#                 final_run_state = 1
#
#     if final_run_state == 0:
#         logging.info("Calculation of genome block info. Done!")
#         return out_mpmat_filename_dict
#
#     elif final_run_state == 1:
#         logging.error("Something wrong with block step!")
#         raise RuntimeError("")
#
#
# def query_site_mpileup_info(site_idx_list, bam_obj, genome_obj, ignore_overlaps=True, min_base_quality=20,
#                             min_mapping_quality=20, dist_cutoff=10000):
#     """
#     INPUT
#         <site_idx_list>
#             list, format like:
#
#                 ["chr1_125178576_CT", "chr1_125178578_CT", "chr1_125178580_CT", "chr1_125178588_CA"]
#
#             The site coordinate is related to the UCSC genome browser index.
#
#             Site index have to be sorted.
#
#         <bam_obj>
#             pysam.AlignmentFile() obj
#
#         <genome_obj>
#             genome obj, which is created by pysam.FastaFile()
#
#     RETURN
#         dict
#             key:
#                 site_index only key chr_name and chr_pos like "chr1_125178576" rather than "chr1_125178576_CT"
#
#             value:
#                 A, T, C, G, N, other count
#     """
#
#     # -------------------------------------------------------->>>>>>>>
#     # check site index list
#     # -------------------------------------------------------->>>>>>>>
#     if len(site_idx_list) == 0:
#         return None
#
#     chr_name = site_idx_list[0].split("_")[0]
#     region_start = region_end = int(site_idx_list[0].split("_")[1])
#     site_dist = 0
#
#     for site_idx in site_idx_list:
#         if site_idx.split("_")[0] != chr_name:
#             raise IOError("<site_idx_list> all site index have to on the same chromosome.")
#
#     if len(site_idx_list) >= 2:
#         region_start = int(site_idx_list[0].split("_")[1])
#         region_end = int(site_idx_list[-1].split("_")[1])
#         site_dist = region_end - region_start
#
#     if site_dist >= dist_cutoff:
#         raise IOError("<site_idx_list> mpileup region is too large, which could take a lot of time!")
#
#     # -------------------------------------------------------->>>>>>>>
#     # make raw dict
#     # -------------------------------------------------------->>>>>>>>
#     site_dict = {}
#     for site_idx in site_idx_list:
#         key = "_".join(site_idx.split("_")[:2])
#         site_dict[key] = {"A": 0, "T": 0, "C": 0, "G": 0, "N": 0, "total": 0}
#
#     # -------------------------------------------------------->>>>>>>>
#     # run mpileup
#     # -------------------------------------------------------->>>>>>>>
#     mpileup_extend_length = 10
#
#     mpileup_iter = bam_obj.pileup(chr_name,
#                                   region_start - mpileup_extend_length,
#                                   region_end + mpileup_extend_length,
#                                   fastafile=genome_obj,
#                                   ignore_overlaps=ignore_overlaps,
#                                   min_base_quality=min_base_quality,
#                                   min_mapping_quality=min_mapping_quality,
#                                   stepper="samtools")
#
#     for pileup in mpileup_iter:
#         run_index = "%s_%s" % (chr_name, pileup.reference_pos + 1)
#
#         if (pileup.reference_pos + 1) > region_end:
#             break
#
#         if site_dict.get(run_index) is not None:
#             base_list = list(map(str.upper, pileup.get_query_sequences()))
#             site_dict[run_index]["A"] = base_list.count("A")
#             site_dict[run_index]["T"] = base_list.count("T")
#             site_dict[run_index]["C"] = base_list.count("C")
#             site_dict[run_index]["G"] = base_list.count("G")
#             site_dict[run_index]["N"] = base_list.count("N")
#             site_dict[run_index]["total"] = len(base_list)
#
#     return site_dict
#
#
# def find_block_info_and_highest_signal(mpmat_info, ctrl_bam_obj, treat_bam_obj, ref_genome_obj, block_mut_num_cutoff=2):
#     """
#     INPUT:
#         <ctrl_bam_obj> and <treat_bam_obj>
#             pysam.AlignmentFile obj
#
#         <ref_genome_obj>
#             pysam.Fasta obj
#
#         <mpmat_info>
#             mpmatLine obj
#
#         <block_mut_num_cutoff>
#             int, site in ctrl sample contain mutation number >= this cutoff will be blocked in the following steps
#
#     RETURN
#         <back_mpmat_line>
#             mpmatLine obj with additional features
#                 1. mpmatLine.block_info_list
#                     list, contain True or False, True means need to block in the following steps.
#
#                 2. mpmatLine.highest_site_dict
#                     dict, format like
#                         {'A': 0,
#                          'T': 6,
#                          'C': 96,
#                          'G': 0,
#                          'N': 0,
#                          'total': 102,
#                          'site_index': 'chr1_121885716',
#                          'full_site_index': 'chr1_121885716_CT',
#                          'mut_num': 6,
#                          'mut_ratio': 0.058823529411764705}
#
#                 3. block_site_num
#                     int, count block site number
#
#     """
#
#     # init params
#     if not hasattr(mpmat_info, "block_info_list"):
#         mpmat_info.block_info_list = [False] * mpmat_info.site_num
#         mpmat_info.block_site_num = 0
#
#     mut_type_base = mpmat_info.mut_type[1]
#
#     highest_treat_mut_site_dict = {
#         "A": 0,
#         "T": 0,
#         "G": 0,
#         "C": 0,
#         "N": 0,
#         "total": 0,
#         "site_index": "",
#         "full_site_index": "",
#         "mut_num": 0,
#         "mut_ratio": 0,
#         "run_state": False
#     }
#
#     # query ctrl and treat pileup info
#     collect_ctrl_site_cover_dict = query_site_mpileup_info(site_idx_list=mpmat_info.site_index_list,
#                                                            bam_obj=ctrl_bam_obj,
#                                                            genome_obj=ref_genome_obj,
#                                                            ignore_overlaps=True,
#                                                            min_base_quality=20,
#                                                            min_mapping_quality=20,
#                                                            dist_cutoff=10000)
#
#     collect_treat_site_cover_dict = query_site_mpileup_info(site_idx_list=mpmat_info.site_index_list,
#                                                             bam_obj=treat_bam_obj,
#                                                             genome_obj=ref_genome_obj,
#                                                             ignore_overlaps=True,
#                                                             min_base_quality=20,
#                                                             min_mapping_quality=20,
#                                                             dist_cutoff=10000)
#
#     # add block info and find highest signal
#     for run_idx, full_site_index in enumerate(mpmat_info.site_index_list):
#         if mpmat_info.SNP_ann_list[run_idx]:
#             mpmat_info.block_info_list[run_idx] = True
#             continue
#
#         # fix site index
#         site_index = "_".join(full_site_index.split("_")[:2])
#
#         # query site cover in ctrl sample
#         ctrl_site_cover_dict = collect_ctrl_site_cover_dict.get(site_index)
#
#         if ctrl_site_cover_dict is not None:
#             ctrl_site_mut_num = ctrl_site_cover_dict.get(mut_type_base)
#         else:
#             continue
#
#         if ctrl_site_mut_num is None:
#             continue
#
#         if ctrl_site_mut_num >= block_mut_num_cutoff:
#             mpmat_info.block_info_list[run_idx] = True
#
#         if not mpmat_info.block_info_list[run_idx]:
#             treat_site_cover_dict = collect_treat_site_cover_dict.get(site_index)
#             if not highest_treat_mut_site_dict["run_state"]:
#                 highest_treat_mut_site_dict = treat_site_cover_dict.copy()
#                 highest_treat_mut_site_dict["site_index"] = site_index
#                 highest_treat_mut_site_dict["full_site_index"] = full_site_index
#                 highest_treat_mut_site_dict["run_state"] = True
#
#             if treat_site_cover_dict[mut_type_base] > highest_treat_mut_site_dict[mut_type_base]:
#                 highest_treat_mut_site_dict = treat_site_cover_dict.copy()
#                 highest_treat_mut_site_dict["site_index"] = site_index
#                 highest_treat_mut_site_dict["full_site_index"] = full_site_index
#                 highest_treat_mut_site_dict["run_state"] = True
#
#     # calculate highest info and signal
#     if highest_treat_mut_site_dict["total"] > 0:
#         highest_treat_mut_site_dict["mut_num"] = highest_treat_mut_site_dict[mut_type_base]
#         highest_treat_mut_site_dict["mut_ratio"] = highest_treat_mut_site_dict["mut_num"] / 1.0 / \
#                                                    highest_treat_mut_site_dict["total"]
#     else:
#         highest_treat_mut_site_dict["mut_num"] = 0
#         highest_treat_mut_site_dict["mut_ratio"] = 0.0
#         highest_treat_mut_site_dict["run_state"] = False
#
#     # add dict into mpmat info
#     mpmat_info.highest_site_dict = highest_treat_mut_site_dict.copy()
#
#     # fix mut_key_list
#     for index, site_index in enumerate(mpmat_info.site_index_list):
#         if mpmat_info.block_info_list[index]:
#             mpmat_info.mut_key_list[index] = "B"
#
#         if mpmat_info.SNP_ann_list[index]:
#             mpmat_info.mut_key_list[index] = "S"
#
#     # fix mut_key
#     mpmat_info.mut_key = "-".join(mpmat_info.mut_key_list)
#
#     # add block_site_num
#     mpmat_info.block_site_num = mpmat_info.block_info_list.count(True)
#
#     return mpmat_info
#
#
# def load_reference_fasta_as_dict(ref_fasta_path, ref_name_list="All", log_verbose=30):
#     """
#     INPUT:
#         <ref_fasta_path>
#             Reference fasta file path
#
#         <ref_name_list>
#             If set All, load all seq info in reference, else only try to load seq_id in the list
#
#     RETURN
#         <ref_seq_dict>
#             A dict, key is seq_id and value is sequence with  .upper()
#
#         None
#             If occur error, return None.
#     """
#
#     # ---------------------------------------------------------->>>>>>>>>>>>>>>>>>>>
#     # log setting
#     # ---------------------------------------------------------->>>>>>>>>>>>>>>>>>>>
#     logging.basicConfig(level=(4 - log_verbose) * 10,
#                         format='%(levelname)-5s @ %(asctime)s: %(message)s ',
#                         datefmt='%Y-%m-%d %H:%M:%S',
#                         stream=sys.stderr,
#                         filemode="w")
#
#     # ---------------------------------------------------------->>>>>>>>>>>>>>>>>>>>
#     # load genome as dict
#     # ---------------------------------------------------------->>>>>>>>>>>>>>>>>>>>
#     try:
#         genome_fa = SeqIO.parse(handle=ref_fasta_path, format="fasta")
#     except:
#         raise IOError("Load file error! %s" % ref_fasta_path)
#
#     # init var
#     ref_seq_dict = {}
#     ref_name_set = set(ref_name_list)
#
#     logging.info("Starting to load the reference genome...")
#
#     for ref in genome_fa:
#         if ref_name_list == "All":
#             ref_seq_dict[ref.id] = ref.seq.upper()
#             logging.debug("Loading genome...\t" + ref.id)
#
#         elif ref.id in ref_name_list:
#             ref_seq_dict[ref.id] = ref.seq.upper()
#             logging.debug("Loading genome...\t" + ref.id)
#
#             # remove already loaded seq
#             ref_name_set.remove(ref.id)
#
#             # load all info
#             if len(ref_name_set) == 0:
#                 break
#
#     logging.info("Loading genome done!")
#
#     return ref_seq_dict
#
#
# def check_input_bam_file(input_bam_filename):
#     """
#     HELP
#
#         1. check exist
#
#         2. check sort state
#
#         3. check .bai index file
#
#     RETURN
#         0 alright
#
#         Not 0:
#             1 not exist
#             2 not sorted
#             3 index not exist
#     """
#     # check exist
#     if not os.path.exists(input_bam_filename):
#         return 1
#
#     # check sort
#     bam_file = pysam.AlignmentFile(input_bam_filename, "rb")
#     if bam_file.header.get("HD").get("SO") != "coordinate":
#         return 2
#
#     # check bai index
#     bai_filename_1 = input_bam_filename + ".bai"
#     bai_filename_2 = os.path.splitext(input_bam_filename)[0] + ".bai"
#     if (not os.path.exists(bai_filename_1)) and (not os.path.exists(bai_filename_2)):
#         return 3
#
#     return 0
#
#
# def get_BAM_ref_name(bam_filename, ref_count_cutoff=0):
#     """
#     INPUT:
#         <bam_filename>
#             BAM file path with .bai index file in a same dir
#
#     RETURN
#         <ref_name_list>
#             return ref name list, which contain more than ref_count_cutoff mapped reads
#
#     """
#     bam_file = pysam.AlignmentFile(bam_filename, "rb")
#     ref_name_list = []
#
#     for bam_index_info in bam_file.get_index_statistics():
#         ref_name = bam_index_info[0]
#         map_count = bam_index_info[1]
#         unmap_count = bam_index_info[2]
#
#         if map_count > ref_count_cutoff:
#             if ref_name != "*":
#                 ref_name_list.append(ref_name)
#
#     bam_file.close()
#
#     return ref_name_list
#
#
# def count_effective_genome(genome_fa_filename, genome_json_out=None, log_verbose=3):
#     """
#     INPUT:
#         <genome_fa_filename>
#             str, Genome FASTA filename
#
#         <genome_json_out>
#             str, Output filename with JSON format, default=stdout means print info to screen.
#
#         <log_verbose>
#             int, log output info level, bigger means more log info
#
#     OUTPUT:
#         Count each chromosome A,T,C,G,N count, and out info with JSON format.
#
#         {
#             "chr1": {
#                 "A": 67070277,
#                 "C": 48055043,
#                 "G": 48111528,
#                 "N": 18475410,
#                 "T": 67244164
#             }
#         }
#     """
#
#     # --------------------------------------------------->>>>>>
#     # log setting
#     # --------------------------------------------------->>>>>>
#     logging.basicConfig(level=(4 - log_verbose) * 10,
#                         format='%(levelname)-5s @ %(asctime)s: %(message)s ',
#                         datefmt='%Y-%m-%d %H:%M:%S',
#                         stream=sys.stderr,
#                         filemode="w")
#
#     logging.info("-" * 80)
#     logging.info("Counting reference genome background info...")
#     logging.info("Processing FASTA file...")
#
#     # --------------------------------------------------->>>>>>
#     # load genome file
#     # --------------------------------------------------->>>>>>
#     genome_fa = SeqIO.parse(handle=genome_fa_filename, format="fasta")
#
#     genome_base_count_dict = {
#         "total": {"A": 0, "T": 0, "C": 0, "G": 0, "N": 0}
#     }
#
#     for ref in genome_fa:
#         logging.debug("Counting with %s" % ref.id)
#
#         chr_seq = ref.seq.upper()
#
#         if genome_base_count_dict.get(ref.id) is None:
#             genome_base_count_dict[ref.id] = {}
#
#         # add into dict
#         genome_base_count_dict[ref.id]["A"] = chr_seq.count("A")
#         genome_base_count_dict[ref.id]["T"] = chr_seq.count("T")
#         genome_base_count_dict[ref.id]["C"] = chr_seq.count("C")
#         genome_base_count_dict[ref.id]["G"] = chr_seq.count("G")
#         genome_base_count_dict[ref.id]["N"] = chr_seq.count("N")
#
#         # add to total
#         genome_base_count_dict["total"]["A"] += genome_base_count_dict[ref.id]["A"]
#         genome_base_count_dict["total"]["T"] += genome_base_count_dict[ref.id]["T"]
#         genome_base_count_dict["total"]["C"] += genome_base_count_dict[ref.id]["C"]
#         genome_base_count_dict["total"]["G"] += genome_base_count_dict[ref.id]["G"]
#         genome_base_count_dict["total"]["N"] += genome_base_count_dict[ref.id]["N"]
#
#     # log
#     logging.info("Processing FASTA file. Done!")
#     logging.info("-" * 80)
#
#     # --------------------------------------------------->>>>>>
#     # output json
#     # --------------------------------------------------->>>>>>
#     genome_base_count_dict["reference"] = os.path.abspath(genome_fa_filename)
#
#     if genome_json_out == "stdout":
#         out_file = sys.stdout
#
#     elif genome_json_out is None:
#         out_file = None
#
#     else:
#         out_file = open(genome_json_out, "w")
#
#     if out_file is not None:
#         json_obj = json.dumps(genome_base_count_dict, encoding="utf-8", sort_keys=True, indent=4,
#                               separators=(', ', ': '))
#         out_file.write(json_obj)
#
#     # close file
#     if (genome_json_out != "stdout") and (out_file is not None):
#         out_file.close()
#
#     genome_fa.close()
#
#     # --------------------------------------------------->>>>>>
#     # return part
#     # --------------------------------------------------->>>>>>
#     return genome_base_count_dict
#
#
# def count_split_chrom_mut_bg(bmat_filename, ref_name, query_mut_type_list=["CT", "GA"], json_out_filename=None,
#                              log_verbose=3):
#     """
#     INPUT
#         <bmat_filename>
#             str, .bmat file
#
#         <ref_name>
#             str, ref_name like chr1
#
#     OUTPUT
#         <json_out_filename>
#             str, default=stdout, output calculation reault in JSON format.
#
#         <log_out_filename>
#             str, default=stderr, output run log.
#
#     RETURN
#         dict, count_dict, with the same info to <json_out_filename>
#
#         return dict like:
#
#         {
#             "A": 523999,
#             "AC": 136,
#             "AG": 553,
#             "AT": 145,
#             "C": 566244,
#             "CA": 593,
#             "CG": 190,
#             "CT": 3845,
#             "G": 540117,
#             "GA": 3405,
#             "GC": 202,
#             "GT": 551,
#             "T": 494954,
#             "TA": 168,
#             "TC": 777,
#             "TG": 154,
#             "query_mut_base_count": 7250,
#             "query_mut_bg_pval": 0.006553014793543879,
#             "query_total_base_count": 1106361,
#             "total_base_count": 2125314,
#             "total_mut_base_count": 10719,
#             "total_mut_bg_pval": 0.005043490044294632
#         }
#
#     """
#     # --------------------------------------------------->>>>>>
#     # log setting
#     # --------------------------------------------------->>>>>>
#     logging.basicConfig(level=(4 - log_verbose) * 10,
#                         format='%(levelname)-5s @ %(asctime)s: %(message)s ',
#                         datefmt='%Y-%m-%d %H:%M:%S',
#                         stream=sys.stderr,
#                         filemode="w")
#
#     # ------------------------------------------------------>>>>>>
#     # open file
#     # ------------------------------------------------------>>>>>>
#     # open json out
#     if json_out_filename == "stdout":
#         json_out = sys.stdout
#
#     elif json_out_filename is None:
#         json_out = None
#
#     else:
#         json_out = open(json_out_filename, "w")
#
#     # open input bmat file
#     if (bmat_filename[-3:] == ".gz") or (bmat_filename[-5:] == ".gzip"):
#         input_file = gzip.open(bmat_filename, "r")
#     else:
#         input_file = open(bmat_filename, "r")
#
#     # log
#     logging.info("-" * 80)
#     logging.info("Start to process bmat file. Set ref_name as: %s" % ref_name)
#
#     # ------------------------------------------------------>>>>>>
#     # check header
#     # ------------------------------------------------------>>>>>>
#     header = input_file.readline()
#
#     if not ("chr_name" in header):
#         input_file.seek(0)
#         logging.warning("Input .bmat file does not contain header line.")
#
#     # ------------------------------------------------------>>>>>>
#     # init count dict
#     # ------------------------------------------------------>>>>>>
#     base_count_dict = {
#         "total_base_count": 0,
#         "total_mut_base_count": 0,
#         "total_mut_bg_pval": 0.0001,
#         "query_total_base_count": 0,
#         "query_mut_base_count": 0,
#         "query_mut_bg_pval": 0.0001
#     }
#
#     for from_base in "ATCGN":
#         for to_base in "ATCGN":
#             if from_base == to_base:
#                 key = from_base
#             else:
#                 key = from_base + to_base
#
#             base_count_dict[key] = 0
#
#     # ------------------------------------------------------>>>>>>
#     # count mutation and cover info
#     # ------------------------------------------------------>>>>>>
#     for run_index, line in enumerate(input_file):
#         if run_index % 1000000 == 0:
#             logging.debug("Processing %s bmat file... %s" % (ref_name, run_index))
#
#         line_list = line.strip().split("\t")
#
#         chr_name = line_list[0]
#         ref_base = line_list[2]
#
#         # check chr name
#         if chr_name != ref_name:
#             logging.error(".bmat file chr_name does not match <ref_name>!")
#             raise IOError(".bmat file chr_name does not match <ref_name>!")
#
#         # init record cover num
#         cover_base_count = 0
#
#         for index, to_base in enumerate("AGCTN"):
#             base_count = int(line_list[3 + index])
#             cover_base_count += base_count
#
#             # record mutation info
#             if ref_base != to_base:
#                 key = ref_base + to_base
#                 base_count_dict[key] += base_count
#
#         # accumulation total cover
#         base_count_dict[ref_base] += cover_base_count
#
#     # ------------------------------------------------------>>>>>>
#     # count query and total, calculate bg p.val
#     # ------------------------------------------------------>>>>>>
#     for query_mut_type in query_mut_type_list:
#         query_ref_base = query_mut_type[0]
#
#         # accumulation query number
#         base_count_dict["query_total_base_count"] += base_count_dict[query_ref_base]
#         base_count_dict["query_mut_base_count"] += base_count_dict[query_mut_type]
#
#     for ref_base in "ATCGN":
#         for to_base in "ATCGN":
#             if ref_base == to_base:
#                 key = ref_base
#                 base_count_dict["total_base_count"] += base_count_dict[key]
#
#             else:
#                 key = ref_base + to_base
#                 base_count_dict["total_mut_base_count"] += base_count_dict[key]
#
#     # calculate background pval
#     base_count_dict["total_mut_bg_pval"] = base_count_dict["total_mut_base_count"] / 1.0 / base_count_dict[
#         "total_base_count"]
#     base_count_dict["query_mut_bg_pval"] = base_count_dict["query_mut_base_count"] / 1.0 / base_count_dict[
#         "query_total_base_count"]
#
#     # ------------------------------------------------------>>>>>>
#     # return and output
#     # ------------------------------------------------------>>>>>>
#     # output JSON result
#     if json_out is not None:
#         json_obj = json.dumps(base_count_dict, encoding="utf-8", sort_keys=True, indent=4, separators=(', ', ': '))
#         json_out.write(json_obj)
#
#         # out log files
#     logging.info("Processing bmat file. DONE!")
#
#     # close files
#     if (json_out_filename != "stdout") and (json_out is not None):
#         json_out.close()
#
#     input_file.close()
#
#     # return
#     return base_count_dict
#
#
# def split_bmat_by_chr_name(input_bmat_filename, temp_dir=None, force_temp_dir=True, log_verbose=3, out_gzip=False):
#     """
#     INPUT
#         <input_bmat_filename>
#             str, .bmat filename, support .gz OR .gzip
#
#         <temp_dir>
#             str, a dir to store temp files, None means the same dir with <input_bmat_filename>
#
#         <log_verbose>
#             int, log output info level, bigger means more log info
#
#     OUTPUT
#         Split files by chr_name
#
#     RETURN
#         A dict, structure like:
#
#         dict = {
#             "chr1" : "file.chr1.saFSDfjsj91.bmat",
#             "chr2" : "file.chr2.saFSDasjfj2.bmat"
#             ... ...
#         }
#
#     """
#     # --------------------------------------------------->>>>>>
#     # log setting
#     # --------------------------------------------------->>>>>>
#     logging.basicConfig(level=(4 - log_verbose) * 10,
#                         format='%(levelname)-5s @ %(asctime)s: %(message)s ',
#                         datefmt='%Y-%m-%d %H:%M:%S',
#                         stream=sys.stderr,
#                         filemode="w")
#
#     # --------------------------------------------------->>>>>>
#     # set temp dir
#     # --------------------------------------------------->>>>>>
#     if temp_dir is None:
#         temp_dir = os.path.abspath(os.path.dirname(input_bmat_filename))
#     else:
#         temp_dir = os.path.abspath(temp_dir)
#
#     if temp_dir[-1] != "/":
#         temp_dir += "/"
#
#     # temp dir check and create
#     if not os.path.exists(temp_dir):
#         if force_temp_dir:
#             logging.warning("<temp_dir> setting is not exist \t %s " % temp_dir)
#             logging.warning("<force_temp_dir> set as True, try to create temp dir \t %s" % temp_dir)
#
#             try:
#                 os.makedirs(os.path.abspath(temp_dir))
#             except:
#                 logging.warning("Temp dir creating error: \t %s" % temp_dir)
#                 logging.warning("set <temp_dir> as the same dir with <input_bmat_filename>")
#                 temp_dir = os.path.abspath(os.path.dirname(input_bmat_filename))
#
#         else:
#             temp_dir = os.path.abspath(os.path.dirname(input_bmat_filename))
#             logging.warning("<temp_dir> setting is not exist, set <temp_dir> as %s" % temp_dir)
#     else:
#         temp_dir = os.path.abspath(temp_dir)
#
#     # --------------------------------------------------->>>>>>
#     # get input basename
#     # --------------------------------------------------->>>>>>
#     input_file_basename = os.path.basename(input_bmat_filename)
#
#     # --------------------------------------------------->>>>>>
#     # make record dict
#     # --------------------------------------------------->>>>>>
#     record_dict = {
#         "chr_name_order": [],
#     }
#
#     # --------------------------------------------------->>>>>>
#     # split file
#     # --------------------------------------------------->>>>>>
#     logging.info("Try to split bmat file...")
#     logging.info("Output dir is %s" % temp_dir)
#
#     # open input bmat file
#     if (input_bmat_filename[-3:] == ".gz") or (input_bmat_filename[-5:] == ".gzip"):
#         input_file = gzip.open(input_bmat_filename, "r")
#     else:
#         input_file = open(input_bmat_filename, "r")
#
#     # set init
#     cur_chr_name = None
#     cur_out_file = None
#
#     for line in input_file:
#         line_list = line.strip().split("\t")
#         chr_name = line_list[0]
#
#         if chr_name == "chr_name":
#             continue
#
#         if cur_chr_name is None:
#             # log
#             logging.info("Processing %s .bmat file" % chr_name)
#
#             cur_chr_name = chr_name
#
#             # make temp filename
#             temp_file_basename = "temp_" + input_file_basename + "." + cur_chr_name + "." + "".join(
#                 random.sample(string.ascii_letters + string.digits, 16))
#             temp_file_name = os.path.join(temp_dir, temp_file_basename)
#
#             if out_gzip:
#                 temp_file_name += ".gz"
#
#             # record info into dict
#             record_dict["chr_name_order"].append(cur_chr_name)
#             record_dict[cur_chr_name] = temp_file_name
#
#             # write
#             if not out_gzip:
#                 cur_out_file = open(temp_file_name, "w")
#             else:
#                 cur_out_file = gzip.open(temp_file_name, "w")
#
#             cur_out_file.write(line)
#
#         elif cur_chr_name == chr_name:
#             cur_out_file.write(line)
#
#         elif cur_chr_name != chr_name:
#             cur_out_file.close()
#
#             # log
#             logging.info("Processing %s .bmat file" % chr_name)
#
#             # set next chr_name
#             cur_chr_name = chr_name
#
#             # make temp filename
#             temp_file_basename = "temp_" + input_file_basename + "." + cur_chr_name + "." + "".join(
#                 random.sample(string.ascii_letters + string.digits, 16))
#             temp_file_name = os.path.join(temp_dir, temp_file_basename)
#
#             if out_gzip:
#                 temp_file_name += ".gz"
#
#             # record info into dict
#             record_dict["chr_name_order"].append(cur_chr_name)
#             record_dict[cur_chr_name] = temp_file_name
#
#             # write
#             if not out_gzip:
#                 cur_out_file = open(temp_file_name, "w")
#             else:
#                 cur_out_file = gzip.open(temp_file_name, "w")
#
#             cur_out_file.write(line)
#
#     # close all files
#     try:
#         cur_out_file.close()
#         input_file.close()
#     except:
#         logging.error("Error occurs at close file step.")
#         raise IOError("Error occurs at close file step.")
#
#     # log
#     logging.info("Try to split bmat file. DONE!")
#
#     # --------------------------------------------------->>>>>>
#     # return
#     # --------------------------------------------------->>>>>>
#     return record_dict
#
#
# def multi_calculate_SNP_bg_pval(split_bmat_dict, threads=1, query_mut_type_list=["CT", "GA"], log_verbose=3):
#     """
#     INPUT:
#         <split_bmat_dict>
#             dict, contain split temp files, format like:
#
#             {
#                 "chr1": "./temp_test.merge_chr1_chr2_chr3.1K.bmat.gz.chr1.k1LKEjDWSagudoQI",
#                 "chr2": "./temp_test.merge_chr1_chr2_chr3.1K.bmat.gz.chr2.KUrM1AL2dCoHYRlE",
#                 "chr3": "./temp_test.merge_chr1_chr2_chr3.1K.bmat.gz.chr3.0TJKF4pHAPErWDMi",
#                 "chr_name_order": [
#                     "chr1",
#                     "chr2",
#                     "chr3"
#                 ]
#             }
#
#         <threads>
#             int, set threads number.
#
#         <query_mut_type_list>
#             list, inherit from FUN: count_split_chrom_mut_bg
#
#     RETURN
#         <mutation_bg_dict>
#             dict, contain genome mutation background info, key is ref_name, val is mutation info
#
#             with format like:
#             {
#                 "chr1":{
#                         'A': 523999,
#                         'AC': 136,
#                         'AG': 553,
#                         'AT': 145,
#                         'C': 566244,
#                         'CA': 593,
#                         'CG': 190,
#                         'CT': 3845,
#                         'G': 540117,
#                         'GA': 3405,
#                         'GC': 202,
#                         'GT': 551,
#                         'T': 494954,
#                         'TA': 168,
#                         'TC': 777,
#                         'TG': 154,
#                         'query_mut_base_count': 7250,
#                         'query_mut_bg_pval': 0.006553014793543879,
#                         'query_total_base_count': 1106361,
#                         'total_base_count': 2125314,
#                         'total_mut_base_count': 10719,
#                         'total_mut_bg_pval': 0.005043490044294632
#                 }, ...
#
#             }
#
#     """
#
#     # --------------------------------------------------->>>>>>
#     # log setting
#     # --------------------------------------------------->>>>>>
#     logging.basicConfig(level=(4 - log_verbose) * 10,
#                         format='%(levelname)-5s @ %(asctime)s: %(message)s ',
#                         datefmt='%Y-%m-%d %H:%M:%S',
#                         stream=sys.stderr,
#                         filemode="w")
#
#     # --------------------------------------------------->>>>>>
#     # make record dict
#     # --------------------------------------------------->>>>>>
#     mutation_bg_dict = {
#         "chr_name_order": []
#     }
#
#     # --------------------------------------------------->>>>>>
#     # run part
#     # --------------------------------------------------->>>>>>
#     logging.info("-" * 80)
#     logging.info("Start to calculate chromosome mutation background...")
#
#     pool = multiprocessing.Pool(processes=threads)
#
#     count_split_result_list = []
#
#     for chr_ref_name in split_bmat_dict["chr_name_order"]:
#         split_bmat_filename = split_bmat_dict[chr_ref_name]
#
#         # pool and results
#         count_split_result_list.append(
#             pool.apply_async(
#                 func=count_split_chrom_mut_bg,
#                 args=(
#                     split_bmat_filename,
#                     chr_ref_name,
#                     query_mut_type_list,
#                     None,
#                     0,
#                 )
#             )
#         )
#
#     pool.close()
#     pool.join()
#
#     logging.info("Calculating chromosome mutation background done!")
#     logging.info("-" * 80)
#     # --------------------------------------------------->>>>>>
#     # get result and return part
#     # --------------------------------------------------->>>>>>
#     for index, res in enumerate(count_split_result_list):
#         chr_count_res = res.get().copy()
#         chr_name = split_bmat_dict["chr_name_order"][index]
#
#         # add into dict
#         mutation_bg_dict[chr_name] = chr_count_res
#         mutation_bg_dict["chr_name_order"].append(chr_name)
#
#     return mutation_bg_dict
#
#
# def split_mpmat_by_chr_name(input_mpmat_filename, temp_dir=None, force_temp_dir=True, log_verbose=3):
#     """
#     INPUT
#         <input_mpmat_filename>
#             str, .mpmat filename, support .gz OR .gzip
#
#         <temp_dir>
#             str, a dir to store temp files, None means the same dir with <input_mpmat_filename>
#
#         <log_verbose>
#             int, log output info level, bigger means more log info
#
#     OUTPUT
#         Split files by chr_name
#
#     RETURN
#         A dict, structure like:
#
#         dict = {
#             "chr1" : "file.chr1.saFSDfjsj91.mpmat",
#             "chr2" : "file.chr2.saFSDasjfj2.mpmat"
#             ... ...
#         }
#
#     """
#     # --------------------------------------------------->>>>>>
#     # log setting
#     # --------------------------------------------------->>>>>>
#     logging.basicConfig(level=(4 - log_verbose) * 10,
#                         format='%(levelname)-5s @ %(asctime)s: %(message)s ',
#                         datefmt='%Y-%m-%d %H:%M:%S',
#                         stream=sys.stderr,
#                         filemode="w")
#
#     # --------------------------------------------------->>>>>>
#     # set temp dir
#     # --------------------------------------------------->>>>>>
#     if temp_dir is None:
#         temp_dir = os.path.abspath(os.path.dirname(input_mpmat_filename))
#     else:
#         temp_dir = os.path.abspath(temp_dir)
#
#     # add back slash
#     if temp_dir[-1] != "/":
#         temp_dir += "/"
#
#     # temp dir check and create
#     if not os.path.exists(temp_dir):
#         if force_temp_dir:
#             logging.warning("<temp_dir> setting is not exist \t %s " % temp_dir)
#             logging.warning("<force_temp_dir> set as True, try to create temp dir \t %s" % temp_dir)
#
#             try:
#                 os.makedirs(os.path.abspath(temp_dir))
#             except:
#                 logging.warning("Temp dir creating error: \t %s" % temp_dir)
#                 logging.warning("set <temp_dir> as the same dir with <input_mpmat_filename>")
#                 temp_dir = os.path.abspath(os.path.dirname(input_mpmat_filename))
#
#         else:
#             temp_dir = os.path.abspath(os.path.dirname(input_mpmat_filename))
#             logging.warning("<temp_dir> setting is not exist, set <temp_dir> as %s" % temp_dir)
#     else:
#         temp_dir = os.path.abspath(temp_dir)
#
#     # --------------------------------------------------->>>>>>
#     # get input basename
#     # --------------------------------------------------->>>>>>
#     input_file_basename = os.path.basename(input_mpmat_filename)
#
#     # --------------------------------------------------->>>>>>
#     # make record dict
#     # --------------------------------------------------->>>>>>
#     record_dict = {
#         "chr_name_order": []
#     }
#
#     # --------------------------------------------------->>>>>>
#     # split file
#     # --------------------------------------------------->>>>>>
#     logging.info("Try to split mpmat file...")
#     logging.info("Output dir is %s" % temp_dir)
#
#     # open input mpmat file
#     if (input_mpmat_filename[-3:] == ".gz") or (input_mpmat_filename[-5:] == ".gzip"):
#         input_file = gzip.open(input_mpmat_filename, "rt")
#     else:
#         input_file = open(input_mpmat_filename, "rt")
#
#     # set init
#     cur_chr_name = None
#     cur_out_file = None
#
#     for line in input_file:
#         line_list = line.strip().split("\t")
#         chr_name = line_list[0]
#
#         if chr_name == "chr_name":
#             continue
#
#         if cur_chr_name is None:
#             # log
#             logging.info("Processing %s .mpmat file" % chr_name)
#
#             cur_chr_name = chr_name
#
#             # make temp filename
#             temp_file_basename = "temp_" + input_file_basename + "." + cur_chr_name + "." + "".join(
#                 random.sample(string.ascii_letters + string.digits, 16))
#             temp_file_name = os.path.join(temp_dir, temp_file_basename)
#
#             # record info into dict
#             record_dict["chr_name_order"].append(cur_chr_name)
#             record_dict[cur_chr_name] = temp_file_name
#
#             # write
#             cur_out_file = open(temp_file_name, "w")
#             cur_out_file.write(line)
#
#         elif cur_chr_name == chr_name:
#             cur_out_file.write(line)
#
#         elif cur_chr_name != chr_name:
#             cur_out_file.close()
#
#             # log
#             logging.info("Processing %s .mpmat file" % chr_name)
#
#             # set next chr_name
#             cur_chr_name = chr_name
#
#             # make temp filename
#             temp_file_basename = "temp_" + input_file_basename + "." + cur_chr_name + "." + "".join(
#                 random.sample(string.ascii_letters + string.digits, 16))
#             temp_file_name = os.path.join(temp_dir, temp_file_basename)
#
#             # record info into dict
#             record_dict["chr_name_order"].append(cur_chr_name)
#             record_dict[cur_chr_name] = temp_file_name
#
#             # write
#             cur_out_file = open(temp_file_name, "w")
#             cur_out_file.write(line)
#
#     # close all files
#     try:
#         cur_out_file.close()
#         input_file.close()
#     except:
#         logging.error("Error occurs at close file step.")
#         raise IOError("Error occurs at close file step.")
#
#     # log
#     logging.info("Try to split .mpmat file. DONE!")
#
#     # --------------------------------------------------->>>>>>
#     # return
#     # --------------------------------------------------->>>>>>
#     return record_dict
#
#
# def _count_mismatch_num(align_mismatch_pairs, region_mut_type_list=["CT", "GA"]):
#     """
#     INPUT:
#         <align_mismatch_pairs>
#             list, generate from
#
#     RETURN:
#         <mismatch_count_dict>
#
#     """
#     total_mut_num = 0
#     query_mut_num = 0
#     other_mut_num = 0
#
#     for ref_idx, align_idx, ref_base, align_base in align_mismatch_pairs:
#         total_mut_num += 1
#         mut_type = ref_base + align_base
#
#         if mut_type in region_mut_type_list:
#             query_mut_num += 1
#         else:
#             other_mut_num += 1
#
#     return query_mut_num, other_mut_num, total_mut_num
#
#
# def count_chr_bam_mut_count(input_bam_filename, select_chr_name, ref_genome_filename, query_mut_type_list=["CT", "GA"],
#                             query_mut_min_cutoff=1, query_mut_max_cutoff=16, total_mut_max_cutoff=20,
#                             other_mut_max_cutoff=16, log_verbose=3):
#     """
#     INPUT:
#         <input_bam_filename>
#             str, bam filename
#
#         <select_chr_name>
#             str, like chr1, chr2, chr3 ...
#
#         <ref_genome_filename>
#             str, reference genome FASTA filename
#
#         <query_mut_type_list>
#             list, like ["CT","GA"]
#
#         <query_mut_min_cutoff>
#             int, if mutation number >= query_mut_min_cutoff in the mpmat region,
#                 will be counted as 'region_mut_count'.
#
#         <query_mut_max_cutoff>
#             int, larger than this will be marked.
#
#         <total_mut_max_cutoff>
#             int, larger than this will be marked.
#
#         <other_mut_max_cutoff>
#             int, larger than this will be marked.
#
#     RETURN:
#         <align_count_dict>
#             dict, key and value like:
#                 {
#                     "all_align_count": 0,
#                     "region_mut_count": 0,
#                     "region_non_mut_count": 0,
#                     "total_high_mismatch_count": 0,
#                     "other_high_mismatch_count": 0,
#                     "query_high_mismatch_count": 0,
#                     "all_filter_count": 0
#                 }
#
#     """
#     # ---------------------------------------------------------->>>>>>>>>>>>>>>>>>>>
#     # log setting
#     # ---------------------------------------------------------->>>>>>>>>>>>>>>>>>>>
#     logging.basicConfig(level=(4 - log_verbose) * 10,
#                         format='%(levelname)-5s @ %(asctime)s: %(message)s ',
#                         datefmt='%Y-%m-%d %H:%M:%S',
#                         stream=sys.stderr,
#                         filemode="w")
#
#     # ---------------------------------------------------------->>>>>>>>>>
#     # load genome and open BAM file
#     # ---------------------------------------------------------->>>>>>>>>>
#     input_bam = pysam.AlignmentFile(input_bam_filename, "rb")
#
#     ref_genome_dict = load_reference_fasta_as_dict(ref_fasta_path=ref_genome_filename,
#                                                    ref_name_list=[select_chr_name])
#     # ---------------------------------------------------------->>>>>>>>>>
#     # init var
#     # ---------------------------------------------------------->>>>>>>>>>
#     logging.info("Starting to count alignment mutation info on %s ..." % select_chr_name)
#
#     # define count dict
#     align_count_dict = {
#         "all_align_count": 0,
#         "region_mut_count": 0,
#         "region_non_mut_count": 0,
#         "total_high_mismatch_count": 0,
#         "other_high_mismatch_count": 0,
#         "query_high_mismatch_count": 0,
#         "all_filter_count": 0
#     }
#
#     # ---------------------------------------------------------->>>>>>>>>>
#     # iter align info
#     # ---------------------------------------------------------->>>>>>>>>>
#     for align in input_bam.fetch(reference=select_chr_name):
#         # count total
#         align_count_dict["all_align_count"] += 1
#
#         # make sure MD state
#         MD_state = True
#         try:
#             MD_str_tag = align.get_tag("MD")
#         except:
#             MD_state = False
#
#         # get mismatch pairs
#         if MD_state:
#             align_mismatch_pairs = get_align_mismatch_pairs(align)
#         else:
#             align_mismatch_pairs = get_No_MD_align_mismatch_pairs(align, ref_genome_dict)
#
#         # analysis of mismatch pairs
#         if align_mismatch_pairs is None:
#             align_count_dict["region_non_mut_count"] += 1
#
#         else:
#             # load mut num
#             align_query_mut_num, align_other_mut_num, align_total_mut_num = _count_mismatch_num(
#                 align_mismatch_pairs=align_mismatch_pairs,
#                 region_mut_type_list=query_mut_type_list)
#
#             # filter high mismatch reads
#             if align_total_mut_num >= total_mut_max_cutoff:
#                 align_count_dict["total_high_mismatch_count"] += 1
#                 align_count_dict["all_filter_count"] += 1
#
#             elif align_other_mut_num >= other_mut_max_cutoff:
#                 align_count_dict["other_high_mismatch_count"] += 1
#                 align_count_dict["all_filter_count"] += 1
#
#             elif align_query_mut_num >= query_mut_max_cutoff:
#                 align_count_dict["query_high_mismatch_count"] += 1
#                 align_count_dict["all_filter_count"] += 1
#
#             else:
#                 if align_total_mut_num == 0:
#                     align_count_dict["region_non_mut_count"] += 1
#
#                 elif align_query_mut_num >= query_mut_min_cutoff:
#                     align_count_dict["region_mut_count"] += 1
#
#     # close file
#     input_bam.close()
#
#     logging.info("Counting Poisson background on %s ... Done!" % select_chr_name)
#
#     return align_count_dict
#
#
# def multi_thread_chr_bam_mut_count(input_bam_filename, select_chr_name_list, ref_genome_filename, thread_num=1,
#                                    query_mut_type_list=["CT", "GA"], query_mut_min_cutoff=1, query_mut_max_cutoff=16,
#                                    total_mut_max_cutoff=20, other_mut_max_cutoff=16, log_verbose=3):
#     """
#     INPUT:
#         <input_bam_filename>
#             str, bam filename
#
#         <select_chr_name_list>
#             str, like ["chr1", "chr2", "chr3"]
#
#         <ref_genome_filename>
#             str, reference genome FASTA filename
#
#         <thread_num>
#             int, thread number to use
#
#         <query_mut_type_list>
#             list, like ["CT","GA"]
#
#         <query_mut_min_cutoff>
#             int, if mutation number >= query_mut_min_cutoff in the mpmat region,
#                 will be counted as 'region_mut_count'.
#
#         <query_mut_max_cutoff>
#             int, larger than this will be marked.
#
#         <total_mut_max_cutoff>
#             int, larger than this will be marked.
#
#         <other_mut_max_cutoff>
#             int, larger than this will be marked.
#
#     RETURN:
#         <meta_count_dict>
#             meta_count_dict = {
#                 "meta_data" : {
#                     "create_date" : time,
#                     "input_BAM" : in_bam_filename,
#                     "reference" : ref_genome_filename,
#                     "threads" : thread_num,
#                     "query_mutation_type" : val,
#                     "query_mutation_min_cutoff" : val,
#                     "query_mutation_max_cutoff" : val,
#                     "other_mismatch_max_cutoff" : val,
#                     "total_mismatch_max_cutoff" : val
#                 },
#                 "count_dict" : {
#                     "chr1" : {
#                         "all_align_count": 0,
#                         "region_mut_count": 0,
#                         "region_non_mut_count": 0,
#                         "total_high_mismatch_count": 0,
#                         "other_high_mismatch_count": 0,
#                         "query_high_mismatch_count": 0,
#                         "all_filter_count": 0
#                     },
#                     "chr2" : {
#                         "all_align_count": 0,
#                         "region_mut_count": 0,
#                         "region_non_mut_count": 0,
#                         "total_high_mismatch_count": 0,
#                         "other_high_mismatch_count": 0,
#                         "query_high_mismatch_count": 0,
#                         "all_filter_count": 0
#                     }....
#                 }
#             }
#     """
#
#     # ---------------------------------------------------------->>>>>>>>>>>>>>>>>>>>
#     # log setting
#     # ---------------------------------------------------------->>>>>>>>>>>>>>>>>>>>
#     logging.basicConfig(level=(4 - log_verbose) * 10,
#                         format='%(levelname)-5s @ %(asctime)s: %(message)s ',
#                         datefmt='%Y-%m-%d %H:%M:%S',
#                         stream=sys.stderr,
#                         filemode="w")
#
#     # ---------------------------------------------------------->>>>>>>>>>>>>>>>>>>>
#     # log setting
#     # ---------------------------------------------------------->>>>>>>>>>>>>>>>>>>>
#     logging.info("Starting to count genome Poisson lambada background...")
#
#     pool = multiprocessing.Pool(processes=thread_num)
#
#     # record info
#     input_BAM_count_result = []
#
#     for ref_name in select_chr_name_list:
#         input_BAM_count_result.append(
#             pool.apply_async(
#                 func=count_chr_bam_mut_count,
#                 args=(
#                     input_bam_filename,
#                     ref_name,
#                     ref_genome_filename,
#                     query_mut_type_list,
#                     query_mut_min_cutoff,
#                     query_mut_max_cutoff,
#                     total_mut_max_cutoff,
#                     other_mut_max_cutoff,
#                     log_verbose,
#                 )
#             )
#         )
#
#     pool.close()
#     pool.join()
#
#     # out pool result
#     input_BAM_count_dict = {}
#     for index, res in enumerate(input_BAM_count_result):
#         run_res = res.get()
#         input_BAM_count_dict[select_chr_name_list[index]] = run_res.copy()
#
#     logging.info("Calculation done!")
#
#     # make meta dict
#     meta_count_dict = {
#         "meta_data": {
#             "create_date": time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()),
#             "input_BAM": os.path.abspath(input_bam_filename),
#             "reference": os.path.abspath(ref_genome_filename),
#             "threads": thread_num,
#             "query_mutation_type": ",".join(query_mut_type_list),
#             "query_mutation_min_cutoff": query_mut_min_cutoff,
#             "query_mutation_max_cutoff": query_mut_max_cutoff,
#             "other_mismatch_max_cutoff": other_mut_max_cutoff,
#             "total_mismatch_max_cutoff": total_mut_max_cutoff
#         },
#         "count_dict": input_BAM_count_dict.copy()
#     }
#
#     return meta_count_dict
#
#
# def count_dict_sum_up(count_dict, key_ref_name_list="All", ignore_key_list=["total"]):
#     """
#     INPUT
#         <count_dict>
#             key = chr1, chr2, chr3 ...
#             value = dict and key is
#
#                   non_mut_align_count
#                   all_align_count
#                   total_high_mismatch_count
#                   query_high_mismatch_count
#                   mut_align_count
#                   other_high_mismatch_coun
#                   all_filter_count
#
#         <key_ref_name_list>
#             use those key to sum up and count
#
#         <ignore_key_list>
#             ignore
#
#     RETURN
#         Dict like <count_dict> and add 'total' key
#
#     INFO
#         Final changes, 2020-11-05
#     """
#
#     if key_ref_name_list == "All":
#         key_ref_name_list = list(count_dict.keys())
#
#     total_dict = {}
#     inner_key_list = list(count_dict[key_ref_name_list[0]].keys())
#
#     for inner_key in inner_key_list:
#         total_dict[inner_key] = 0
#
#     for ref_key in key_ref_name_list:
#         if ref_key not in ignore_key_list:
#             for ref_inner_key in inner_key_list:
#                 total_dict[ref_inner_key] += count_dict[ref_key][ref_inner_key]
#
#     count_dict["total"] = total_dict.copy()
#
#     return count_dict
#
#
# def calculate_bg_scale_dict(meta_data_dict, to_large_state=False, scale_reads_count=None):
#     """
#     INPUT:
#         <meta_data_dict>
#             A dict generated from raw BAM file, and data structure like
#                 meta_count_dict = {
#                     "meta_data" : {
#                         "create_date" : ,
#                         "ctrl_BAM" : ,
#                         "treat_BAM" : ,
#                         "reference" : ,
#                         "threads" : ,
#                         "query_mutation_type" : ,
#                         "query_mutation_min_cutoff" : ,
#                         "query_mutation_max_cutoff" : ,
#                         "other_mismatch_max_cutoff" : ,
#                         "total_mismatch_max_cutoff" : ,
#                         "filter_log" :
#                     },
#                     "ctrl" : ,
#                     "treat" : ,
#                 }
#
#         <to_large_state>
#             Default is False
#                 Means scale the larger sample up to the smaller sample.
#                 Usually, scaling down will bring down background noise, which is good for enrichment detection.
#
#     RETURN:
#         <scale_dict>
#
#
#     INFO:
#         Final changes, 2020-11-05
#     """
#
#     # init vars
#     ctrl_scale_factor = 1
#     treat_scale_factor = 1
#     scale_dict = {"ctrl": {}, "treat": {}}
#
#     select_key_list = [
#         "all_align_count",
#         "region_mut_count",
#         "region_non_mut_count"
#     ]
#
#     if scale_reads_count == 0:
#         raise ValueError("scale_reads_count can't be 0!")
#
#     elif scale_reads_count is not None:
#         ctrl_scale_factor = meta_data_dict["ctrl"]["total"]["all_align_count"] / 1.0 / int(scale_reads_count)
#         treat_scale_factor = meta_data_dict["treat"]["total"]["all_align_count"] / 1.0 / int(scale_reads_count)
#
#     for ref_name in meta_data_dict["ctrl"].keys():
#         # make all count scale factor
#         ctrl_all_count = meta_data_dict["ctrl"][ref_name]["all_align_count"]
#         treat_all_count = meta_data_dict["treat"][ref_name]["all_align_count"]
#
#         if to_large_state:
#             all_scale_count = max(ctrl_all_count, treat_all_count)
#         else:
#             all_scale_count = min(ctrl_all_count, treat_all_count)
#
#         if all_scale_count == 0:
#             raise ValueError("%s all_align_count is 0!!!" % ref_name)
#
#         ctrl_all_scale_factor = ctrl_all_count / 1.0 / all_scale_count
#         treat_all_scale_factor = treat_all_count / 1.0 / all_scale_count
#
#         # make dict
#         if scale_dict["ctrl"].get(ref_name) is None:
#             scale_dict["ctrl"][ref_name] = {}
#
#         if scale_dict["treat"].get(ref_name) is None:
#             scale_dict["treat"][ref_name] = {}
#
#         for key in meta_data_dict["ctrl"][ref_name].keys():
#             ctrl_count = meta_data_dict["ctrl"][ref_name][key]
#             treat_count = meta_data_dict["treat"][ref_name][key]
#
#             if scale_reads_count is None:
#                 if to_large_state:
#                     scale_value = max(ctrl_count, treat_count)
#                 else:
#                     scale_value = min(ctrl_count, treat_count)
#
#                 if scale_value == 0:
#                     if key in select_key_list:
#                         raise ValueError("%s %s scale_value is 0!!!" % (ref_name, key))
#
#                     else:
#                         sys.stderr.write("Warning %s %s scale_value is 0!!!\n" % (ref_name, key))
#                         sys.stderr.write("Use all_align_count scale info instead." + "\n")
#                         ctrl_scale_factor = ctrl_all_scale_factor
#                         treat_scale_factor = treat_all_scale_factor
#
#                 else:
#                     ctrl_scale_factor = ctrl_count / 1.0 / scale_value
#                     treat_scale_factor = treat_count / 1.0 / scale_value
#
#             # make new dict key
#             if key in select_key_list:
#                 scale_dict_key = key + ".scale_factor"
#                 scale_dict["ctrl"][ref_name][scale_dict_key] = ctrl_scale_factor
#                 scale_dict["treat"][ref_name][scale_dict_key] = treat_scale_factor
#
#     # return part
#     return scale_dict
#
#
# def calculate_effective_genome(ref_genome_filename, log_verbose=3):
#     """
#     INPUT:
#         <ref_genome_filename>
#             str, ref FASTA filename
#
#     RETURN:
#         <genome_base_count_dict>
#             dict, key is chr1, chr2 ... and total
#                   value is A,G,C,T,N count
#     """
#     # log setting
#     logging.basicConfig(level=(4 - log_verbose) * 10,
#                         format='%(levelname)-5s @ %(asctime)s: %(message)s ',
#                         datefmt='%Y-%m-%d %H:%M:%S',
#                         stream=sys.stderr,
#                         filemode="w")
#
#     logging.info("Starting to calculate genome effective length...")
#
#     # init vars
#     genome_base_count_dict = {"total": {"A": 0, "T": 0, "C": 0, "G": 0, "N": 0}}
#
#     # open genome
#     genome_fa = SeqIO.parse(handle=ref_genome_filename, format="fasta")
#
#     for ref in genome_fa:
#         logging.debug("Run on %s ..." % ref.id)
#
#         chr_seq = ref.seq.upper()
#
#         if genome_base_count_dict.get(ref.id) is None:
#             genome_base_count_dict[ref.id] = {}
#
#         # add into dict
#         genome_base_count_dict[ref.id]["A"] = chr_seq.count("A")
#         genome_base_count_dict[ref.id]["T"] = chr_seq.count("T")
#         genome_base_count_dict[ref.id]["C"] = chr_seq.count("C")
#         genome_base_count_dict[ref.id]["G"] = chr_seq.count("G")
#         genome_base_count_dict[ref.id]["N"] = chr_seq.count("N")
#
#         # add to total
#         genome_base_count_dict["total"]["A"] += genome_base_count_dict[ref.id]["A"]
#         genome_base_count_dict["total"]["T"] += genome_base_count_dict[ref.id]["T"]
#         genome_base_count_dict["total"]["C"] += genome_base_count_dict[ref.id]["C"]
#         genome_base_count_dict["total"]["G"] += genome_base_count_dict[ref.id]["G"]
#         genome_base_count_dict["total"]["N"] += genome_base_count_dict[ref.id]["N"]
#
#     # log
#     logging.info("Calculating genome effective length, done!")
#
#     return genome_base_count_dict
#
#
# def back_all_normalization_scale_dict(meta_data_dict, genome_base_count_dict, seq_reads_length=150,
#                                       norm_scale_reads_count=1e6):
#     """
#     INPUT:
#         <meta_dict>
#             dict, contain genome wide alignment value info, that is a combination of ctrl and treat result
#                   generated from FUN multi_thread_chr_bam_mut_count
#
#         <genome_base_count_dict>
#             dict, generate from FUN calculate_effective_genome
#
#         <norm_scale_reads_count>
#             int, if set <scale_reads_count> is 1e6, the normalization value equals to CPM.
#
#     RETURN:
#         <scale_factor_dict>
#
#         <normalize_scale_factor_dict>
#
#         <genome_bg_dict>
#
#     INFO
#         Final changes, 2020-11-06
#     """
#     # -------------------------------------------------->>>>>>>>>>
#     # init vars
#     # -------------------------------------------------->>>>>>>>>>
#     genome_bg_dict = {"ctrl": {
#         "scale_mut_bg": {},
#         "norm_scale_mut_bg": {},
#         "scale_all_bg": {},
#         "norm_scale_all_bg": {}
#     }, "treat": {
#         "scale_mut_bg": {},
#         "norm_scale_mut_bg": {},
#         "scale_all_bg": {},
#         "norm_scale_all_bg": {}
#     }}
#
#     # -------------------------------------------------->>>>>>>>>>
#     # make scale factor dict
#     # -------------------------------------------------->>>>>>>>>>
#     scale_factor_dict = calculate_bg_scale_dict(meta_data_dict, to_large_state=False)
#
#     # -------------------------------------------------->>>>>>>>>>
#     # make normalize scale factor
#     # -------------------------------------------------->>>>>>>>>>
#     normalize_scale_factor_dict = calculate_bg_scale_dict(meta_data_dict,
#                                                           to_large_state=False,
#                                                           scale_reads_count=int(norm_scale_reads_count))
#
#     # -------------------------------------------------->>>>>>>>>>
#     # genome background
#     # -------------------------------------------------->>>>>>>>>>
#     genome_effective_len = 0
#     for dna_base in "AGCT":
#         genome_effective_len += genome_base_count_dict["total"][dna_base]
#
#     # genome level background
#     ctrl_scale_factor = scale_factor_dict["ctrl"]["total"]["all_align_count.scale_factor"]
#     treat_scale_factor = scale_factor_dict["ctrl"]["total"]["all_align_count.scale_factor"]
#
#     ctrl_norm_scale_factor = normalize_scale_factor_dict["ctrl"]["total"]["all_align_count.scale_factor"]
#     treat_norm_scale_factor = normalize_scale_factor_dict["treat"]["total"]["all_align_count.scale_factor"]
#
#     # genome level mutation background
#     ctrl_mut_bg_genome_lambda = meta_data_dict["ctrl"]["total"][
#                                     "region_mut_count"] * seq_reads_length / 1.0 / ctrl_scale_factor / genome_effective_len
#     treat_mut_bg_genome_lambda = meta_data_dict["treat"]["total"][
#                                      "region_mut_count"] * seq_reads_length / 1.0 / treat_scale_factor / genome_effective_len
#
#     ctrl_norm_mut_bg_genome_lambda = meta_data_dict["ctrl"]["total"][
#                                          "region_mut_count"] * seq_reads_length / 1.0 / ctrl_norm_scale_factor / genome_effective_len
#     treat_norm_mut_bg_genome_lambda = meta_data_dict["treat"]["total"][
#                                           "region_mut_count"] * seq_reads_length / 1.0 / treat_norm_scale_factor / genome_effective_len
#
#     genome_bg_dict["ctrl"]["scale_mut_bg"]["genome_bg"] = ctrl_mut_bg_genome_lambda
#     genome_bg_dict["treat"]["scale_mut_bg"]["genome_bg"] = treat_mut_bg_genome_lambda
#
#     genome_bg_dict["ctrl"]["norm_scale_mut_bg"]["genome_bg"] = ctrl_norm_mut_bg_genome_lambda
#     genome_bg_dict["treat"]["norm_scale_mut_bg"]["genome_bg"] = treat_norm_mut_bg_genome_lambda
#
#     # genome level all align background
#     ctrl_all_bg_genome_lambda = meta_data_dict["ctrl"]["total"][
#                                     "all_align_count"] * seq_reads_length / 1.0 / ctrl_scale_factor / genome_effective_len
#     treat_all_bg_genome_lambda = meta_data_dict["treat"]["total"][
#                                      "all_align_count"] * seq_reads_length / 1.0 / treat_scale_factor / genome_effective_len
#
#     ctrl_norm_all_bg_genome_lambda = meta_data_dict["ctrl"]["total"][
#                                          "all_align_count"] * seq_reads_length / 1.0 / ctrl_norm_scale_factor / genome_effective_len
#     treat_norm_all_bg_genome_lambda = meta_data_dict["treat"]["total"][
#                                           "all_align_count"] * seq_reads_length / 1.0 / treat_norm_scale_factor / genome_effective_len
#
#     genome_bg_dict["ctrl"]["scale_all_bg"]["genome_bg"] = ctrl_all_bg_genome_lambda
#     genome_bg_dict["treat"]["scale_all_bg"]["genome_bg"] = treat_all_bg_genome_lambda
#
#     genome_bg_dict["ctrl"]["norm_scale_all_bg"]["genome_bg"] = ctrl_norm_all_bg_genome_lambda
#     genome_bg_dict["treat"]["norm_scale_all_bg"]["genome_bg"] = treat_norm_all_bg_genome_lambda
#
#     # -------------------------------------------------->>>>>>>>>>
#     # chromosome background
#     # -------------------------------------------------->>>>>>>>>>
#     for ref_name in meta_data_dict["ctrl"].keys():
#         # chr effective length
#         chr_effective_len = 0
#         for dna_base in "AGCT":
#             chr_effective_len += genome_base_count_dict[ref_name][dna_base]
#
#         # chr scale factor
#         chr_ctrl_scale_factor = scale_factor_dict["ctrl"][ref_name]["all_align_count.scale_factor"]
#         chr_treat_scale_factor = scale_factor_dict["treat"][ref_name]["all_align_count.scale_factor"]
#
#         chr_ctrl_norm_scale_factor = normalize_scale_factor_dict["ctrl"][ref_name]["all_align_count.scale_factor"]
#         chr_treat_norm_scale_factor = normalize_scale_factor_dict["treat"][ref_name]["all_align_count.scale_factor"]
#
#         # chr mutation background
#         ctrl_mut_bg_chr_lambda = meta_data_dict["ctrl"][ref_name][
#                                      "region_mut_count"] * seq_reads_length / 1.0 / chr_ctrl_scale_factor / chr_effective_len
#         treat_mut_bg_chr_lambda = meta_data_dict["treat"][ref_name][
#                                       "region_mut_count"] * seq_reads_length / 1.0 / chr_treat_scale_factor / chr_effective_len
#
#         norm_ctrl_mut_bg_chr_lambda = meta_data_dict["ctrl"][ref_name][
#                                           "region_mut_count"] * seq_reads_length / 1.0 / chr_ctrl_norm_scale_factor / chr_effective_len
#         norm_treat_mut_bg_chr_lambda = meta_data_dict["treat"][ref_name][
#                                            "region_mut_count"] * seq_reads_length / 1.0 / chr_treat_norm_scale_factor / chr_effective_len
#
#         # store in dict
#         genome_bg_dict["ctrl"]["scale_mut_bg"][ref_name] = ctrl_mut_bg_chr_lambda
#         genome_bg_dict["treat"]["scale_mut_bg"][ref_name] = treat_mut_bg_chr_lambda
#
#         genome_bg_dict["ctrl"]["norm_scale_mut_bg"][ref_name] = norm_ctrl_mut_bg_chr_lambda
#         genome_bg_dict["treat"]["norm_scale_mut_bg"][ref_name] = norm_treat_mut_bg_chr_lambda
#
#         # chr all background
#         ctrl_all_bg_chr_lambda = meta_data_dict["ctrl"][ref_name][
#                                      "all_align_count"] * seq_reads_length / 1.0 / chr_ctrl_scale_factor / chr_effective_len
#         treat_all_bg_chr_lambda = meta_data_dict["treat"][ref_name][
#                                       "all_align_count"] * seq_reads_length / 1.0 / chr_treat_scale_factor / chr_effective_len
#
#         norm_ctrl_all_bg_chr_lambda = meta_data_dict["ctrl"][ref_name][
#                                           "all_align_count"] * seq_reads_length / 1.0 / chr_ctrl_norm_scale_factor / chr_effective_len
#         norm_treat_all_bg_chr_lambda = meta_data_dict["treat"][ref_name][
#                                            "all_align_count"] * seq_reads_length / 1.0 / chr_treat_norm_scale_factor / chr_effective_len
#
#         # store in dict
#         genome_bg_dict["ctrl"]["scale_all_bg"][ref_name] = ctrl_all_bg_chr_lambda
#         genome_bg_dict["treat"]["scale_all_bg"][ref_name] = treat_all_bg_chr_lambda
#
#         genome_bg_dict["ctrl"]["norm_scale_all_bg"][ref_name] = norm_ctrl_all_bg_chr_lambda
#         genome_bg_dict["treat"]["norm_scale_all_bg"][ref_name] = norm_treat_all_bg_chr_lambda
#
#     return scale_factor_dict, normalize_scale_factor_dict, genome_bg_dict
#
#
#
# def _log_cmd_str(args):
#     full_cmd_str = f"""python find-significant-mpmat.py
#     --mpmat_table {args["mpmat_table"]}
#     --output {args["output"]}
#     --ctrl_BAM {args["ctrl_BAM"]}
#     --treat_BAM {args["treat_BAM"]}
#     --reference {args["reference"]}
#     --thread {args["thread"]}
#     --query_mutation_type {args["query_mutation_type"]}
#     --verbose {args["verbose"]}
#     --keep_temp_file {args["keep_temp_file"]}
#     --mpmat_filter_info_col_index {args["mpmat_filter_info_col_index"]}
#     --mpmat_block_info_col_index {args["mpmat_block_info_col_index"]}
#     --region_block_mut_num_cutoff {args["region_block_mut_num_cutoff"]}
#     --query_mut_min_cutoff {args["query_mut_min_cutoff"]}
#     --query_mut_max_cutoff {args["query_mut_max_cutoff"]}
#     --total_mut_max_cutoff {args["total_mut_max_cutoff"]}
#     --other_mut_max_cutoff {args["other_mut_max_cutoff"]}
#     --seq_reads_length {args["seq_reads_length"]}
#     --scale_reads_count {args["scale_reads_count"]}
#     --lambda_method {args["lambda_method"]}
#     --poisson_method {args["poisson_method"]}
#     """
#
#     return full_cmd_str
#
#
# def _check_file_exist(args):
#     """
#     INPUT:
#         <args>
#             obj, argparse obj
#
#     RETURN:
#         <check_all_state>
#             bool, True means all files are exist, False means at least one of files can't pass file check step.
#
#         <check_exist_dict>
#             dict, each item contain 3 elements:
#                 1.input filename
#                 2.check state
#                 3. check reason
#     """
#     # init list
#     check_all_state = True
#     check_exist_dict = {}
#
#     # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#     # 1. check ctrl bam file
#     # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#     ctrl_check_res = check_input_bam_file(os.path.abspath(args.ctrl_BAM))
#
#     if ctrl_check_res == 1:
#         check_exist_dict["ctrl_BAM"] = [args['ctrl_BAM'], False, "BAM not exist!"]
#         check_all_state = False
#
#     elif ctrl_check_res == 2:
#         check_exist_dict["ctrl_BAM"] = [args['ctrl_BAM'], False, "BAM not sorted by coordinate!"]
#         check_all_state = False
#
#     elif ctrl_check_res == 3:
#         check_exist_dict["ctrl_BAM"] = [args['ctrl_BAM'], False, "BAM doesn't contain index file!"]
#         check_all_state = False
#
#     else:
#         check_exist_dict["ctrl_BAM"] = [args['ctrl_BAM'], True, None]
#
#     # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#     # 2. check PD bam file
#     # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#     treat_check_res = check_input_bam_file(os.path.abspath(args.treat_BAM))
#
#     if treat_check_res == 1:
#         check_exist_dict["treat_BAM"] = [args['treat_BAM'], False, "BAM not exist!"]
#         check_all_state = False
#
#     elif treat_check_res == 2:
#         check_exist_dict["treat_BAM"] = [args['treat_BAM'], False, "BAM not sorted by coordinate!"]
#         check_all_state = False
#
#     elif treat_check_res == 3:
#         check_exist_dict["treat_BAM"] = [args['treat_BAM'], False, "BAM doesn't contain index file!"]
#         check_all_state = False
#
#     else:
#         check_exist_dict["treat_BAM"] = [args['treat_BAM'], True, None]
#
#     # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#     # 3. check other file exist
#     # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#     # mpmat_table
#     if not os.path.exists(os.path.abspath(args['mpmat_table'])):
#         check_exist_dict["mpmat_table"] = [args['mpmat_table'], False, "--mpmat_table does not exist!"]
#         check_all_state = False
#     else:
#         check_exist_dict["mpmat_table"] = [args['mpmat_table'], True, None]
#
#     # reference
#     if not os.path.exists(os.path.abspath(args['reference'])):
#         check_exist_dict["reference"] = [args['reference'], False, "--reference does not exist!"]
#         check_all_state = False
#     else:
#         check_exist_dict["reference"] = [args['reference'], True, None]
#
#     # make log
#     logging.basicConfig(level=10,
#                         format='%(levelname)-5s @ %(asctime)s: %(message)s ',
#                         datefmt='%Y-%m-%d %H:%M:%S',
#                         stream=sys.stderr,
#                         filemode="w")
#
#     for file_type in ["mpmat_table", "ctrl_BAM", "treat_BAM", "reference"]:
#         if check_exist_dict[file_type][1]:
#             logging.info("Check exist \n \t--%s %s" % (file_type, check_exist_dict[file_type][0]))
#             logging.info("Yes!")
#         else:
#             logging.info("Check exist \n \t--%s %s" % (file_type, check_exist_dict[file_type][0]))
#             logging.info("No! Reason: %s" % (check_exist_dict[file_type][2]))
#
#         sys.stderr.write("-" * 80 + "\n")
#
#     # return part
#     return check_all_state, check_exist_dict
#
#
# header_list = [
#     "chr_name",
#     "region_start",
#     "region_end",
#     "mpmat_index",
#     "region_site_num",
#     "region_block_site_num",
#     "region_mut_site_num",
#     "region_site_index",
#     "region_block_state",
#     "region_highest_site_index",
#     "region_highest_site_mut_num",
#     "region_highest_site_cover_num",
#     "region_highest_site_mut_ratio",
#     "ctrl_count",
#     "treat_count",
#     "ctrl_mut_count",
#     "treat_mut_count",
#     "ctrl_count.norm",
#     "treat_count.norm",
#     "ctrl_mut_count.norm",
#     "treat_mut_count.norm",
#     "count_info",
#     "log2_FC",
#     "log2_FC_mut",
#     "test_state",
#     "p_value",
#     "FDR"
# ]
