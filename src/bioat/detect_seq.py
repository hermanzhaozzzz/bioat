# import os
# import sys
# import logging
# from bioat.logger import set_logging_level
# from bioat.lib.libdetect_seq import (
#     check_input_bam_file,
#     clear_temp_files_by_dict,
#     split_mpmat_by_chr_name,
#     make_qvalue_with_BH_method,
#     multi_thread_chr_bam_mut_count,
#     calculate_effective_genome,
#     count_dict_sum_up,
#     back_all_normalization_scale_dict,
#     multi_run_mpmat_poisson_test,
#     merge_split_files,
#     _log_cmd_str,
#     _check_file_exist
# )
#
#
# class DetectSeq:
#     def __init__(self):
#         pass
#
#     def find_significant_mpmat(
#             self,
#             mpmat_table: str,
#             ctrl_BAM: str,
#             treat_BAM: str,
#             reference: str,
#             output: str = "./poisson_output.tsv",
#             thread: int = 1,
#             query_mutation_type: str = "CT,GA",
#             keep_temp_file=False,
#             mpmat_filter_info_col_index: int = -1,
#             mpmat_block_info_col_index: int = -1,
#             region_block_mut_num_cutoff: int = 2,
#             query_mut_min_cutoff: int = 1,
#             query_mut_max_cutoff: int = 16,
#             total_mut_max_cutoff: int = 20,
#             other_mut_max_cutoff: int = 22,
#             seq_reads_length: int = 150,
#             scale_reads_count: int = 1e6,
#             lambda_method: str = "ctrl_max",
#             poisson_method: str = "mutation",
#             log_level: str = 'INFO'
#     ):
#         """
#         A tool to find significant enrichment mpmat region. You can set a lot of parameters with this program,
#         usually default parameters can work well.
#
#         :param mpmat_table: .mpmat table, which can be generated from <pmat-merge.py> code
#         :param ctrl_BAM: Control BAM file
#         :param treat_BAM: Treatment BAM file
#         :param reference: Genome FASTA file
#         :param output: Output Poisson test result, default=poisson_output.tsv
#
#         :param thread: Number of threads used to process data
#         :param query_mutation_type: Query mutation type, which will be considered as mutation signals, default=CT,GA
#         :param keep_temp_file: If keep temp files, default=False
#         :param mpmat_filter_info_col_index: Column index for filter info, -1 means no such info. Default=-1. Usually can set as 13
#         :param mpmat_block_info_col_index: Column index for region block info, -1 means no such info. Default=-1.
#         :param region_block_mut_num_cutoff: Site filter cutoff, if a site has a mutation signal >= this cutoff in ctrl sample, the site will be blocked in the downstream analysis. Default=2
#         :param query_mut_min_cutoff: An alignment contains query mutation count lower than this cutoff in mpmat region will considered as 'non_mut_align', default=1
#         :param query_mut_max_cutoff: An alignment contain mutation number higher than this, considered as 'query_high_mismatch', which often caused by Bismark mapping error, default=16
#         :param total_mut_max_cutoff: An alignment contain mutation number higher than this, considered as 'total_high_mismatch', default=20
#         :param other_mut_max_cutoff: An alignment contain mutation number higher than this, considered as 'other_high_mismatch', default=12
#         :param seq_reads_length: Sequencing reads length, default=150
#         :param scale_reads_count: Scaled final output region signal, default=1000000 Usually, default=1e6, means CPM
#         :param lambda_method: Can be set as 'ctrl_max', 'treat_max', 'max', 'raw', default=ctrl_max
#         :param poisson_method: Can be set as 'mutation' OR 'all', default=mutation. 'mutation' means only use mutation alignments to run Poisson test,'all' means use all alignments to run Poisson, which similar to MACS2
#         :param log_level: 'CRITICAL', 'ERROR', 'WARNING', 'INFO', 'DEBUG', 'NOTSET'
#         :return:
#         """
#         set_logging_level(level=log_level)
#         # set logger
#         lib_name = __name__
#         function_name = sys._getframe().f_code.co_name
#         logger = logging.getLogger(f'{lib_name}.{function_name} ==> ')
#
#         ARGS = dict(
#             mpmat_table=mpmat_table,
#             output=output,
#             ctrl_BAM=ctrl_BAM,
#             treat_BAM=treat_BAM,
#             reference=reference,
#             thread=thread,
#             query_mutation_type=query_mutation_type,
#             verbose=log_level,
#             keep_temp_file=keep_temp_file,
#             mpmat_filter_info_col_index=mpmat_filter_info_col_index,
#             mpmat_block_info_col_index=mpmat_block_info_col_index,
#             region_block_mut_num_cutoff=region_block_mut_num_cutoff,
#             query_mut_min_cutoff=query_mut_min_cutoff,
#             query_mut_max_cutoff=query_mut_max_cutoff,
#             total_mut_max_cutoff=total_mut_max_cutoff,
#             other_mut_max_cutoff=other_mut_max_cutoff,
#             seq_reads_length=seq_reads_length,
#             scale_reads_count=scale_reads_count,
#             lambda_method=lambda_method,
#             poisson_method=poisson_method
#         )
#
#         # output full cmd
#         sys.stderr.write("\n" + "-" * 80 + "\n")
#         sys.stderr.write(_log_cmd_str(args=ARGS))
#
#         # check file exist
#         sys.stderr.write("\n" + "-" * 80 + "\n")
#         check_res = _check_file_exist(args=ARGS)
#         if not check_res[0]:
#             raise IOError("Please make sure each file is exist!")
#
#         # check output file dir
#         if output == "stdout":
#             temp_dir = os.getcwd()
#         else:
#             temp_dir = os.path.dirname(output)
#
#         # fix ARGS
#         query_mutation_type = query_mutation_type.split(",")
#
#         # ---------------------------------------------------------------------------->>>>>>>>>>
#         # Part I split mpmat
#         # ---------------------------------------------------------------------------->>>>>>>>>>
#         # split
#         sys.stderr.write("\n" + "-" * 80 + "\n")
#         split_mpmat_filename_dict = split_mpmat_by_chr_name(
#             input_mpmat_filename=mpmat_table,
#             temp_dir=temp_dir,
#             force_temp_dir=True)
#
#         # meta dict
#         meta_dict = {"filename": {}}
#         meta_dict["filename"]["mpmat_region_AddBlockInfo"] = mpmat_table
#         meta_dict["filename"]["mpmat_region_AddBlockInfo_split"] = split_mpmat_filename_dict.copy()
#
#         # make chr order
#         ref_order_dict = {}
#         ref_order_list = []
#
#         for index, chr_name in enumerate(meta_dict["filename"]["mpmat_region_AddBlockInfo_split"]["chr_name_order"]):
#             ref_order_dict[chr_name] = index
#             ref_order_list.append(chr_name)
#
#         # ---------------------------------------------------------------------------->>>>>>>>>>
#         # Part II calculation of normalization scale factor
#         # ---------------------------------------------------------------------------->>>>>>>>>>
#         sys.stderr.write("-" * 80 + "\n")
#         logger.info("Calculating ctrl sample Poisson background...")
#
#         ctrl_poisson_bg_meta_dict = multi_thread_chr_bam_mut_count(
#             input_bam_filename=ARGS.ctrl_BAM,
#             select_chr_name_list=ref_order_list,
#             ref_genome_filename=ARGS.reference,
#             thread_num=ARGS.thread,
#             query_mut_type_list=query_mutation_type,
#             query_mut_min_cutoff=ARGS.query_mut_min_cutoff,
#             query_mut_max_cutoff=ARGS.query_mut_max_cutoff,
#             total_mut_max_cutoff=ARGS.total_mut_max_cutoff,
#             other_mut_max_cutoff=ARGS.other_mut_max_cutoff,
#             log_verbose=ARGS.verbose)
#
#         sys.stderr.write("-" * 80 + "\n")
#         logger.info("Calculating treat sample Poisson background...")
#
#         treat_poisson_bg_meta_dict = multi_thread_chr_bam_mut_count(
#             input_bam_filename=ARGS.treat_BAM,
#             select_chr_name_list=ref_order_list,
#             ref_genome_filename=ARGS.reference,
#             thread_num=ARGS.thread,
#             query_mut_type_list=query_mutation_type,
#             query_mut_min_cutoff=ARGS.query_mut_min_cutoff,
#             query_mut_max_cutoff=ARGS.query_mut_max_cutoff,
#             total_mut_max_cutoff=ARGS.total_mut_max_cutoff,
#             other_mut_max_cutoff=ARGS.other_mut_max_cutoff,
#             log_verbose=ARGS.verbose)
#
#         # genome effective length
#         sys.stderr.write("-" * 80 + "\n")
#         genome_base_count_dict = calculate_effective_genome(ref_genome_filename=ARGS.reference,
#                                                             log_verbose=ARGS.verbose)
#
#         # make meta data dict
#         meta_data_dict = {"ctrl": count_dict_sum_up(ctrl_poisson_bg_meta_dict["count_dict"]),
#                           "treat": count_dict_sum_up(treat_poisson_bg_meta_dict["count_dict"])}
#
#         # make scale dict
#         scale_factor_dict, normalize_scale_factor_dict, genome_bg_dict = back_all_normalization_scale_dict(
#             meta_data_dict=meta_data_dict,
#             genome_base_count_dict=genome_base_count_dict,
#             seq_reads_length=ARGS.seq_reads_length,
#             norm_scale_reads_count=int(ARGS.scale_reads_count))
#
#         # ---------------------------------------------------------------------------->>>>>>>>>>
#         # Part III Poisson test on mpmat region
#         # ---------------------------------------------------------------------------->>>>>>>>>>
#         poisson_out_split_dict = multi_run_mpmat_poisson_test(
#             mpmat_block_split_dict=meta_dict["filename"]["mpmat_region_AddBlockInfo_split"],
#             ctrl_bam_filename=ARGS.ctrl_BAM,
#             treat_bam_filename=ARGS.treat_BAM,
#             ref_genome_fa_filename=ARGS.reference,
#             scale_factor_dict=scale_factor_dict,
#             normalize_scale_factor_dict=normalize_scale_factor_dict,
#             genome_bg_dict=genome_bg_dict,
#             lambda_bg_method=ARGS.lambda_method,
#             poisson_method=ARGS.poisson_method,
#             log_verbose=ARGS.verbose,
#             thread=ARGS.thread,
#             region_block_mut_num_cutoff=ARGS.region_block_mut_num_cutoff,
#             reads_query_mut_min_cutoff=ARGS.query_mut_min_cutoff,
#             reads_query_mut_max_cutoff=ARGS.query_mut_max_cutoff,
#             reads_total_mut_max_cutoff=ARGS.total_mut_max_cutoff,
#             reads_other_mut_max_cutoff=ARGS.other_mut_max_cutoff,
#             mpmat_filter_info_col_index=ARGS.mpmat_filter_info_col_index,
#             mpmat_block_info_col_index=ARGS.mpmat_block_info_col_index
#         )
#
#         meta_dict["filename"]["Poisson_out_split"] = poisson_out_split_dict.copy()
#         meta_dict["filename"]["Poisson_out_merge_NoFDR"] = ARGS.output + ".NoFDR"
#
#         # ---------------------------------------------------------------------------->>>>>>>>>>
#         # Part VI output Poisson test result
#         # ---------------------------------------------------------------------------->>>>>>>>>>
#         sys.stderr.write("-" * 80 + "\n")
#         merge_state, pval_list_str = merge_split_files(split_file_dict=meta_dict["filename"]["Poisson_out_split"],
#                                                        key_order_list=ref_order_list,
#                                                        out_filename=ARGS.output + ".NoFDR",
#                                                        header_list=None,
#                                                        in_sep="\t", out_sep="\t",
#                                                        log_verbose=ARGS.verbose,
#                                                        return_col_index=-1)
#
#         sys.stderr.write("-" * 80 + "\n")
#         logger.info("Calculating FDR value and make final output table...")
#
#         # change pval into float
#         pval_list = []
#         for pval_str in pval_list_str:
#             if pval_str != "NA":
#                 pval_list.append(eval(pval_str))
#             else:
#                 pval_list.append("NA")
#
#         # make qval with BH method
#         qval_list = make_qvalue_with_BH_method(pval_list)
#
#         # make final output
#         final_out_file = open(ARGS.output, "wt")
#         final_out_file.write("\t".join(header_list) + "\n")
#
#         with open(meta_dict["filename"]["Poisson_out_merge_NoFDR"], "rt") as no_fdr_file:
#             for index, line in enumerate(no_fdr_file):
#                 line_list = line.strip().split("\t")
#                 line_list.append(str(qval_list[index]))
#                 final_out_file.write("\t".join(line_list) + "\n")
#
#         final_out_file.close()
#
#         logger.info("Done!")
#
#         # ---------------------------------------------------------------------------->>>>>>>>>>
#         # Part VII remove temp files
#         # ---------------------------------------------------------------------------->>>>>>>>>>
#         if ARGS.keep_temp_file == "False":
#             sys.stderr.write("-" * 80 + "\n")
#             logger.info("Set <keep_temp_file> as False... Removing temp files...")
#
#             sys.stderr.write("-" * 80 + "\n")
#             logger.info("Removing mpmat split files...")
#             clear_temp_files_by_dict(temp_file_dict=meta_dict["filename"]["mpmat_region_AddBlockInfo_split"],
#                                      log_verbose=ARGS.verbose)
#
#             sys.stderr.write("-" * 80 + "\n")
#             logger.info("Removing Poisson test split files...")
#             clear_temp_files_by_dict(temp_file_dict=meta_dict["filename"]["Poisson_out_split"],
#                                      log_verbose=ARGS.verbose)
#
#             sys.stderr.write("-" * 80 + "\n")
#             logger.info("Removing Poisson test raw file...")
#             try:
#                 os.remove(meta_dict["filename"]["Poisson_out_merge_NoFDR"])
#             except Warning as w:
#                 logger.warning("Removing file error: \n\t%s" % meta_dict["filename"]["Poisson_out_merge_NoFDR"])
#
#         logger.info("Everything done!")
#
#
# # ---------------------------------------------------------------------------->>>>>>>>>>
# #  main part
# # ---------------------------------------------------------------------------->>>>>>>>>>
# if __name__ == '__main__':
#     main()
