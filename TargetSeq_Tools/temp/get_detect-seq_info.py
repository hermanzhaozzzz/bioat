import argparse
import logging
import sys

import numpy as np
import pandas as pd
import pysam

from loader import load_everybase_from_bowtie_table, load_reference_fasta_as_dict
from parse_pysam import get_query_site_mpileup_info

"""initialize"""
pd.set_option("max_colwidth", 60)  # column最大宽度
pd.set_option("display.width", 150)  # dataframe宽度
pd.set_option("display.max_columns", None)  # column最大显示数
pd.set_option("display.max_rows", 100)  # row最大显示数


def logging_setting():
    # log setting
    logging.basicConfig(
        level=logging.DEBUG,
        format="%(levelname)-5s @ %(asctime)s: %(message)s ",
        datefmt="%Y-%m-%d %H:%M:%S",
        stream=sys.stderr,
        filemode="w",
        force=True,
    )


def argparse_setting():
    # argparse setting
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--bowtie_table",
        type=str,
        required=True,
        help="""
            the result of: 
                bowtie \
                    path/genome.fa.bowtie1_index \
                    target-region.fa -f \
                    > target-region.align.tsv"
        """,
    )
    parser.add_argument(
        "--reference",
        type=str,
        required=True,
        help="""
            fasta, the same genome used to mapping bam file.
        """,
    )
    parser.add_argument(
        "--bam",
        type=str,
        required=True,
        help="""
            Such as Detect-seq Pull-down bam file.
            Must be sorted by position.
            Must have .bai index
        """,
    )
    parser.add_argument(
        "--out",
        type=str,
        default=None,
        help="""
            The path to put out csv file.
        """,
    )
    parser.add_argument(
        "--treatment",
        type=str,
        default=None,
        help="""
            Treatment of this sample
        """,
    )
    parser.add_argument(
        "--rep",
        type=str,
        default=None,
        help="""
            Replication of this sample: rep1, rep2
        """,
    )
    return parser.parse_args()


if __name__ == "__main__":
    # use logging
    logging_setting()
    # load params
    args = argparse_setting()
    pt_aligninfo = args.bowtie_table
    pt_ref = args.reference
    pt_bam = args.bam
    pt_out = args.out
    treatment = args.treatment
    rep = args.rep

    # load bam
    bam_file = pysam.AlignmentFile(pt_bam, "rb")
    # load ref obj
    genome_ref = pysam.FastaFile(filename=pt_ref)
    # load ref dict
    dt_ref = load_reference_fasta_as_dict(
        ref_fasta_path=pt_ref,
        ref_name="All",
        log_verbose=logging.DEBUG,
    )
    # load region info
    df_region = load_everybase_from_bowtie_table(
        bowtie_table=pt_aligninfo, narrow_cutoff=10, near_seq_extend=10
    )
    region_count = len(df_region.groupby("region_id"))
    now_count = 0

    dt_all_region = {}
    for idx_df_one_region, df_one_region in df_region.groupby("region_id"):
        now_count += 1
        ls_site_idx_all = df_one_region.site_idx.tolist()
        for i in range(0, len(ls_site_idx_all), 100):
            ls_site_idx = ls_site_idx_all[i : i + 100]

            dt_one_region = get_query_site_mpileup_info(
                site_idx_list=ls_site_idx,
                bam_obj=bam_file,
                genome_obj=genome_ref,
            )
            dt_all_region.update(dt_one_region)

        logging.info("Parsed %s/%s aim regions." % (now_count, region_count))

    df_detect = pd.DataFrame(dt_all_region).T
    df_detect.reset_index(inplace=True)
    df_detect.rename(columns={"index": "site_idx"}, inplace=True)
    df_detect_single_base = df_region.merge(df_detect, how="left", on=["site_idx"])

    if df_detect_single_base.isna().sum().sum() != 0:
        logging.fatal("merging file fail!!")
        raise ValueError("merging file fail!!")
    
    # 有BUG，不能用，先放着回头改！
    for col_base in df_detect_single_base[["A", "G", "C", "T"]]:
        df_detect_single_base["%s_ratio" % col_base] = (
            df_detect_single_base[col_base] / df_detect_single_base["total"]
        )
    if treatment is not None:
        treatment = treatment.strip()
        df_detect_single_base.insert(0, "treatment", treatment, allow_duplicates=False)
    else:
        df_detect_single_base.insert(0, "treatment", np.NaN, allow_duplicates=False)
    if rep is not None:
        df_detect_single_base.insert(1, "rep", rep, allow_duplicates=False)
    else:
        df_detect_single_base.insert(1, "rep", np.NaN, allow_duplicates=False)

    df_detect_single_base.to_csv(pt_out, header=True, index=False, sep=",")
