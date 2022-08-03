import argparse
import logging
import os
import sys

import numpy as np
import pandas as pd

"""initialize"""
pd.set_option("max_colwidth", 60)  # column最大宽度
pd.set_option("display.width", 150)  # dataframe宽度
pd.set_option("display.max_columns", None)  # column最大显示数
pd.set_option("display.max_rows", 100)  # row最大显示数

dt_type = {
    "chr_name": str,
    "chr_index": str,
    "ref_base": str,
    "A": int,
    "G": int,
    "C": int,
    "T": int,
    "del_count": int,
    "insert_count": int,
    "ambiguous_count": int,
    "deletion": str,
    "insertion": str,
    "ambiguous": str,
    "mut_num": int,
}


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
        "--bmat_folder",
        type=str,
        required=True,
        help="""
            All the bmat file should be formated as:
            TargetSeq-TREATMENT_x_REP_x_REGIONID_x_CUTOFF_x.bmat.gz
            and put in this folder.
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
    return parser.parse_args()


def return_base_ratio(x):
    if x["ref_base"] == "C":
        return x["T"] / x["depth"]
    elif x["ref_base"] == "G":
        return x["A"] / x["depth"]
    else:
        return np.NaN


if __name__ == "__main__":
    # logging_setting
    logging_setting()
    # argparse_setting
    args = argparse_setting()
    pt_bmat = args.bmat_folder
    pt_out = args.out
    ls_bmat = [os.path.join(pt_bmat, i) for i in os.listdir(pt_bmat) if "bmat" in i]
    # # TargetSeq-TREATMENT_x_REP_x_REGIONID_x_CUTOFF_x.bmat.gz

    # """ load bmat and merge """
    df_merge = pd.DataFrame()
    for path in ls_bmat:
        df = pd.read_csv(path, header=0, index_col=False, sep="\t", dtype=dt_type)
        df["bmat_name"] = path.replace("bmat/", "")
        df_merge = pd.concat([df_merge, df], ignore_index=True)

    # fix region_id for ND4 ND5.1 ND6
    # df_merge.chr_name = df_merge.chr_name.map(lambda x: x.replace("SHARE", "share"))
    # df_merge.chr_name = df_merge.chr_name.map(lambda x: x.replace("ONLY", "only"))
    # """ fix rep info """
    # remove no need cols
    df = df_merge[
        [
            "chr_name",
            "chr_index",
            "ref_base",
            "A",
            "G",
            "C",
            "T",
            "mut_num",
            "bmat_name",
        ]
    ].copy()

    # form rep info
    df["rep"] = df.bmat_name.map(
        lambda x: "rep" + x.split("_REP_")[1].split("_CUTOFF_")[0]
    )
    # form treatment info
    df["treatment"] = df.bmat_name.map(
        lambda x: x.split("-TREATMENT_")[1].split("_REP_")[0]
    )
    # form cutoff info
    df["cutoff"] = df.bmat_name.map(
        lambda x: int(x.split("_CUTOFF_")[1].split("_REGIONID_")[0])
    )

    if df.isna().sum().sum() != 0:
        raise ValueError("bmat file name should be formated")

    # select
    df_select = df[
        [
            "chr_name",
            "chr_index",
            "ref_base",
            "A",
            "G",
            "C",
            "T",
            "bmat_name",
            "rep",
            "treatment",
            "cutoff",
        ]
    ].copy()

    # df_select
    df_select["depth"] = df_select[["A", "G", "C", "T"]].apply(sum, axis=1)
    # 防止分母为0
    df_select.depth = df_select.depth.map(lambda x: x if x != 0 else 1)

    df_bmat_melt = df_select.melt(
        id_vars=[
            "treatment",
            "rep",
            "cutoff",
            "chr_name",
            "ref_base",
            "chr_index",
            "bmat_name",
            "depth",
        ],
        value_vars=["A", "G", "C", "T"],
        var_name="mut_base",
        value_name="mut_count",
    )

    df_bmat_melt.columns = [
        "treatment",
        "rep",
        "cutoff",
        "region_id",
        "ref_base",
        "relative_pos",
        "bmat_name",
        "total_count",
        "mut_base",
        "mut_count",
    ]
    df_bmat_melt.to_csv(pt_out, index=False, header=True)
    logging.info("script done.")
