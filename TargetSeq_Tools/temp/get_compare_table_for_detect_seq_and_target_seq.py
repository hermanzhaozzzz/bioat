import logging
import math
import os
import sys

import numpy as np
import pandas as pd

"""initialize"""
pd.set_option("max_colwidth", 60)  # column最大宽度
pd.set_option("display.width", 200)  # dataframe宽度
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


if __name__ == "__main__":
    # logging_setting
    logging_setting()
    # file list path
    pt_target = "../additional_file/20220104_DdCBE-target-seq-info.csv,"
    pt_detect = "../additional_file/20220104_DdCBE-detect-seq-info_ND4-Det_rep1.csv, ../additional_file/20220104_DdCBE-detect-seq-info_ND4-Det_rep2.csv,../additional_file/20220104_DdCBE-detect-seq-info_ND5-1-Det_rep1.csv,../additional_file/20220104_DdCBE-detect-seq-info_ND5-1-Det_rep2.csv,../additional_file/20220104_DdCBE-detect-seq-info_ND6-Det_rep1.csv,../additional_file/20220104_DdCBE-detect-seq-info_ND6-Det_rep2.csv"
    pt_out = "../additional_file/20220104_DdCBE-Merged-info.csv"
    # deal with file list
    ls_target = [i.strip() for i in pt_target.split(",") if i.strip() != ""]
    ls_detect = [i.strip() for i in pt_detect.split(",") if i.strip() != ""]

    df_target = pd.DataFrame()

    for pt_df in ls_target:
        df = pd.read_csv(pt_df, sep=",", header=0, index_col=False)
        df_target = pd.concat([df_target, df], axis=0, ignore_index=True)

    df_detect = pd.DataFrame()

    for pt_df in ls_detect:
        df = pd.read_csv(pt_df, sep=",", header=0, index_col=False, dtype={'A_ratio': np.float64, 'T_ratio': np.float64, 'C_ratio': np.float64, 'G_ratio': np.float64})
        df_detect = pd.concat([df_detect, df], axis=0, ignore_index=True)

    # for test!!
    # "".join(
    #     df_target[
    #         (df_target.treatment == "N6-2000-12")
    #         & (df_target.rep == "rep1")
    #         & (df_target.cutoff == 3)
    #         & (df_target.region_id == "ND6-only-1")
    #         & (df_target.relative_pos >= 15)
    #         & (df_target.mut_base == "C")
    #     ].ref_base
    # )

    # "".join(
    #     df_detect[
    #         (df_detect.treatment == "ND6-Det")
    #         & (df_detect.rep == "rep1")
    #         & (df_detect.region_id == "ND6-only-1")
    #         & (df_detect.relative_pos >= 15)
    #     ].genome_base
    # )
    # fix target rep3 to detect rep1
    df_target["rep_back_target"] = pd.Series(dtype=object)
    df_target.rep_back_target = df_target.rep
    df_target.loc[df_target.rep == "rep3", "rep"] = "rep1"
    df_target[df_target.rep_back_target == "rep3"]
    # fix na
    df_detect.fillna(0)
    for col in df_detect[df_detect.A_ratio.map(math.isnan)][
        ["A_ratio", "G_ratio", "C_ratio", "T_ratio"]
    ]:
        df_detect.loc[:, col] = 0
    df_detect.isna().sum()
    # merge detect-seq and target-seq
    df = pd.merge(
        left=df_target,
        right=df_detect,
        how="inner",
        on=["rep", "region_id", "relative_pos"],
        suffixes=("_TargetSeq", "_DetectSeq"),
        indicator=True,
    )
    # check na
    ifna = df.isna().sum().sum()

    if ifna != 0:
        logging.warn("{} NA values were found in merged table!".format(ifna))

    logging.info("Merging done!")
    df
    df.to_csv(pt_out, sep=",", index=False, header=True)
