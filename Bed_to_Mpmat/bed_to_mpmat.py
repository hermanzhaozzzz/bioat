import argparse
import gzip
import logging
import os
import sys

import numpy as np
import pandas as pd
import pysam
from Bio import SeqIO

logging.basicConfig(
    level=(4 - 30) * 10,
    format="%(levelname)-5s @ %(asctime)s: %(message)s ",
    datefmt="%Y-%m-%d %H:%M:%S",
    stream=sys.stderr,
    filemode="w",
)


class siteIndex(object):
    """
    INPUT
        <site_index>
            str, like chr1_10000_CT

    RETURN
        siteIndex obj
    """

    def __init__(self, _site_index):
        site_index_split = _site_index.split("_")

        self.chr_name = site_index_split[0]
        self.site_pos = int(site_index_split[1])
        self.site_index = site_index_split[0] + "_" + site_index_split[1]
        self.site_index_raw = _site_index

        if len(site_index_split) > 2:
            self.mut_type = site_index_split[2]


class bmatLine(object):
    """
    INPUT:
        <bmat_line_list>
            list, info like
            [
               'chr1',
               '10013',
               'T',
               '0',
               '0',
               '0',
               '6',
               '0',
               '0',
               '0',
               '.',
               '.',
               '.',
               '0'
            ]

    RETURN:
        A bmatLine obj
    """

    def __init__(self, bmat_line_list):
        self.chr_name = bmat_line_list[0]
        self.chr_index = int(bmat_line_list[1])
        self.ref_base = bmat_line_list[2]

        self.count_del = int(bmat_line_list[7])
        self.count_insert = int(bmat_line_list[8])
        self.count_ambis = int(bmat_line_list[9])

        self.deletion = bmat_line_list[10]
        self.ambiguous = bmat_line_list[11]
        self.insertion = bmat_line_list[12]
        self.mut_num = int(bmat_line_list[13])

        # make dict
        self.count_dict = {
            "A": int(bmat_line_list[3]),
            "G": int(bmat_line_list[4]),
            "C": int(bmat_line_list[5]),
            "T": int(bmat_line_list[6]),
        }


def parser():
    args = argparse.ArgumentParser(
        description="A tools for formating MPMAT file from BED6/BED12 file with necessary fake information."
    )
    args.add_argument(
        "--bedHasHeader",
        type=str,
        required=True,
        help="True if bed has header, else False",
    )
    args.add_argument(
        "--bmatHasHeader",
        type=str,
        required=True,
        help="True if bmat has header, else False",
    )
    args.add_argument("-b", "--bed", required=True, help="Path of the bed file")
    args.add_argument(
        "-m",
        "--bmat",
        required=True,
        help="Path of the bmat file, gzip file is accepted",
    )
    args.add_argument(
        "-g", "--genome", required=True, help="Path of the genome fasta file"
    )
    args.add_argument(
        "-o",
        "--output_mpmat",
        default="./output.mpmat",
        help="Path of the fake mpmat file, default: ./output.mpmat",
    )
    args.add_argument(
        "--extend_length_PAM_term",
        type=int,
        default=0,
        help="Pextend length at PAM term, default: int 0",
    )
    args.add_argument(
        "--extend_length_far_end_term_from_PAM",
        type=int,
        default=0,
        help="Pextend length at far end term from PAM, default: int 0",
    )
    return args.parse_args()


def cmp_site_index(site_index_a, site_index_b, ref_order_dict):
    """
    INPUT:
        <site_index_a>, <site_index_b>
            str, like chr1_629627_CT

        <ref_order>
            dict, format like:

            ref_order_dict = {
                'chr1': 0,
                'chr19': 1,
                'chr20': 2
            }

    RETURN:
        0, site_index_a == site_index_b at position level, ignore mutation info;
        -1, site_a at upstream of site_b;
        1, site_a at downstream of site_b
    """
    # print(site_index_a)
    # print(site_index_b)
    site_A = siteIndex(site_index_a)
    site_B = siteIndex(site_index_b)

    if site_A.chr_name == site_B.chr_name:
        if site_A.site_pos == site_B.site_pos:
            return 0

        elif site_A.site_pos < site_B.site_pos:
            return -1

        elif site_A.site_pos > site_B.site_pos:
            return 1

    else:
        # print(site_A)
        site_A_chr_order = ref_order_dict.get(site_A.chr_name)
        site_B_chr_order = ref_order_dict.get(site_B.chr_name)

        if (site_A_chr_order is not None) and (site_B_chr_order is not None):
            if site_A_chr_order < site_B_chr_order:
                return -1

            elif site_A_chr_order > site_B_chr_order:
                return 1

        else:
            raise TypeError("Site index not in your reference!")


def load_reference_fasta_as_dict(ref_fasta_path, ref_name_list="All", log_verbose=30):
    """
    INPUT:
        <ref_fasta_path>
            Reference fasta file path

        <ref_name_list>
            If set All, load all seq info in reference, else only try to load seq_id in the list

    RETURN
        <ref_seq_dict>
            A dict, key is seq_id and value is sequence with  .upper()

        None
            If occur error, return None.
    """

    # ---------------------------------------------------------->>>>>>>>>>>>>>>>>>>>
    # log setting
    # ---------------------------------------------------------->>>>>>>>>>>>>>>>>>>>
    logging.basicConfig(
        level=(4 - log_verbose) * 10,
        format="%(levelname)-5s @ %(asctime)s: %(message)s ",
        datefmt="%Y-%m-%d %H:%M:%S",
        stream=sys.stderr,
        filemode="w",
    )

    # ---------------------------------------------------------->>>>>>>>>>>>>>>>>>>>
    # load genome as dict
    # ---------------------------------------------------------->>>>>>>>>>>>>>>>>>>>
    try:
        genome_fa = SeqIO.parse(handle=ref_fasta_path, format="fasta")
    except:
        raise IOError("Load file error! %s" % ref_fasta_path)

    # init var
    ref_seq_dict = {}
    ref_name_set = set(ref_name_list)

    logging.info("Starting to load the reference genome...")

    for ref in genome_fa:
        if ref_name_list == "All":
            ref_seq_dict[ref.id] = ref.seq.upper()
            logging.debug("Loading genome...\t" + ref.id)

        elif ref.id in ref_name_list:
            ref_seq_dict[ref.id] = ref.seq.upper()
            logging.debug("Loading genome...\t" + ref.id)

            # remove already loaded seq
            ref_name_set.remove(ref.id)

            # load all info
            if len(ref_name_set) == 0:
                break

    logging.info("Loading genome done!")

    return ref_seq_dict


def query_region_bmat_info(bmat_file, site_index_list, genome_order_dict):
    """
    INPUT:
        <bmat_file>
            file.obj, bmat file handle

        <site_index_list>
            list, like [chr1_20452_CT, chr1_20467_C., chr1_20474_CT]

        <genome_order_dict>
            dict, for FUN <cmp_site_index>

    RETURN
        <site_dict>
            dict, key is site_index, value bmat line list

    """

    # init
    site_index = siteIndex(site_index_list[0])
    bmat_line = bmat_file.readline()
    if type(bmat_line) == str:
        pass
    else:
        bmat_line = bmat_line.decode(encoding="utf8")

    # define dict
    query_total_num = len(site_index_list)
    query_site_num = 0
    query_res_dict = {"site_index_list": []}

    while bmat_line != "":

        bmat_line_list = bmat_line.strip().split("\t")
        bmat_site_index = "_".join(bmat_line_list[0:3])

        cmp_res = cmp_site_index(
            site_index.site_index_raw, bmat_site_index, genome_order_dict
        )

        if cmp_res == 0:
            query_site_num += 1

            # add info into dict
            query_res_dict["site_index_list"].append(site_index.site_index)
            query_res_dict[site_index.site_index] = bmatLine(bmat_line_list)

            if query_site_num >= query_total_num:
                break

            else:
                # read new
                site_index = siteIndex(site_index_list[query_site_num])
                bmat_line = bmat_file.readline()
                if type(bmat_line) == str:
                    pass
                else:
                    bmat_line = bmat_line.decode(encoding="utf8")

        elif cmp_res == 1:
            bmat_line = bmat_file.readline()
            if type(bmat_line) == str:
                pass
            else:
                bmat_line = bmat_line.decode(encoding="utf8")

        elif cmp_res == -1:
            query_site_num += 1

            # add info into dict
            query_res_dict["site_index_list"].append(site_index.site_index)
            query_res_dict[site_index.site_index] = None

            if query_site_num >= query_total_num:
                break

            else:
                # read new
                site_index = siteIndex(site_index_list[query_site_num])
                bmat_line = bmat_file.readline()
                if type(bmat_line) == str:
                    pass
                else:
                    bmat_line = bmat_line.decode(encoding="utf8")

    return query_res_dict


if __name__ == "__main__":
    # BED_PATH = "Bed2MpmatFaker/test/20200611-293T-EMX1-Detect-seq_pRBS.bed"
    # # BMAT_PATH = "/Users/mac/mac_data2/FilesForTest/bmat/293T-bat_VEGFA-All-PD_rep1_hg38.MAPQ20.bmat"
    # BMAT_PATH = "/Users/mac/mac_data2/FilesForTest/bmat/293T-bat_VEGFA-All-PD_rep1_hg38.select.C.bmat.gz"
    # REF_PATH = "/Users/mac/Nutstore/Coding/github/snakepipes_bioinformatics_hermanzhaozzzz/genome_fa/genome_ucsc_hg38.fa"
    # bmatHasHeader = False
    # bedHasHeader = True
    # OUT_MPMAT = "/Users/mac/Downloads/test.txt"
    # farTermExtend = 0
    # pamTermExtend = 0

    # parse params
    parser = parser()
    BED_PATH = parser.bed
    BMAT_PATH = parser.bmat
    REF_PATH = parser.genome
    bmatHasHeader = parser.bmatHasHeader
    bedHasHeader = parser.bedHasHeader
    OUT_MPMAT = parser.output_mpmat
    farTermExtend = parser.extend_length_far_end_term_from_PAM
    pamTermExtend = parser.extend_length_PAM_term

    # load ref_fasta as dict
    REF_FA = load_reference_fasta_as_dict(ref_fasta_path=REF_PATH, log_verbose=30)
    # load bmat file
    if BMAT_PATH[-3:] == ".gz" or BMAT_PATH[-5:] == ".gzip":
        bmat_file = gzip.open(BMAT_PATH, "r")
    else:
        bmat_file = open(BMAT_PATH, "r")
    # create output mpmat
    out_mpmat = open(OUT_MPMAT, "a")

    # str to bool
    if bedHasHeader == "False":
        bedHasHeader = False
    elif bedHasHeader == "True":
        bedHasHeader = True
    else:
        raise ValueError("bedHasHeader must be True or False")
    if bmatHasHeader == "False":
        bmatHasHeader = False
    elif bmatHasHeader == "True":
        bmatHasHeader = True
    else:
        raise ValueError("bmatHasHeader must be True or False")
    if bmatHasHeader:
        fake_var = next(bmat_file)

    # 判断bed是否是空的,如果是空的，直接返回空mpmat，如果不是，继续
    with open(BED_PATH, "r") as f:
        if f.readlines() == []:
            bl_bed_empty = True
        else:
            bl_bed_empty = False

    # parse bed file
    if bl_bed_empty == True:
        out_mpmat.write("the bed file is empty!")
    else:
        # 将bed读作Dataframe
        if bedHasHeader:
            df = pd.read_csv(BED_PATH, sep="\t")
        else:
            df = pd.read_csv(BED_PATH, sep="\t", header=None)
        # 只取bed的前6列信息：BED6
        df = df.iloc[:, 0:7]
        df.columns = ["chrom", "start", "end", "name", "score", "strand"]
        ls_chrom = ["chr%s" % i for i in list(range(1, 23)) + ["X", "Y", "M"]]

        # 将每个 染色体 内的 bed信息 提取至词典中 供后续 按照 染色体 调用
        dt_chrom_bed_info = {}

        for chrom in ls_chrom:
            df_chrom = df[df["chrom"] == chrom]
            df_chrom = df_chrom.sort_values(by=["start"])
            df_chrom.index = range(df_chrom.shape[0])
            dt_chrom_bed_info[chrom] = df_chrom

        for chrom in ls_chrom:
            if dt_chrom_bed_info[chrom].empty:
                # 如果 词典中，对应染色体的那部分bed是空的，就跳过，继续分析下一个染色体中的bed
                continue
            else:
                # 对于某特定 染色体 中的bed
                for row in range(dt_chrom_bed_info[chrom].shape[0]):
                    # 遍历该Dataframe，每次取一行bed信息
                    # 并指定 start、end、name、score、strand
                    # 已知这一个bed区块的chrom信息
                    # 则指定了chrom、start、end、name、score、strand
                    start = dt_chrom_bed_info[chrom].loc[row, "start"]
                    end = dt_chrom_bed_info[chrom].loc[row, "end"]
                    name = dt_chrom_bed_info[chrom].loc[row, "name"]
                    score = dt_chrom_bed_info[chrom].loc[row, "score"]
                    strand = dt_chrom_bed_info[chrom].loc[row, "strand"]
                    # seq and abs_index of bases
                    # 从REF_FA的dict中将对应这一行bed信息的region的seq切片出来
                    seq = REF_FA[chrom][start - 1 : end]
                    seq = seq.upper()
                    # 生成seq对应的染色体位置的绝对index
                    seq_idx = np.arange(start, end + 1)
                    # form <site_index_list>
                    # It's a list
                    # like [chr1_20452_CT, chr1_20467_C., chr1_20474_CT]
                    # C don't mut to T --> C.
                    # C mut to T --> CT

                    site_index_list = []
                    for idx, base in enumerate(seq):
                        # 调整start和end的值
                        # region为正向，就看正链REF中C的突变情况
                        # region为反向，就看正链REF中G的突变情况，对应到反向链，其实就还是看C
                        if strand == "+" and base == "C":
                            base2base = "CT"
                            start -= farTermExtend
                            end += pamTermExtend
                        elif strand == "-" and base == "G":
                            base2base = "GA"
                            start -= pamTermExtend
                            end += farTermExtend
                        else:
                            continue
                        site_index_list.append(
                            "{chrom}_{abs_idx}_{base2base}".format(
                                chrom=chrom,
                                abs_idx=seq_idx[idx],
                                base2base=base2base,
                            )
                        )
                    if len(site_index_list) != 0 and site_index_list != [""]:
                        query_mut_info = query_region_bmat_info(
                            bmat_file=bmat_file,
                            site_index_list=site_index_list,
                            genome_order_dict=REF_FA,
                        )
                        logging.debug("=" * 20)
                        seqstr = "this seq: ", site_index_list
                        logging.debug(seqstr)
                        site_index_list_mut_count = []
                        site_index_list_coverage = []
                        for base in site_index_list:
                            basestr = "this base:" + base
                            logging.debug(basestr)
                            try:
                                dt_this_base = query_mut_info[base[:-3]].count_dict
                                logging.debug(dt_this_base)
                                site_index_list_mut_count.append(
                                    int(dt_this_base[base[-1]])
                                )
                                coverage = 0
                                for key in dt_this_base:
                                    coverage += int(dt_this_base[key])
                                site_index_list_coverage.append(int(coverage))
                            except:
                                site_index_list_mut_count.append(0)
                                site_index_list_coverage.append(0)

                        # logging.debug("=" * 20)
                        ls_ratio = []
                        for idx in range(len(site_index_list)):
                            if site_index_list_coverage[idx] == 0:
                                ls_ratio.append(0.0)
                            else:
                                ls_ratio.append(
                                    site_index_list_mut_count[idx]
                                    / (site_index_list_coverage[idx])
                                )
                        # logging.debug(ls_ratio)

                        count_mut_site_in_tandom = len(site_index_list)
                        for i in range(len(site_index_list_mut_count)):
                            if site_index_list_mut_count[i] == 0:
                                site_index_list[i] = site_index_list[i][:-1] + "."
                                count_mut_site_in_tandom -= 1

                        mpmat_line = (
                            "{chrom}\t{start_tandom}\t{end_tandom}\t{count_tandom_site}\t".format(
                                chrom=chrom,
                                start_tandom=site_index_list[0].split("_")[1],
                                end_tandom=site_index_list[-1].split("_")[1],
                                count_tandom_site=len(site_index_list),
                            )
                            + "{count_mut_site_in_tandom}\t{count_SNP_site_in_tandom}\t".format(
                                count_mut_site_in_tandom=count_mut_site_in_tandom,
                                count_SNP_site_in_tandom=0,
                            )
                            + "{mut_site_index}\t{mut_count_this_site}\t{mut_coverage_this_site}\t".format(
                                mut_site_index=",".join(site_index_list),
                                mut_count_this_site=",".join(
                                    [str(i) for i in site_index_list_mut_count]
                                ),
                                mut_coverage_this_site=",".join(
                                    [str(i) for i in site_index_list_coverage]
                                ),
                            )
                            + "{mut_ratio_this_site}\t".format(
                                mut_ratio_this_site=",".join([str(i) for i in ls_ratio])
                            )
                            + "{isSNP}\t{zeros}\t{passTest}\n".format(
                                isSNP=",".join(["False"] * len(site_index_list)),
                                zeros=",".join(["0"] * len(site_index_list)),
                                passTest=",".join(["Pass"] * len(site_index_list)),
                            )
                        )
                        logging.debug("DEBUG: " + mpmat_line)
                        out_mpmat.write(mpmat_line)
                    else:
                        # 如果seq中没有要check的base，比如+没C或-没G就pass这个点
                        pass

logging.debug("The program done!")

bmat_file.close()
out_mpmat.close()
