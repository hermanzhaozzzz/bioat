import argparse
import gzip
import logging
import os
import sys

import numpy as np
import pandas as pd
from Bio import SeqIO

logging.basicConfig(
    level=(4 - 30) * 10,
    format="%(levelname)-5s @ %(asctime)s: %(message)s ",
    datefmt="%Y-%m-%d %H:%M:%S",
    stream=sys.stderr,
    filemode="w",
)


class SiteIndex(object):
    """SiteIndex.
    SiteIndex object can parse and return base mutation info.

    Args
        site_index (str): Like chr1_10000_CT or chr1_10000.

    Methods
        no methods

    Attributes
        chr_name (str): chr_name.
        chr_index (int): chr_index.
        site_index (str): like chr1_10000.
        site_index_raw (str): like chr1_10000_CT.

    Raises
        no raises.
    """

    def __init__(self, site_index):
        site_index_split = site_index.split("_")

        self.chr_name = site_index_split[0]
        self.chr_index = int(site_index_split[1])
        self.site_index = site_index_split[0] + "_" + site_index_split[1]
        self.site_index_raw = site_index

        if len(site_index_split) > 2:
            self.mut_type = site_index_split[2]


class BmatLine(object):
    """BmatLine.
    SiteIndex object can parse and return base mutation info.

    Args
        bmat_line_list (list): Like this:
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

    Methods
        no methods

    Attributes
        chr_name (str): chr_name.
        chr_index (int): chr_index.
        ref_base (str): ref_base.
        count_del (int): count_del.
        count_insert (int): count_insert.
        count_ambis (int): count_ambis.
        deletion (str): deletion.
        ambiguous (str): ambiguous.
        insertion (str): insertion.
        mut_num (int): mut_num.
        count_dict (dict): like this:
            {
                'A': 0,
                'G': 0,
                'C': 10,
                'T': 25,
            }.

    Raises
        no raises.
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
    """CLI interface.

    Args
        no args.

    Returns
        argparse.Namespace
    Raises
        no raises.
    """
    args = argparse.ArgumentParser(
        description="A tools for formating MPMAT file from BED6/BED12 file with necessary fake information."
    )
    args.add_argument(
        "--bedHasHeader",
        type=bool,
        required=True,
        help="True if bed has header, " "else False",
    )
    args.add_argument(
        "--bmatHasHeader",
        type=bool,
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
        "-g", "--genome", required=True, help="Path of the genome fastx file"
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
    parser = args.parse_args()

    return parser


def cmp_site_index(site_index_a, site_index_b, ref_order_dict):
    """cmp_site_index.

    Args
        site_index_a (str): like chr1_629627_CT.
        site_index_b (str): like chr1_629627_CT.
        ref_order_dict (dict): like this:
            {
                'chr1': 0,
                'chr19': 1,
                'chr20': 2
            }

    Returns
        status (int): -1|0|1, details:
            -1, site_a at upstream of site_b
            0, site_index_a == site_index_b at position level, ignore mutation info
            1, site_a at downstream of site_b

    Raises
        TypeError: see codes for more.
    """
    # 这里可能有个问题就是有ref 没chr， 有ref 有chr
    # 比如 chr1_939943_CT 1_10006_C这两个去比较
    # 策略是加上chr，肯定不够通用，以后再说改的事儿
    if not site_index_a.startswith("chr"):
        site_index_a = f"chr{site_index_a}"

    if not site_index_b.startswith("chr"):
        site_index_b = f"chr{site_index_b}"

    site_A = SiteIndex(site_index_a)
    site_B = SiteIndex(site_index_b)

    if site_A.chr_name == site_B.chr_name:
        if site_A.chr_index == site_B.chr_index:
            return 0

        elif site_A.chr_index < site_B.chr_index:
            return -1

        elif site_A.chr_index > site_B.chr_index:
            return 1

    else:
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
    """load_reference_fasta_as_dict.

    Args
        ref_fasta_path (str): Reference fastx file path.
        ref_name_list (str): If set All, load all seq info in reference,
            else only try to load seq_id in the list.
        ref_order_dict (dict): like this:
            {
                'chr1': 0,
                'chr19': 1,
                'chr20': 2
            }

    Returns
        ref_seq_dict (dict): A dict, key is seq_id and value is sequence with .upper()
        or
        None (NoneType): If occur error, return None.

    Raises
        IOError: see codes for more.
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
        if ref_fasta_path.endswith(".gz"):
            genome_fa = SeqIO.parse(
                handle=gzip.open(ref_fasta_path, "rt"), format="fastx"
            )
        else:
            genome_fa = SeqIO.parse(handle=open(ref_fasta_path, "rt"), format="fastx")
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
        # for debug
        # break

    del genome_fa
    logging.info("Loading genome done!")

    return ref_seq_dict


def query_region_bmat_info(bmat_file, site_index_list, genome_order_dict):
    """query_region_bmat_info.

    Args
        bmat_file (BmatLine object): Reference fastx file path.
        site_index_list (list): like [chr1_20452_CT, chr1_20467_C., chr1_20474_CT].
        genome_order_dict (dict): or FUN <cmp_site_index>, like this:
            {
                'chr1': 0,
                'chr19': 1,
                'chr20': 2
            }

    Returns
        site_dict (dict): ey is site_index, value bmat line list.

    Raises
        IOError: see codes for more.
    """
    site_index = SiteIndex(site_index_list[0])
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
            query_res_dict[site_index.site_index] = BmatLine(bmat_line_list)

            if query_site_num >= query_total_num:
                break

            else:
                # read new
                site_index = SiteIndex(site_index_list[query_site_num])
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
                site_index = SiteIndex(site_index_list[query_site_num])
                bmat_line = bmat_file.readline()
                if type(bmat_line) == str:
                    pass
                else:
                    bmat_line = bmat_line.decode(encoding="utf8")

    return query_res_dict


if __name__ == "__main__":
    DEBUG = True
    if DEBUG:
        BED_PATH = "../Test_Resources/test.bed.gz"
        # BMAT_PATH = "../Test_Resources/test.bmat"
        BMAT_PATH = "../Test_Resources/test.bmat.gz"
        REF_PATH = "/Users/zhaohuanan/igv/genomes/seq/hg38.fa.gz"
        # mm10.fa.gz
        bmatHasHeader = True
        bedHasHeader = True
        OUT_MPMAT = "../Test_Outfiles/test.mpmat.gz"
        farTermExtend = 0
        pamTermExtend = 0
    else:
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

    REF_FA = load_reference_fasta_as_dict(ref_fasta_path=REF_PATH, log_verbose=30)
    bmat_file = (
        gzip.open(BMAT_PATH, "rt")
        if BMAT_PATH.endswith(".gz")
        else open(BMAT_PATH, "rt")
    )
    # create output mpmat
    out_mpmat = (
        gzip.open(OUT_MPMAT, "at")
        if OUT_MPMAT.endswith(".gz")
        else open(OUT_MPMAT, "at")
    )

    if bmatHasHeader:
        fake_var = next(bmat_file)

    # 判断bed是否是空的, 空的抛出错误
    assert os.path.getsize(BED_PATH)

    # 将bed读作Dataframe
    if bedHasHeader:
        df = pd.read_csv(BED_PATH, sep="\t")
    else:
        df = pd.read_csv(BED_PATH, sep="\t", header=None)
    # 只取bed的前6列信息：BED6
    df = df.iloc[:, 0:7]
    df.columns = ["chrom", "start", "end", "name", "score", "strand"]
    # 将每个 染色体 内的 bed信息 提取至词典中 供后续 按照 染色体 调用
    dt_chrom_bed_info = {}

    for chrom, df_ in df.groupby("chrom"):
        df_ = df_.sort_values(by=["start"])
        df_.index = range(df_.shape[0])
        dt_chrom_bed_info[chrom] = df_

    for chrom, df_ in df.groupby("chrom"):
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
                try:
                    seq = REF_FA[chrom][(start - 1) : end]
                except KeyError:
                    logging.warning(f"Skip Bed Coordinate: {chrom}\t{start}\t{end}")
                    continue

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
                    # logging.debug("=" * 20)
                    # seqstr = "this seq: ", site_index_list
                    # logging.debug(seqstr)
                    site_index_list_mut_count = []
                    site_index_list_coverage = []
                    for base in site_index_list:
                        # basestr = "this base:" + base
                        # logging.debug(basestr)
                        try:
                            dt_this_base = query_mut_info[base[:-3]].count_dict
                            # logging.debug(dt_this_base)
                            site_index_list_mut_count.append(
                                int(dt_this_base[base[-1]])
                            )
                            coverage = 0
                            for key in dt_this_base:
                                coverage += int(dt_this_base[key])
                            site_index_list_coverage.append(int(coverage))
                        except KeyError:
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
                    # logging.debug("DEBUG: " + mpmat_line)
                    out_mpmat.write(mpmat_line)
                else:
                    # 如果seq中没有要check的base，比如+没C或-没G就pass这个点
                    pass

logging.info("The program done!")

bmat_file.close()
out_mpmat.close()
