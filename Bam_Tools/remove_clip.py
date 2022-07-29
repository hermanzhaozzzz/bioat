# ——————————————————>>>>>>>>>>
# Project Name: remove_clip.py
# Author: Herman ZHAO
# E-mail: hermanzhaozzzz@gmail.com
# Update log:
#     2022-07-29: start project
#     YYYY-MM-DD: fix project #
# ——————————————————>>>>>>>>>>
import argparse
import os
import sys

import pysam


def compared_version(ver1, ver2):
    """Compared_version number.
    Compared_version number like: "9.1.1">"10.1.2", or "10.12.2.6.5">"10.12.2.6".

    Args:
        ver1: version 1
        ver2: version 2
    Returns:
        int, -1/0/1, ver1<ver2, return -1; ver1==ver2; return 0; ver1 > ver2, return 1.
    """
    list1 = str(ver1).split(".")
    list2 = str(ver2).split(".")
    # 循环次数为短的列表的len
    for i in range(len(list1)) if len(
            list1) < len(list2) else range(len(list2)):
        if int(list1[i]) == int(list2[i]):
            pass
        elif int(list1[i]) < int(list2[i]):
            return -1
        else:
            return 1
    # 循环结束，哪个列表长哪个版本号高
    if len(list1) == len(list2):
        return 0
    elif len(list1) < len(list2):
        return -1
    else:
        return 1


def check_version():
    if compared_version(pysam.__version__, "0.19.1") == -1:
        raise ValueError(
            f"pysam version is {pysam.__version__} and it should >=0.19.1")
    if compared_version(sys.version.split()[0], "3.9.1") == -1:
        raise ValueError(
            f"Python version is {sys.version.split()[0]} and it should >= 3.9.1:)")


def get_args():
    parser = argparse.ArgumentParser(
        description="remove softclip reads in BAM file")

    parser.add_argument(
        "-i",
        "--input",
        help="BAM file sorted by coordinate with soft/hard clip reads",
        type=str,
        nargs='?',
        default=sys.stdin)

    parser.add_argument(
        "-o",
        "--output",
        help="BAM file sorted by coordinate without soft/hard clip reads",
        type=str,
        nargs='?',
        default=sys.stdout)
    parser.add_argument(
        "-t",
        "--threads",
        help="threads used by pysam and samtools core, default=os.cpu_count() - 1",
        type=str,
        default=os.cpu_count() -
                1)
    parser.add_argument(
        "-b",
        "--output_format",
        help="BAM, SAM(default)",
        type=str,
        default="SAM"
    )
    parser.add_argument(
        "-p",
        "--remove_as_paired",
        help="True, False. remove single/paired-reads if any clip in read1 or read2(default: True)",
        type=bool,
        default=True)
    parser.add_argument(
        "-c",
        "--max_clip",
        help="the maximum clips allowed per read, default=0",
        type=int,
        default=0
    )
    return parser.parse_args()


if __name__ == "__main__":
    args = get_args()

    # [print(i) for i in dir(pysam)]

    # input_bam = "../Test_Resources/test_sorted.bam"
    # # input_bam = "../Test_Resources/test_sorted_n.bam"
    # output_bam = "../Test_Outfiles/test_sorted_rmclip.bam"
    check_version()

    if sys.stdin.isatty():
        try:
            bam_in = pysam.AlignmentFile(args.input, "r", threads=args.threads, check_sq=False)
        except ValueError:
            bam_in = pysam.AlignmentFile(args.input, "rb", threads=args.threads, check_sq=False)
    else:
        try:
            bam_in = pysam.AlignmentFile("-", "r", threads=args.threads, check_sq=False)
        except ValueError:
            bam_in = pysam.AlignmentFile("-", "rb", threads=args.threads, check_sq=False)

    try:
        so = bam_in.header['HD']['SO']
    except KeyError:
        so = None

    if so != "queryname" or so == "coordinate" or so is None:
        raise ValueError(f"the input BAM|SAM must be sorted by name and has header SO:queryname!\nyour header: {so}")

    write_mode = 'wb' if args.output_format == "BAM" else 'w'
    bam_out = pysam.AlignmentFile(args.output, "w", template=bam_in)

    # Iterate through reads.
    read1, read2 = None, None

    for read in bam_in:
        if not read.is_paired or read.is_unmapped or read.mate_is_unmapped:
            # or read.is_duplicate:
            # only use mapped paired reads
            continue
        # ------------------------------------------------------------------->>>>>>>>>>
        # 如果先看到read1，将read1赋给read， 判断是否read1、read2都存在
        #    不都存在：
        #        即read2不存在，进入下个循环
        #    都存在，判断query name是否一致
        # ------------------------------------------------------------------->>>>>>>>>>
        if read.is_read2:
            # if first is read2, then use it to find read2
            read2 = read
        else:
            # if first is read1, then use it to find read1, and next query
            read1 = read
            read2 = None
            continue

        if None not in [
            read1,
            read2] and read1.query_name == read2.query_name:
            # print("found a pair!")
            # print(f'Cigar: Read1-{read1.cigarstring}, Read2-{read1.cigarstring}')
            # print(read.get_cigar_stats())
            softclip1 = read1.get_cigar_stats()[0][4]
            softclip2 = read2.get_cigar_stats()[0][4]

            if not args.remove_as_paired:
                if softclip1 <= args.max_clip:
                    bam_out.write(read1)
                if softclip2 <= args.max_clip:
                    bam_out.write(read2)
            elif args.remove_as_paired:
                if softclip1 <= args.max_clip and softclip2 <= args.max_clip:
                    bam_out.write(read1)
                    bam_out.write(read2)
                else:
                    # print("Filtering pairs")
                    # print(f'Cigar: Read1-{read1.cigarstring}, Read2-{read2.cigarstring}')
                    # print(read1.get_cigar_stats()[0][4])
                    # print(read2.get_cigar_stats()[0][4])
                    pass
            else:
                pass
                # print(f'Cigar: Read1-{read1.cigarstring}, Read2-{read2.cigarstring}')
                # raise ValueError(
                #     f"remove_se_or_pe can not be {args.remove_se_or_pe}")

    bam_in.close()
    bam_out.close()
