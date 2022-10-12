import argparse
import gzip
import os
import random
import string
import sys
import time
from multiprocessing import Process

import numpy as np


# ------------------------------------------------------------------->>>>>>>>>>
# define params
# ------------------------------------------------------------------->>>>>>>>>>
def get_parser():
    """get params from cmd python bmatFilter.py -<item> [parma]

    Returns:
        object: Object instance of argparse.ArgumentParser
    """
    parser = argparse.ArgumentParser(
        description="filter of bmat/bmat.gz, filter param全是与或非中的与关系，也就是条件同时成立！"
    )
    parser.add_argument(
        "-i", "--Input", required=True, type=str, help="输入bmat或者bmat.gz"
    )
    parser.add_argument(
        "-o", "--Output", required=True, type=str, help="输出bmat或者bmat.gz"
    )
    parser.add_argument(
        "-p",
        "--Threads",
        default=1,
        type=int,
        help="Multiple threads number, default=1",
    )
    parser.add_argument(
        "--TempDir",
        help="Where to keep temp files, default is the same dir with --Input",
        default=None,
    )
    parser.add_argument(
        "--MinDepth",
        default=0,
        type=int,
        help="能够接受的单碱基上覆盖的最小read count数，default=0不做过滤",
    )
    parser.add_argument(
        "--MinMutNum",
        default=-1,
        type=int,
        help=""""能够接受的单碱基最小突变数，
        小于MinMutNum的位点会被过滤；
        如果MinMutNum=0，则输出所有位点；
        如果MinMutNum=1，则不含突变的位点会被过滤；
        如果MinMutNum=2，则含有0个、1个突变的位点会被过滤。default=-1输出所有位点（和0作用一致）"
        """,
    )
    parser.add_argument(
        "--MaxMutNum",
        default=-1,
        type=int,
        help="""能够接受的单碱基最大突变数, 大于MaxMutNum的位点会被过滤；
        如果MaxMutNum=1，则只输出无突变的位点；
        如果MaxMutNum=2，则输出含有0个、1个突变的位点；
        如果MaxMutNum=-1，则不进行过滤。default=-1。
        default=-1，则不使用该过滤参数
        """,
    )
    parser.add_argument(
        "--MinMutRatio",
        default=-1,
        type=float,
        help="""能够接受的单碱基最小突变比率，default=-1，则不使用该过滤参数""",
    )
    parser.add_argument(
        "--MaxMutRatio",
        default=-1,
        type=float,
        help="""能够接受的单碱基最小突变比率，default=-1，则不使用该过滤参数""",
    )
    return parser


# ------------------------------------------------------------------->>>>>>>>>>
# split file for multiprocess
# ------------------------------------------------------------------->>>>>>>>>>
def split_file_and_make_temp(input_filename, threads_num=1, temp_dir=None):
    """split input file for multi threads calculating

    Args:
        input_filename (str): path of input file
        threads_num (int, optional): number of cores to use. Defaults to 1.
        temp_dir (str, optional): path of splited temp files. Defaults to None.

    Returns:
        list: path list of split files
    """
    # --------------------------------------------------->>>>>>
    # set temp dir
    # --------------------------------------------------->>>>>>
    if temp_dir is None:
        temp_dir = os.path.dirname(input_filename)

    # --------------------------------------------------->>>>>>
    # get input basename
    # --------------------------------------------------->>>>>>
    input_file_basename = os.path.basename(input_filename)

    # --------------------------------------------------->>>>>>
    # creat temp file name and open file
    # --------------------------------------------------->>>>>>
    # >>>>>> LOG
    sys.stderr.write(
        "Making temp files... \t %s \n"
        % (time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()))
    )

    temp_filename_list = []
    temp_file_list = []

    for index in range(threads_num):
        # make filename
        temp_file_basename = (
            "temp_"
            + input_file_basename
            + "."
            + str(index)
            + "."
            + "".join(random.sample(string.ascii_letters + string.digits, 16))
        )
        temp_file_name = os.path.join(temp_dir, temp_file_basename)
        temp_filename_list.append(temp_file_name)

        # open file
        if IfGzipInput:
            temp_file = gzip.open(temp_file_name, "w")
        elif not IfGzipInput:
            temp_file = open(temp_file_name, "w")
        temp_file_list.append(temp_file)

    # --------------------------------------------------->>>>>>
    # counting input file line number
    # --------------------------------------------------->>>>>>
    # >>>>>> LOG
    sys.stderr.write(
        "Counting input file... \t %s \n"
        % (time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()))
    )

    if IfGzipInput:
        input_file = gzip.open(input_filename, "r")
    elif not IfGzipInput:
        input_file = open(input_filename, "r")

    total_input_line_num = 0
    for line in input_file:
        total_input_line_num += 1

    if total_input_line_num % threads_num == 0:
        each_file_line_num = total_input_line_num // threads_num
    else:
        each_file_line_num = (total_input_line_num // threads_num) + 1

    input_file.close()

    # >>>>>> LOG
    sys.stderr.write(
        "Done! \t %s \n" % (time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()))
    )

    # --------------------------------------------------->>>>>>
    # split temp files
    # --------------------------------------------------->>>>>>
    # write into output filename
    if IfGzipInput:
        input_file = gzip.open(input_filename, "r")
    else:
        input_file = open(input_filename, "r")

    for index, line in enumerate(input_file):
        file_index = index // each_file_line_num
        temp_file_list[file_index].write(line)

    # close output filename
    input_file.close()
    [temp_file.close() for temp_file in temp_file_list]

    # >>>>>> LOG
    sys.stderr.write(
        "Make temp files done! \t %s \n"
        % (time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()))
    )

    return temp_filename_list


# ------------------------------------------------------------------->>>>>>>>>>
# merge temp to output
# ------------------------------------------------------------------->>>>>>>>>>
def merge_out_bmat(temp_filename_list, output_filename, remove_temp=True):
    """merge temp files to output file

    Args:
        temp_filename_list (list): the returned list of func: split_file_and_make_temp
        output_filename (str): path of output file
        remove_temp (bool, optional): Whether to delete the temp files after the program runs. Defaults to True.
    """

    # write header
    header = [
        "chr_name",
        "chr_index",
        "ref_base",
        "A",
        "G",
        "C",
        "T",
        "del_count",
        "insert_count",
        "ambiguous_count",
        "deletion",
        "insertion",
        "ambiguous",
        "mut_num",
    ]

    # open output file
    if output_filename == "Stdout":
        output_file = sys.stdout
    elif IfGzipOutput:
        output_file = gzip.open(output_filename, "w")
    elif not IfGzipOutput:
        output_file = open(output_filename, "w")

    # write header
    if IfGzipOutput:
        output_file.write(("\t".join(header) + "\n").encode("utf-8"))
    elif not IfGzipOutput:
        output_file.write("\t".join(header) + "\n")

    # open temp bmat file and write to merge file
    for temp_filename in temp_filename_list:
        # write file
        if IfGzipInput:
            with gzip.open(temp_filename + ".bmat", "r") as temp_file:
                for line in temp_file:
                    output_file.write(line)
        elif not IfGzipInput:
            with open(temp_filename + ".bmat", "r") as temp_file:
                for line in temp_file:
                    output_file.write(line)
    # remove temp files
    if remove_temp:
        for temp_filename in temp_filename_list:
            os.remove(temp_filename)
            os.remove(temp_filename + ".bmat")


# ------------------------------------------------------------------->>>>>>>>>>
# filter core function
# ------------------------------------------------------------------->>>>>>>>>>
def core_filter(
    Input,
    Output="Stdout",
    IfGzipInput=False,
    IfGzipOutput=False,
    MinDepth=0,
    MinMutNum=-1,
    MaxMutNum=-1,
    MinMutRatio=-1,
    MaxMutRatio=-1,
):
    """the filter function to screen bmat lines

    Args:
        Input (str): one splited path name from the returned list of function: split_file_and_make_temp
        Output (str, optional): path of output file. Defaults to "Stdout".
        MinDepth (int, optional): min read count of one base. Defaults to 0.
        MinMutNum (int, optional): . Defaults to 0.
    """
    # open input file
    if IfGzipInput:
        InputFile = gzip.open(Input, "r")
    elif not IfGzipInput:
        InputFile = open(Input, "r")
    # open output file
    if Output == "Stdout":
        OutputFile = sys.stdout
    elif IfGzipOutput:
        OutputFile = gzip.open(Output, "a")
    elif not IfGzipOutput:
        OutputFile = open(Output, "a")

    # iterator
    for line in InputFile:
        if IfGzipInput:
            line = line.decode("utf-8")
        elif not IfGzipInput:
            pass

        ls_line = line.strip().split("\t")

        # 不分析header行
        if ls_line[0] == "chr_name":
            continue

        array_line = np.array([int(i) for i in ls_line[3:7]])
        LineDepth = np.sum(array_line)

        # 不满足最小深度的，过滤掉
        if LineDepth < MinDepth:
            continue

        MutCount = int(ls_line[-1])
        # 小于可接受最小突变数的，过滤掉
        if (MinMutNum != -1) and (MutCount < MinMutNum):
            continue
        # 大于可接受最大突变数的，过滤掉
        if (MaxMutNum != -1) and (MutCount > MaxMutNum):
            continue

        # ratio和Num一起作用
        ThisBaseMutRatio = MutCount / LineDepth if LineDepth > 0.0 else 0.0
        # 小于可接受最小突变率的，过滤掉
        if (MinMutRatio != -1.0) and (ThisBaseMutRatio < MinMutRatio):
            continue
        # 大于可接受最大突变率的，过滤掉
        if (MaxMutRatio != -1.0) and (ThisBaseMutRatio > MaxMutRatio):
            continue

        # ------------------------------------------------------------------->>>>>>>>>>
        # 临时条件
        # if($total <= 100 && $total >= 50){
        #     next if $mutations < ($total-10);#mismatch < depth - 10 # == if True: continue
        #     print OUT "$_\n";
        # }
        # elsif($total > 100){
        #     next if $mutations < ($total-20);#mismatch < depth - 20 # == if True: continue
        #     print OUT "$_\n";
        # ------------------------------------------------------------------->>>>>>>>>>
        # if LineDepth < 50:
        #     continue
        # if (LineDepth <= 100) and (LineDepth >= 50):
        #     if MutCount < LineDepth - 10:
        #         continue
        # elif LineDepth > 100:
        #     if MutCount < LineDepth - 20:
        #         continue
        # ------------------------------------------------------------------->>>>>>>>>>

        # 通过 filter 的行被写入输出文件
        str_out = "\t".join(ls_line)
        str_out = str_out + "\n"
        if IfGzipOutput:
            OutputFile.write(str_out.encode("utf-8"))
        elif not IfGzipOutput:
            OutputFile.write(str_out)
    InputFile.close()
    OutputFile.close()


if __name__ == "__main__":
    # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # load the paramters
    # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    parser = get_parser()
    args = parser.parse_args()
    InputFilePath = args.Input
    OutputFilePath = args.Output
    Threads = args.Threads
    TempDir = args.TempDir
    MinDepth = args.MinDepth
    MinMutNum = args.MinMutNum
    MaxMutNum = args.MaxMutNum
    MinMutRatio = args.MinMutRatio
    MaxMutRatio = args.MaxMutRatio
    print("*" * 30)
    print("MinDepth\t", MinDepth)
    print("MinMutNum\t", MinMutNum)
    print("MaxMutNum\t", MaxMutNum)
    print("MinMutRatio\t", MinMutRatio)
    print("MaxMutRatio\t", MaxMutRatio)
    print("*" * 30)
    IfGzipInput = True if InputFilePath[-3:] == ".gz" else False
    IfGzipOutput = True if OutputFilePath[-3:] == ".gz" else False
    # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # remove output file
    # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    if os.path.exists(OutputFilePath):
        print("file %s exists & remove it" % OutputFilePath)
        os.remove(OutputFilePath)
    # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # make temp files
    # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    temp_filename_list = split_file_and_make_temp(
        input_filename=InputFilePath, threads_num=Threads, temp_dir=TempDir
    )
    # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # multiple processing part
    # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    sys.stderr.write(
        "Parsing files... \t %s \n"
        % (time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()))
    )

    # multiprocess
    procs_list = []

    for index, temp_filename in enumerate(temp_filename_list):
        out_filename = temp_filename + ".bmat"

        # add multiprocess
        # params

        sub_proc = Process(
            target=core_filter,
            args=(
                temp_filename,  # Input,
                out_filename,  # Output="Stdout",
                IfGzipInput,  # IfGzipInput=False,
                IfGzipOutput,  # IfGzipOutput=False,
                MinDepth,  # MinDepth=0,
                MinMutNum,  # MinMutNum=-1,
                MaxMutNum,  # MaxMutNum=-1,
                MinMutRatio,  # MinMutRatio=-1,
                MaxMutRatio,  # MaxMutRatio=-1,
            ),
        )
        procs_list.append(sub_proc)
        sub_proc.start()

    for sub_proc in procs_list:
        sub_proc.join()

    sys.stderr.write(
        "Done! \t %s \n" % (time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()))
    )
    # >>>>>> LOG
    sys.stderr.write(
        "Merging files... \t %s \n"
        % (time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()))
    )

    merge_out_bmat(
        temp_filename_list=temp_filename_list,
        output_filename=OutputFilePath,
        remove_temp=True,
    )

    # >>>>>> LOG
    sys.stderr.write(
        "Done!... \t %s \n" % (time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()))
    )
