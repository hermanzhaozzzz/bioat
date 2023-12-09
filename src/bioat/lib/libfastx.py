import sys
import re
import os
import subprocess
from Bio import SeqIO
from bioat.logger import get_logger

__module_name__ = "bioat.lib.libfastx"


def casfinder(input: str, output: str, lmin: int | None = None,
              lmax: int | None = None, log_level='DEBUG') -> None:
    # set logger
    logger = get_logger(level=log_level, module_name=__module_name__, func_name=sys._getframe().f_code.co_name)
    # filter length of assembly contigs
    assembly_input = SeqIO.parse(input, 'fasta')

    # define filter_func
    if lmin is None and lmax is not None:
        filter_func = lambda contig: len(contig) <= lmax
    elif lmin is not None and lmax is None:
        filter_func = lambda contig: len(contig) >= lmin
    else:
        filter_func = lambda contig: lmin <= len(contig) <= lmax

    contigs_input = (contig for contig in assembly_input)
    contigs_output = (contig for contig in contigs_input if filter_func(contig=contig))
    SeqIO.write(contigs_output, output, 'fasta')
    # percent = count_out / count_in * 100
    # logger.info(f"loaded {count_in}, filtered {count_in - count_out}, {percent:.1%} writen out")

    # for test
    prodigal = "/Users/zhaohuanan/micromamba/envs/snakepipes_Cas-mining/bin/prodigal"
    pilercr = "/Users/zhaohuanan/micromamba/envs/snakepipes_Cas-mining/bin/pilercr"
    work_path = '/Users/zhaohuanan/Downloads/test'

    # 定义一些常量和变量
    file = f"{work_path}/202155.assembled.fna.filtered.fasta"
    crisper_scaffold_file = f"{file}.crisper.scaffold.fa"
    pep_file = f"{file}.pep.fasta"
    pep_cas_file = f"{file}.pep.cas.fasta"
    crispr_loc_file = f"{file}.crispr.loc"
    crispr_spacer_file = f"{file}.crispr.spacer"
    temp_dir = f"{work_path}/"  # 用于存储临时文件的目录

    # filtered_fasta = SeqIO.parse(file, "fasta")

    subprocess.check_call(
        [prodigal, '-a', f'{work_path}/temp.pep', '-i', file, '-p', 'single', '-f', 'gff', '-o',
         f'{work_path}/temp.gff'],
    )
    subprocess.check_call(
        [pilercr, '-in', file, '-out', f'{work_path}/temp.spacer']
    )

    # 2.get cas locs
    with open(f'{work_path}/temp.spacer', 'rt') as spacer_file, open('temp.spacer.loc', 'wt') as loc_file:
        s = 0
        number = ""
        start = 0
        end = 0
        start1 = 0
        end1 = 0
        pre = ""
        for line in spacer_file:
            line = line.rstrip()
            if line.startswith("SUMMARY BY POSITION"):
                s = 1
            elif line.startswith("Help on reading this report"):
                s = 0
            if s > 0:
                bb = line.split()
                len_bb = len(bb)
                if line.startswith(">"):
                    ar = line.split()
                    scaffold = ar[0].replace(">", "")
                if len_bb > 7 and int(bb[1]) > 0:
                    if start1 < 0:
                        start1 = 1
                    if pre == "=":
                        start = int(bb[-6])
                        end = int(bb[-6]) + int(bb[-5]) + 1
                        start1 = start - 10000
                        end1 = end + 10000
                        if start1 < 1:
                            start1 = 1
                        print(
                            f"{number}_{scaffold}\t{start1}\t{end1}\t{start}\t{end}\t{bb[-1]}\t{bb[-4]}\t{bb[-3]}\t{bb[-2]}\t{int(bb[1])}\n",
                            file=loc_file)
                    else:
                        start = int(bb[-7])
                        end = int(bb[-7]) + int(bb[-6]) + 1
                        start1 = start - 10000
                        end1 = end + 10000
                        if start1 < 1:
                            start1 = 1
                        print(
                            f"{number}_{scaffold}\t{start1}\t{end1}\t{start}\t{end}\t{bb[-1]}\t{bb[-5]}\t{bb[-4]}\t{bb[-3]}\t{int(bb[1])}\n",
                            file=loc_file)
                pre_char = line[0]
                pre = pre_char
    # 3.get protein locs
    with open('temp.gff', 'r') as A, open('temp.gff.loc', 'w') as B:
        loc_dict = {}
        for line in A:
            line = line.rstrip()
            bb = line.split('\t')
            if bb[2] == "CDS":
                scaffold = f"{number}_{bb[0]}"
                ll = bb[8]
                ll = ll.replace(";", "\t")
                ll = ll.replace("_", "\t")
                b = ll.split("\t")
                gene = f"{scaffold}_{b[1]}"
                print(f"{scaffold}\t{bb[3]}\t{bb[4]}\t{bb[6]}\t{gene}", file=B)
                loc_dict[gene] = f"{scaffold}\t{bb[3]}\t{bb[4]}\t{bb[6]}"

    # 4.cas locs vs protein locs
    # 设置文件路径
    temp_gff_loc = "temp.gff.loc"
    temp_spacer_loc = "temp.spacer.loc"
    temp_bed = "temp.bed"

    # 使用bedtools进行交集操作
    os.system(f"bedtools intersect -wo -a {temp_gff_loc} -b {temp_spacer_loc} > {temp_bed}")

    # 打开并读取temp.bed文件
    with open(temp_bed, "r") as A:
        for line in A:
            line = line.rstrip()
            bb = line.split()
            line = line.replace("\_metacontig\_\_", "\t")
            ar = line.split()
            good[bb[4]] = 1
            cpscaffold[ar[1]] = 1
    # 5.get cas locs protein
    with open('temp.pep', 'r') as A, open('temp.pep.fasta', 'w') as B, open('temp.pep.filtered.fasta', 'w') as C:
        s = 0
        for line in A:
            line = line.rstrip()
            line = re.sub(r'^>', f'>{number}_', line)
            if line.startswith(">"):
                ll = line.replace(">", "")
                bb = ll.split()
                print(f">{bb[0]}\t{loc[bb[0]]}", file=B)
                if bb[0] in good:
                    s = 1
                    print(f">{bb[0]}\t{loc[bb[0]]}", file=C)
                else:
                    s = 0
            else:
                print(line, file=B)
                if s > 0:
                    print(line, file=C)
    # 6.get cas locas scaffold
    with open('temp', 'r') as A, open('temp-crisper-scaffold.fa', 'w') as B:
        for line in A:
            ll = line.rstrip()
            ll = ll.replace(">", "")
            bb = ll.split()
            try:
                line2 = next(A)
                if bb[0] in cpscaffold:
                    print(f"{line.rstrip()}{line2}", file=B)
            except StopIteration:
                break

    # 7.save result files
    os.system(f"cat temp-crisper-scaffold.fa >> {file}.crisper.scaffold.fa")
    os.system(f"cat temp.pep.fasta >> {file}.pep.fasta")
    os.system(f"cat temp.pep.filtered.fasta >> {file}.pep.cas.fasta")
    os.system(f"cat temp.spacer.loc >> {file}.crispr.loc")
    os.system(f"cat temp.spacer >> {file}.crispr.spacer")
    os.system("rm temp*")
    # os.system(f"rm {file}")

    counter = 0
    with open('temp', 'w') as AAA:
        pass


if __name__ == '__main__':
    casfinder(
        input="/Users/zhaohuanan/Downloads/test/202155.assembled.fna",
        output="/Users/zhaohuanan/Downloads/test/202155.assembled.fna.filtered.fasta",
        lmin=3001,  # 3001 in Nature Methods paper
        lmax=None
    )
