import sys
import re
import os
import subprocess
from Bio import SeqIO
from bioat.logger import get_logger

__module_name__ = "bioat.lib.libfastx"


def casfinder(input_fa: str, output_faa: str | None = None, lmin: int | None = None,
              lmax: int | None = None, extend: int = 10_000, log_level='DEBUG') -> None:
    # set logger
    logger = get_logger(level=log_level, module_name=__module_name__, func_name=sys._getframe().f_code.co_name)
    # -----------------------------
    # FOR TEST
    # -----------------------------
    prodigal = "/Users/zhaohuanan/micromamba/envs/snakepipes_Cas-mining/bin/prodigal"
    pilercr = "/Users/zhaohuanan/micromamba/envs/snakepipes_Cas-mining/bin/pilercr"
    work_path = '/Users/zhaohuanan/Downloads/test'

    # 定义一些常量和变量
    # input
    f_fa_input = input_fa  # input_fa = '202155.assembled.fna'
    # output
    f_faa_pep_cas = f"{input_fa}.pep.cas.faa" if output_faa is None else output_faa
    # others
    f_fa_filtered = f"{input_fa}.filtered.fa"
    f_pilercr_crispr_spacer = f"{input_fa}.crispr.spacer.pilercr"
    f_faa_pep = f"{input_fa}.pep.faa"
    f_gff = f"{input_fa}.gff"
    f_tsv_pep_loc = f"{input_fa}.pep.loc.tsv"
    f_tsv_crispr_loc = f"{input_fa}.crispr.loc.tsv"
    f_fa_crisper_scaffold = f"{input_fa}.crisper.scaffold.fa"

    tests = {
        0: False,  # 0. filter contigs
        1: False,  # 1. cas & protein annotation
        2: False,  # 2. get cas locs  # TODO 这里和原始代码的计算得出的坐标位置不同，仔细check每行计算哪里出了问题
        3: True,  # 3. get protein locs
        4: False,  # 4. cas locs vs protein locs
        5: False,  # 5. get cas locs protein
        6: False,  # 6. get cas locas scaffold
        7: False,  # 7. save result files
    }
    # -----------------------------
    # 0. filter contigs
    # -----------------------------
    logger.info('0. filter contigs')
    if tests[0]:
        # filter length of assembly contigs
        logger.debug(f"filter length of contigs for {f_fa_input}")
        assembly_input = SeqIO.parse(f_fa_input, 'fasta')
        # define filter_func
        if lmin is None and lmax is not None:
            filter_func = lambda contig: len(contig) <= lmax
        elif lmin is not None and lmax is None:
            filter_func = lambda contig: len(contig) >= lmin
        else:
            filter_func = lambda contig: lmin <= len(contig) <= lmax

        contigs_input = (contig for contig in assembly_input)
        contigs_output = (contig for contig in contigs_input if filter_func(contig=contig))
        logger.debug(f'writing filtered contigs to {f_fa_filtered}')
        SeqIO.write(contigs_output, f_fa_filtered, 'fasta')
    # -----------------------------
    # 1. cas & protein annotation
    # -----------------------------
    logger.info('1. cas & protein annotation')
    if tests[1]:
        logger.debug("subprocess call for prodigal")
        subprocess.check_call(
            [prodigal, '-a', f_faa_pep, '-i', f_fa_filtered, '-p', 'single', '-f', 'gff', '-o', f_gff],
            stderr=open('/dev/null', 'wt'),
        )
        logger.debug("call has returned")
        logger.debug("subprocess call for pilercr")
        subprocess.check_call(
            [pilercr, '-in', f_fa_filtered, '-out', f_pilercr_crispr_spacer],
            stderr=open('/dev/null', 'wt')
        )
        logger.debug("call has returned")
    # -----------------------------
    # 2.get cas locs
    # -----------------------------
    logger.info('2.get cas locs')
    if tests[2]:
        with open(f_pilercr_crispr_spacer, 'rt') as f_spacer, open(f_tsv_crispr_loc, 'wt') as f_loc:
            lines = f_spacer.readlines()
            lines = [line.rstrip() for line in lines if len(line.rstrip()) > 0]  # skip empty line, drop "\n"

            parse_status = False

            for line in lines:
                if line.startswith("SUMMARY BY POSITION"):  # the first line to parse
                    logger.debug('find `SUMMARY BY POSITION part`')
                    parse_status = True
                    continue
                if not parse_status:
                    continue

                # when parse_status = True
                # skip empty line when loading file (see beginning of the `loading`)
                # then, the format to parse, eg.
                """
                # in
                >Ga0307431_1000033
                Array          Sequence    Position      Length  # Copies  Repeat  Spacer    Distance  Consensus
                =====  ================  ==========  ==========  ========  ======  ======  ==========  =========
                    1  Ga0307431_100003       60530        1523        21      37      37              GTTT...AC
                >Ga0307431_1000385
                Array          Sequence    Position      Length  # Copies  Repeat  Spacer    Distance  Consensus
                =====  ================  ==========  ==========  ========  ======  ======  ==========  =========
                    5  Ga0307431_100038          56        1025        16      30      36              CTTCCAATCCTACCTATGAGGAATTGAAAT
                    6  Ga0307431_100038        1250         823        13      29      37         133  CTTCCAATCCTACCTATGAGGAATTGAAA
                >Ga0307431_1000754
                Array          Sequence    Position      Length  # Copies  Repeat  Spacer    Distance  Consensus
                =====  ================  ==========  ==========  ========  ======  ======  ==========  =========
                   11  Ga0307431_100075       18230         366         6      36      30              GCTGTGATAGACCTCGATTTGTGGGGTAGTAACAGC
                   12  Ga0307431_100075       20063         179         3      18      62        1437  ATGCCTAAGTATCTTAGG
                # out
                202155.assembled.fna_metacontig__Ga0307431_1000033	50828	72054	60828	62054	GTTT...AC	17	37	37	1
                """
                if line.startswith(">"):
                    contig = line[1:]
                    logger.debug(f'find contig line, the contig = {contig}')
                    continue
                elif line.startswith("Array") or line.startswith("====="):
                    logger.debug('find annotation line, skip')
                    continue
                else:
                    logger.debug('find a crispr array, try to parse it')
                # when find a crispr array
                # try to parse
                info = line.split()  # split all symbol than can not see '\t', ' ', '\n', et al.
                logger.debug(f'splitted line info is: info = {info}')

                if not info[0].isnumeric():
                    raise ValueError(f'info = {info}')
                logger.debug(f'contig = {contig}, info[1] = {info[1]}')
                # assert contig == info[1]  # check contig name?  # pilercr返回的info[1]可能是contig的name的截短
                array_index, _, crispr_start, crispr_length = info[:4]
                copies, repeat, spacer = info[4:7]
                seq = info[-1]  # 因为Distance有可能为空，所以不放在一起解包
                array_index = int(array_index)
                crispr_start = int(crispr_start)
                crispr_length = int(crispr_length)
                copies = int(copies)
                repeat = int(repeat)
                spacer = int(spacer)
                logger.debug(
                    f'array_index = {array_index}, contig = {contig}, '
                    f'crispr_start = {crispr_start}, crispr_length = {crispr_length}'
                )

                crispr_end = crispr_start + crispr_length + 1
                # extending
                locus_start = crispr_start - extend if crispr_start > extend else 1
                locus_end = crispr_end + extend
                # TODO 这里和原始代码的计算得出的坐标位置不同，仔细check每行计算哪里出了问题
                # write to loc_file
                line_out = (f"{os.path.basename(f_fa_input)}_metacontig_{contig}\t"  # 202155.assembled.fna_metacontig__Ga0307431_1000033
                     f"{locus_start}\t{locus_end}\t"  # 50828	72054
                     f"{crispr_start}\t{crispr_end}\t"  # 60828	62054
                     f"{seq}\t"  # GTTT...AC
                     f"{copies}\t{repeat}\t"  # 17	37  # Copies Repeat
                     f"{spacer}\t{array_index}\n"  # 37	1  # Spacer  Array
                     )
                f_loc.write(line_out)
    # -----------------------------
    # 3.get protein locs
    # -----------------------------
    logger.info('3.get protein locs')
    if tests[3]:
        with open(f_gff, 'rt') as f_gff_raw, open(f_tsv_pep_loc, 'wt') as f_loc:
            lines = f_gff_raw.readlines()
            lines = [line.rstrip() for line in lines if len(line.rstrip()) > 0]  # skip empty line, drop "\n"
            dt_loc = dict()

            for line in lines:
                info = line.split('\t')
                if info[2] == "CDS":
                    scaffold = f"{number}_{info[0]}"
                    ll = info[8]
                    ll = ll.replace(";", "\t")
                    ll = ll.replace("_", "\t")
                    b = ll.split("\t")
                    gene = f"{scaffold}_{b[1]}"
                    print(f"{scaffold}\t{info[3]}\t{info[4]}\t{info[6]}\t{gene}", file=B)
                    dt_loc[gene] = f"{scaffold}\t{info[3]}\t{info[4]}\t{info[6]}"
    # -----------------------------
    # 4.cas locs vs protein locs
    # -----------------------------
    logger.info('4.cas locs vs protein locs')
    if tests[4]:
        # 设置文件路径
        f_tsv_pep_loc = "temp.gff.loc"
        f_tsv_crispr_loc = "temp.spacer.loc"
        temp_bed = "temp.bed"

        # 使用bedtools进行交集操作
        os.system(f"bedtools intersect -wo -a {f_tsv_pep_loc} -b {f_tsv_crispr_loc} > {temp_bed}")

        # 打开并读取temp.bed文件
        with open(temp_bed, "r") as A:
            for line in A:
                line = line.rstrip()
                info = line.split()
                line = line.replace("\_metacontig\_\_", "\t")
                ar = line.split()
                good[info[4]] = 1
                cpscaffold[ar[1]] = 1
    # -----------------------------
    # 5.get cas locs protein
    # -----------------------------
    logger.info('5.get cas locs protein')
    if tests[5]:
        with open('temp.pep', 'r') as A, open('temp.pep.fasta', 'w') as B, open('temp.pep.filtered.fasta', 'w') as C:
            parse = 0
            for line in A:
                line = line.rstrip()
                line = re.sub(r'^>', f'>{number}_', line)
                if line.startswith(">"):
                    ll = line.replace(">", "")
                    info = ll.split()
                    print(f">{info[0]}\t{loc[info[0]]}", file=B)
                    if info[0] in good:
                        parse = 1
                        print(f">{info[0]}\t{loc[info[0]]}", file=C)
                    else:
                        parse = 0
                else:
                    print(line, file=B)
                    if parse > 0:
                        print(line, file=C)
    # -----------------------------
    # 6.get cas locas scaffold
    # -----------------------------
    logger.info('6.get cas locas scaffold')
    if tests[6]:
        with open('temp', 'r') as A, open('temp-crisper-scaffold.fa', 'w') as B:
            for line in A:
                ll = line.rstrip()
                ll = ll.replace(">", "")
                info = ll.split()
                try:
                    line2 = next(A)
                    if info[0] in cpscaffold:
                        print(f"{line.rstrip()}{line2}", file=B)
                except StopIteration:
                    break
    # -----------------------------
    # 7.save result files
    # -----------------------------
    logger.info('7.save result files')
    if tests[7]:
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

    logger.info('End, exit.')


if __name__ == '__main__':
    # casfinder(
    #     input="/Users/zhaohuanan/Downloads/test/202155.assembled.fna",
    #     output="/Users/zhaohuanan/Downloads/test/202155.assembled.fna.filtered.fasta",
    #     lmin=3001,  # 3001 in Nature Methods paper
    #     lmax=None
    # )
    pass
