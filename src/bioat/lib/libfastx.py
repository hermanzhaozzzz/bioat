import sys
import re
import subprocess
from Bio import SeqIO
from bioat.logger import get_logger

__module_name__ = "bioat.lib.libfastx"


def filter_fasta_length(input: str, output: str, lmin: int | None = None,
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


def casfinder():
    # 定义一些常量和变量
    path = '/Users/zhaohuanan/Downloads/test'
    file = f"{path}/202155.assembled.fna.filtered.fasta"
    crisper_scaffold_file = f"{file}.crisper.scaffold.fa"
    pep_file = f"{file}.pep.fasta"
    pep_cas_file = f"{file}.pep.cas.fasta"
    crispr_loc_file = f"{file}.crispr.loc"
    crispr_spacer_file = f"{file}.crispr.spacer"
    temp_dir = f"{path}/"  # 用于存储临时文件的目录

    # filtered_fasta = SeqIO.parse(file, "fasta")

    subprocess.check_call(
        ['/Users/zhaohuanan/micromamba/envs/snakepipes_Cas-mining/bin/prodigal', '-a', 'temp.pep', '-i', file, '-p', 'single', '-f', 'gff', '-o', 'temp.gff'],

    )



if __name__ == '__main__':
    # filter_fasta_length(
    #     input="/Users/zhaohuanan/Downloads/test/202155.assembled.fna",
    #     output="/Users/zhaohuanan/Downloads/test/202155.assembled.fna.filtered.fasta",
    #     lmin=3001,  # 3001 in Nature Methods paper
    #     lmax=None
    # )
    casfinder()