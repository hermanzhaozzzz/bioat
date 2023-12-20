"""Doc.
TODO
"""
import os
import subprocess
import sys

from Bio import SeqIO
from pybedtools import BedTool

from bioat.logger import get_logger

__module_name__ = "bioat.lib.libfastx"


def casfinder(
    input_fa: str,
    output_faa: str | None = None,
    lmin: int | None = None,
    lmax: int | None = None,
    extend: int = 10_000,
    temp_dir: str | None = None,
    prodigal: str | None = None,
    pilercr: str | None = None,
    log_level="DEBUG",
) -> None:
    # set logger
    logger = get_logger(
        level=log_level,
        module_name=__module_name__,
        func_name=sys._getframe().f_code.co_name,
    )
    workspace = os.getcwd()  # where I am

    if temp_dir is None:
        dirname = os.path.dirname(input_fa)
        if dirname in ("", "./"):
            temp_dir = workspace
        else:
            temp_dir = dirname
    else:
        pass

    fa_input = input_fa  # input_fa = '202155.assembled.fna'
    fa_pep_cas = os.path.join(
        workspace,
        f"{input_fa}.pep.cas.faa" if output_faa is not None else str(output_faa),
    )
    fa_filtered = os.path.join(temp_dir, f"{input_fa}.filtered.fa")
    f_pilercr = os.path.join(temp_dir, f"{input_fa}.crispr.spacer.pilercr")
    fa_pep = os.path.join(temp_dir, f"{input_fa}.pep.faa")
    gff = os.path.join(temp_dir, f"{input_fa}.gff")
    bed_pep = os.path.join(temp_dir, f"{input_fa}.pep.loc.bed")
    bed_crispr = os.path.join(temp_dir, f"{input_fa}.crispr.loc.bed")
    fa_crisper_scaffold = os.path.join(temp_dir, f"{input_fa}.crisper.scaffold.fa")
    bed_cas = os.path.join(temp_dir, f"{input_fa}.cas.loc.bed")
    assembly_id = f"{os.path.basename(fa_input)}"

    tests = {
        0: False,  # PASS # 0. filter contigs
        1: False,  # PASS # 1. cas & protein annotation
        2: False,  # PASS # 2. get crispr loci
        3: False,  # PASS # 3. get protein locs
        4: True,  # PASS # 4. cas locs vs protein locs
        5: False,  # PASS # 5. get cas locs protein
        6: False,  # PASS # 6. get cas loc as scaffold
        7: False,  # PASS # 7. save result files
    }
    # -----------------------------
    # 0. filter contigs
    # -----------------------------
    logger.info("0. filter contigs")
    if tests[0]:
        # filter length of assembly contigs
        logger.debug(f"filter length of contigs from file @ {fa_input}")
        assembly_input = SeqIO.parse(fa_input, "fasta")

        # define filter_func
        def filter_func(contig, lmin, lmax) -> bool:
            if lmin is None and lmax is not None:
                return len(contig) <= lmax
            elif lmin is not None and lmax is None:
                return len(contig) >= lmin
            else:
                return lmin <= len(contig) <= lmax

        contigs_input = (contig for contig in assembly_input)
        contigs_output = (
            contig for contig in contigs_input if filter_func(contig, lmin, lmax)
        )
        logger.debug(f"writing filtered contigs to file @ {fa_filtered}")
        SeqIO.write(contigs_output, fa_filtered, "fasta")
    # -----------------------------
    # 1. cas & protein annotation
    # -----------------------------
    logger.info("1. cas & protein annotation")
    if tests[1]:
        logger.debug("subprocess call for prodigal")
        subprocess.check_call(
            [
                prodigal,
                "-a",
                fa_pep,
                "-i",
                fa_filtered,
                "-p",
                "single",
                "-f",
                "gff",
                "-o",
                gff,
            ],
            stderr=open("/dev/null", "wt"),
        )
        logger.debug(f"call has returned, check output @ {fa_pep}, {gff}")
        logger.debug("subprocess call for pilercr")
        subprocess.check_call(
            [pilercr, "-in", fa_filtered, "-out", f_pilercr],
            stderr=open("/dev/null", "wt"),
        )
        logger.debug(f"call has returned, check output @ {f_pilercr}")
    # -----------------------------
    # 2.get crispr loci
    # -----------------------------
    logger.info("2.get crispr loci")
    if tests[2]:
        with open(f_pilercr, "rt") as wrapper_i, open(bed_crispr, "wt") as wrapper_o:
            lines = wrapper_i.readlines()
            lines = [
                line.rstrip() for line in lines if len(line.rstrip()) > 0
            ]  # skip empty line, drop "\n"

            parse_status = False

            for line in lines:
                if line.startswith("SUMMARY BY POSITION"):  # the first line to parse
                    logger.debug("find `SUMMARY BY POSITION part`")
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
                    logger.debug(f"find contig line, the contig = {contig}")
                    continue
                elif line.startswith("Array") or line.startswith("====="):
                    logger.debug("find annotation line, skip")
                    continue
                else:
                    logger.debug("find a crispr array, try to parse it")
                # when find a crispr array
                # try to parse
                info = (
                    line.split()
                )  # split all symbol than can not see '\t', ' ', '\n', et al.
                logger.debug(f"splitted line info is: info = {info}")

                if not info[0].isnumeric():
                    raise ValueError(f"info = {info}")  # TODO
                logger.debug(f"contig = {contig}, info[1] = {info[1]}")
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
                    f"array_index = {array_index}, contig = {contig}, "
                    f"crispr_start = {crispr_start}, crispr_length = {crispr_length}"
                )

                crispr_end = crispr_start + crispr_length + 1
                # extending
                locus_start = crispr_start - extend if crispr_start > extend else 1
                locus_end = crispr_end + extend
                # write to loc_file
                line_out = (
                    f"assembly_id@{assembly_id};contig_id@{contig}\t"  # chrom
                    f"{locus_start}\t"  # chromStart
                    f"{locus_end}\t"  # chromEnd
                    f"crispr_start@{crispr_start};"
                    f"crispr_end@{crispr_end};"
                    f"crispr_extend_start@{locus_start};"
                    f"crispr_extend_end@{locus_end};"
                    f"repeat_seq@{seq};"
                    f"n_copies@{copies};"
                    f"n_repeats@{repeat};"
                    f"n_spacers@{spacer};"
                    f"crispr_loci_id@{array_index}\t"  # name
                    "0\t"  # score
                    ".\n"  # strand
                )
                wrapper_o.write(line_out)
    # -----------------------------
    # 3.get protein locs
    # -----------------------------
    logger.info("3.get protein locs")
    if tests[3]:
        with open(gff, "rt") as wrapper_i, open(bed_pep, "wt") as wrapper_o:
            lines = wrapper_i.readlines()
            # skip annotation lines wiht # and blank lines
            lines = (
                line.rstrip()
                for line in lines
                if len(line.rstrip()) > 0 and not line.startswith("#")
            )

            for line in lines:
                info = line.split("\t")
                # TODO 这里的start和stop可能需要+-1？
                (
                    contig,
                    source,
                    feature,
                    start,
                    stop,
                    score,
                    strand,
                    phase,
                    atributes,
                ) = info
                _temp = [i.split("=") for i in atributes.strip().split(";") if i != ""]
                atributes = {k.strip(): v.strip() for k, v in _temp}
                atributes_in_bed = "|".join([f"{k}={v}" for k, v in atributes.items()])

                if feature != "CDS":
                    continue

                idx = atributes["ID"].split("_")[-1]
                start, stop, strand = info[3], info[4], info[6]
                # write to loc_file
                line_out = (
                    f"assembly_id@{assembly_id};contig_id@{contig}\t"  # chrom
                    f"{start}\t"  # chromStart
                    f"{stop}\t"  # chromEnd
                    f"CDS_id@{assembly_id}_{contig}_{idx};"
                    f"CDS_start@{start};"
                    f"CDS_end@{stop};"
                    f"CDS_predict_score@{score};"
                    f"CDS_strand@{strand};"
                    f"CDS_phase@{phase};"
                    f"Predigal_info@{atributes_in_bed}\t"  # name
                    f"{score}\t"
                    f"{strand}\n"
                )
                wrapper_o.write(line_out)
    # -----------------------------
    # 4.cas locs vs protein locs
    # -----------------------------
    logger.info("4.cas locs vs protein locs")
    if tests[4]:
        # 使用bedtools进行交集操作
        """
        -wo
        Write the original A and B entries plus the number of base pairs of overlap
            between the two features.
            - Overlaps restricted by -f and -r.
            - Only A features with overlap are reported.
        """
        bed_a = BedTool(fn=bed_pep)
        bed_b = BedTool(fn=bed_crispr)
        bed_intersect = bed_a.intersect(bed_b, wo=True)
        bed_intersect.moveto(bed_cas)
        logger.debug(f"generate cas.bed file, check output @ {bed_cas}")
        # bed cas is now ready for fetch fasta of pep cas!
        exit()  # TODO
        pick_this_crispr_scaffold = dict()
        pick_this_cas = dict()
        with open(f_bed_cas_loc, "rt") as bed:
            bedlines = bed.readlines()
            bedlines = [i.rstrip() for i in bedlines]
            for line in bedlines:
                info = line.split("\t")
                print(f"info = {info}")
                scaffold, start, stop, strand, gene_name = info[:5]
                this_gene_name = gene_name.split("_metacontig__")[-1]
                this_contig = scaffold.split("_metacontig__")[-1]
                pick_this_cas[this_gene_name] = {
                    "contig": this_contig,
                    "gene": this_gene_name,
                    "start": start,
                    "stop": stop,
                    "strand": strand,
                }
                # Ga0307431_1000001_1
                pick_this_crispr_scaffold[this_contig] = True
                # Ga0307431_1000033
    # -----------------------------
    # 5.get cas locs protein
    # -----------------------------
    logger.info("5.get cas locs protein")
    if tests[5]:
        assembly_input = SeqIO.parse(f_faa_pep, "fasta")
        contigs_input = (contig for contig in assembly_input)
        contigs_output = []

        for contig in contigs_input:
            print(contig.__dir__())
            header = contig.id  # header = Ga0307431_1014098_3

            if header in pick_this_cas.keys():
                info = pick_this_cas[header]
                contig.id = f"{info['gene']}\t{info['contig']}\t{info['start']}\t{info['stop']}\t{info['strand']}"
                contigs_output.append(contig)
        SeqIO.write(contigs_output, f_faa_pep_cas, "fasta")
    # -----------------------------
    # 6.get cas locas scaffold
    # -----------------------------
    logger.info("6.get cas locas scaffold")
    if tests[6]:
        with open("temp", "r") as A, open("temp-crisper-scaffold.fa", "w") as B:
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
    logger.info("7.save result files")
    if tests[7]:
        os.system(f"cat temp-crisper-scaffold.fa >> {file}.crisper.scaffold.fa")
        os.system(f"cat temp.pep.fasta >> {file}.pep.fasta")
        os.system(f"cat temp.pep.filtered.fasta >> {file}.pep.cas.fasta")
        os.system(f"cat temp.spacer.loc >> {file}.crispr.loc")
        os.system(f"cat temp.spacer >> {file}.crispr.spacer")
        os.system("rm temp*")
        counter = 0

    logger.info("End, exit.")


if __name__ == "__main__":
    # casfinder(
    #     input="/Users/zhaohuanan/Downloads/test/202155.assembled.fna",
    #     output="/Users/zhaohuanan/Downloads/test/202155.assembled.fna.filtered.fasta",
    #     lmin=3001,  # 3001 in Nature Methods paper
    #     lmax=None
    # )
    pass
