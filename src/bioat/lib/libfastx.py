"""Doc.
TODO
"""
import os
import subprocess
import sys
import re
import gzip

import pandas as pd
import numpy as np
from Bio import SeqIO, bgzf
from pybedtools import BedTool
from bioat.logger import get_logger

__module_name__ = "bioat.lib.libfastx"


def filter_fasta_length(contig, lmin, lmax) -> bool:
    if lmin is None and lmax is not None:
        return len(contig) <= lmax
    elif lmin is not None and lmax is None:
        return len(contig) >= lmin
    else:
        return lmin <= len(contig) <= lmax


def cas_finder(
        input_fa: str,
        output_faa: str | None = None,
        output_crispr_info_tab: str | None = None,
        lmin: int | None = None,
        lmax: int | None = None,
        extend: int = 10_000,
        temp_dir: str | None = None,
        prodigal: str | None = None,
        pilercr: str | None = None,
        rm_temp: bool = True,
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
    fa_pep_cas = os.path.join(  # output_faa, 不指定 | 指定
        workspace,
        f"{input_fa}.pep.cas.faa" if output_faa is None else str(output_faa),
    )
    dirname = os.path.dirname(fa_pep_cas)
    os.makedirs(dirname, exist_ok=True)
    # temp files
    fa_crisper_scaffold = os.path.join(temp_dir, f"{fa_pep_cas}.crisper.scaffold.fa")
    f_pilercr = os.path.join(temp_dir, f"{fa_pep_cas}.crispr.spacer.pilercr")
    fa_pep = os.path.join(temp_dir, f"{input_fa}.pep.faa")
    fa_filtered = os.path.join(temp_dir, f"{input_fa}.filtered.fa")
    gff = os.path.join(temp_dir, f"{input_fa}.gff")
    bed_pep = os.path.join(temp_dir, f"{input_fa}.pep.loc.bed")
    bed_crispr = os.path.join(temp_dir, f"{input_fa}.crispr.loc.bed")
    bed_cas = os.path.join(temp_dir, f"{input_fa}.cas.loc.bed")
    # info
    assembly_id = f"{os.path.basename(fa_input)}"

    tests = {
        0: True,  # PASS # 0. filter contigs
        1: True,  # PASS # 1. cas & protein annotation
        2: True,  # PASS # 2. get crispr loci
        3: True,  # PASS # 3. get protein cds
        4: True,  # PASS # 4. cas loci vs protein cds
        5: True,  # PASS # 5. get cas faa
        6: True,  # PASS # 6. get crispr loci as scaffold
    }
    # -----------------------------
    # 0. filter contigs
    # -----------------------------
    if tests[0]:
        logger.info("0.filter contigs")
        # filter length of assembly contigs
        logger.debug(f"filter length of contigs from file @ {fa_input}")
        handler = gzip.open(fa_input, "rt") if fa_input.endswith(".gz") else open(fa_input, "rt")
        contigs_input = SeqIO.parse(handler, "fasta")
        contigs_output = (
            contig
            for contig in contigs_input
            if filter_fasta_length(contig, lmin, lmax)
        )
        logger.debug(f"writing filtered contigs to file @ {fa_filtered}")
        # don't consider gz file. because it is a temp file.
        SeqIO.write(contigs_output, fa_filtered, "fasta")
        handler.close()

        # if fa_filtered is empty, just return empty fa_pep_cas
        logger.debug(f"checking if @ {fa_filtered} is empty.")

        if os.path.getsize(fa_filtered) == 0:
            # don't consider gz file. because it is a temp file.
            logger.warning(f"fa_filtered @ {fa_filtered} is empty! will touch an empty output file @ {fa_pep_cas}")
            with open(fa_pep_cas, 'wt') as f:
                f.write('')
            if output_crispr_info_tab:
                with open(output_crispr_info_tab, 'wt') as f:
                    f.write('')
            logger.info("End, exit.")
            return  # just return output file as an empty file
    # -----------------------------
    # 1. cas & protein annotation
    # -----------------------------
    if tests[1]:
        logger.info("1.cas & protein annotation")
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
            stderr=open("/dev/null", "wt") if log_level != "DEBUG" else sys.stderr,
        )
        logger.debug(f"call has returned, check output @ {fa_pep}, {gff}")
        # the output .pep.faa's header: >name # start # stop # 1/-1 # prodigalInfo
        # scaffold_1025_6
        # the output .gff
        # scaffold_1025 ID=1_1
        # scaffold_1025 ID=1_2
        # scaffold_1303 ID=2_1
        # ...
        logger.debug("subprocess call for pilercr")
        subprocess.check_call(  # don't export help doc for output crispr info
            [pilercr, "-noinfo", "-in", fa_filtered, "-out", f_pilercr],
            stderr=open("/dev/null", "wt") if log_level != "DEBUG" else sys.stderr,
        )
        logger.debug(f"call has returned, check output @ {f_pilercr}")
        # if fa_filtered is empty, just return empty fa_pep_cas
        logger.debug(f"checking if @ {fa_filtered} is empty.")
    # -----------------------------
    # 2.get crispr loci
    # -----------------------------
    if tests[2]:
        logger.info("2.get crispr loci")
        with open(f_pilercr, "rt") as wrapper_i, open(bed_crispr, "wt") as wrapper_o:
            lines = wrapper_i.readlines()
            lines = [
                line.rstrip() for line in lines if len(line.rstrip()) > 0
            ]  # skip empty line, drop "\n"

            parse_status = False
            parse_detail = False
            ls_detail = []

            for line in lines:
                # parse detail info
                if line.startswith("DETAIL REPORT"):  # the first line to parse detail info
                    logger.debug("find `DETAIL REPORT part`")
                    parse_detail = True
                if line.startswith("SUMMARY BY SIMILARITY"):  # the end of detail info
                    parse_detail = False
                if parse_detail:
                    ls_detail.append(line)
                # parse repeat info
                if line.startswith("SUMMARY BY POSITION"):  # the first line to parse repeat info
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
                    # 2024-01-12 fix bug, if contig with space
                    # contig = contig.replace(" ", "_")
                    contig = contig.split(" ")[0]
                    logger.debug(f"find contig line, the contig = {contig}")
                    continue
                elif line.startswith("Array") or line.startswith("====="):
                    logger.debug("find annotation line, skip")
                    continue
                else:
                    logger.debug("find a crispr array, try to parse it")
                # when find a crispr array
                # try to parse
                # 2024-01-12 fix Sequence have space bug: invalid literal for int() with base 10: 'T'
                line = line[:7] + line[7:23].replace(" ", "_") + line[23:]
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
                    f"crispr_id@{assembly_id}+{contig}+{array_index};"
                    f"crispr_start@{crispr_start};"
                    f"crispr_end@{crispr_end};"
                    f"crispr_extend_start@{locus_start};"
                    f"crispr_extend_end@{locus_end};"
                    f"repeat_seq@{seq};"
                    f"n_copies@{copies};"
                    f"n_repeats@{repeat};"
                    f"n_spacers@{spacer}\t"
                    "0\t"  # score
                    ".\n"  # strand
                )
                wrapper_o.write(line_out)

            if output_crispr_info_tab:
                logger.info(f'2.additional job: writing crispr info table to @ {output_crispr_info_tab}')
                if not len(ls_detail) >= 5:
                    # is empty
                    with open(output_crispr_info_tab, 'wt') as f:
                        f.write('')

                ls_detail = ls_detail[1:]
                ls_all = []
                ls_temp = None
                ls_rep_info = []

                for info in ls_detail:
                    info = info.strip()
                    if info.startswith('Array'):
                        array_id = info.replace('Array ', '')
                        if ls_temp:
                            ls_all.append(ls_temp)
                        continue
                    elif info.startswith('>'):
                        contig_id = info[1:]
                        ls_temp = []
                        continue
                    elif info.startswith('Pos'):
                        continue
                    elif info.startswith('='):
                        continue
                    elif len(info.split()) == 4:
                        repeat_seq_rep = info.split()[3]
                        repeat_length_rep = info.split()[1]
                        spacer_length_rep = info.split()[2]
                        ls_rep_info.append([repeat_seq_rep, repeat_length_rep, spacer_length_rep])
                        continue
                    ls_temp.append([f'{assembly_id}+{contig_id}+{array_id}', info])
                ls_all.append(ls_temp)

                ls_df = []

                if not ls_all or ls_all == [None]:
                    # don't consider gz file. because it is a temp file.
                    logger.warning(
                        f'find no crispr info in @ {fa_filtered}! will touch an empty output file @ {fa_pep_cas}')
                    with open(fa_pep_cas, 'wt') as f:
                        f.write('')
                    if output_crispr_info_tab:
                        with open(output_crispr_info_tab, 'wt') as f:
                            f.write('')
                    logger.info("End, exit.")
                    return  # just return output file as an empty file

                for idx, one in enumerate(ls_all):
                    ls_one_info = [i[:1] + ls_rep_info[idx] + i[1].split() for i in one]
                    if len(ls_one_info[0]) == len(ls_one_info[-1]):
                        pass
                    else:
                        ls_one_info[-1].insert(7, np.NaN)
                    df = pd.DataFrame(ls_one_info, columns=[
                        'crispr_id',
                        'representative_repeat_seq',
                        'representative_repeat_length',
                        'representative_spacer_length',
                        'position',
                        'repeat_length',
                        '%identity',
                        'spacer_length',
                        'left_flank',
                        'repeat_mismatch',
                        'spacer_seq'
                    ])
                    ls_df.append(df)

                df_all = pd.concat(ls_df)
                df_all = df_all[[
                    'crispr_id',
                    'representative_repeat_seq',
                    'repeat_mismatch',
                    'representative_repeat_length',
                    'repeat_length',
                    'spacer_seq',
                    'representative_spacer_length',
                    'spacer_length',
                    'position',
                    '%identity',
                    'left_flank',
                ]]
                df_all.to_csv(output_crispr_info_tab, index=False)
    # -----------------------------
    # 3.get protein locs
    # -----------------------------
    if tests[3]:
        logger.info("3.get protein loci")
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
                # DEBUG 这里的start和stop可能需要+-1？
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
                idx = atributes["ID"].split("_")[-1]
                atributes["ID"] = idx
                atributes_in_bed = "|".join([f"{k}={v}" for k, v in atributes.items()])

                if feature != "CDS":
                    continue

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
    # 4.cas locus vs protein locus
    # -----------------------------
    if tests[4]:
        logger.info("4.cas locus vs protein locus")
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
    # -----------------------------
    # 5.get cas candidate proteins
    # -----------------------------
    if tests[5]:
        logger.info("5.get cas candidate proteins")
        # -----------------------------
        # special part for a dict info
        # -----------------------------
        # load keys from cas bed
        with open(bed_cas, "rt") as f:
            contigs_cas = f.readlines()
        cas_loc_info = {}

        for loc_info in contigs_cas:
            dt_this = {}
            this_loc_info = loc_info.rstrip().split("\t")
            assert this_loc_info[0] == this_loc_info[6]
            for temp in this_loc_info[0].split(";"):
                k, v = temp.split("@")
                dt_this[k] = v

            for temp in this_loc_info[3].split(";"):
                k, v = temp.split("@")
                dt_this[k] = v

            for temp in this_loc_info[9].split(";"):
                k, v = temp.split("@")
                dt_this[k] = v

            dt_this["overlapped_base_count"] = this_loc_info[12]
            # cas_loc_info[dt_this["CDS_id"].split(".fna_")[-1]] = dt_this

            for temp in dt_this['Predigal_info'].split('|'):
                k, v = temp.split("=")
                if k == 'ID':
                    CDS_id = v
                    break
            cas_loc_info[f'{dt_this["contig_id"]}_{CDS_id}'] = dt_this

        crispr_ids = ["_".join(i.split("_")[:-1]) for i in cas_loc_info.keys()]

        # logger.debug(
        #     f'cas_loc_info = \n{cas_loc_info}\n'
        # )

        contigs_input = SeqIO.parse(fa_pep, "fasta")
        contigs_output = []

        for contig in contigs_input:
            header = contig.id  # header = Ga0307431_1014098_3
            contig.description = contig.description.split('# ')[-1]

            if header in cas_loc_info.keys():
                info = cas_loc_info[header]
                contig.id = (
                        info["CDS_id"]
                        + " "
                        + ";".join([f"{k}@{v}" for k, v in info.items()])
                )
                contigs_output.append(contig)
        handler = gzip.open(fa_pep_cas, 'wt') if fa_pep_cas.endswith(".gz") else open(fa_pep_cas, 'wt')
        SeqIO.write(contigs_output, handler, "fasta")
        handler.close()
        logger.debug(f"generate cas.faa file, check output @ {fa_pep_cas}")
    # -----------------------------
    # 6.get cas candidate scaffolds
    # -----------------------------
    if tests[6]:
        logger.info("6.get cas candidate scaffolds")
        handler = gzip.open(fa_input, "rt") if fa_input.endswith(".gz") else open(fa_input, "rt")
        assembly_input = SeqIO.parse(handler, "fasta")
        contigs_output = []

        for contig in assembly_input:
            if contig.id in crispr_ids:
                contigs_output.append(contig)
        SeqIO.write(contigs_output, fa_crisper_scaffold, "fasta")
        handler.close()
        logger.debug(
            f"generate crispr_scaffold.fa file, check output @ {fa_crisper_scaffold}"
        )

    if rm_temp:
        logger.info(f"removing temp files @ {temp_dir}")
        os.remove(fa_crisper_scaffold)
        os.remove(f_pilercr)
        os.remove(fa_pep)
        os.remove(fa_filtered)
        os.remove(gff)
        os.remove(bed_pep)
        os.remove(bed_crispr)
        os.remove(bed_cas)

    logger.info("End, exit.")


def cas13_finder(
        input_faa: str,
        output_faa: str | None = None,
        lmin: int | None = None,
        lmax: int | None = None,
        log_level="DEBUG",
) -> None:
    # set logger
    logger = get_logger(
        level=log_level,
        module_name=__module_name__,
        func_name=sys._getframe().f_code.co_name,
    )
    workspace = os.getcwd()  # where I am

    fa_input = input_faa  # input_fa = '202155.assembled.fna'
    fa_pep_cas13 = os.path.join(
        workspace,
        f"{input_faa}.with.HEPNs.faa" if output_faa is None else str(output_faa),
    )
    dirname = os.path.dirname(fa_pep_cas13)
    os.makedirs(dirname, exist_ok=True)

    # start to parse HEPN pattern
    handler = gzip.open(fa_input, "rt") if fa_input.endswith(".gz") else open(fa_input, "rt")
    fa_cases = SeqIO.parse(handler, format="fasta")
    fa_cases_filter_length = [
        cas for cas in fa_cases if filter_fasta_length(cas, lmin, lmax)
    ]
    handler.close()

    # pattern = re.compile(r"R[NHD][A-Z]{3,5}H")
    # pattern = re.compile(r"R[NHQD][A-Z]{3,5}H")
    pattern = re.compile(r"(R[NHQD][A-Z]{3,5}H|H[A-Z]{3,5}[NHQD]R)")
    # R: Arginine
    # N: Asparagine
    # H: Histidine
    # Q: Glutarnine
    # D: Asparticacid
    # the pattern is based on three papers:
    # 1.Anantharaman, V., Makarova, K. S., Burroughs, A. M., Koonin, E. V. & Aravind, L. Comprehensive analysis of the HEPN superfamily: identification of novel roles in intra-genomic conflicts, defense, pathogenesis and RNA processing. Biol. Direct 8, 15 (2013).
    # see ![](http://_pic.zhaohuanan.cc:7777/images/2024/01/08/20240108002122.png)
    # 2.Xu, C. et al. Programmable RNA editing with compact CRISPR–Cas13 systems from uncultivated microbes. Nat. Methods 18, 499–506 (2021).
    # see ![](http://_pic.zhaohuanan.cc:7777/images/2024/01/08/20240108004055.png)
    # 3. what? for figure below
    # see ![](http://_pic.zhaohuanan.cc:7777/images/2023/12/24/20231224215827.png)

    ls_fa_cases_out = []

    for fa_cas in fa_cases_filter_length:
        # matches = pattern.findall(string=str(fa_cas.seq))
        matches = pattern.finditer(string=str(fa_cas.seq))

        ls_span = []

        for match in matches:
            dt_span = dict()
            dt_span["span"] = match.span()
            dt_span["target"] = match.string[match.start(): match.end()]
            ls_span.append(dt_span)

        if len(ls_span) > 0:
            logger.debug("↓" * 20)
            logger.debug(f"fa_cas.id = {fa_cas.id}")
            logger.debug(f"fa_cas.description = {fa_cas.description}")
            logger.debug(f"fa_cas.seq = {fa_cas.seq}")
            logger.debug(f"fa_cas.length = {len(fa_cas)}")
            logger.debug(f"find pattern = r'{match.re.pattern}'")
            logger.debug(f"HEPN info = {ls_span}")

            for i, dt in enumerate(ls_span):
                span = dt["span"]
                target = dt["target"]
                s = f";HEPN-{i + 1}@[span={span},target={target}]"
                fa_cas.description += s
            fa_cas.description += f";HEPN-count@{len(ls_span)}"
            logger.debug(f"new fa_cas.description = {fa_cas.description}")
            ls_fa_cases_out.append(fa_cas)
    handler = gzip.open(fa_pep_cas13, 'wt') if fa_pep_cas13.endswith(".gz") else open(fa_pep_cas13, 'wt')
    SeqIO.write(ls_fa_cases_out, handler, "fasta")
    logger.info("End, exit.")


def format_this_fastx(file: str, new_file: str | None = None, log_level: str = "DEBUG"):
    logger = get_logger(
        level=log_level,
        module_name=__module_name__,
        func_name=sys._getframe().f_code.co_name,
    )
    logger.debug("start to format fastx file")
    f_input = file
    temp_file = f".bioat_temp_{f_input}"
    new_file = f_input if new_file is None else new_file

    handler = gzip.open(f_input, 'wt') if f_input.endswith(".gz") else open(f_input, 'wt')
    fasta = SeqIO.parse(handler, "fasta")

    iter_out = (i for i in fasta)
    handler_out = gzip.open(temp_file, 'wt') if temp_file.endswith(".gz") else open(temp_file, 'wt')
    SeqIO.write(iter_out, handler_out, "fasta")
    handler.close()
    handler_out.close()
    logger.debug(f"rename {temp_file} to {new_file}")
    os.rename(temp_file, new_file)
    logger.debug("Done. Exit")


if __name__ == "__main__":
    pass
