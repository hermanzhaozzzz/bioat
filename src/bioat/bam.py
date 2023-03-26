from __future__ import absolute_import

import logging
import gzip
import os
import sys
import string
import random
from multiprocessing import Process
from bioat.logger import set_logging_level


class Bam:
    def __init__(self):
        pass

    def mpileup_to_table(
            self,
            mpileup: str,
            output: str = sys.stdout,
            threads: int = os.cpu_count() - 1,
            mutation_number_threshold: int = None,
            temp_dir: str = f"/tmp/bioat_{''.join(random.sample(string.ascii_letters + string.digits, 16))}",
            remove_temp: bool = True,
            log_level:str = 'INFO'
    ):
        """Convert mpileup file to info file.

        :param mpileup: samtools mpileup format file
        :param output: Output parsed file
        :param threads: Multiple threads number
        :param mutation_number_threshold: Only contain mutation info go to the output, set 0 mean output all site
        :param temp_dir: Where to keep temp files, default is /tmp
        :param remove_temp: Where to keep temp files, default is /tmp
        :param log_level: 'CRITICAL', 'ERROR', 'WARNING', 'INFO', 'DEBUG', 'NOTSET'
        """
        set_logging_level(level=log_level)

        def _temp_split_bam(mpileup, threads, temp_dir):
            # creat temp file name and open file
            logging.debug("Making temp files...")

            input_file_basename = os.path.basename(mpileup)

            temp_filename_list = []
            temp_file_list = []

            os.makedirs(temp_dir, exist_ok=True)

            for index in range(threads):
                # make filename
                temp_file_basename = "temp_" + input_file_basename + "." + str(index) + "." + "".join(
                    random.sample(string.ascii_letters + string.digits, 16)) + '.gz'
                temp_file_name = os.path.join(temp_dir, temp_file_basename)
                temp_filename_list.append(temp_file_name)

                # open file
                temp_file = gzip.open(temp_file_name, "wt")
                temp_file_list.append(temp_file)

            # counting input file line number
            logging.debug("Counting input file...")

            total_input_line_num = 0
            input_file = gzip.open(mpileup, "rt") if mpileup.endswith('.gz') else open(mpileup, "rt")

            # count total_input_line_num
            for _ in input_file:
                total_input_line_num += 1

            if total_input_line_num % threads == 0:
                each_file_line_num = (total_input_line_num // threads)
            else:
                each_file_line_num = (total_input_line_num // threads) + 1

            input_file.close()
            logging.debug("Done!")

            # split temp files
            # write into output filename
            input_file = gzip.open(mpileup, "rt") if mpileup.endswith('.gz') else open(mpileup, "rt")

            for index, line in enumerate(input_file):
                file_index = index // each_file_line_num
                temp_file_list[file_index].write(line)

            # close output filename
            input_file.close()
            [temp_file.close() for temp_file in temp_file_list]
            logging.debug("Make temp files done!")

            return temp_filename_list

        def _merge_out_bmat(temp_filename_list, output, remove_temp=True):
            if type(output) is str:
                output = gzip.open(output, "wt") if output.endswith('.gz') else open(output, "wt")
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
                "mut_num"
            ]
            output.write("\t".join(header) + "\n")

            # open temp bmat file
            for temp_filename in temp_filename_list:
                # write file
                with gzip.open(temp_filename + ".bmat.gz", "rt") as temp_file:
                    for line in temp_file:
                        output.write(line)

            # remove temp files
            if remove_temp:
                for temp_filename in temp_filename_list:
                    os.remove(temp_filename)
                    os.remove(temp_filename + ".bmat.gz")
                os.rmdir(temp_dir)

        # make temp files
        temp_filename_list = _temp_split_bam(mpileup=mpileup, threads=threads, temp_dir=temp_dir)
        # multiple processing part
        logging.debug("Parsing files...")

        procs_list = []

        for index, temp_mpileup in enumerate(temp_filename_list):
            temp_output = temp_mpileup + ".bmat.gz"
            # add multiprocess
            sub_proc = Process(target=_run_mpileup_to_table,
                               args=(temp_mpileup, temp_output, mutation_number_threshold,))
            procs_list.append(sub_proc)
            sub_proc.start()

        for sub_proc in procs_list:
            sub_proc.join()
        logging.debug("Done!")
        # merge output files
        logging.debug("Merging files...")
        _merge_out_bmat(temp_filename_list=temp_filename_list, output=output, remove_temp=remove_temp)
        logging.debug("Done!...")


def _run_mpileup_to_table(temp_mpileup, temp_output, mutation_number_threshold):
    # open input file
    temp_mpileup = gzip.open(temp_mpileup, "rt") if temp_mpileup.endswith('.gz') else open(temp_mpileup, "rt")
    # open output file
    temp_output = gzip.open(temp_output, "wt") if temp_output.endswith('.gz') else open(temp_output, "wt")

    # parse mpileup file
    for line in temp_mpileup:
        data = line.strip().split('\t')
        chr_name = data[0]
        bp = data[1]
        bases = data[4].upper()
        ref = data[2].upper()
        types = {'A': 0, 'G': 0, 'C': 0, 'T': 0, '-': [], '+': [], 'Not': []}

        # coverage > 0
        i = 0
        while i < len(bases):
            base = bases[i]

            if base == '^':
                i += 2

            elif base == '$':
                i += 1

            elif base == '-':
                i += 1
                del_str_num = ""

                while bases[i].isdigit():
                    del_str_num += bases[i]
                    i += 1

                delNum = int(del_str_num)
                delSeq = ""
                for a in range(delNum):
                    delSeq += bases[i]
                    i += 1
                types['-'].append(delSeq)

            elif base == '*':
                types['-'].append(bases[i])
                i += 1

            elif base == '+':
                i += 1
                add_str_num = ""

                while bases[i].isdigit():
                    add_str_num += bases[i]
                    i += 1

                addNum = int(add_str_num)
                addSeq = ""
                for a in range(addNum):
                    addSeq += bases[i]
                    i += 1
                types['+'].append(addSeq)

            elif (base == '.') or (base == ','):
                types[ref] += 1
                i += 1

            else:
                if base in types:
                    types[base] += 1
                else:
                    types['Not'].append(base)
                i += 1

        adds = '.'
        adds_count = 0
        if len(types['+']) > 0:
            adds = ','.join(types['+'])
            adds_count = len(types['+'])

        dels = "."
        dels_count = 0
        if len(types['-']) > 0:
            dels = ",".join(types['-'])
            dels_count = len(types['-'])

        amb = '.'
        amb_count = 0
        if len(types['Not']) > 0:
            amb = ','.join(types['Not'])
            amb_count = len(types['Not'])

        # get other mutation number
        if ref in types:
            mut_num = types["A"] + types["T"] + types["G"] + types["C"] - types[ref]
        else:
            mut_num = 0

        out_list = [
            chr_name,
            bp,
            ref,
            types['A'],
            types['G'],
            types['C'],
            types['T'],
            dels_count,
            adds_count,
            amb_count,
            dels,
            adds,
            amb,
            mut_num
        ]

        # make a filter
        if mutation_number_threshold:
            if mut_num >= mutation_number_threshold:
                temp_output.write("\t".join(map(str, out_list)) + "\n")
        else:
            temp_output.write("\t".join(map(str, out_list)) + "\n")

    temp_mpileup.close()
    temp_output.close()


if __name__ == '__main__':
    pass
