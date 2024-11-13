import gzip
import sys

import pandas as pd

from bioat.exceptions import BioatFileFormatError, BioatFileNotCompleteError
from bioat.lib.libfastx import calculate_length_distribution, format_this_fastx
from bioat.logger import get_logger

__module_name__ = "bioat.fastxtools"


class FastxTools:
    """FASTA & FASTQ toolbox."""

    def __init__(self):
        self.fastx = None
        pass

    def fmt_this(
        self, file: str, new_file: str | None = None, force=False, log_level="WARNING"
    ):
        """Formats a FASTA file to improve readability.

        Args:
            file (str): The input filename for the FASTA file.
            new_file (str | None): The output filename. If None, the `file` will be replaced. Default is None.
            force (bool): If True, forces the formatting even if the output file exists. Default is False.
            log_level (str): The logging level for messages. Default is "WARNING".

        This function calls 'format_this_fastx' to perform the actual formatting on the specified FASTA file.
        """
        format_this_fastx(
            old_file=file, new_file=new_file, force=force, log_level=log_level
        )

    def plot_length_distribution(
        self,
        file: str,
        table: str | None = None,
        image: str | None = None,
        plt_show: bool = False,
        log_level="WARNING",
    ):
        """Plots the length distribution of a FASTA file.

        Args:
            file (str): The input filename for the FASTA file.
            table (str | None, optional): The output filename for the length distribution table. Default is <file>.lengths.
            image (str | None, optional): The output filename for the length distribution figure. Default is <file>.lengths.pdf.
            plt_show (bool, optional): If True, shows the plot. Default is False.
            log_level (str): The logging level for messages. Default is "WARNING".
        """
        calculate_length_distribution(
            file=file, table=table, image=image, plt_show=plt_show, log_level=log_level
        )

    def mgi_parse_md5(self, file: str, log_level="WARNING"):
        """Converts a mgi-like MD5 file into a standard MD5 file.

        Args:
            file (str): The name of the mgi-like MD5 file to read.
            log_level (str, optional): The logging level to use.
                It can be INFO, DEBUG, WARNING, or ERROR. The default is WARNING.

        """

        logger = get_logger(
            level=log_level,
            module_name=__module_name__,
            func_name="mgi_parse_md5",
        )
        try:
            df = pd.read_csv(
                file,
                header=None,
                index_col=False,
                sep="  ",
                engine="python",
            )
            df.columns = ["md5", "filename"]
            df.md5 = df.md5.map(lambda x: x.lower())
        except ValueError as e:
            logger.error(BioatFileFormatError(e))
            exit(1)
        df.filename = df.filename + ".gz"
        to_path = file.replace(".txt", "") + ".fix.md5"
        logger.info(f"write to {to_path}")
        df.to_csv(to_path, header=False, index=False, sep="\t")

    @staticmethod
    def _load_fastx_generator(file, log_level="WARNING"):
        """

        :param str file: path of input <fastq | fastq.gz | fastx | fastx.gz>
        :return: a generator for all reads:
            i.e. print(next(obj)) -> ['header', 'seq', 'info', 'quality']
        :rtype: generator
        """
        logger = get_logger(
            level=log_level,
            module_name=__module_name__,
            func_name="_load_fastx_generator",
        )
        f = open(file, "rt") if not file.endswith(".gz") else gzip.open(file, "rt")
        # FASTQ @  | FASTA >
        symbol = f.read(1)
        f.close()
        f = open(file, "rt") if not file.endswith(".gz") else gzip.open(file, "rt")

        if symbol == "@":
            logger.debug("detect FASTQ file")
            read = []
            line = f.readline().rstrip()

            while True:
                if not line:
                    break
                else:
                    if line.startswith("@"):
                        # header or not complete
                        n = len(read)
                        if n == 0:
                            # header, 第一次循环开始
                            read.append(line)  # add header!
                            line = f.readline().rstrip()  # next line
                            continue
                        elif n == 4:
                            # read == ['header', 'seq', 'info', 'quality']
                            yield read
                            read = []
                            read.append(line)
                            line = f.readline().rstrip()
                        elif n == 3:
                            """
                            @@IIEIBCE>IC<IBIIIIEAIEIEB<IDECCD6 # line! 期望它是 header！现在它是 quality
                            [
                                '@Beta12AdemL1C001R00100001768/1',
                                'ATCCCCGTATCTTCACCCCACCACAAACTATTAG',
                                '+',
                            ]'@@IIEIBCE>IC<IBIIIIEAIEIEB<IDECCD6'

                            """
                            read.append(line)
                            line = f.readline().rstrip()

                            if not line.startswith("@"):
                                raise ValueError("The file may be incomplete!")
                        else:
                            f.close()
                            raise ValueError("The file may be incomplete!")

                        # header line!
                        # read.append(line)  # add header !
                    else:
                        # not header line!
                        read.append(line)
                        line = f.readline().rstrip()

            # ls.append(read)
            yield read
            f.close()
            # return ls

        elif symbol == ">":
            logger.debug("detect FASTA file")
            read = []
            seq = ""
            line = f.readline().rstrip()

            while True:
                if not line:
                    f.close()
                    break
                else:
                    if line.startswith(">"):
                        # 读取 header！
                        n = len(read)

                        if n == 0:
                            # 第一次循环
                            read.append(line)
                            line = f.readline().rstrip()
                        elif n == 1:
                            # 已经有一个 header 了！现在缺 seq
                            read.append(seq)  # add seq line
                            # ls.append(read)
                            yield read
                            read = []  # 重置 read 这个 list
                            read.append(line)
                            line = f.readline().rstrip()
                            seq = ""  # 重置 seq 这个 str
                        else:
                            f.close()
                            logger.error(
                                BioatFileNotCompleteError("The file may be incomplete!")
                            )
                            exit(1)
                    else:
                        # 读取并添加 seq
                        seq += line
                        line = f.readline().rstrip()

            read.append(seq)
            # ls.append(read)
            yield read
            # return ls
        else:
            f.close()
            logger.error(
                BioatFileFormatError(
                    "Input line one must starts with `@` for FASTQ or `>` for FASTA!"
                )
            )
            exit(1)
        f.close()

    # CLI subcommand for filtering sequences
    def filter_read_contains_n(self, file: str, output=sys.stdout.name):
        """Filter reads that contain the 'N' base in FASTA or FASTQ formats.

        This function processes a given FASTA or FASTQ file and
        filters out reads that contain the 'N' base.
        The result is directed to an output file, or to stdout
        if no output file is specified.

        Args:
            file (str): The name of the FASTA or FASTQ file to be processed.
            output (str): The name of the output file. Defaults to stdout
                          if not specified, and the format matches the input.

        Returns:
            None
        """

        if self.fastx is None:
            self.fastx = self._load_fastx_generator(file)

        self.file = file
        output = sys.stdout if output == sys.stdout.name else output

        if isinstance(output, str):
            f = (
                gzip.open(output, "wt")
                if output.endswith(".gz")
                else open(output, "wt")
            )
        else:
            f = output

        fx = self._load_fastx_generator(file)

        for read in fx:
            seq = read[1]

            if seq.upper().__contains__("N"):
                continue
            else:
                f.write("\n".join(read) + "\n")
        f.close()


if __name__ == "__main__":
    pass
