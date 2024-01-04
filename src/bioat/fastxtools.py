import gzip
import sys
import pandas as pd
from bioat.exceptions import BioatFileFormatError, BioatFileNotCompleteError
from bioat.logger import get_logger
from bioat.lib.libfastx import format_this_fastx

__module_name__ = "bioat.fastxtools"


class FastxTools:
    """FASTA & FASTQ toolbox."""

    def __init__(self):
        self.fastx = None
        pass

    def fmt_this(self, file: str, new_file: str | None = None, log_level="WARNING"):
        """format fasta file and make it easier to see.

        :param file: input filename, fasta
        :param new_file: output filename, default is None, and when it is None, `file` will be replaced
        :param log_level: like others
        """
        format_this_fastx(file=file, new_file=new_file, log_level=log_level)

    def mgi_parse_md5(self, file: str, log_level="WARNING"):
        """Read mgi-like md5 file and convert to a normal md5 file.

        :param file: file name of a mgi-like md5 file.
        :param log_level: INFO/DEBUG/WARNING/ERROR, default is WARNING.
        """
        logger = get_logger(
            level=log_level,
            module_name=__module_name__,
            func_name=sys._getframe().f_code.co_name,
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

    def check_completeness(self, file: str, fmt="FASTQ", log_level="WARNING"):
        """正在开发中

        :param file: :param str file: path of input <fastq | fastq.gz | fastx | fastx.gz>
        :param fmt: FASTQ | FASTA
        :param log_level:
        :return: str, PASS | FILENAME_LOG_FAIL
        """
        logger = get_logger(
            level=log_level,
            module_name=__module_name__,
            func_name=sys._getframe().f_code.co_name,
        )
        pass

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
            func_name=sys._getframe().f_code.co_name,
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

    def get_headers(self):
        # if self.fastx is None:
        self.fastx = self._load_fastx_generator()

        return [i[0][1:] for i in self.fastx]

    def get_sequence(self):
        # if self.fastx is None:
        self.fastx = self._load_fastx_generator()

        return [i[1] for i in self.fastx]

    def get_record_dict(self):
        # if self.fastx is None:
        self.fastx = self._load_fastx_generator()

        return {i[0][1:]: i[1] for i in self.fastx}

    # cli subcmd
    def filter_read_contains_N(self, file: str, output=sys.stdout.name):
        """Filter read contains N base in FASTA or FASTQ



        :param file: FASTA | FASTQ file name.
        :param output: FASTA | FASTQ file name file, stdout if not assigned, format as input.
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
    # fastx = FastxTools(file='../../data/random_l-180_n-100_hg38_minus-test.fa')
    pass
