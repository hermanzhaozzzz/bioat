from __future__ import absolute_import
import gzip
import sys
from bioat.logger import get_logger

__module_name__ = "bioat.fastxtools"


class FastxTools:
    def __init__(self):
        pass

    @staticmethod
    def _load_fastx_generator(fx, log_level='WARNING'):
        """

        :param str file: path of input <fastq | fastq.gz | fastx | fastx.gz>
        :return: a generator for all reads:
            i.e. print(next(obj)) -> ['header', 'seq', 'info', 'quality']
        :rtype: generator
        """
        logger = get_logger(level=log_level, module_name=__module_name__, func_name=sys._getframe().f_code.co_name)
        f = open(fx, 'rt') if not fx.endswith('.gz') else gzip.open(fx, 'rt')
        # FASTQ @  | FASTA >
        symbol = f.read(1)
        f.close()
        f = open(fx, 'rt') if not fx.endswith('.gz') in fx else gzip.open(fx, 'rt')

        if symbol == '@':
            # FASTQ
            logger.debug('is fastq!')
            read = []
            line = f.readline().rstrip()
            # TODO

            while True:
                if not line:
                    break
                else:
                    if line.startswith('@'):
                        n = len(read)
                        if n == 0:
                            # read == []
                            # 第一次循环开始
                            read.append(line)  # add header!
                            line = f.readline().rstrip()
                        elif n == 4:
                            # read == ['header', 'seq', 'info', 'quality']
                            # ls.append(read)
                            yield read
                            read = []
                            read.append(line)
                            line = f.readline().rstrip()
                        elif n == 3:
                            # TODO
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

                            if not line.startswith('@'):
                                raise ValueError('The file may be incomplete!')
                        else:
                            f.close()
                            raise ValueError('The file may be incomplete!')

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
        elif symbol == '>':
            # FASTA
            # print('is fastx!')
            ls = []
            read = []
            seq = ''
            # raw_info = [i.rstrip() for i in f.readlines()]
            line = f.readline().rstrip()
            # print(raw_info)

            while True:
                if not line:
                    break
                else:
                    if line.startswith('>'):
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
                            seq = ''  # 重置 seq 这个 str
                        else:
                            f.close()
                            raise ValueError('The file may be incomplete!')
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
            raise ValueError('Input line one must starts with @ for FQ or > for FA!')
        f.close()

    def _load_fastx(self):
        if self.iterable:
            return self.__load_fastx_generator()
        else:
            return self.__load_fastx_list()

    def get_headers(self):
        # if self.fastx is None:
        self.fastx = self._load_fastx()

        return [i[0][1:] for i in self.fastx]

    def get_sequence(self):
        # if self.fastx is None:
        self.fastx = self._load_fastx()

        return [i[1] for i in self.fastx]

    def get_record_dict(self):
        # if self.fastx is None:
        self.fastx = self._load_fastx()

        return {i[0][1:]: i[1] for i in self.fastx}

    # cli subcmd
    def filter_read_contains_N(self, file: str, output=sys.stdout.name):
        """Filter read contains N base in FASTA or FASTQ



        :param file: FASTA | FASTQ file name.
        :param output: FASTA | FASTQ file name file, stdout if not assigned, format as input.
        """
        if self.fastx is None:
            self.fastx = self._load_fastx()

        self.file = file
        output = sys.stdout if output == sys.stdout.name else output

        if isinstance(output, str):
            f = gzip.open(output, 'wt') if output.endswith('.gz') else open(output, 'wt')
        else:
            f = output

        fx = self._load_fastx()

        for read in fx:
            seq = read[1]

            if seq.upper().__contains__('N'):
                continue
            else:
                f.write('\n'.join(read) + '\n')
        f.close()


if __name__ == '__main__':
    # fastx = FastxTools(file='../../data/random_l-180_n-100_hg38_minus-test.fa')
    pass
