"""
FASTQ 读取
FASTA 读取
...

load_fastx
load_fastx_generator

"""
import gzip


def load_fastx(file: str) -> list:
    """

    :param str file: path of input <fastq | fastq.gz | fasta | fasta.gz>
    :return: a list for all reads:
        i.e. [['header', 'seq', 'info', 'quality'], [...],[...]]
    :rtype: list
    """
    f = open(file, 'rt') if not '.gz' in file else gzip.open(file, 'rt')
    # FASTQ @  FASTA >
    symbol = f.read(1)
    # print(symbol)
    f.close()
    # 全部载入内存就好，当文件很小的时候
    f = open(file, 'rt') if not '.gz' in file else gzip.open(file, 'rt')

    if symbol == '@':
        # FASTQ
        # print('is fastq')
        ls = []
        read = []
        raw_info = [i.rstrip() for i in f.readlines()]
        line_fix = None  # 为了 fix quality 中开头为@的情况

        # for line in tqdm(raw_info):
        for line in raw_info:
            if line.startswith('@'):
                n = len(read)
                if n == 0:
                    # 第一次循环开始
                    read.append(line)  # add header
                elif n == 4:
                    # 后面的循环开始
                    ls.append(read)
                    read = []
                    read.append(line)  # add header
                elif n == 3:
                    # try to add quality line
                    """Fix a bug.
                    如果第 4 行是@开头,先正常操作一行 read 的 append,再判断下一行是否是@开头,如果不是,文件损坏,如果是,没问题继续
                    ['@Beta12AdemL1C001R00100001768/1', 
                    'ATCCCCGTATCTTCACCCCACCACAAACTATTAGCTTTAGA', 
                    '+']
                    '@@IIEIBCE>IC<IBIIIIEAIEIEB<IDECCD6ICBCED<'
                    """
                    # print(read)
                    # print(line)
                    if line_fix:
                        # 遇到需要fix的下一次循环
                        read.append(line_fix)
                        ls.append(read)
                        read = []
                        read.append(line)
                        line_fix = None
                    else:
                        # 第一次遇到需要 fix 的 line
                        line_fix = line
                else:
                    f.close()
                    raise ValueError('The file may be incomplete!')
            else:
                read.append(line)  # add other three lines: seq, info, quality

        ls.append(read)
        f.close()
        return ls


    elif symbol == '>':
        # FASTA
        # print('is fasta!')
        ls = []
        read = []
        seq = ''
        raw_info = [i.rstrip() for i in f.readlines()]
        # print(raw_info)

        for line in raw_info:
            if line.startswith('>'):
                # 读取 header！
                n = len(read)

                if n == 0:
                    # 第一次循环
                    read.append(line)
                elif n == 1:
                    # 已经有一个 header 了！现在缺 seq
                    read.append(seq)  # add seq line
                    ls.append(read)
                    read = []  # 重置 read 这个 list
                    read.append(line)
                    seq = ''  # 重置 seq 这个 str
                else:
                    f.close()
                    raise ValueError('The file may be incomplete!')
            else:
                # 读取并添加 seq
                seq += line

        read.append(seq)
        ls.append(read)
        return ls
    else:
        raise ValueError('Input line one must starts with @ for FQ or > for FA!')
    # f.close()

# 当文件很大的时候，用生成器函数，产生一个迭代器来进行文件读取
def load_fastx_generator(file):
    """

    :param str file: path of input <fastq | fastq.gz | fasta | fasta.gz>
    :return: a generator for all reads:
        i.e. print(next(obj)) -> ['header', 'seq', 'info', 'quality']
    :rtype: generator
    """
    f = open(file, 'rt') if not '.gz' in file else gzip.open(file, 'rt')
    # FASTQ @  FASTA >
    symbol = f.read(1)
    # print(symbol)
    f.close()
    # 全部载入内存就好，当文件很小的时候
    f = open(file, 'rt') if not '.gz' in file else gzip.open(file, 'rt')

    if symbol == '@':
        # FASTQ
        # print('is fastq!')
        # ls = []
        read = []
        line = f.readline().rstrip()

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
        # print('is fasta!')
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

if __name__ == '__main__':
    pass