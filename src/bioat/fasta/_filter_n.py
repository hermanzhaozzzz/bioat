import os
import gzip
from sys import stdout
from bioat._loader import load_fastx, load_fastx_generator


def main(input, output=stdout.name):
    # fix params
    output = stdout if output == stdout.name else output

    if isinstance(output, str):
        f = open(output, 'wt') if not output.endswith('.gz') else gzip.open(output, 'wt')
    else:
        f = output

    fa = load_fastx_generator(input)


    for read in fa:
        header, seq = read

        if seq.upper().__contains__('N'):
            continue
        else:
            f.write(f'{header}\n{seq}\n')

    f.close()


if __name__ == "__main__":
    main(input='../../../data/random_l-180_n-100_hg38_minus-test.fa', output='/Users/zhaohuanan/Downloads/testout.fa.gz')
