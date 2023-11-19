import sys
import time
import pandas as pd
from bioat import get_logger
from pandarallel import pandarallel
from bioat.lib._dev_tools import profile
import gc

__module_name__ = 'bioat.hictools'

try:
    pandarallel.initialize(
        progress_bar=False,
        use_memory_fs=True,
        verbose=1
    )
except SystemError:
    pandarallel.initialize(
        progress_bar=False,
        use_memory_fs=False,
        verbose=1
    )


class HiCTools:
    """Hi-C toolbox."""

    def __init__(self):
        pass

    @profile
    def get_effective_resolutions(
            self,
            genome_index,
            valid_pairs,
            output,
            log_level: str = 'INFO'
    ):
        """Get deepest resolution of cis-interation.

        输入为result文件，输出为matrix的大致resolution.

        :param genome_index: genome.fa.fai
        :param valid_pairs: rm_dup_pairs.allValidPairs by HiCTools-Pro
        :param output: table_output, TSV file
        :param log_level: 'CRITICAL', 'ERROR', 'WARNING', 'INFO', 'DEBUG', 'NOTSET'
        """
        logger = get_logger(level=log_level, module_name=__module_name__, func_name=sys._getframe().f_code.co_name)
        logger.debug('Develop mode is on')
        # load ref genome lengths
        df_chromosome = pd.read_csv(
            genome_index,
            sep='\t',
            usecols=[0, 1],
            names=['chromosome', 'length'],
            index_col='chromosome',
        )
        # load sample valid pairs by HiCTools-Pro after remove duplication and filter.
        reader = pd.read_csv(
            valid_pairs,
            sep='\t',
            header=None,
            names=[
                'read_name', 'chr_A', 'start_A', 'strand_A', 'chr_B', 'start_B', 'strand_B', 'insert_size',
                'fragment_name_A', 'fragment_name_B', 'mapQ_A', 'mapQ_B', 'allele_specific_info'
            ],
            usecols=['chr_A', 'start_A', 'chr_B', 'start_B'],
            dtype={
                'chr_A': str,
                'start_A': int,
                'chr_B': str,
                'start_B': int,
            },
            iterator=True
        )
        chunks = []

        while True:
            try:
                chunk = reader.get_chunk(10_000_000)
                chunks.append(chunk)
            except StopIteration:
                logger.debug("pandas reader iteration is done.")
                break

        df_valid_pairs = pd.concat(chunks, ignore_index=True)
        logger.debug(df_valid_pairs.info(memory_usage='deep'))
        # 思路：二分法确认new_bin下线
        output = open(output, 'wt')
        write_lines = 0
        output.write('test_bin(kb)\ttotal_bins\tthreshold(%)\tis_effective\n')
        write_lines += 1
        # 思路：二分法确认new_bin下线
        circle = 1
        new_bin = 512_000  # 512kb

        def get_bin_index(x):
            cumsum = df_chromosome['bins_cumsum'][x['chr_A']]
            bin_count_in_this_chrom_a = int(x['start_A'] / new_bin) + 1
            bin_count_in_this_chrom_b = int(x['start_B'] / new_bin) + 1
            bin_index_a = cumsum + bin_count_in_this_chrom_a - 1
            bin_index_b = cumsum + bin_count_in_this_chrom_b - 1
            return bin_index_a, bin_index_b

        while True:
            # 每次迭代都把new_bin / 2
            new_bin = int(new_bin / circle)
            circle = 2

            # calculate total bins
            df_chromosome['bins'] = df_chromosome['length'].map(lambda x: int(x / new_bin) + 1)
            df_chromosome['bins_cumsum'] = df_chromosome['bins'].cumsum()
            total_bins = df_chromosome.iloc[-1, 2]

            # calculate mapped bin index
            df_valid_pairs[['bin_index_A', 'bin_index_B']] = df_valid_pairs.parallel_apply(
                # df_valid_pairs[['bin_index_A', 'bin_index_B']] = df_valid_pairs.apply(
                get_bin_index,
                axis=1,
                result_type='expand'
            )

            # print(df_valid_pairs)
            # threshold calculation
            s_bin_idx_count1 = df_valid_pairs['bin_index_A'].value_counts()
            s_bin_idx_count2 = df_valid_pairs['bin_index_B'].value_counts()
            df_bin_interaction_count = pd.concat([s_bin_idx_count1, s_bin_idx_count2], axis=1)
            df_bin_interaction_count.fillna(0, inplace=True)
            df_bin_interaction_count['total'] = df_bin_interaction_count['bin_index_A'] + df_bin_interaction_count[
                'bin_index_B'
            ]
            threshold = round((df_bin_interaction_count['total'] >= 1000).sum() / total_bins * 100, 1)

            if threshold > 80:
                # 输出结果为有效,进入下一个循环
                output.write(f'{int(new_bin / 1000)}\t{total_bins}\t{threshold}\tYES\n')
                write_lines += 1
            else:
                # if write_lines != 1:
                #     # 第1+n次计算无效,直接终止程序
                #     output.write(f'{int(new_bin / 1000)}\t{total_bins}\t{threshold}\tNO\n')
                # else:
                #     # 第一次计算就无效,直接终止程序
                output.write(f'{int(new_bin / 1000)}\t{total_bins}\t{threshold}\tNO\n')
                break
            del s_bin_idx_count1, s_bin_idx_count2, df_bin_interaction_count
            gc.collect()
        output.close()


if __name__ == '__main__':
    a = time.time()
    hic = HiCTools()
    hic.get_effective_resolutions(
        genome_index="/Volumes/zhaohn_HD/Bio/1.database/db_genomes/genome_fa/genome_ucsc_hg38/genome_ucsc_hg38.fa.fai",
        valid_pairs="../../data/hic/test.rm_dup_pairs.allValidPairs",
        output="/Users/zhaohuanan/Downloads/testout"
    )
    a = time.time() - a
    # print(f'\n{a:.3f}s')
