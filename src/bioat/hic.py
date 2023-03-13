from genome import Genome


# 常量定义
# HG19_LEN = [249250621, 243199373, 198022430, 191154276, 180915260, 171115067, 159138663, 146364022, 141213431,
#             135534747,
#             135006516, 133851895, 115169878, 107349540, 102531392, 90354753, 81195210, 78077248, 59128983, 63025520,
#             48129895,
#             51304566, 155270560, 59373566, 42999]
#
# HG38_LEN = [248956422, 242193529, 198295559, 190214555, 181538259, 170805979, 159345973, 145138636, 138394717,
#             133797422,
#             135086622, 133275309, 114364328, 107043718, 101991189, 90338345, 83257441, 80373285, 58617616, 64444167,
#             46709983,
#             50818468, 156040895, 57227415, 42999]


class HiC():
    def __init__(self):
        pass

    def get_deepest_resolution(self, matrix=None, genome=None, index=None, pattern=None, ):
        """Get deepest resolution of cis-interation.

        输入为result文件，输出为matrix的大致resolution.

        :param matrix: 输入需要的result文件
        :param genome:
        :param index: 记录chr的index
        :param pattern: 输入hg19_N_pattern.bed, block_file, 需要屏蔽的region
        :return:
        """

        def _return_loci_count(hic_result, binsize, chr_len, region=None):
            # 生成存储矩阵的list
            if not region:
                loci_count = chr_len // binsize + 1
            else:
                region_min = region[0]
                region_max = region[1]
                loci_count = (region_max - region_min) // binsize + 1

            loci = [0] * loci_count

            # 生成loci
            f_in = open(hic_result, "rt")
            for hic_inter in f_in:
                hic_inter = hic_inter.split(" ")
                position_1 = int(hic_inter[3])
                position_2 = int(hic_inter[6])

                if not region:
                    bin_index_1 = position_1 // binsize
                    bin_index_2 = position_2 // binsize

                    if bin_index_1 == bin_index_2:
                        loci[bin_index_1] += 1
                    else:
                        loci[bin_index_1] += 1
                        loci[bin_index_2] += 1

                else:
                    if (region_min <= position_1 <= region_max) and (region_min <= position_2 <= region_max):
                        bin_index_1 = (position_1 - region_min) // binsize
                        bin_index_2 = (position_2 - region_max) // binsize

                        if bin_index_1 == bin_index_2:
                            loci[bin_index_1] += 1
                        else:
                            loci[bin_index_1] += 1
                            loci[bin_index_2] += 1
            f_in.close()
            return loci

        def _count_block(block_file, chr_index, binsize, region=None):
            """
            返回对应矩阵中N的列数
            :param block_file: hg19_N_pattern.bed
            :param chr_index: 1..23
            :param binsize: int num
            :param region: 记录需要提取的region 默认全部提取，[0,1000000]
            :return: N block的数量
            """
            block_count = 0
            count_state = False

            with open(block_file, "r") as f_in:
                for block_line in f_in:
                    block_line = block_line.split(" ")
                    block_index_1 = int(block_line[1])
                    block_index_2 = int(block_line[2])

                    if block_line[0] == "chr{0}".format(chr_index):
                        if not count_state:
                            count_state = True

                        if not region:  # 没有region参数的情况
                            loci_index_1 = block_index_1 // binsize
                            loci_index_2 = block_index_2 // binsize
                            block_count += loci_index_2 - loci_index_1

                        else:  # 设置region参数的情况
                            if block_index_1 <= region[0] <= block_index_2:
                                # 当region的左边界在block区域中
                                loci_index_1 = region[0] // binsize
                                loci_index_2 = min(region[1], block_index_2) // binsize
                                block_count += loci_index_2 - loci_index_1

                            elif block_index_1 <= region[0] <= block_index_2:
                                # 当region的右边界在block区域中
                                loci_index_1 = block_index_1 // binsize
                                loci_index_2 = region[1] // binsize
                                block_count += loci_index_2 - loci_index_1

                    elif count_state:  # 判断是否读取完毕
                        break

                return block_count

        def _quantile(num_list, quantile, start_index=0):
            """
            提取到
            :param num_list: 输入的原始数字列表
            :param quantile: 需要取得的分位数
            :param start_index: 需要从哪个位置开始计算quantile
            :return: 对应的分位数
            """
            # 按从小到大排序
            sort_list = sorted(num_list, key=lambda x: x)
            # 忽略掉多个之前的数据
            sort_list = sort_list[start_index:]
            quantile_index = int(len(sort_list) * quantile // 1) - 1
            return sort_list[quantile_index]

        # genome = '/Volumes/zhaohn_HD/Bio/1.database/db_genomes/genome_fa/genome_ucsc_hg38/genome_ucsc_hg38.fa'
        # genome = '/Volumes/zhaohn_HD/Bio/1.database/db_genomes/genome_fa/genome_ucsc_hg38/genome_ucsc_hg38.fa2.gz'
        genome = '../../data/test_genome.fa'
        matrix = '../../data/hic/SRR400264_01_1000_raw.matrix'
        chr_lengths = Genome(file=genome).chr_lengths
        binsize_list = [1000, 2000, 3000, 4000, 5000, 8000, 10000, 15000, 20000, 25000, 30000, 50000, 100000,
                        1000000]
        PRINT_STATE = False

        for matrix_binsize in binsize_list:
            # 返回loci
            loci = _return_loci_count(hic_result=input, chr_len=chr_lengths['chr1'], binsize=matrix_binsize)
            # 计算在该binsize的情况下N block对应的bin的数量
            loci_n_count = _count_block(block_file=pattern, chr_index=index, binsize=matrix_binsize)
            # 计算matrix的quantile
            loci_quantile = _quantile(loci, 0.8, start_index=loci_n_count)

            if loci_quantile >= 1000:
                print(f"chr{index} {matrix_binsize} {loci_quantile}")
                PRINT_STATE = True
                break

        if not PRINT_STATE:
            print(f"chr{index} None None")

    # elif args.command == "ResolutionAll":
    #     # 计算cis矩阵的分辨率，方便计算所有的matrix的resolution
    #     # 输入为result文件所在的dir，输出为matrix的大致resolution
    #     # -i input 参数输入需要的result文件的目录
    #     if not FilePath(args.input).check_is_dir():
    #         raise IOError("Your -i --input is not right!")
    #     # -o output 参数输出resolution的信息
    #     try:
    #         f_out = open(args.output, "w")
    #         f_out.close()
    #     except:
    #         raise IOError("Your -o --output is not right!")
    #     # --pattern file 输入hg19_N_pattern.bed
    #     if not FilePath(args.pattern).check_path():
    #         raise IOError("Your --pattern is not right!")
    #     # --MAPQ
    #     if args.MAPQ is None:
    #         raise IOError("We need --MAPQ argue!")
    #     else:
    #         args.MAPQ = int(args.MAPQ)
    #
    #     # 定义常用的变量
    #     binsize_list = [3000, 4000, 5000, 8000, 10000, 15000, 20000, 25000, 30000, 50000, 100000, 1000000]
    #     result_file_str = "chr_{0}_MAPQ{1}_result.txt"
    #
    #     # 为了节约内存，需要通过循环进行求resolution，每次1个chr
    #     for chr_index in range(1, 25):
    #         # 打开需要写入的文件，使用mode a形式打开
    #         f_out = open(args.output, "a")
    #         # 记录是否有符合要求的resolution
    #         print_state = False
    #         # 生成result文件path
    #         result_file = args.input + result_file_str.format(chr_index, args.MAPQ)
    #         # 计算相应染色体的长度
    #         chr_len = Chromosome(chr_index, version="hg19").get_length()
    #
    #         for matrix_binsize in binsize_list:
    #             # 返回loci
    #             loci = _return_loci_count(hic_result=result_file, chr_len=chr_len, binsize=matrix_binsize)
    #             # 计算在该binsize的情况下N block对应的bin的数量
    #             loci_n_count = _count_block(block_file=args.pattern, chr_index=chr_index, binsize=matrix_binsize)
    #             # 计算matrix的quantile
    #             loci_quantile = _quantile(loci, 0.8, start_index=loci_n_count)
    #
    #             if loci_quantile >= 1000:
    #                 f_out.write("chr{0} {1} {2}".format(chr_index, matrix_binsize, loci_quantile) + "\n")
    #                 print_state = True
    #                 break
    #
    #             # 释放空间
    #             loci = None
    #
    #         # 如果没有合适的resolution则输出None
    #         if not print_state:
    #             f_out.write("chr{0} {1} {2}".format(chr_index, "None", "None") + "\n")
    #
    #         f_out.close()  # 关闭文件
    #


# class Genome(object):
#     """
#     定义Genome类，可以直接调用长度酶切位点等信息
#     """
#
#     def __init__(self, file, version="hg19"):
#         self.version = version
#         self.file = file
#         self.genome_list = []
#         self.genome_seq = []
#
#     # def get_genome(self):
#     #     self.genome_list = FastaFile(self.file).get_sequence()
#     #     self.genome_seq = self.genome_list[1]
#     #     return self.genome_seq


# class Matrix(object):
#     """
#     输入为二维列表，当做矩阵来处理
#     """
#
#     def __init__(self, input_matrix, **keys):
#         self.matrix = input_matrix  # matrix itself
#         self.nrow = len(input_matrix)
#         self.ncol = len(input_matrix[0])
#         self.dim = [len(input_matrix), len(input_matrix[0])]
#
#     def t_matrix(self):
#         # return transposition matrix
#         return zip(*self.matrix)
#
#     def col_sum(self):
#         # return column sum of raw matrix
#         col_sum = []
#         trans_matrix = zip(*self.matrix)
#         for col in trans_matrix:
#             col_sum.append(sum(col))
#         return col_sum
#
#     def row_sum(self):
#         # return row sum of raw matrix
#         row_sum = []
#         for row in self.matrix:
#             row_sum.append(sum(row))
#         return row_sum


# class HicResultFile(object):
#     """
#     定义result文件类型，从而快速生成matrix文件
#     """
#
#     def __init__(self, result_file, version):
#         self.path = result_file
#         self.version = version
#
#         # 定义初始值
#         self.ncol = 1
#         self.nrow = 1
#         self.matrix = []
#         self.matrix_type = ""
#
#     def get_dim(self, binsize):
#         """
#         计算result文件中对应的矩阵应该的row与col大小
#         :param binsize: 需要生成的binsize大小
#         :return:返回hic矩阵的行，列
#         """
#         hic_result_file = open(self.path, "r")
#         hic_inter = hic_result_file.readline().strip().split(" ")
#         hic_result_file.close()
#         hic_inter_obj = HicInteraction(hic_inter)
#
#         if hic_inter_obj.chr_index_1 == hic_inter_obj.chr_index_2:
#             self.matrix_type = "cis"
#             chr_len = Chromosome(hic_inter_obj.chr_index_1,
#                                  version=self.version
#                                  ).get_length()
#
#             self.ncol = chr_len // binsize + 1
#             self.nrow = self.ncol
#
#         elif hic_inter_obj.chr_index_1 != hic_inter_obj.chr_index_2:
#             self.matrix_type = "trans"
#             chr_index_min = min(hic_inter_obj.chr_index_1,
#                                 hic_inter_obj.chr_index_2)
#
#             chr_index_max = max(hic_inter_obj.chr_index_1,
#                                 hic_inter_obj.chr_index_2)
#
#             chr_len_1 = Chromosome(chr_index_min,
#                                    version=self.version
#                                    ).get_length()
#
#             chr_len_2 = Chromosome(chr_index_max,
#                                    version=self.version
#                                    ).get_length()
#
#             self.nrow = chr_len_1 // binsize + 1
#             self.ncol = chr_len_2 // binsize + 1
#
#         return [self.nrow, self.ncol]
#
#     def get_matrix(self, binsize):
#         """
#         生成hic矩阵，需要定义binsize
#         :return:返回的为1个二维矩阵
#         """
#         # 通过get_dim的函数来获得行列
#         # 当cis的时候，返回值相等；当trans的时候nrow > ncol
#         self.nrow, self.ncol = self.get_dim(binsize)
#
#         # 生成2维列表
#         matrix = [[0] * self.ncol for index in range(self.nrow)]
#
#         # 打开result 文件
#         with open(self.path, "r") as hic_result_file:
#
#             for hic_inter in hic_result_file:
#                 hic_inter = hic_inter.strip().split(" ")
#                 hic_inter_obj = HicInteraction(hic_inter)
#                 matrix_index_1 = hic_inter_obj.chr_position_1 // binsize
#                 matrix_index_2 = hic_inter_obj.chr_position_2 // binsize
#
#                 if self.nrow == self.ncol:  # cis
#                     if matrix_index_1 == matrix_index_2:
#                         matrix[matrix_index_1][matrix_index_2] += 1
#                     else:
#                         matrix[matrix_index_1][matrix_index_2] += 1
#                         matrix[matrix_index_2][matrix_index_1] += 1
#
#                 elif self.nrow != self.ncol:  # trans
#                     if hic_inter_obj.chr_index_1 > hic_inter_obj.chr_index_2:
#                         matrix[matrix_index_2][matrix_index_1] += 1
#                     elif hic_inter_obj.chr_index_1 < hic_inter_obj.chr_index_2:
#                         matrix[matrix_index_1][matrix_index_2] += 1
#
#         return matrix
#
#     def get_part_matrix(self, binsize, region=None):
#         """
#         生成hic矩阵，需要定义binsize
#         return: 1个二维矩阵
#         """
#         # region 参数
#         if region is None:
#             raise IOError("No region input!")
#         # 通过region生成matrix的row和col
#         region_min = region[0]
#         region_max = region[1]
#         matrix_nrow = (region_max - region_min) // binsize
#         matrix_ncol = matrix_nrow
#         # 生成2维列表
#         matrix = [[0] * matrix_ncol for index in range(matrix_nrow)]
#
#         # 打开result 文件
#         hic_result_file = open(self.path, "r")
#
#         for hic_inter in hic_result_file:
#             hic_inter = hic_inter.strip().split(" ")
#             hic_inter_obj = HicInteraction(hic_inter)
#             if hic_inter_obj.in_region([region_min, region_max]):
#                 matrix_index_1 = (hic_inter_obj.chr_position_1 - region_min) // binsize
#                 matrix_index_2 = (hic_inter_obj.chr_position_2 - region_min) // binsize
#
#                 if matrix_index_1 == matrix_index_2:
#                     matrix[matrix_index_1][matrix_index_2] += 1
#                 else:
#                     matrix[matrix_index_1][matrix_index_2] += 1
#                     matrix[matrix_index_2][matrix_index_1] += 1
#
#         hic_result_file.close()
#
#         return matrix
#
#     def get_rDNA_matrix(self, chr_index, chr_binsize, rDNA_binsize):
#         """
#         生成hic矩阵，需要定义binsize
#         :return:返回的为1个二维矩阵
#         """
#         # 通过get_dim的函数来获得行列
#         # 当cis的时候，返回值相等；当trans的时候nrow > ncol
#         if chr_index == 25:
#             matrix_ncol = matrix_nrow = Chromosome(chr_index=25, version="hg19").get_length() // rDNA_binsize + 1
#         else:
#             matrix_nrow = Chromosome(chr_index, version="hg19").get_length() // chr_binsize + 1
#             matrix_ncol = Chromosome(chr_index=25, version="hg19").get_length() // rDNA_binsize + 1
#
#         # 生成2维列表
#         matrix = [[0] * matrix_ncol for index in range(matrix_nrow)]
#
#         # 打开result 文件
#         if chr_index == 25:
#             with open(self.path, "r") as hic_result_file:
#                 for hic_inter in hic_result_file:
#                     hic_inter = hic_inter.strip().split(" ")
#                     chr_index_1 = int(hic_inter[2])
#                     chr_index_2 = int(hic_inter[5])
#                     matrix_index_1 = chr_index_1 // rDNA_binsize
#                     matrix_index_2 = chr_index_2 // rDNA_binsize
#
#                     if matrix_index_1 == matrix_index_2:
#                         matrix[matrix_index_1][matrix_index_2] += 1
#                     else:
#                         matrix[matrix_index_1][matrix_index_2] += 1
#                         matrix[matrix_index_2][matrix_index_1] += 1
#
#         else:
#             with open(self.path, "r") as hic_result_file:
#                 for hic_inter in hic_result_file:
#                     hic_inter = hic_inter.strip().split(" ")
#                     chr_index_1 = int(hic_inter[2])
#                     chr_index_2 = int(hic_inter[5])
#                     if chr_index_1 == 25:
#                         matrix_index_col = chr_index_1 // rDNA_binsize
#                         matrix_index_row = chr_index_2 // chr_binsize
#                         matrix[matrix_index_row][matrix_index_col] += 1
#
#                     elif chr_index_2 == 25:
#                         matrix_index_col = chr_index_2 // rDNA_binsize
#                         matrix_index_row = chr_index_1 // chr_binsize
#                         matrix[matrix_index_row][matrix_index_col] += 1
#
#         return matrix


# class HicInteraction(object):
# """
# 方便处理Hic interaction的每一行
# """
#
# def __init__(self, interaction, read_len=36):
#
#     self.interaction = interaction
#     self.read_length = read_len
#
#     # 输入的为1个interaction列表，总共7个元素分别进行如下赋值
#     self.header = interaction[0]
#
#     self.chr_strand_1 = int(interaction[1])
#     self.chr_index_1 = int(interaction[2])
#     self.chr_position_1 = int(interaction[3])
#
#     self.chr_strand_2 = int(interaction[4])
#     self.chr_index_2 = int(interaction[5])
#     self.chr_position_2 = int(interaction[6])
#
#     # 初始化数据
#     self.loci_index_1 = 0
#     self.loci_index_2 = 0
#     self.fragment_dist = 0
#     self.type = "0"
#     self.enzyme = ""
#
# def show(self):
#     return " ".join(self.interaction)
#
# def in_region(self, region):
#     region_min = region[0]
#     region_max = region[1]
#     if (region_min <= self.chr_position_1, self.chr_position_2 <= region_max) == (True, True):
#         return True
#     else:
#         return False
#
# def loci_index(self, binsize):
#     # 根据self.chr_index和binsize计算bin的index 返回对应的lociIndex
#     self.loci_index_1 = self.chr_position_1 // binsize
#     self.loci_index_2 = self.chr_position_2 // binsize
#     return [self.loci_index_1, self.loci_index_2]
#
# def get_type(self):
#     # 根据reads的上下游关系以及正负链关系，判断为inward，outward，samestrand
#
#     if self.chr_strand_1 == self.chr_strand_2:
#         self.type = "3"  # samestrand
#
#     elif self.chr_position_1 < self.chr_position_2:
#         if self.chr_strand_1 == "0":
#             self.type = "1"  # inward
#         else:
#             self.type = "2"  # outward
#
#     elif self.chr_position_1 > self.chr_position_2:
#         if self.chr_position_1 == "0":
#             self.type = "2"  # outward
#         else:
#             self.type = "1"  # inward
#     return self.type
#
# def get_frag_distance(self, enzyme, chr_seq=""):
#     # 初始化赋值
#     self.enzyme = enzyme
#
#     # 计算两条reads的所在fragment之间的距离
#     chr_position_min = min(self.chr_position_1, self.chr_position_2)
#     chr_position_max = max(self.chr_position_1, self.chr_position_2)
#
#     add_length = 10000
#
#     fragment_1 = chr_seq[chr_position_min:chr_position_min + add_length]
#     fragment_2 = chr_seq[chr_position_max - add_length:chr_position_max][::-1]
#
#     fragment_1_border = fragment_1.find(self.enzyme)
#     fragment_2_border = fragment_2.find(self.enzyme[::-1])
#
#     if fragment_1_border == -1 or fragment_2_border == -1:
#         self.fragment_dist = 0
#     else:
#         self.fragment_dist = chr_position_max \
#                              - chr_position_min \
#                              - fragment_1_border \
#                              - fragment_2_border \
#                              - self.read_length
#
#         if self.fragment_dist < 0:  # 在同一个fragment上
#             self.fragment_dist = 0
#
#     return self.fragment_dist


# class HicMapRead(object):
#     """
#     处理 reads_name strand_info chr_index chr_position
#     可以使用 if_DownStreamRE() 来判断是否有酶切位点
#     """
#
#     def __init__(self, hic_read, read_len=36):
#
#         if len(hic_read) < 4:
#             raise IOError("Your input hic_map_read is illegal.")
#
#         self.header = hic_read[0]
#
#         # 储存FLAG信息，可以修正
#         self.chr_strand = hic_read[1]
#         # if self.chr_strand == "0":
#         #     self.chr_strand = "0"  # 正链
#         #
#         # elif self.chr_strand == "16" or self.chr_strand == "1" :
#         #     self.chr_strand = "1"  # 负链
#         #
#         # else:
#         #     self.chr_strand = "2"  # 其他情况
#
#         # 储存chr_index信息
#         self.chr_index = int(hic_read[2])
#         # if hic_read[2] == "rDNA":
#         #     self.chr_index = 25
#         # elif hic_read[2] == "chrY":
#         #     self.chr_index = 24
#         # elif hic_read[2] == "chrX":
#         #     self.chr_index = 23
#         # elif len(hic_read[2]) <= 5:
#         #     self.chr_index = int(hic_read[2].replace("chr",""))
#
#         self.chr_position = int(hic_read[3])
#
#         self.read_len = 36  # read_len 是map时候选择trim的reads的长度我们用36bp
#         self.fragment = ""
#         self.if_enzyme = False
#
#     def downstream_enzyme(self, chr_seq, enzyme, distance):
#         # 用于计算hic read在一定distance内是否有酶切位点
#         # 如果没有返回True，有返回False
#
#         if self.chr_strand == "0":
#             # 如果map到了chromosome的正链上
#             self.fragment = chr_seq[self.chr_position:self.chr_position + distance + self.read_len]
#
#         elif self.chr_strand == "1":
#             # 如果map到了chromosome的负链上
#             self.fragment = chr_seq[self.chr_position - distance:self.chr_position]
#
#         if self.fragment.find(enzyme) != -1:
#             self.if_enzyme = True
#         else:
#             self.if_enzyme = False
#
#         return self.if_enzyme
#
#     def show(self):
#         show_str = " ".join(map(str, [self.header, self.chr_strand, self.chr_index, self.chr_position]))
#         return show_str

#
# # 常用函数
# def change_flag(hic_read):
#     if len(hic_read) < 4:
#         raise IOError("Your input hic_map_read is illegal.")
#
#     # 储存FLAG信息，可以修正
#     chr_strand = hic_read[1]
#     if chr_strand == "0":
#         chr_strand = "0"  # 正链
#     elif chr_strand == "16" or chr_strand == "1":
#         chr_strand = "1"  # 负链
#     else:
#         chr_strand = "2"  # 其他情况
#
#     # 储存chr_index信息
#     chr_index = ""
#     if hic_read[2] == "rDNA":
#         chr_index = "25"
#     elif hic_read[2] == "chrY":
#         chr_index = "24"
#     elif hic_read[2] == "chrX":
#         chr_index = "23"
#     elif len(hic_read[2]) <= 5:
#         try:
#             chr_index = hic_read[2][3:]
#         except:
#             chr_index = "0"
#
#     # 返回新的list
#     return [hic_read[0], chr_strand, chr_index, hic_read[3]]


# def downstream_enzyme(hic_read, chr_seq, enzyme, distance, read_len=36):
#     # 用于计算hic read在一定distance内是否有酶切位点
#     # 如果没有返回True，有返回False
#     fragment = ""
#     chr_strand = hic_read[1]
#     chr_position = int(hic_read[3])
#
#     if chr_strand == "0":
#         # 如果map到了chromosome的正链上
#         fragment = chr_seq[chr_position: chr_position + distance + read_len]
#     elif chr_strand == "1":
#         # 如果map到了chromosome的负链上
#         fragment = chr_seq[chr_position - distance: chr_position]
#
#     if fragment.find(enzyme) != -1:
#         if_enzyme = True
#     else:
#         if_enzyme = False
#
#     return if_enzyme
#
#
# def cis_add_info(hic_inter, chr_seq, enzyme, read_len=36):
#     # 根据reads的上下游关系以及正负链关系，判断为inward，outward，samestrand
#     inter_type = ""
#
#     chr_strand_1 = hic_inter[1]
#     chr_strand_2 = hic_inter[4]
#     chr_position_1 = int(hic_inter[3])
#     chr_position_2 = int(hic_inter[6])
#
#     if chr_strand_1 == chr_strand_2:
#         inter_type = "3"  # samestrand
#
#     elif chr_position_1 < chr_position_2:
#         if chr_strand_1 == "0":
#             inter_type = "1"  # inward
#         else:
#             inter_type = "2"  # outward
#
#     elif chr_position_1 > chr_position_2:
#         if chr_strand_1 == "0":
#             inter_type = "2"  # outward
#         else:
#             inter_type = "1"  # inward
#
#     # 计算两条reads的所在fragment之间的距离
#     chr_position_min = min(chr_position_1, chr_position_2)
#     chr_position_max = max(chr_position_1, chr_position_2)
#
#     add_length = 20000
#
#     fragment_1 = chr_seq[chr_position_min:chr_position_min + add_length]
#     fragment_2 = chr_seq[chr_position_max - add_length:chr_position_max][::-1]
#
#     fragment_1_border = fragment_1.find(enzyme)
#     fragment_2_border = fragment_2.find(enzyme[::-1])
#
#     if fragment_1_border == -1 or fragment_2_border == -1:
#         fragment_dist = 0
#     else:
#         fragment_dist = chr_position_max - chr_position_min \
#                         - fragment_1_border - fragment_2_border - read_len
#
#         if fragment_dist < 0:
#             fragment_dist = 0  # 在同一个fragment上
#
#     return [inter_type, fragment_dist]


# class Chromosome(object):
#     '''Define chromosome obj'''
#
#     def __init__(self, chr_index=1, version="hg19"):
#         self.index = chr_index - 1
#         self.version = version
#         self.length = 0
#
#     def get_length(self):
#         # 返回染色体的长度
#         if self.version == "hg19":
#             self.length = HG19_LEN[self.index]
#
#         elif self.version == "hg38":
#             self.length = HG38_LEN[self.index]
#
#         return self.length
#
#     def get_loci(self, bin_size=1000000):
#         # 返回1个列表，分别是start，end，正链，负链reads数目
#         self.get_length()
#         loci_count = self.length // bin_size + 1
#         chr_loci = []
#         for index in range(loci_count):
#             chr_loci.append([index * bin_size, (index + 1) * bin_size, 0])
#         return chr_loci
#
#     def get_hic_matrix(self, bin_size=1000000):
#         # 返回HicMatrix
#         self.get_length()
#         loci_count = self.length // bin_size + 1
#         chr_matrix = []
#         for i in range(loci_count):
#             chr_matrix.append([0] * loci_count)
#         return chr_matrix
#

# class BedRead(object):
#     '''Define the bed read object'''
#
#     def __init__(self, bed_read):
#         self.bed_list = bed_read
#         self.start = int(bed_read[1])
#         self.end = int(bed_read[2])
#
#         if bed_read[0] == "rDNA":
#             self.chr_index = 25
#         elif bed_read[0] == "chrY":
#             self.chr_index = 24
#         elif bed_read[0] == "chrX":
#             self.chr_index = 23
#         else:
#             self.chr_index = int(bed_read[0][3:])
#
#         self.bed_list[0] = self.chr_index
#
#         if len(bed_read) >= 4:
#             self.name = bed_read[3]
#             self.bed_list[3] = self.name
#
#         if len(bed_read) >= 5:
#             self.score = int(bed_read[4])
#             self.bed_list[4] = self.score
#
#         if len(bed_read) >= 6:
#             self.strand = bed_read[5]
#             self.bed_list[5] = self.strand
#
#     def show(self):
#         write_str = "\t".join(map(str, self.bed_list))
#         return write_str
#
#     def get_loci_index(self, bin_size=1000000):
#         loci_index = self.start // bin_size
#         return loci_index
#
#
# def _return_cis_matrix(hic_result, chr_len, binsize, region=None):
#     """
#     生成cis matrix
#     :param hic_result:hic matrix的result文件，chr1_MAPQ1_result.txt
#     :param chr_index: 染色体编号，1
#     :param binsize: 需要生成矩阵的resolution大小 10000
#     :param chr_len: 对应染色体的长度
#     :param region: 记录需要提取的region 默认全部提取，[0,1000000]
#     :return: 返回2维列表
#     """
#
#     # 生成存储矩阵的list
#     if not region:
#         bin_count = chr_len // binsize + 1
#     else:
#         region_min = region[0]
#         region_max = region[1]
#         bin_count = (region_max - region_min) // binsize + 1
#
#     matrix = [([0] * bin_count) for i in range(bin_count)]
#
#     # 生成矩阵
#     with open(hic_result, "r") as f_in:
#         for hic_inter in f_in:
#             hic_inter = hic_inter.split(" ")
#             position_1 = int(hic_inter[3])
#             position_2 = int(hic_inter[6])
#
#             if not region:
#                 bin_index_1 = position_1 // binsize
#                 bin_index_2 = position_2 // binsize
#
#                 if bin_index_1 == bin_index_2:
#                     matrix[bin_index_1][bin_index_2] += 1
#                 else:
#                     matrix[bin_index_1][bin_index_2] += 1
#                     matrix[bin_index_2][bin_index_1] += 1
#
#             else:
#                 if (region_min <= position_1 <= region_max) and (region_min <= position_2 <= region_max):
#                     bin_index_1 = (position_1 - region_min) // binsize
#                     bin_index_2 = (position_2 - region_max) // binsize
#
#                     if bin_index_1 == bin_index_2:
#                         matrix[bin_index_1][bin_index_2] += 1
#                     else:
#                         matrix[bin_index_1][bin_index_2] += 1
#                         matrix[bin_index_2][bin_index_1] += 1
#
#         return matrix
#
#

#
#

#
#
# def _write_matrix(matrix_file, matrix, sep=","):
#     """
#     :param matrix_file: 输出文件名
#     :param matrix: 二维列表
#     :return: 成功写入返回True 否则返回False
#     """
#     with open(matrix_file, "w") as f_matrix:
#         for row in matrix:
#             f_matrix.write(sep.join(map(str, row)) + "\n")
#         return True
# class FastaFile(object):
#     '''FASTA object'''
#
#     def __init__(self, fasta_file_path):
#         self.path = fasta_file_path
#
#     def get_length(self):  # get fasta file length
#         fasta_file = open(self.path, "r")
#         fasta_header_list = []  # 用来储存fasta所有的head信息
#         fasta_seq_len_list = []  # 用来储存fasta所有的seq信息
#         fasta_seq_list = []
#
#         fasta_header = fasta_file.readline().strip()
#         for fasta_line in fasta_file:
#             if fasta_line[0] != ">":
#                 # 将每一行的长度存入到列表中
#                 fasta_seq_list.append(len(fasta_line.strip()))
#             else:
#                 fasta_header_list.append(fasta_header)
#                 # 对列表进行求和，计算总体的长度
#                 fasta_seq_len_list.append(sum(fasta_seq_list))
#                 fasta_seq_list = []
#                 fasta_header = fasta_line.strip()
#
#         fasta_header_list.append(fasta_header)
#         fasta_seq_len_list.append(len(fasta_seq_list))
#         fasta_file.close()
#
#         return [fasta_header_list, fasta_seq_len_list]
#
#     def get_sequence(self, upper=True):  # get fasta seq
#
#         fasta_file = open(self.path, "r")
#         fasta_header_list = []
#         fasta_seq_combine = []
#         fasta_seq_list = []
#         fasta_header = fasta_file.readline().strip()
#
#         for fasta_line in fasta_file:
#             if fasta_line[0] != ">":
#                 if upper == True:
#                     fasta_seq_list.append(fasta_line.strip().upper())
#                 else:
#                     fasta_seq_list.append(fasta_line.strip())
#             else:
#                 fasta_header_list.append(fasta_header)
#                 fasta_seq_combine.append("".join(fasta_seq_list))
#                 fasta_seq_list = []
#                 fasta_header = fasta_line.strip()
#
#         fasta_header_list.append(fasta_header)
#         fasta_seq_combine.append("".join(fasta_seq_list))
#         fasta_file.close()
#
#         return [fasta_header_list, fasta_seq_combine]
#
#
# class FastqRead(object):
#     '''Define fastq read object'''
#
#     def __init__(self, fastqRead, phred=33):
#         self.head = fastqRead[0]
#         self.sequence = fastqRead[1]
#         self.info = fastqRead[2]
#         self.quality = fastqRead[3]
#         self.phred = phred
#         self.length = len(self.sequence)
#
#     def aver_quality(self):
#         # 计算序列的平均质量
#         sum_quality = 0
#         for char in self.sequence:
#             score = ord(char) - self.phred
#             sum_quality += score
#         aver_quality = sum_quality * 1.0 / self.length
#         return aver_quality
#
#     def trimmer(self, start=0, end=None):
#         # 返回trim的fastq read对象
#         if (start < 0) or (start >= self.length): start = 0
#         if (end >= self.length) or (end <= start): end = None
#
#         trim_sequence = self.sequence[start:end]
#         trim_quality = self.quality[start:end]
#         trim_info = self.info + " Trimmed {0},{1}".format(start, end)
#         trim_read = [self.head, trim_sequence, trim_info, trim_quality]
#         return FastqRead(trim_read)
#
#     def write_format(self):
#         # 返回一个str可以用于直接写入文件
#         write_str = "{0}\n{1}\n{2}\n{3}".format(self.head, self.sequence, self.info, self.quality)
#         return write_str


if __name__ == '__main__':
    # fastx = Fastx(file='../../data/test_genome.fa')
    hic = HiC()
    hic.get_deepest_resolution()
    # parser = argparse.ArgumentParser(description="This tool is for creating Hic matrix. <meng_howard@yahoo.com>")
    #
    # parser.add_argument("command",
    #                     help="command list: <MakeMatrix> <MakeAllMatrix> <Resolution> <MakeResult> <MakePartMatrix>")
    # parser.add_argument("-i", "--input",
    #                     help="input filtered Cis file OR result file", required=True)
    # parser.add_argument("-o", "--output",
    #                     help="output matrix file OR result file")
    # parser.add_argument("-s", "--binsize",
    #                     help="matrix binsize")
    #
    # parser.add_argument("--index",
    #                     help="chromosome index, like 1")
    # parser.add_argument("--region",
    #                     help="output part of matrix, like [0,10000]")
    # parser.add_argument("--transfile",
    #                     help="input transfile when you use <MakeResult>")
    # parser.add_argument("--pattern",
    #                     help="input N pattern file when you use <Resolution>")
    # parser.add_argument("--MAPQ",
    #                     help="input MAPQ information")
    #
    # args = parser.parse_args()
    #
    # if args.command == "MakeMatrix":
    #     # 输入为result文件，输出为matrix
    #     # -i input 参数 / -o output参数
    #     # input为result文件，output为matrix文件
    #     try:
    #         f_result = open(args.input, "r")
    #         f_out = open(args.output, "w")
    #     except:
    #         raise IOError("Your -i OR -o is not right!")
    #     # -s binsize 参数
    #     if args.binsize is None:
    #         raise IOError("We need -s --binsize argue!")
    #     else:
    #         args.binsize = int(args.binsize)
    #     # --region
    #     if args.region is not None:
    #         try:
    #             args.region = eval(args.region)
    #         except:
    #             raise IOError("We need --region with a right format!")
    #
    #     hic_result_obj = HicResultFile(args.input, version="hg19")
    #
    #     # 获得hic matrix
    #     if args.region is not None:
    #         hic_matrix = hic_result_obj.get_part_matrix(binsize=args.binsize, region=args.region)
    #     else:
    #         hic_matrix = hic_result_obj.get_matrix(binsize=args.binsize)
    #
    #     # 输出matrix
    #     for matrix_row in hic_matrix:
    #         try:
    #             f_out.write(",".join(map(str, matrix_row)) + "\n")
    #         except:
    #             raise IOError("Error in {0}".format(args.output))
    #
    # elif args.command == "MakeAllMatrix":
    #     # 输入为result文件目录，输出为全部matrix
    #     # -i input 参数
    #     if not FilePath(args.input).check_is_dir():
    #         raise IOError("Your -i --input is not right!")
    #     # -o output 参数
    #     if not FilePath(args.output).check_is_dir():
    #         raise IOError("Your -o --output is not right!")
    #     # -s binsize 参数
    #     if args.binsize is None:
    #         raise IOError("We need -s --binsize argue!")
    #     else:
    #         args.binsize = int(args.binsize)
    #     # MAPQ参数
    #     if args.MAPQ is None:
    #         raise IOError("We need --MAPQ argue!")
    #     else:
    #         args.MAPQ = int(args.MAPQ)
    #
    #     result_file_cis = args.input + "chr_{0}_MAPQ{1}_result.txt"
    #     result_file_trans = args.input + "chr_{0}_{1}_MAPQ{2}_result.txt"
    #     matrix_file_cis = args.output + "chr_{0}_{1}_MAPQ{2}.txt"
    #     matrix_file_trans = args.output + "chr_{0}_{1}_{2}_MAPQ{3}.txt"
    #
    #     # 生成并写入矩阵
    #     for chr_index_1 in range(1, 25):
    #         for chr_index_2 in range(chr_index_1, 25):
    #             hic_matrix = None
    #             if chr_index_1 == chr_index_2:
    #                 result_file = result_file_cis.format(chr_index_1, args.MAPQ)
    #                 matrix_file = matrix_file_cis.format(chr_index_1, args.binsize, args.MAPQ)
    #
    #                 hic_result_obj = HicResultFile(result_file, version="hg19")
    #                 hic_matrix = hic_result_obj.get_matrix(binsize=args.binsize)
    #                 _write_matrix(matrix_file, hic_matrix, sep=",")
    #
    #             elif chr_index_1 < chr_index_2:
    #                 result_file = result_file_trans.format(chr_index_1, chr_index_2, args.MAPQ)
    #                 matrix_file = matrix_file_trans.format(chr_index_1, chr_index_2, args.binsize, args.MAPQ)
    #
    #                 hic_result_obj = HicResultFile(result_file, version="hg19")
    #                 hic_matrix = hic_result_obj.get_matrix(binsize=args.binsize)
    #                 _write_matrix(matrix_file, hic_matrix, sep=",")
    #
    #     # 生成rDNA矩阵
    #     for chr_index_1 in range(1, 26):
    #         hic_matrix = None
    #
    #         if chr_index_1 == 25:
    #             result_file = result_file_cis.format(chr_index_1, args.MAPQ)
    #             matrix_file = matrix_file_cis.format(chr_index_1, 1000, args.MAPQ)
    #
    #             hic_result_obj = HicResultFile(result_file, version="hg19")
    #             hic_matrix = hic_result_obj.get_rDNA_matrix(chr_index=25, chr_binsize=args.binsize, rDNA_binsize=1000)
    #             _write_matrix(matrix_file, hic_matrix, sep=",")
    #
    #         elif chr_index_1 < 25:
    #             result_file = result_file_trans.format(chr_index_1, 25, args.MAPQ)
    #             matrix_file = matrix_file_trans.format(chr_index_1, 25, args.binsize, args.MAPQ)
    #
    #             hic_result_obj = HicResultFile(result_file, version="hg19")
    #             hic_matrix = hic_result_obj.get_rDNA_matrix(chr_index=chr_index_1, chr_binsize=args.binsize,
    #                                                         rDNA_binsize=1000)
    #             _write_matrix(matrix_file, hic_matrix, sep=",")
    #
    # elif args.command == "LociCount":
    #     # result 文件的目录
    #     if not FilePath(args.input).check_is_dir():
    #         raise IOError("Your -i --input is not right!")
    #     # -o output 参数输出resolution的信息
    #     if not FilePath(args.output).check_is_dir():
    #         raise IOError("Your -o --output is not right!")
    #     # --MAPQ
    #     if args.MAPQ is None:
    #         raise IOError("We need --MAPQ argue!")
    #     else:
    #         args.MAPQ = int(args.MAPQ)
    #     # binsize 设置
    #     if not args.binsize:
    #         raise IOError("We need -s --binsize argue!")
    #     else:
    #         args.binsize = int(args.binsize)
    #
    #     result_file_str = "chr_{0}_MAPQ{1}_result.txt"
    #     loci_file_str = "chr_{0}_{1}_MAPQ{2}_loci.txt"
    #     loci_str = "chr{0}\t{1}\t{2}\t{3}"
    #
    #     for chr_index in range(1, 25):
    #         result_file = args.input + result_file_str.format(chr_index, args.MAPQ)
    #         loci_file = args.output + loci_file_str.format(chr_index, args.binsize, args.MAPQ)
    #         chr_len = Chromosome(chr_index, version="hg19").get_length()
    #
    #         loci = _return_loci_count(hic_result=result_file, chr_len=chr_len, binsize=args.binsize)
    #
    #         with open(loci_file, "w") as f_out:
    #             for loci_index, loci_count in enumerate(loci):
    #                 loci_min = loci_index * args.binsize
    #                 loci_max = loci_min + args.binsize
    #                 loci_write = loci_str.format(chr_index, loci_min, loci_max, loci_count) + "\n"
    #                 f_out.write(loci_write)
    #
    # elif args.command == "MakeResult":
    #     # -i input 参数
    #     if not FilePath(args.input).check_path():
    #         raise IOError("Your -i --input is not right!")
    #     # --transfile 参数
    #     if not FilePath(args.transfile).check_path():
    #         raise IOError("Your --transfile is not right!")
    #     # -o output 参数
    #     if not FilePath(args.output).check_is_dir():
    #         raise IOError("Your -o --output is not right!")
    #     # --index
    #     if args.index is None:
    #         raise IOError("We need --index argue!")
    #     else:
    #         args.index = int(args.index)
    #     # --MAPQ
    #     if args.MAPQ is None:
    #         raise IOError("We need --MAPQ argue!")
    #     else:
    #         args.MAPQ = int(args.MAPQ)
    #
    #     chr_cis_result = []
    #     # 生成长度为25的列表
    #     chr_trans_result = [[] for index in range(25)]
    #
    #     # 读取cis的过滤后文件，并分染色体放到矩阵中
    #     with open(args.input, 'r') as f_cis:
    #         for hic_cis_inter in f_cis:
    #             try:
    #                 chr_index = int(hic_cis_inter.split(" ")[2])
    #             except:
    #                 chr_index = 26
    #             if chr_index == args.index:
    #                 chr_cis_result.append(hic_cis_inter)
    #     # 写入文件并释放内存
    #     cis_file = args.output + "chr_{0}_MAPQ{1}_result.txt".format(args.index, args.MAPQ)
    #
    #     with open(cis_file, "w") as f_cis_out:
    #         f_cis_out.write("".join(chr_cis_result))
    #         chr_cis_result = None
    #
    #     # 读取trans文件，与上部分相同
    #     with open(args.transfile, 'r') as f_trans:
    #         for hic_trans_inter in f_trans:
    #             try:
    #                 chr_index_1 = int(hic_trans_inter.split(" ")[2])
    #                 chr_index_2 = int(hic_trans_inter.split(" ")[5])
    #             except:
    #                 chr_index_1 = chr_index_2 = 26
    #
    #             # 如果与args.index相同，则写入到相应的trans list 位置
    #             if (chr_index_1 == args.index) and (chr_index_2 > args.index):
    #                 chr_trans_result[chr_index_2 - 1].append(hic_trans_inter)
    #
    #             elif (chr_index_2 == args.index) and (chr_index_1 > args.index):
    #                 chr_trans_result[chr_index_1 - 1].append(hic_trans_inter)
    #
    #     # 写入trans文件
    #     trans_file_tmp = "chr_{0}_{1}_MAPQ{2}_result.txt"
    #     for index, trans_result in enumerate(chr_trans_result):
    #         if index + 1 > args.index:
    #             trans_file = args.output + trans_file_tmp.format(args.index, (index + 1), args.MAPQ)
    #
    #             with open(trans_file, "w") as f_out:
    #                 f_out.write("".join(trans_result))
    #                 chr_trans_result[index] = None
    #

    # else:
    #     print("Your <command> option is illegal!")
