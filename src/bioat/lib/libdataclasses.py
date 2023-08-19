from dataclasses import dataclass


@dataclass
class Bed:
    chromosome: str
    start: int
    end: int
    name: str
    score: int
    strand: str


@dataclass
class Bam:
    pass


@dataclass
class Fasta:
    header: str
    sequence: str

    def __init__(self):
        pass
    # def format_fasta(self, length: 80):
    #     # 格式化输出fasta的sequence
    #     pass


@dataclass
class Fastq:
    # https://upload-images.jianshu.io/upload_images/7976641-18ceaafbce3d93d7.png?imageMogr2/auto-orient/strip|imageView2/2/w/1199/format/webp
    header: str
    sequence: str
    extras: str
    quality: str

    # def get_phread_info(self):
    #     # 判断是Phred33还是Phred64
    #     # https://www.jianshu.com/p/248308513e2e
    #     pass


@dataclass
class VCF:
    pass