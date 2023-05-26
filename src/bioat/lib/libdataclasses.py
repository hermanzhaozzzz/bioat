import dataclasses


@dataclasses.dataclass
class Bed:
    chromosome: str
    start: int
    end: int
    name: str
    score: int | str
    strand: str


@dataclasses.dataclass
class Bam:
    pass


@dataclasses.dataclass
class Fasta:
    header: str
    sequence: str

    def format_fasta(self, length: 80):
        # 格式化输出fasta的sequence
        pass


@dataclasses.dataclass
class Fastq:
    # https://upload-images.jianshu.io/upload_images/7976641-18ceaafbce3d93d7.png?imageMogr2/auto-orient/strip|imageView2/2/w/1199/format/webp
    header: str
    sequence: str
    extras: str
    quality: str

    def get_phread_info(self):
        # 判断是Phred33还是Phred64
        # https://www.jianshu.com/p/248308513e2e
        pass


@dataclasses.dataclass
class VCF:
    pass