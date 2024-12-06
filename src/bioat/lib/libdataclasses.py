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


@dataclass(frozen=True)
class Fasta:
    """
    usage:
        fa = Fasta('this_header', 'AGCTACGTCCCCTGA')
        fa
        print(fa)
    """
    header: str
    sequence: str

    def __eq__(self, other):
        if self.header == other.header and self.sequence == other.sequence:
            return True
        else:
            return False

    def __len__(self):
        return len(self.sequence)

    @property
    def length(self):
        return len(self.sequence)

    def __str__(self, line_width: int = 20):
        ls = []

        for i in range(self.length // line_width + 1):
            ls.append(self.sequence[i * line_width:(i + 1) * line_width])
        s = '\n'.join(ls)
        s = f'>{self.header}\n{s}'
        return s

    def __repr__(self, show_all: bool = False, show_length: int = 10):
        if show_all:
            return f"Fasta(header='{self.header}', sequence='{self.sequence}')"
        else:
            if self.length > show_length:
                s = self.sequence[: show_length // 2] + '...' + self.sequence[-show_length // 2:]
                return f"Fasta(header='{self.header}', sequence='{s}')"
            else:
                return f"Fasta(header='{self.header}', sequence='{self.sequence}')"


@dataclass(frozen=True)
class Fastq(Fasta):
    header: str
    sequence: str
    extras: str
    quality: str

    def __eq__(self, other, check_quality: bool = False, check_extras: bool = False):
        # rewrite from class Fasta
        checker = [self.sequence == other.sequence]

        if check_extras:
            checker.append(self.extras == other.extras)
        if check_quality:
            checker.append(self.quality == other.quality)

        if all(checker):
            return True
        else:
            return False

    def __str__(self, line_width: int = 20):
        # rewrite from class Fasta
        return f'@{self.header}\n{self.sequence}\n{self.extras}\n{self.quality}\n'

    def __repr__(self):
        # rewrite from class Fasta
        return (f"Fastq(header='{self.header}', sequence='{self.sequence}', "
                f"extras='{self.extras}', quality='{self.quality}')")

    def is_valid(self):
        if self.length == len(self.quality):
            return True
        else:
            return False

    # def get_phread_info(self):
    #     # 判断是Phred33还是Phred64
    #     # https://www.jianshu.com/p/248308513e2e
    #     pass


@dataclass(frozen=True)
class Assembly(object):
    contigs: list[Fasta]
    length: int
    path: str

    def __init__(self, path: str):
        pass


@dataclass
class VCF:
    pass
