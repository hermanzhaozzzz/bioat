"""Doc.
dataclasses:
    - https://docs.python.org/zh-cn/3/library/dataclasses.html
    - python version >= 3.7.0, PEP 557, PEP 526
    - 这个模块提供了一个装饰器和一些函数，用于自动为用户自定义的类添加生成的 special method 例如 __init__() 和 __repr__()。

---

demo:
    from dataclasses import dataclass

    @dataclass
    class InventoryItem:
        '''annotates'''
        name: str
        unit_price: float
        quantity_on_hand: int = 0

        def total_cost(self) -> float:
            return self.unit_price * self.quantity_on_hand

    # 除其他内容以外，还将添加如下所示的 __init__():
    def __init__(self, name: str, unit_price: float, quantity_on_hand: int = 0):
        self.name = name
        self.unit_price = unit_price
        self.quantity_on_hand = quantity_on_hand

    # 请注意，此方法会自动添加到类中：而不是在如上所示的 InventoryItem 定义中被直接指定。


    # 下面这三种 dataclass() 用法是等价的
    @dataclass
    class C:
        ...

    @dataclass()
    class C:
        ...

    @dataclass(init=True, repr=True, eq=True, order=False, unsafe_hash=False, frozen=False,
               match_args=True, kw_only=False, slots=False, weakref_slot=False)
    class C:
        ...

    具体说明查看官方文档
"""

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
