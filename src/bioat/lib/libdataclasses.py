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
