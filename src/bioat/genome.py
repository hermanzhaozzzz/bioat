import gzip
from Bio import SeqIO
from tqdm import tqdm


class Genome():
    def __init__(self, file):
        self.file = file
        self.genome = gzip.open(self.file, 'rt') if self.file.endswith('.gz') else open(self.file, 'rt')
        # 耗时!
        self.chr_lengths = \
            {record.id: len(record.seq) for record in tqdm(SeqIO.parse(self.genome, "fasta"))}
        self.genome.close()
    # def get_lengths(self):
    #     pass
    #
    # def get_predicted_centromere(self):
    #     pass


if __name__ == '__main__':
    # genome = '/Volumes/zhaohn_HD/Bio/1.database/db_genomes/genome_fa/genome_ucsc_hg38/genome_ucsc_hg38.fa'
    genome = '/Volumes/zhaohn_HD/Bio/1.database/db_genomes/genome_fa/genome_ucsc_hg38/genome_ucsc_hg38.fa2.gz'
    genome = Genome(genome)
    print(genome.chr_lengths)
