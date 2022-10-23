from sys import stdout
from fire import Fire
from bioat import set_logging_level
from bioat.fasta._filter_n import main as filter_n


class Cli():
    """Fasta tools
    """

    def filter_n(self, input, output=stdout.name):
        """filter reads contains N

        Params
        :param str input: path of input <fasta|fasta.gz>
        :param str output: path of output <fasta|fasta.gz>
        """
        return filter_n(input, output)


def main() -> int:
    cli = Cli()
    Fire(cli, name='bioat_fastatools')


if __name__ == '__main__':
    main()
