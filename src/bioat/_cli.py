import fire
from bioat import __name__ as name
from bioat import __version__ as version
from bioat import __author__ as author
from bioat import __email__ as email
from bioat import __upgrade_date__ as upgrade_date
from bioat import __pysam_version__
from bioat import __pysamstats_version__


about = f"""
    - bioat version:
        {version}
    - repository address:
        https://github.com/hermanzhaozzzz/BioinformaticAnalysisTools
    - author:
        - name: {author}
        - email: {email}
    - last update:
        - {upgrade_date}
    - key dependencies:
        - pysam == {__pysam_version__}
        - pysamstats == {__pysamstats_version__}
"""

class Cli(object):
    """Cli interface of python package <bioat>

    - bioat is a commandline tool and a python package, which is short for
        Bioinformatic Analysis Tools.
    - bioat has many subcommand to deal with different bio-format:
        BED, BAM, FASTA, FASTQ, VCF, et al.
    """

    def add(self, x, y):
        return x + y

    def multiply(self, x, y):
        return x * y
    def about(self):
        """Print information about <bioat>
        """
        return about


def main() -> int:
    calculator = Cli()
    fire.Fire(calculator, name=name)


if __name__ == "__main__":
    main()
