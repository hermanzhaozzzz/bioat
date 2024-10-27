"""
crisprtools.py

This module provides a toolbox for mining CRISPR-related sequences in metagenomic data.
It includes functionality for identifying Cas candidates associated with CRISPR loci
and for annotating Cas13 candidates from protein sequences. The main class, `CrisprTools`,
provides methods for calling external executables such as Prodigal and Pilercr for protein
prediction and CRISPR identification.

Classes:
    CrisprTools: A class containing methods for identifying Cas candidates and Cas13
                 candidates from genomic data.

Methods:
    cas_finder: A method for de novo annotation of Cas candidates from CRISPR loci,
                utilizing input fasta files and producing various output files.
    cas13_finder: A method for the annotation of Cas13 candidates from protein
                  fasta files.

Usage:
    To use this module, create an instance of `CrisprTools` and call the methods
    `cas_finder` or `cas13_finder` with the appropriate parameters to perform
    CRISPR analysis on your dataset.
"""

from bioat.lib.libfastx import cas13_finder, cas_finder
from bioat.lib.libpath import check_executable


class CrisprTools:
    """CRISPR mining toolbox.

    This class provides methods for performing CRISPR analysis on datasets, including
    finding Cas proteins and CRISPR sequences.

    Attributes:
        None
    """

    def __init__(self):
        pass

    def cas_finder(
        self,
        input_fa,
        output_faa=None,
        output_contig_fa=None,
        output_crispr_info_tab=None,
        lmin=3000,  # 3001 in Nature Methods paper
        lmax=None,
        extend=10_000,
        temp_dir=None,
        prodigal=None,
        prodigal_mode="meta",
        pilercr=None,
        rm_temp=True,
        log_level="INFO",
    ):
        """De novo annotation for Cas candidates from neighbors of CRISPR loci.

        Args:
            input_fa (str): Path to the input metagenome fasta file containing many contigs.
            output_faa (str, optional): Path to save the de novo annotated Cas candidates.
            output_contig_fa (str, optional): Path to save the whole contigs of de novo annotated Cas candidates.
            output_crispr_info_tab (str, optional): Path to save the de novo annotated CRISPR info table (CSV format).
            lmin (int, optional): Minimum length for a contig. Defaults to 3000.
            lmax (int, optional): Maximum length for a contig. Defaults to None.
            extend (int, optional): Distance over which proteins are considered from the start/end of the CRISPR loci. Defaults to 10000.
            temp_dir (str, optional): Directory to store temporary files. Defaults to None.
            prodigal (str, optional): Path to the Prodigal executable. Defaults to None.
            prodigal_mode (str, optional): Mode for Prodigal annotation. Can be "meta" or "single". Defaults to "meta".
            pilercr (str, optional): Path to the Pilercr executable. Defaults to None.
            rm_temp (bool, optional): If False, temporary files will be kept. Defaults to True.
            log_level (str, optional): Logging level; set to "DEBUG" to see logs from prodigal and pilercr. Defaults to "INFO".
        """
        check_executable(name=prodigal)
        check_executable(name=pilercr)
        cas_finder(
            input_fa=input_fa,
            output_faa=output_faa,
            output_contig_fa=output_contig_fa,
            output_crispr_info_tab=output_crispr_info_tab,
            lmin=lmin,
            lmax=lmax,
            extend=extend,
            temp_dir=temp_dir,
            prodigal=prodigal,
            prodigal_mode=prodigal_mode,
            pilercr=pilercr,
            rm_temp=rm_temp,
            log_level=log_level,
        )

    # def cas12_finder(self):
    #     """Cas12 mining toolbox.
    #
    #     Under development!
    #     """
    #     pass

    def cas13_finder(
        self, input_faa, output_faa=None, lmin=200, lmax=1500, log_level="INFO"
    ):
        """De novo annotation for Cas13 candidates from proteins.faa.

        Args:
            input_faa (str): The input file containing Cas candidates in .faa format.
            output_faa (str, optional): The output file for Cas13 candidates in .faa format. Defaults to None.
            lmin (int, optional): Minimum length for a Cas candidate. Defaults to 200.
            lmax (int, optional): Maximum length for a Cas candidate. Defaults to 1500.
            log_level (str, optional): The logging level. Options are 'CRITICAL', 'ERROR', 'WARNING', 'INFO', 'DEBUG', 'NOTSET'. Defaults to 'INFO'.
        """
        cas13_finder(
            input_faa=input_faa,
            output_faa=output_faa,
            lmin=lmin,
            lmax=lmax,
            log_level=log_level,
        )


if __name__ == "__main__":
    pass
