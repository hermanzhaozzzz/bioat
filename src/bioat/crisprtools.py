from bioat.lib.libpath import HOME, exists_in_PATH
from bioat.lib.libfastx import cas_finder, cas13_finder


class CrisprTools:
    """CRISPR mining toolbox."""

    # for development test
    _prodigal = (
        "prodigal"
        if exists_in_PATH("prodigal")
        else f"{HOME}/micromamba/envs/snakepipes_Cas-mining/bin/prodigal"
    )
    _pilercr = (
        "pilercr"
        if exists_in_PATH("pilercr")
        else f"{HOME}/micromamba/envs/snakepipes_Cas-mining/bin/pilercr"
    )
    # /for development test

    def cas_finder(
        self,
        input_fa,
        output_faa=None,
        lmin=3000,  # 3001 in Nature Methods paper
        lmax=None,
        extend=10_000,
        temp_dir=None,
        prodigal=_prodigal,
        pilercr=_pilercr,
        rm_temp=True,
        log_level="INFO",
    ):
        """De novo annotation for Cas candidates from neighbor of CRISPR loci.

        :param input_fa: metagenome.fa (with many contigs in it)
        :param output_faa: De novo annotated Cas candidates
        :param lmin: min length for a contig
        :param lmax: max length for a contig
        :param extend: proteins are considered over how far from the start/end of the CRISPR loci
        :param temp_dir: folder to put temp files in
        :param prodigal: the executable prodigal path
        :param pilercr: the executable pilercr path
        :param rm_temp: set False to keep the temp files
        :param log_level: set DEBUG to see output from prodigal and pilercr
        """
        cas_finder(
            input_fa=input_fa,
            output_faa=output_faa,
            lmin=lmin,
            lmax=lmax,
            extend=extend,
            temp_dir=temp_dir,
            prodigal=prodigal,
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

        :param input_faa: cas_candidates.faa
        :param output_faa: cas13_candidates.faa
        :param lmin: min length for a cas candidate
        :param lmax: max length for a cas candidate
        :param log_level: 'CRITICAL', 'ERROR', 'WARNING', 'INFO', 'DEBUG', 'NOTSET'
        """
        """Cas13 mining toolbox."""
        # start to test HEPN filter
        cas13_finder(
            input_faa=input_faa,
            output_faa=output_faa,
            lmin=lmin,
            lmax=lmax,
            log_level=log_level,
        )


if __name__ == "__main__":
    pass
