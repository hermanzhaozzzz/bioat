from bioat.lib.libpath import HOME
from bioat.lib.libfastx import cas_finder, cas13_finder


class CrisprTools:
    """CRISPR mining toolbox."""

    _prodigal = f"{HOME}/micromamba/envs/snakepipes_Cas-mining/bin/prodigal"
    _pilercr = f"{HOME}/micromamba/envs/snakepipes_Cas-mining/bin/pilercr"

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
