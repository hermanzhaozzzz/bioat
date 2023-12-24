from bioat.lib.libpath import HOME
from bioat.lib.libfastx import cas_finder, cas13_finder


class CrisprTools:
    """CRISPR mining toolbox."""

    prodigal = f"{HOME}/micromamba/envs/snakepipes_Cas-mining/bin/prodigal"
    pilercr = f"{HOME}/micromamba/envs/snakepipes_Cas-mining/bin/pilercr"

    def casfinder(
        self,
    ):
        cas_finder(
            input_fa="/Users/zhaohuanan/Downloads/test/202155.assembled.fna",
            output_faa="/Users/zhaohuanan/Downloads/test/final.cas.faa",
            lmin=3000,  # 3001 in Nature Methods paper
            lmax=None,
            extend=10_000,
            temp_dir=None,
            prodigal=self.prodigal,
            pilercr=self.pilercr,
            rm_temp=True,
            log_level="INFO",
        )

    # def cas12_finder(self):
    #     """Cas12 mining toolbox.
    #
    #     Under development!
    #     """
    #     pass
    def cas13_finder(
        self,
    ):
        """Cas13 mining toolbox
        Under development!
        """
        # casfinder(
        #     input_fa="/Users/zhaohuanan/Downloads/test/202155.assembled.fna",
        #     output_faa="/Users/zhaohuanan/Downloads/test/final.cas.faa",
        #     lmin=3000,  # 3001 in Nature Methods paper
        #     lmax=None,
        #     extend=10_000,
        #     temp_dir=None,
        #     prodigal=self.prodigal,
        #     pilercr=self.pilercr,
        #     rm_temp=True,
        #     log_level="INFO",
        # )

        # start to test HEPN filter
        cas13_finder(
            input_faa=f"{HOME}/Downloads/test/final.cas.faa",
            # output_faa=f"{HOME}/Downloads/test/final.cas.HEPN0.faa",
            # output_faa=f"{HOME}/Downloads/test/final.cas.HEPN1.faa",
            output_faa=f"{HOME}/Downloads/test/final.cas.HEPN2.faa",
            log_level="DEBUG",
        )


if __name__ == "__main__":
    pass
