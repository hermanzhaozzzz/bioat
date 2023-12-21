from bioat.lib.libfastx import casfinder


class CrisprTools:
    """CRISPR mining toolbox."""

    prodigal = "/Users/zhaohuanan/micromamba/envs/snakepipes_Cas-mining/bin/prodigal"
    pilercr = "/Users/zhaohuanan/micromamba/envs/snakepipes_Cas-mining/bin/pilercr"

    def casfinder(
        self,
    ):
        casfinder(
            input_fa="/Users/zhaohuanan/Downloads/test/202155.assembled.fna",
            output_faa="/Users/zhaohuanan/Downloads/test/final.cas.faa",
            lmin=3000,  # 3001 in Nature Methods paper
            lmax=None,
            extend=10_000,
            temp_dir=None,
            prodigal=self.prodigal,
            pilercr=self.pilercr,
            rm_temp=True,
            log_level="DEBUG",
        )

    # def cas12_finder(self):
    #     """Cas12 mining toolbox.
    #
    #     Under development!
    #     """
    #     pass
    def cas13_finder(
        self,
        input_fa,
    ):
        """Cas13 mining toolbox
        Under development!
        """
        casfinder(
            input_fa="/Users/zhaohuanan/Downloads/test/202155.assembled.fna",
            output_faa="/Users/zhaohuanan/Downloads/test/final.cas.faa",
            lmin=3000,  # 3001 in Nature Methods paper
            lmax=None,
            extend=10_000,
            temp_dir=None,
            prodigal=self.prodigal,
            pilercr=self.pilercr,
            rm_temp=True,
            log_level="DEBUG",
        )


if __name__ == "__main__":
    pass
