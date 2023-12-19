from bioat.lib.libfastx import casfinder


class CrisprTools:
    """CRISPR mining toolbox."""

    def casfinder(
        self,
    ):
        prodigal = (
            "/Users/zhaohuanan/micromamba/envs/snakepipes_Cas-mining/bin/prodigal"
        )
        pilercr = "/Users/zhaohuanan/micromamba/envs/snakepipes_Cas-mining/bin/pilercr"

        casfinder(
            input_fa="/Users/zhaohuanan/Downloads/test/202155.assembled.fna",
            output_faa="/Users/zhaohuanan/Downloads/test/final.cas.faa",
            lmin=3001,  # 3001 in Nature Methods paper
            lmax=None,
            extend=10_000,
            temp_dir=None,
            prodigal=prodigal,
            pilercr=pilercr,
            log_level="DEBUG",
        )

    # def cas12_finder(self):
    #     """Cas12 mining toolbox.
    #
    #     Under development!
    #     """
    #     pass
    # def cas13_finder(self):
    #     """Cas13 mining toolbox
    #     Under development!
    #     """
    #     pass


if __name__ == "__main__":
    pass
