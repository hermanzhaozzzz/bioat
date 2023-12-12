from bioat.lib.libfastx import casfinder


class CrisprTools:
    """CRISPR mining toolbox."""

    def casfinder(
            self,

    ):
        casfinder(
            input_fa="/Users/zhaohuanan/Downloads/test/202155.assembled.fna",
            output_faa=None,
            lmin=3001,  # 3001 in Nature Methods paper
            lmax=None
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


# """
# alias lipercr=/Users/zhaohuanan/micromamba/envs/snakepipes_Cas-mining/bin/pilercr
# alias prodigal=/Users/zhaohuanan/micromamba/envs/snakepipes_Cas-mining/bin/prodigal
# alias bedtools=/Users/zhaohuanan/micromamba/envs/snakepipes_Cas-mining/bin/bedtools
# """

if __name__ == "__main__":
    pass
