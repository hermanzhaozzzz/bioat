import pandas as pd


class MgiTools():
    """Raw MGI sequencing data toolbox."""

    def __init__(self):
        pass

    def parse_md5(self, file: str):
        """Read mgi-like md5 file and convert to a normal md5 file.

        :param file: file name of a mgi-like md5 file.
        """
        df = pd.read_csv(
            file,
            header=None,
            index_col=False,
            sep="  ",
            engine="python",
        )
        df.columns = ["md5", "filename"]
        df.md5 = df.md5.map(lambda x: x.lower())
        df.filename = df.filename + ".gz"
        df.to_csv(file.replace(".txt", "") + ".fix.md5", header=False, index=False, sep="\t")


if __name__ == "__main__":
    pass
