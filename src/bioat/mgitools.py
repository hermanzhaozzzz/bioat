from __future__ import absolute_import
from bioat.logger import get_logger
import pandas as pd

__module_name__ = 'bioat.mgitools'


class MgiTools():
    """Raw MGI sequencing data toolbox."""

    def __init__(self):
        pass

    def parse_md5(self, file: str, log_level='WARNING'):
        """Read mgi-like md5 file and convert to a normal md5 file.

        :param file: file name of a mgi-like md5 file.
        :param log_level: INFO/DEBUG/WARNING/ERROR, default is WARNING.
        """
        logger = get_logger(level=log_level, module_name=__module_name__, func_name=sys._getframe().f_code.co_name)
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
