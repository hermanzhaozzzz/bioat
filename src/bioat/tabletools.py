"""_summary_.

author: Herman Huanan Zhao
email: hermanzhaozzzz@gmail.com
homepage: https://github.com/hermanzhaozzzz

_description_

example 1:
    bioat list
        <in shell>:
            $ bioat list
        <in python consolo>:
            >>> from bioat.cli import Cli
            >>> bioat = Cli()
            >>> bioat.list()
            >>> print(bioat.list())

example 2:
    _example_
"""

import sys

import pandas as pd

from bioat.logger import LoggerManager

lm = LoggerManager(mod_name="bioat.tabletools")


class TableTools:
    """To integrate tables."""

    lm.set_names(cls_name="TableTools")

    def __init__(self):
        pass

    def merge(
        self,
        inputs,
        tags,
        output,
        input_fmt="tsv",
        output_fmt="tsv",
        input_header=False,
        output_header=False,
        log_level="WARNING",
    ):
        """A simple tool to merge same formatted tables from different sample.

        Params
        :param inputs: table files
        :param tags: tags for each table
        :param output: merged file
        :param input_fmt: tsv | csv
        :param output_fmt: tsv | csv
        :param input_header: True | False, input table has header or not
        :param output_header: True | False, output table has header or not
        :param log_level: log status
        """
        lm.set_names(func_name="merge")
        lm.set_level(log_level)

        # fix params
        inputs = list(inputs) if isinstance(inputs, tuple) else inputs.split(",")
        tags = list(tags) if isinstance(tags, tuple) else tags.split(",")

        lm.logger.debug(f"input_fmt={input_fmt}")
        lm.logger.debug(f"output_fmt={output_fmt}")
        lm.logger.debug(f"input_header={input_header}")
        lm.logger.debug(f"output_header={output_header}")

        input_header = 0 if input_header else None
        output = sys.stdout if output == sys.stdout.name else output
        dt_sep = {"csv": ",", "tsv": "\t"}

        df = pd.DataFrame()

        for file, tag in zip(inputs, tags, strict=False):
            _df = pd.read_csv(
                file, sep=dt_sep[input_fmt], header=input_header, index_col=None
            )
            _df.insert(0, column="<sample>", value=tag)
            df = pd.concat([df, _df], axis=0, ignore_index=True)

        df.to_csv(output, sep=dt_sep[output_fmt], header=output_header, index=None)

    def split(
        self,
        input: str,
        n: int,
        output_prefix=None,
        input_fmt="tsv",
        output_fmt="tsv",
        input_header=False,
        output_header=False,
        compress=False,
        log_level="WARNING",
    ):
        """A simple tool to split table into parts.

        Params
        :param input: table to split
        :param n: split table into n parts
        :param output_prefix: name prefix for splitted parts, the same with input if not defined
        :param input_fmt: tsv | csv
        :param output_fmt: tsv | csv
        :param input_header: True | False, input table has header or not
        :param output_header: True | False, output table has header or not
        :param compress: True | False, gzip the output table or not
        :param log_level: log status
        """
        lm.set_names(func_name="split")
        lm.set_level(log_level)

        lm.logger.debug(f"input_fmt={input_fmt}")
        lm.logger.debug(f"output_fmt={output_fmt}")
        lm.logger.debug(f"input_header={input_header}")
        lm.logger.debug(f"output_header={output_header}")
        # fix params
        input_header = 0 if input_header else None
        output_prefix = output_prefix if output_prefix else f"{input}_"
        compress = ".gz" if compress else ""

        dt_sep = {"csv": ",", "tsv": "\t"}
        # load table
        df = pd.read_csv(input, header=input_header, sep=dt_sep[input_fmt])

        # max lines for each part table: chunk_size
        chunk_size = df.shape[0] // n if df.shape[0] % n == 0 else df.shape[0] // n + 1

        # write out
        for idx in range(n):
            df_write_out = df.iloc[idx * chunk_size : (idx + 1) * chunk_size,]
            df_write_out.to_csv(
                f"{output_prefix}{idx}.%s{compress}" % output_fmt,
                index=None,
                header=output_header,
                sep=dt_sep[output_fmt],
            )


if __name__ == "__main__":
    pass
