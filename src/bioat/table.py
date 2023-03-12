import sys
import pandas as pd


class Table():
    def __init__(self):
        pass

    def merge(
            self,
            inputs,
            tags,
            output,
            input_fmt,
            output_fmt='csv',
            input_header=True,
            output_header=True
    ):
        """A simple tool to merge same formatted tables from different sample.

        Params
        :param inputs: csv files joint by comma
        :param tags: tags for each file
        :param output: output file name
        :param output_fmt: csv | tsv, define the output format
        :param input_header: True or False
        :param output_header: True of False
        """
        # fix params
        inputs = list(inputs) if isinstance(inputs, tuple) else inputs.split(',')
        tags = list(tags) if isinstance(tags, tuple) else tags.split(',')
        input_header = 0 if input_header else None
        output = sys.stdout if output == sys.stdout.name else output
        dt_sep = {'csv': ',', 'tsv': '\t'}

        df = pd.DataFrame()

        for input, tag in zip(inputs, tags):
            _df = pd.read_csv(
                input,
                sep=dt_sep[input_fmt],
                header=input_header,
                index_col=None)
            _df.insert(0, column='<sample>', value=tag)
            df = pd.concat([df, _df], axis=0, ignore_index=True)

        df.to_csv(output, sep=dt_sep[output_fmt], header=output_header, index=None)


if __name__ == "__main__":
    pass
