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
            input_fmt='tsv',
            output_fmt='tsv',
            input_header=False,
            output_header=False
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
        """
        # fix params
        inputs = list(inputs) if isinstance(inputs, tuple) else inputs.split(',')
        tags = list(tags) if isinstance(tags, tuple) else tags.split(',')
        input_header = 0 if input_header else None
        output = sys.stdout if output == sys.stdout.name else output
        dt_sep = {'csv': ',', 'tsv': '\t'}

        df = pd.DataFrame()

        for file, tag in zip(inputs, tags):
            _df = pd.read_csv(
                file,
                sep=dt_sep[input_fmt],
                header=input_header,
                index_col=None)
            _df.insert(0, column='<sample>', value=tag)
            df = pd.concat([df, _df], axis=0, ignore_index=True)

        df.to_csv(output, sep=dt_sep[output_fmt], header=output_header, index=None)


if __name__ == "__main__":
    pass
