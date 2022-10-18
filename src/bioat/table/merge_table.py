import pandas as pd
from sys import stdout
from fire import Fire


class MergeTableCli():
    """A simple tool to merge tables from different sample.

    """

    def from_csv(
            self,
            inputs: str,
            tags: str,
            output=stdout.name,
            output_fmt='csv',
            input_header=True,
            output_header=True):
        """Merge from csv files

        Params
        :param inputs: csv files joint by comma
        :param tags: tags for each file
        :param output: output file name
        :param output_fmt: define the output format
        :param input_header: True or False
        :param output_header: True of False
        """
        return merge_table(
            inputs=inputs,
            tags=tags,
            output=output,
            input_fmt='csv',
            output_fmt=output_fmt,
            input_header=input_header,
            output_header=output_header)

    def from_tsv(
            self,
            inputs: str,
            tags: str,
            output=stdout.name,
            output_fmt='tsv',
            input_header=True,
            output_header=True):
        """Merge from tsv files
        Params
        :param inputs: csv files joint by comma
        :param tags: tags for each file
        :param output: output file name
        :param output_fmt: define the output format
        :param input_header: True or False
        :param output_header: True of False
        """
        return merge_table(
            inputs=inputs,
            tags=tags,
            output=output,
            input_fmt='tsv',
            output_fmt=output_fmt,
            input_header=input_header,
            output_header=output_header)


def deal_with_comma_params(params) -> list:
    params = list(params) if isinstance(params, tuple) else params.split(',')
    return [i.strip() for i in params]


def merge_table(
        inputs,
        tags,
        output,
        input_fmt,
        output_fmt,
        input_header=True,
        output_header=True):
    # fix params
    inputs = deal_with_comma_params(inputs)
    tags = deal_with_comma_params(tags)
    input_header = 0 if input_header else None
    output = stdout if output == stdout.name else output
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

    df.to_csv(output, sep=dt_sep[output_fmt], header=output_header, index=None, )


def main():
    Fire(MergeTableCli, name='merge-table')


if __name__ == "__main__":
    main()
