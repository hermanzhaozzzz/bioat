#!/usr/bin/env python
# encoding: utf-8
import argparse

import pandas as pd


def get_parser():
    """Get param and return a parser object."""
    parser = argparse.ArgumentParser(
        description="convert MGI md5 file to md5sum supported file"
    )
    parser.add_argument("filename", type=str, help="input MGI md5 file")
    return parser


def main(filename):
    """Read mgi-like md5 file and convert to a normal md5 file."""
    df = pd.read_csv(
        filename,
        header=None,
        index_col=False,
        sep="  ",
        engine="python",
    )
    df.columns = ["md5", "filename"]
    df.md5 = df.md5.map(lambda x: x.lower())
    df.filename = df.filename + ".gz"
    df.to_csv(
        filename.replace(".txt", "") + ".fix.md5", header=False, index=False, sep="\t"
    )


if __name__ == "__main__":
    parser = get_parser()
    args = parser.parse_args()
    main(filename=args.filename)
