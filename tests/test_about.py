"""Test for `about` subcommand

author: Herman Huanan Zhao
email: hermanzhaozzzz@gmail.com
homepage: https://github.com/hermanzhaozzzz

>>> example 1:
    bioat about
"""

import subprocess

import bioat
from bioat.cli import Cli

cli = Cli()


def test_about():
    """test about function"""
    assert hasattr(cli, "about")
    assert isinstance(cli.about(), str)
    # test about text attribute
    assert isinstance(cli.about(), str)


def test_cli_about():
    """test about cmd"""
    args = ["bioat", "about"]
    assert subprocess.run(args, check=True)
