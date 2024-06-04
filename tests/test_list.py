"""Tests for `list` subcommand.

author: Herman Huanan Zhao
email: hermanzhaozzzz@gmail.com
homepage: https://github.com/hermanzhaozzzz

_description_

>>> example 1:
    bioat list
"""

import subprocess

from bioat.cli import Cli

cli = Cli()


def test_list():
    assert hasattr(cli, 'list')
    out = cli.list()
    assert isinstance(out, str)


def test_cli_list():
    args = ['bioat', 'list']
    subprocess.run(args, check=True)
