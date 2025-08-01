"""Tests for `version` module."""

import subprocess

from bioat.cli import Cli

cli = Cli()


def test_version():
    assert hasattr(cli, "version")
    out = cli.version()
    assert isinstance(out, str)


def test_cli_version():
    args = ["bioat", "version"]
    assert subprocess.run(args, check=False).returncode == 0  # noqa
