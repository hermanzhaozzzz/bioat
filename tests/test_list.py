import subprocess
from bioat.cli import Cli

cli = Cli()


def test_list():
    assert hasattr(cli, 'list')
    out = cli.list()
    assert isinstance(out, str)


def test_cli_list():
    args = ['bioat', 'list']
    assert subprocess.run(args).returncode == 0
