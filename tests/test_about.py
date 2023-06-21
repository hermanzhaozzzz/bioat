import subprocess
import bioat
from bioat.cli import Cli

cli = Cli()


def test_about():
    # test about function
    assert hasattr(cli, 'about')
    out = cli.about()
    assert isinstance(out, str)
    # test about text attribute
    assert bioat.about is not None
    assert isinstance(bioat.about, str)


def test_cli_about():
    args = ['bioat', 'about']
    assert subprocess.check_call(args) == 0
