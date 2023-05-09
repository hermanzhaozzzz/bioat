from bioat.cli import Cli


def test_list():
    cli = Cli()
    assert hasattr(cli, 'version')
    out = cli.list()
    assert isinstance(out, str)