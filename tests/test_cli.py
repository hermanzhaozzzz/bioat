from bioat.cli import Cli

bioat_cli = Cli()


def test_about():
    assert isinstance(bioat_cli.about(), str)
    assert len(bioat_cli.about()) >= 10


def test_list_happy_path():
    assert bioat_cli.list() != ""  # list should not be empty even with no subcommands
    result = bioat_cli.list()
    sub_cmds = [
        # "bam",
        "bed",
        "crispr",
        "fastx",
        "meta",
        "search",
        # "system",
        "table",
        "target_seq",
    ]
    for sub_cmd in sub_cmds:
        assert sub_cmd in result
    assert "_private" not in result  # ensure private attributes are not included


def test_version():
    assert isinstance(bioat_cli.version(), str)
    assert bioat_cli.version() >= "0.12.15"
