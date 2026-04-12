from types import SimpleNamespace

import bioat
import bioat.cli as cli_module
from bioat.cli import Cli


def test_package_exports_are_lazy(monkeypatch):
    bioat.__dict__.pop("BamTools", None)
    calls = []

    class DummyBamTools:
        pass

    def fake_import_module(name):
        calls.append(name)
        return SimpleNamespace(BamTools=DummyBamTools)

    monkeypatch.setattr(bioat, "import_module", fake_import_module)

    assert "BamTools" not in bioat.__dict__
    assert bioat.BamTools is DummyBamTools
    assert calls == ["bioat.bamtools"]
    assert bioat.BamTools is DummyBamTools
    assert calls == ["bioat.bamtools"]


def test_cli_tools_are_loaded_on_demand(monkeypatch):
    cli = Cli()
    calls = []

    class DummyBedTools:
        pass

    def fake_import_module(name):
        calls.append(name)
        return SimpleNamespace(BedTools=DummyBedTools)

    monkeypatch.setattr(cli_module, "import_module", fake_import_module)

    assert cli._tool_cache == {}
    assert cli.version() >= "0.12.15"
    assert calls == []

    tool = cli.bed
    assert isinstance(tool, DummyBedTools)
    assert calls == ["bioat.bedtools"]
    assert cli.bed is tool
    assert calls == ["bioat.bedtools"]


def test_cli_list_does_not_load_tools(monkeypatch):
    cli = Cli()
    calls = []

    def fake_import_module(name):
        calls.append(name)
        raise AssertionError(name)

    monkeypatch.setattr(cli_module, "import_module", fake_import_module)

    output = cli.list()
    assert "bam" in output
    assert "mpileup2table" in output
    assert "table" in output
    assert "merge" in output
    assert calls == []
