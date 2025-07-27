from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import pytest

from bioat.searchtools import SearchTools

DEBUG = False


def test_google_scholar_dataframe_structure():
    if not DEBUG:
        pytest.skip("Skipping...")
    tool = SearchTools()
    df = tool.google_scholar(
        keyword="prime editing",
        n_results=10,
        plot=False,
        output=None,
        sort_by="CitePerYear",
        _for_test=True,
    )
    assert isinstance(df, pd.DataFrame)
    assert set(df.columns) == {
        "Author",
        "Citations",
        "CitePerYear",
        "Year",
        "Venue",
        "Title",
        "Publisher",
        "Source",
    }
    assert len(df) > 0
    assert df["CitePerYear"].dtype in [int, float]


@pytest.mark.parametrize(
    "ext, save_fn",
    [
        (".csv", pd.read_csv),
        (".tsv", lambda f: pd.read_csv(f, sep="\t")),
        (".xlsx", pd.read_excel),
    ],
)
def test_save_table(ext, save_fn, tmp_path):
    if not DEBUG:
        pytest.skip("Skipping...")
    tool = SearchTools()
    filepath = tmp_path / Path(f"output{ext}")
    tool.google_scholar(keyword="CRISPR", n_results=10, output=filepath)
    assert Path.exists(filepath)
    df_loaded = save_fn(filepath)
    assert isinstance(df_loaded, pd.DataFrame)
    assert "Title" in df_loaded.columns


def test_plot_show_called(monkeypatch):
    if not DEBUG:
        pytest.skip("Skipping...")
    called = {"show": False}

    def fake_show():
        called["show"] = True

    monkeypatch.setattr(plt, "show", fake_show)
    tool = SearchTools()
    _ = tool.google_scholar(keyword="CRISPR", n_results=10, plot=True)
    assert called["show"] is True


def test_plot_show_not_called(monkeypatch):
    if not DEBUG:
        pytest.skip("Skipping...")
    called = {"show": False}

    def fake_show():
        called["show"] = True

    monkeypatch.setattr(plt, "show", fake_show)
    tool = SearchTools()
    _ = tool.google_scholar(keyword="CRISPR", n_results=10, plot=False)
    assert called["show"] is False


def test_invalid_output_format():
    if not DEBUG:
        pytest.skip("Skipping...")
    tool = SearchTools()
    with pytest.raises(Exception):  # noqa
        tool.google_scholar(keyword="CRISPR", n_results=10, output="result.txt")
