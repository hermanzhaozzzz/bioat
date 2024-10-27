"""test_searchtools

author: Herman Huanan Zhao
email: hermanzhaozzzz@gmail.com
homepage: https://github.com/hermanzhaozzzz

_description_

>>> example 1:
    _example_

>>> example 2:
    _example_
"""


import pytest

from bioat.cli import Cli

# from .settings import DATA_PATH, TESTOUT_PATH

bioat_cli = Cli()

DEBUG = False


def test_google_scholar():
    if not DEBUG:
        pytest.skip('skip test for google scholar')
    bioat_cli.search.google_scholar(keyword="test", save_table=False)


def test_google_scholar_no_keyword():
    if not DEBUG:
        pytest.skip('skip test for google scholar')
    bioat_cli.search.google_scholar(save_table=False)


def test_google_scholar_xlsx():
    if not DEBUG:
        pytest.skip('skip test for google scholar')
    bioat_cli.search.google_scholar(
        keyword="test", output="/dev/null/xlsx.xlsx", save_table=False
    )
