"""
test_searchtools
"""
import os
import pytest
import subprocess

from bioat.cli import Cli

# from .settings import DATA_PATH, TESTOUT_PATH

cli = Cli()

DEBUG = False


def test_google_scholar():
    if not DEBUG:
        pytest.skip('skip test for google scholar')
    cli.search.google_scholar(keyword='test', save_table=False)


def test_google_scholar_no_keyword():
    if not DEBUG:
        pytest.skip('skip test for google scholar')
    cli.search.google_scholar(save_table=False)


def test_google_scholar_xlsx():
    if not DEBUG:
        pytest.skip('skip test for google scholar')
    cli.search.google_scholar(keyword='test', output='/dev/null/xlsx.xlsx', save_table=False)
