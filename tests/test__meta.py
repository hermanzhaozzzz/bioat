from importlib.metadata import metadata

import bioat

META = metadata("bioat")  # __META__


def test_metadata():
    assert META["Name"] == "bioat"  # __PKG_NAME__
    assert "Author" in META  # __AUTHOR__
    assert "Author-email" in META  # __AUTHOR_EMAIL__
    assert "Summary" in META  # __DESCRIPTION__
    assert "License" in META  # __LICENSE__
    assert "Version" in META  # __VERSION__


def test_doc_format():
    assert bioat.__DOC_FORMAT__ == "google"


def test_doc_page():
    assert bioat.__DOC_PAGE__ == "https://bioat.readthedocs.io/en/latest/"


def test_issue_page():
    assert bioat.__ISSUE_PAGE__ == "https://github.com/hermanzhaozzzz/bioat/issues"


def test_license():
    assert bioat.__LICENSE__ == "Apache-2.0"


def test_version():
    assert bioat.__VERSION__ >= "0.12.15"
