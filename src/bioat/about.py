"""Cli about

author: Herman Huanan Zhao
email: hermanzhaozzzz@gmail.com
homepage: https://github.com/hermanzhaozzzz

>>> example 1:
    bioat about
"""

from bioat.logger import __PKG_NAME__
from bioat.version import (
    __author__,
    __doc_address__,
    __email__,
    __issue_address__,
    __repo_address__,
    __upgrade_date__,
    __version__,
)

__ABOUT__ = f"""
BioAT ({__PKG_NAME__})
    - {__PKG_NAME__} version:
        {__version__}
    - last update:
        - {__upgrade_date__}
    - repository page:
        {__repo_address__}
    - doc page:
        {__doc_address__}
    - issue / new feature page:
        {__issue_address__}
    - author:
        - name: {__author__}
        - email: {__email__}
    ---
"""
