"""
about.py

This module provides information about the BioAT package, including
version, author details, and links to related resources.

Author:
    Herman Huanan Zhao (hermanzhaozzzz AT gmail.com)

Example:
    To display the information about the package, run the command:
        bioat about

Attributes:
    __ABOUT__ (str): A formatted string containing information about
                     the BioAT package, including version, repository
                     page, documentation page, issue tracking page,
                     and author details.
"""

from bioat import (
    __AUTHOR__,
    __AUTHOR_EMAIL__,
    __DOC_PAGE__,
    __HOME_PAGE__,
    __ISSUE_PAGE__,
    __PKG_NAME__,
    __VERSION__,
)

__ABOUT__ = f"""\
BioAT ({__PKG_NAME__})
    - {__PKG_NAME__} version:
        {__VERSION__}
    - Home page:
        {__HOME_PAGE__}
    - Doc page:
        {__DOC_PAGE__}
    - Issue & feature request page:
        {__ISSUE_PAGE__}
    - Author:
        - name: {__AUTHOR__}
        - email: {__AUTHOR_EMAIL__}
"""
