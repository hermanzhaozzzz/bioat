from __future__ import absolute_import
from bioat.version import (__version__, __author__, __upgrade_date__, __email__, __doc_format__)
from bioat import __name__ as name
about = f"""
BioinformaticAnalysisTools ({name})
    - bioat version:
        {__version__}
    - repository address:
        https://github.com/hermanzhaozzzz/bioat

    - author:
        - name: {__author__}
        - email: {__email__}
    - doc:
        - format: {__doc_format__}
    - last update:
        - {__upgrade_date__}
"""
