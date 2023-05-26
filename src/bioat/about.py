from __future__ import absolute_import
from bioat.version import (__version__, __author__, __upgrade_date__, __email__, __doc_address__, __repo_address__)
from bioat import __name__ as name

about = f"""
BioinformaticAnalysisTools ({name})
    - bioat version:
        {__version__}
    - last update:
        - {__upgrade_date__}
    - doc:
        - address: {__doc_address__}
    - repository address:
        {__repo_address__}
    - author:
        - name: {__author__}
        - email: {__email__}
"""
