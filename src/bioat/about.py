from bioat.version import (__version__, __author__, __upgrade_date__, __email__, __doc_address__, __repo_address__,
                           __issue_address__)
from bioat import __name__ as name

about = f"""
BioinformaticAnalysisTools ({name})
    - bioat version:
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
        
    **Copyright**:
        For researchers: the authors appreciate the citations to this tool (citations are not mandatory),
            please cite my papers:
                1. aaaaaaa
                2. bbbbbbb
        For commercial use:
            NOT PERMITTED unless permission is obtained from the author
"""
