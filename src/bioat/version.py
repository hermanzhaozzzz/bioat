with open('../../pyproject.toml', 'rt') as f:
    __version__ = [i.rstrip().split('= ')[-1]
                   for i in f.readlines() if 'version =' in i][0].replace('"', '')

__upgrade_date__ = '2023-03-31'
__author__ = 'Hua-nan ZHAO @ Tsinghua University'
__email__ = 'hermanzhaozzzz@gmail.com'
__doc_format__ = "restructuredtext"
