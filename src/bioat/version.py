from __future__ import absolute_import
import bioat
import datetime
import os

project_toml: str = os.path.join(bioat.__path__[0], '../../pyproject.toml')
sec = os.path.getmtime(project_toml)
__upgrade_date__ = datetime.date.fromtimestamp(sec)

with open(project_toml, 'rt') as f:
    __version__ = [i.rstrip().split('= ')[-1]
                   for i in f.readlines() if 'version =' in i][0].replace('"', '')

__author__ = 'Hua-nan ZHAO @ Tsinghua University'
__email__ = 'hermanzhaozzzz@gmail.com'
__doc_format__ = "restructuredtext"
