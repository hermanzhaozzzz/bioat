import os
import sys
import toml
from datetime import datetime

# 将项目路径添加到 sys.path，以便 Sphinx 能找到代码
sys.path.insert(0, os.path.abspath('../src'))

# 解析 pyproject.toml 文件
pyproject_path = os.path.abspath("../pyproject.toml")
if not os.path.exists(pyproject_path):
    raise FileNotFoundError(f"pyproject.toml not found at {pyproject_path}")

with open(pyproject_path, "r") as f:
    pyproject_data = toml.load(f)

# 从 pyproject.toml 中提取元信息
project_metadata = pyproject_data.get("tool", {}).get("poetry", {})
project = project_metadata.get("name", "Unknown Project")
author = ", ".join(project_metadata.get("authors", ["Unknown Author"]))
release = project_metadata.get("version", "0.1.0")
description = project_metadata.get("description", "Project Description")
year = datetime.now().year

# Sphinx 配置
extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.napoleon',
    'sphinx_rtd_theme',
    'sphinx.ext.viewcode',
    'myst_parser',
]

# 模板路径和主题设置
templates_path = ['_templates']
html_theme = 'sphinx_rtd_theme'

# 支持的文档格式
source_suffix = ['.rst', '.md']

# 项目语言
language = 'en'

# 静态文件路径
html_static_path = ['_static']
static_path = os.path.join(os.path.dirname(__file__), "_static")
if not os.path.exists(static_path):
    os.makedirs(static_path)

# Napoleon 配置
napoleon_google_docstring = True
napoleon_numpy_docstring = False
napoleon_include_init_with_doc = True
napoleon_include_private_with_doc = False
napoleon_include_special_with_doc = True
napoleon_use_ivar = True
napoleon_use_param = True
napoleon_use_rtype = True

# Myst 配置
myst_enable_extensions = [
    "deflist",
    "dollarmath",
    "colon_fence",
    "fieldlist",
    "substitution",
]
myst_heading_anchors = 3

# 项目信息自动插入
rst_epilog = f"""
.. |project_name| replace:: {project}
.. |project_description| replace:: {description}
"""

# 动态生成版本号
version = release.rsplit('.', 1)[0]  # 主版本号

# 输出设置
html_title = f"{project} Documentation"
html_short_title = project
pygments_style = 'sphinx'
highlight_language = 'python'