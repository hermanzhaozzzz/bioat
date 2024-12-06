import os
import sys
import toml  # 用于解析 pyproject.toml

# 将项目路径添加到 sys.path，以便 Sphinx 能找到代码
sys.path.insert(0, os.path.abspath('../src'))

# 解析 pyproject.toml 文件
pyproject_path = os.path.abspath("../pyproject.toml")
with open(pyproject_path, "r") as f:
    pyproject_data = toml.load(f)

# 从 pyproject.toml 中提取元信息
project_metadata = pyproject_data.get("tool", {}).get("poetry", {})
project = project_metadata.get("name", "Unknown Project")
author = project_metadata.get("authors", ["Unknown Author"])[0]
release = project_metadata.get("version", "0.1.0")
description = project_metadata.get("description", "Project Description")

# Sphinx 配置
extensions = [
    'sphinx.ext.autodoc',        # 自动生成文档
    'sphinx.ext.napoleon',       # 支持 Google 和 NumPy 风格的 docstring
    'sphinx.ext.viewcode',       # 为代码添加链接
    'myst_parser',               # 支持 Markdown 文件
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

# Napoleon 配置（针对 Google 风格 docstring）
napoleon_google_docstring = True  # 启用 Google 风格
napoleon_numpy_docstring = False  # 关闭 NumPy 风格，避免冲突
napoleon_include_init_with_doc = True  # 初始化方法包含在文档中
napoleon_include_private_with_doc = False  # 排除私有方法
napoleon_include_special_with_doc = True  # 包括特殊方法（如 __str__）
napoleon_use_ivar = True  # 使用 `ivar` 标签支持类变量
napoleon_use_param = True  # 使用 `param` 标签解析参数
napoleon_use_rtype = True  # 使用 `rtype` 标签解析返回值类型

# 项目信息自动插入
rst_epilog = f"""
.. |project_name| replace:: {project}
.. |project_description| replace:: {description}
"""

# 输出设置
html_title = f"{project} Documentation"
html_short_title = project