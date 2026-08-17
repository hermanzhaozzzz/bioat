import sys
from datetime import datetime
from pathlib import Path

import tomllib

# 将项目路径添加到 sys.path, 以便 Sphinx 能找到代码
project_root = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(project_root / "src"))

# 解析 pyproject.toml 文件
pyproject_path = project_root / "pyproject.toml"
if not pyproject_path.exists():
    msg = f"pyproject.toml not found at {pyproject_path}"
    raise FileNotFoundError(msg)

with pyproject_path.open("rb") as f:
    pyproject_data = tomllib.load(f)

# 从 pyproject.toml 中提取元信息
project_metadata = pyproject_data.get("project", {})
project = project_metadata.get("name", "Unknown Project")
author = ", ".join(
    item.get("name", "Unknown Author") for item in project_metadata.get("authors", [{}])
)
release = project_metadata.get("version", "0.1.0")
description = project_metadata.get("description", "Project Description")
year = datetime.now().year

# Sphinx 配置
extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon",
    "sphinx.ext.viewcode",
    "myst_parser",
]

# 模板路径和主题设置
templates_path = ["_templates"]
html_theme = "alabaster"

# 支持的文档格式
source_suffix = [".rst", ".md"]

# 项目语言
language = "en"

# 静态文件路径
html_static_path = ["_static"]
static_path = Path(__file__).parent / "_static"
static_path.mkdir(exist_ok=True)

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
version = release.rsplit(".", 1)[0]  # 主版本号

# 输出设置
master_doc = "index"
html_title = f"{project} Documentation"
html_short_title = project
pygments_style = "sphinx"
highlight_language = "python"
