# scripts/toml2yaml.py

import os
import sys

import toml


def parse_pyproject(pyproject_path):
    """解析 pyproject.toml 文件"""
    with open(pyproject_path, "r") as f:
        pyproject = toml.load(f)

    # 提取基础信息
    poetry = pyproject.get("tool", {}).get("poetry", {})
    package_info = {
        "name": poetry.get("name", "unknown"),
        "version": poetry.get("version", "0.0.1"),
        "description": poetry.get("description", "No description available."),
        "license": poetry.get("license", "No license provided"),
        "authors": poetry.get("authors", []),
        "maintainers": poetry.get("maintainers", []),
        "repository": poetry.get("repository", ""),
        "homepage": poetry.get("homepage", ""),
        "readme": poetry.get("readme", "README.md"),
        "documentation": poetry.get("documentation", ""),
    }

    # 提取依赖项和开发依赖项
    dependencies = poetry.get("dependencies", {})
    dev_dependencies = poetry.get("group", {}).get("dev", {}).get("dependencies", {})

    # 处理依赖项中的格式问题
    fixed_dependencies = fix_dependency_format(dependencies)

    # 替换 matplotlib 为 matplotlib-base
    if "matplotlib" in fixed_dependencies:
        fixed_dependencies["matplotlib-base"] = fixed_dependencies.pop("matplotlib")

    return package_info, fixed_dependencies, dev_dependencies


def fix_dependency_format(dependencies):
    """修正依赖项中的格式问题，确保符合 conda 要求"""
    need_fix_bioconda = ("pysam", "dna-features-viewer")
    fixed_dependencies = {}
    for package, version in dependencies.items():
        if isinstance(version, str):
            # 替换 "^" 和 "≥" 为 ">="，并移除不必要的空格
            version = version.replace("^", ">=").replace("≥", ">=")
        if package in need_fix_bioconda:
            # 对 pysam 等依赖项特殊处理，确保来源为 bioconda
            fixed_dependencies[package] = {"version": version, "channel": "bioconda"}
        else:
            fixed_dependencies[package] = version
    return fixed_dependencies


def create_meta_yaml(package_info, dependencies, conda_recipe_dir):
    """根据解析的依赖项和项目信息创建 meta.yaml 文件"""
    name = package_info["name"]
    version = package_info["version"]
    description = package_info["description"]
    license_name = package_info["license"]
    repository = package_info["repository"]
    homepage = package_info["homepage"]
    readme = package_info["readme"]
    # maintainers = [
    #     author.split("<")[0].strip() for author in package_info["maintainers"]
    # ]
    # maintainer_str = "\n    - ".join(maintainers)
    maintainer = package_info["repository"].split("/")[-2]
    conda_requirements = []
    channels = set()

    # 处理依赖项并生成格式化的 requirements
    for pkg, info in dependencies.items():
        if isinstance(info, dict) and "channel" in info:
            # 对 bioconda 来源的包做特殊处理
            conda_requirements.append(f"{pkg} {info['version']}  # [{info['channel']}]")
            channels.add(info["channel"])
        else:
            conda_requirements.append(f"{pkg} {info}")

    meta_yaml = f"""
package:
  name: {name}
  version: "{version}"

source:
  git_url: {repository}
  git_rev: v{version}

build:
  noarch: python
  number: 0
  script: "{{{{ PYTHON }}}} -m pip install ."

requirements:
  host:
    - python >=3.10
    - pip
    - setuptools
    - poetry
  run:
    - {"\n    - ".join(conda_requirements)}

test:
  commands:
    - {name} --help

about:
  home: {homepage}
  license: {license_name}
  license_file: LICENSE
  summary: "{description}"

extra:
  recipe-maintainers:
    - {maintainer}

channels:
    - bioconda
    - conda-forge
    - hcc
    - pkgs/main
"""

    # 确保 Conda 配方目录存在
    if not os.path.exists(conda_recipe_dir):
        os.makedirs(conda_recipe_dir)

    # 写入 meta.yaml 文件
    with open(f"{conda_recipe_dir}/meta.yaml", "w") as f:
        f.write(meta_yaml)

    print(f"Conda recipe created at {conda_recipe_dir}/meta.yaml")


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python toml2yaml.py <PYPROJECT_TOML_PATH> <CONDA_RECIPE_DIR>")
        sys.exit(1)

    pyproject_path = sys.argv[1]
    conda_recipe_dir = sys.argv[2]

    # 解析 pyproject.toml 文件
    package_info, dependencies, dev_dependencies = parse_pyproject(pyproject_path)

    # 生成 meta.yaml 配方
    create_meta_yaml(package_info, dependencies, conda_recipe_dir)
