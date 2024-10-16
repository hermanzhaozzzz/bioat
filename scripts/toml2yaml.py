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
        "maintainer": poetry.get("repository", "").split("/")[-2],
        "repository": poetry.get("repository", ""),
        "homepage": poetry.get("homepage", ""),
        "readme": poetry.get("readme", "README.md"),
        "documentation": poetry.get("documentation", ""),
    }

    # 提取依赖项
    dependencies = poetry.get("dependencies", {})

    # 提取可选依赖项（extras）
    extras = poetry.get("extras", {})

    # 处理依赖项中的格式问题
    fixed_dependencies = fix_dependency_format(dependencies)
    fixed_extras = fix_dependency_format_from_extras(extras, dependencies)

    return package_info, fixed_dependencies, fixed_extras


def fix_dependency_format(dependencies):
    """修正依赖项中的格式问题，确保符合 conda 要求，并替换 matplotlib 为 matplotlib-base"""
    fixed_dependencies = {}
    for pkg, ver in dependencies.items():
        if isinstance(ver, dict):
            # 只提取 'ver' 字段
            ver = ver.get("version", "")
        if isinstance(ver, str):
            # 确保没有空格，并替换 "^" 和 "≥" 为 ">="
            ver = ver.replace("^", ">=").replace("≥", ">=").strip()

        # 替换 matplotlib 为 matplotlib-base
        if pkg == "matplotlib":
            pkg = "matplotlib-base"

        fixed_dependencies[pkg] = ver
    return fixed_dependencies


def fix_dependency_format_from_extras(extras, dependencies):
    """将 extras 中的可选依赖项合并到主依赖项中，并修正格式"""
    combined_dependencies = dependencies.copy()  # 保留原始依赖项
    for extra_deps in extras.values():
        for dep in extra_deps:
            if dep in dependencies:
                version = dependencies.get(dep)
                # 仅提取版本信息
                if isinstance(version, dict):
                    version = version.get("version", "")
                    # 替换 matplotlib 为 matplotlib-base
                if dep == "matplotlib":
                    dep = "matplotlib-base"
                combined_dependencies[dep] = version
    return combined_dependencies


def create_meta_yaml(package_info, dependencies, extras_dependencies, conda_recipe_dir):
    """根据解析的依赖项和项目信息创建 meta.yaml 文件"""
    name = package_info["name"]
    version = package_info["version"]
    description = package_info["description"]
    license_name = package_info["license"]
    repository = package_info["repository"]
    homepage = package_info["homepage"]
    maintainer = package_info["repository"].split("/")[-2]

    # 获取 Python 版本信息，默认使用 >=3.10
    python_version = dependencies.pop("python", ">=3.10")

    conda_requirements = dict()

    # 合并主依赖和可选依赖，并生成格式化的 requirements
    all_dependencies = {**dependencies, **extras_dependencies}
    for pkg, ver in all_dependencies.items():
        if pkg == "matplotlib":
            pkg = "matplotlib-base"
        if pkg in ("playwright", "selenium", "beautifulsoup4"):
            continue
        conda_requirements[pkg] = ver

    conda_requirements_ls = []
    for pkg, ver in conda_requirements.items():
        # if pkg == "python":
        #     # If python is a host requirement, it should be a run requirement.
        #     conda_requirements_ls.append(f"    - {pkg}")
        # else:
        #     conda_requirements_ls.append(f"    - {pkg} {ver}")
        conda_requirements_ls.append(f"    - {pkg} {ver}")

    conda_requirements_ls = "\n".join(conda_requirements_ls)

    # 基础 meta.yaml 配置
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
    - python {python_version} # Non noarch packages should have python requirement without any version constraints
    - pip
    - setuptools
    - poetry
  run:
{conda_requirements_ls}

test:
  commands:
    - {name} --help  #  `{name}` is a command-line tool
  imports:
    - {name}

about:
  home: {homepage}
  license: {license_name}
  license_file: LICENSE
  summary: "{description}"

extra:
  recipe-maintainers:
    - {maintainer}
"""

    # 确保 Conda 配方目录存在
    os.makedirs(conda_recipe_dir, exist_ok=True)

    # 写入 meta.yaml 文件
    with open(os.path.join(conda_recipe_dir, "meta.yaml"), "w") as f:
        f.write(meta_yaml)

    print(f"Conda recipe created at {conda_recipe_dir}/meta.yaml")


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python toml2yaml.py <PYPROJECT_TOML_PATH> <CONDA_RECIPE_DIR>")
        sys.exit(1)

    pyproject_path = sys.argv[1]
    conda_recipe_dir = sys.argv[2]

    # 解析 pyproject.toml 文件
    (
        package_info,
        dependencies,
        extras_dependencies,
    ) = parse_pyproject(pyproject_path)

    # 生成 meta.yaml 配方
    create_meta_yaml(
        package_info,
        dependencies,
        extras_dependencies,
        conda_recipe_dir,
    )
