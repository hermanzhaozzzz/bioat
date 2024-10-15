#!/bin/bash
# --------------------------
# 提交新版本
# chmod +x scripts/publish.sh
# ./scripts/publish.sh
# --------------------------
# [pypi + poetry]设置密钥和授权
# 确保 PyPI Token 被正确设置。
# poetry config pypi-token.pypi pypi-somekeys
# poetry config --list | grep pypi-token
# --------------------------
# [conda forge  github]授权脚本自动创建 Pull Request
# 确保你已登录 GitHub CLI：使用以下命令确保你已经登录到 GitHub CLI，并授权脚本自动创建 Pull Request：
# gh auth login
# cd staged-recipes 
# gh repo set-default conda-forge/staged-recipes
# --------------------------
# 提交前,先tag版本号
# git checkout main  # 或其他主分支
# 如果tag和以前的相同,可以使用 -f覆盖旧标签,慎用,不利于协作开发
# git tag v0.12.13
# git push origin v0.12.13
# --------------------------
# 如果PR提交失败,则需要根据PR Comments的提示
# 对PR进行新的更改, 手动提交PR更新
# cd staged-recipes/recipes/bioat
# git add meta.yaml
# git commit -m "Add python >=3.8 requirement to meta.yaml"
# git push
# push后即进入自动化check流程
# --------------------------
# 初次提交版本时可以直接自动化提交版本
# ./scripts/publish.sh
# --------------------------
# 项目信息
PACKAGE_NAME="bioat"
GITHUB_REPO_URL="https://github.com/hermanzhaozzzz/bioat"
VERSION=$(poetry version --short)
CONDA_RECIPE_DIR="conda-recipe"
FORKED_REPO_URL="https://github.com/hermanzhaozzzz/staged-recipes"  # 你的 fork 仓库地址

# 获取 GitHub 用户信息
MAINTAINER="hermanzhaozzzz"  # 您的 GitHub 用户名

# 检查版本是否已经存在于 PyPI
function check_pypi_version {
  echo "Checking if version $VERSION already exists on PyPI..."
  response=$(curl -s -o /dev/null -w "%{http_code}" https://pypi.org/pypi/$PACKAGE_NAME/$VERSION/json)

  if [ "$response" == "200" ]; then
    echo "Version $VERSION already exists on PyPI, skipping PyPI upload."
    return 1
  else
    echo "Version $VERSION does not exist on PyPI, proceeding with the upload."
    return 0
  fi
}

# 发布到 PyPI
function publish_to_pypi {
  echo "Publishing to PyPI..."
  poetry build
  if [ $? -ne 0 ]; then
    echo "Error: Failed to build the package."
    exit 1
  fi

  poetry publish
  if [ $? -ne 0 ];then
    echo "Error: Failed to publish to PyPI."
    exit 1
  fi
  echo "Successfully published to PyPI!"
}

# 创建 Conda 配方文件并自动修复 linter 错误
function create_conda_recipe {
  echo "Creating Conda recipe..."
  rm -rf $CONDA_RECIPE_DIR
  mkdir -p $CONDA_RECIPE_DIR

  # 确保 LICENSE 文件存在
  if [ ! -f $LICENSE_FILE ];then
    echo "Error: LICENSE file not found. Please add a LICENSE file."
    exit 1
  fi

  # 创建 meta.yaml 文件，包含 setuptool 作为构建后端和 noarch 配置
  cat > $CONDA_RECIPE_DIR/meta.yaml << EOL
package:
  name: $PACKAGE_NAME
  version: "$VERSION"

source:
  git_url: $GITHUB_REPO_URL
  git_rev: v$VERSION

build:
  noarch: python  # 使用 noarch: python 来表示该包与平台无关
  number: 0
  script: "{{ PYTHON }} -m pip install ."

requirements:
  host:
    - python >=3.8  # 在 host 中指定 Python 版本下限
    - pip
    - setuptools  # 明确指定 setuptools 作为构建后端
    - poetry
  run:
    - python >=3.8  # 在 run 中也指定 Python 版本下限

test:
  commands:
    - $PACKAGE_NAME --help

about:
  home: $GITHUB_REPO_URL
  license: AGPL-3.0-only  # 更新为 GNU AGPLv3
  license_file: $LICENSE_FILE
  summary: "bioat, A python package & command line toolkit for Bioinformatics and data science!"

extra:
  recipe-maintainers:
    - $MAINTAINER
EOL

  echo "Conda recipe created and common linting issues fixed!"
}

# 确保设置默认仓库为 myfork
function ensure_default_remote {
  # 提取 OWNER/REPO 格式
  REPO_PATH=conda-forge/staged-recipes
  echo "Ensuring default remote repository is set to $REPO_PATH..."
  gh repo set-default $REPO_PATH || echo "Could not set default repository"
}

# 提交到 Conda Forge 或更新现有 PR
function submit_to_conda_forge {
  ensure_default_remote  # 调用设置默认仓库的函数
  echo "Submitting to Conda Forge..."
  
  # 从conda forge/staged-recipes 仓库克隆
  if [ ! -d "staged-recipes" ]; then
    git clone https://github.com/conda-forge/staged-recipes.git
  else
    cd staged-recipes
    git pull origin main  # 使用main分支代替master
  fi

  # 检查是否有 myfork 远程仓库存在于我的 fork仓库中
  if ! git remote | grep -q myfork; then
    git remote add myfork $FORKED_REPO_URL
  fi

  git checkout add-$PACKAGE_NAME-$VERSION || git checkout -b add-$PACKAGE_NAME-$VERSION
  
  mkdir -p recipes/$PACKAGE_NAME
  cp ../$CONDA_RECIPE_DIR/meta.yaml recipes/$PACKAGE_NAME/
  
  git add .
  git commit -m "Update recipe for $PACKAGE_NAME v$VERSION"
  git push myfork add-$PACKAGE_NAME-$VERSION

  existing_pr=$(gh pr list --base main --head myfork:add-$PACKAGE_NAME-$VERSION --json number --jq '.[0].number')

  if [ -z "$existing_pr" ]; then
    gh pr create --title "Add $PACKAGE_NAME v$VERSION" --body "This PR adds the Conda recipe for $PACKAGE_NAME version $VERSION."
    echo "New pull request created for Conda Forge!"
  else
    echo "Pull request #$existing_pr already exists. Pushing updates..."
  fi

  cd ..
}

# 主函数，发布到 PyPI 和 Conda Forge
function main {
  if ! command -v poetry &> /dev/null; then
    echo "Error: poetry is not installed. Please install poetry first."
    exit 1
  fi
  if ! command -v gh &> /dev/null; then
    echo "Error: GitHub CLI (gh) is not installed. Please install gh first."
    exit 1
  fi

  check_pypi_version
  if [ $? -eq 0 ]; then
    publish_to_pypi
  fi

  create_conda_recipe
  submit_to_conda_forge
}

main