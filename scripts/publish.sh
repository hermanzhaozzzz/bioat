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
# --------------------------
# 项目信息
PACKAGE_NAME="bioat"
GITHUB_REPO_URL="https://github.com/hermanzhaozzzz/bioat"
VERSION=$(poetry version --short)
CONDA_RECIPE_DIR="conda-recipe"
FORKED_REPO_URL="https://github.com/hermanzhaozzzz/staged-recipes"  # 你的 fork 仓库地址

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
  if [ $? -ne 0 ]; then
    echo "Error: Failed to publish to PyPI."
    exit 1
  fi
  echo "Successfully published to PyPI!"
}

# 创建 Conda 配方文件
function create_conda_recipe {
  echo "Creating Conda recipe..."
  mkdir -p $CONDA_RECIPE_DIR
  cat > $CONDA_RECIPE_DIR/meta.yaml << EOL
package:
  name: $PACKAGE_NAME
  version: "$VERSION"

source:
  git_url: $GITHUB_REPO_URL
  git_rev: v$VERSION

build:
  number: 0
  script: "{{ PYTHON }} -m pip install ."

requirements:
  build:
    - python
    - pip

  run:
    - python >=3.7

test:
  commands:
    - $PACKAGE_NAME --help

about:
  home: $GITHUB_REPO_URL
  license: MIT
  summary: "A bioinformatics analysis toolkit"
EOL
  echo "Conda recipe created!"
}

# 提交到 Conda Forge
function submit_to_conda_forge {
  echo "Submitting to Conda Forge..."

  # 克隆 Conda Forge 的 staged-recipes 仓库
  if [ ! -d "staged-recipes" ]; then
    git clone https://github.com/conda-forge/staged-recipes.git
  fi
  cd staged-recipes

  # 检查是否有 myfork 远程仓库
  if ! git remote | grep -q myfork; then
    git remote add myfork $FORKED_REPO_URL
  fi

  # 创建一个分支并复制 Conda recipe
  git checkout -b add-$PACKAGE_NAME-$VERSION
  mkdir -p recipes/$PACKAGE_NAME
  cp ../$CONDA_RECIPE_DIR/meta.yaml recipes/$PACKAGE_NAME/

  # 提交并推送到自己的 GitHub 仓库
  git add .
  git commit -m "Add recipe for $PACKAGE_NAME v$VERSION"
  git push myfork add-$PACKAGE_NAME-$VERSION

  # 使用 GitHub CLI 创建 Pull Request
  gh pr create --title "Add $PACKAGE_NAME v$VERSION" --body "This PR adds the Conda recipe for $PACKAGE_NAME version $VERSION."
  cd ..
  echo "Pull request created for Conda Forge!"
}

# 主函数，发布到 PyPI 和 Conda Forge
function main {
  check_pypi_version
  if [ $? -eq 0 ]; then
    publish_to_pypi
  fi

  create_conda_recipe
  submit_to_conda_forge
}

main