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
# 项目信息从 pyproject.toml 中解析
PYPROJECT_FILE="pyproject.toml"
PACKAGE_NAME=$(grep -E '^name = "' $PYPROJECT_FILE | sed 's/name = "//;s/"//')
VERSION=$(grep -E '^version = "' $PYPROJECT_FILE | sed 's/version = "//;s/"//')
GITHUB_REPO_URL=$(grep -E '^repository = "' $PYPROJECT_FILE | sed 's/repository = "//;s/"//')
LICENSE_FILE="LICENSE"
CONDA_RECIPE_DIR="conda-recipe"
FORKED_REPO_URL="https://github.com/hermanzhaozzzz/staged-recipes"  # 你的 fork 仓库地址

# 获取 GitHub 用户信息
MAINTAINER=$(grep -E '^authors = \["' $PYPROJECT_FILE | sed 's/authors = \["//;s/".*//')

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

# 创建 Conda 配方文件
function create_conda_recipe {
  echo "Creating Conda recipe..."
  rm -rf $CONDA_RECIPE_DIR
  mkdir -p $CONDA_RECIPE_DIR

  # 确保 LICENSE 文件存在
  if [ ! -f $LICENSE_FILE ]; then
    echo "Error: LICENSE file not found. Please add a LICENSE file."
    exit 1
  fi

  # 调用 Python 脚本生成 Conda 配方
  python3 scripts/toml2yaml.py $PYPROJECT_FILE $CONDA_RECIPE_DIR
}

# 提交到 Conda Forge 或更新现有 PR
function submit_to_conda_forge {
  echo "Submitting to Conda Forge..."
  
  # 从conda forge/staged-recipes 仓库克隆
  if [ ! -d "staged-recipes" ]; then
    git clone https://github.com/conda-forge/staged-recipes.git
  else
    cd staged-recipes
    git pull origin main  # 使用main分支代替master
  fi

  # 检查是否有 myfork 远程仓库存在于我的 fork 仓库中
  if ! git remote | grep -q myfork; then
    git remote add myfork $FORKED_REPO_URL
  fi

  git checkout add-$PACKAGE_NAME-$VERSION || git checkout -b add-$PACKAGE_NAME-$VERSION
  
  mkdir -p recipes/$PACKAGE_NAME
  cp ../$CONDA_RECIPE_DIR/meta.yaml recipes/$PACKAGE_NAME/
  
  git add .
  git commit -m "Update recipe for $PACKAGE_NAME v$VERSION"
  git push myfork add-$PACKAGE_NAME-$VERSION
  # 检查是否已存在的 PR
  existing_pr=$(gh pr list --base main --head myfork:add-$PACKAGE_NAME-$VERSION --json number --jq '.[0].number')

  if [ -z "$existing_pr" ]; then
    # 如果没有现有的 PR，则创建新的 PR
    gh pr create --title "Add $PACKAGE_NAME v$VERSION" --body "This PR adds the Conda recipe for $PACKAGE_NAME version $VERSION."
    echo "New pull request created for Conda Forge!"
  else
    # 如果已有 PR，推送更新
    echo "Pull request #$existing_pr already exists. Pushing updates..."
    git push myfork add-$PACKAGE_NAME-$VERSION
  fi

  cd ..
}

# 主函数，发布到 PyPI 和 Conda Forge
function main {
  # 检查是否安装了 poetry
  if ! command -v poetry &> /dev/null; then
    echo "Error: poetry is not installed. Please install poetry first."
    exit 1
  fi
  # 检查是否安装了 gh (GitHub CLI)
  if ! command -v gh &> /dev/null; then
    echo "Error: GitHub CLI (gh) is not installed. Please install gh first."
    exit 1
  fi

  # 检查 PyPI 版本是否存在
  check_pypi_version
  if [ $? -eq 0 ]; then
    # 如果 PyPI 上没有相同版本，则发布到 PyPI
    publish_to_pypi
  fi

  # 生成 Conda 配方并提交到 Conda Forge
  create_conda_recipe
  submit_to_conda_forge
}

main