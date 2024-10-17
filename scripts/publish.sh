#!/bin/bash
# --------------------------
# 提交新版本
# chmod +x scripts/publish.sh
# ./scripts/publish.sh
# --------------------------
# [pypi + poetry]设置密钥和授权
# 确保 PyPI Token 被正确设置。
# poetry config pypi-token.pypi pypi-somekeys
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
AUTHOR=$(echo "$FORKED_REPO_URL" | cut -d'/' -f4)

# 获取 GitHub 用户信息
MAINTAINER=$(grep -E '^authors = \["' $PYPROJECT_FILE | sed 's/authors = \["//;s/".*//')

# 检查版本是否已经存在于 PyPI
function check_pypi_version {
    echo "检查版本 $VERSION 是否已经存在于 PyPI"
    response=$(curl -s -o /dev/null -w "%{http_code}" https://pypi.org/pypi/$PACKAGE_NAME/$VERSION/json)

    if [ "$response" == "200" ]; then
        echo "版本 $VERSION 已经存在于 PyPI, 跳过PyPI发布."
        return 1
    else
        echo "版本 $VERSION 正在进行PyPI发布..."
        return 0
    fi
}

# 发布到 PyPI
function publish_to_pypi {
    echo "发布至PyPI..."
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
    echo "成功发布版本 $VERSION 至PyPI!"
}

# 创建 Conda 配方文件
function create_conda_recipe {
    echo "创建 Conda 配方(recipe)文件..."
    rm -rf $CONDA_RECIPE_DIR
    mkdir -p $CONDA_RECIPE_DIR

    # 确保 LICENSE 文件存在
    if [ ! -f $LICENSE_FILE ]; then
        echo "Error: LICENSE file not found. Please add a LICENSE file."
        exit 1
    fi

    # # 调用 Python 脚本生成 Conda 配方
    # python3 scripts/toml2yaml.py $PYPROJECT_FILE $CONDA_RECIPE_DIR
    grayskull pypi bioat --output staged-recipes/recipes
}

# 提交到 Conda Forge 或更新现有 PR
function submit_to_conda_forge {
    echo "向 conda-forge提交中..."
    
    # 从conda forge/staged-recipes 仓库克隆
    if [ ! -d "staged-recipes" ]; then
        echo "本地缺少staged-recipes仓库,重新从conda-forge官方仓库克隆到本地..."
        git clone https://github.com/conda-forge/staged-recipes.git
    else
        echo "本地存在staged-recipes仓库,尝试拉取conda-forge官方仓库更新..."
        cd staged-recipes
        echo 
        echo "=========== git ===========↓↓↓↓↓↓↓↓"
        git pull origin main --rebase
        echo "=========== git ===========↑↑↑↑↑↑↑↑"
        echo 
    fi

    # 检查是否有 myfork 远程仓库存在于我的 fork 仓库中
    if ! git remote | grep -q myfork; then
        echo "对本地的staged-recipes仓库,添加我fork的远程仓库,即$FORKED_REPO_URL,并命名为 myfork 远程仓库..."
        git remote add myfork $FORKED_REPO_URL
    fi

    echo "签出分支 add-$PACKAGE_NAME-$VERSION"
    echo 
    echo "=========== git ===========↓↓↓↓↓↓↓↓"
    git checkout add-$PACKAGE_NAME-$VERSION || git checkout -b add-$PACKAGE_NAME-$VERSION
    echo "=========== git ===========↑↑↑↑↑↑↑↑"
    echo 
    echo "复制配方文件到 recipes/$PACKAGE_NAME 文件夹下"
    mkdir -p recipes/$PACKAGE_NAME
    cp ../$CONDA_RECIPE_DIR/meta.yaml recipes/$PACKAGE_NAME/
    
    echo "现在在位置 $PWD/staged-recipes 下"
    echo "将本地分支 add-$PACKAGE_NAME-$VERSION 中配方文件的更改提交到我的 myfork 远程仓库中,即$FORKED_REPO_URL"

    echo 
    echo "=========== git ===========↓↓↓↓↓↓↓↓"
    git add .
    # git commit -m "Update recipe for $PACKAGE_NAME v$VERSION"
    # 通过添加 --allow-empty 标志来强制 Git 提交，即使没有实际内容的变更，也会进行一次提交。这在自动化脚本中是非常有用的方式。
    git commit --allow-empty -m "Update recipe for $PACKAGE_NAME v$VERSION @`date`"
    # 如果你不介意覆盖历史提交，你可以在 git push 时使用 -f 标志强制推送更新，确保每次提交都能强制推送
    git push myfork add-$PACKAGE_NAME-$VERSION -f
    # git push myfork add-$PACKAGE_NAME-$VERSION

    echo "=========== git ===========↑↑↑↑↑↑↑↑"
    echo 
    echo "推送成功"

    # 检查是否已存在的 PR
    TITLE="Add recipe for $PACKAGE_NAME v$VERSION"
    BRANCH_NAME="add-$PACKAGE_NAME-$VERSION"
    echo "检查在 conda-forge/staged-recipes 仓库中是否已存在名为 [ $TITLE ] 的 Pull Request(PR)..."

    # 使用 gh pr list 来检查是否已经有对应的 PR
    existing_pr=$(gh pr list --base main -A $AUTHOR --json number --jq '.[0].number')

    if [ -z "$existing_pr" ]; then
        # 如果没有现有的 PR，则创建新的 PR
        echo "未发现现有的 PR, 创建新的 PR, 即标题为 [ $TITLE ] 的PR..."
        pr_url=$(gh pr create --title "$TITLE" --body "This PR adds the Conda recipe for $PACKAGE_NAME version $VERSION.")
        echo "已经在 conda-forge/staged-recipes 仓库中创建了新的 PR, 名为 [ $TITLE ]."
        echo "PR URL: $pr_url"
    else
        # 如果已有 PR，推送更新
        echo "发现现有的 PR #$existing_pr, 直接推送更新至我的 myfork 远程仓库,并请求 PR 自动工作流..."
        git push myfork "$BRANCH_NAME"
        pr_url="https://github.com/conda-forge/staged-recipes/pull/$existing_pr"
        echo "现有的 PR URL (点击查看): $pr_url"
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
    # submit_to_conda_forge
}

main
