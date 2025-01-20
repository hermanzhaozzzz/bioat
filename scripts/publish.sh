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
# 然后正式提交
# ./scripts/publish.sh
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
    # 从conda forge/staged-recipes 仓库克隆
    if [ ! -d "staged-recipes" ]; then
        echo "本地缺少staged-recipes仓库,重新从conda-forge官方仓库克隆到本地..."
        git clone https://github.com/conda-forge/staged-recipes.git
    else
        echo "本地存在staged-recipes仓库,删除后拉取"
        /bin/rm -rf staged-recipes
        git clone https://github.com/conda-forge/staged-recipes.git
    fi

    echo "开始使用grayskull创建本次更新的 Conda 配方(recipe)文件,见下面的grayskull日志(建议等待120s以上以同步bioat的pypi版本更新状态)..."
    rm -rf $CONDA_RECIPE_DIR
    mkdir -p $CONDA_RECIPE_DIR

    # 确保 LICENSE 文件存在
    if [ ! -f $LICENSE_FILE ]; then
        echo "Error: LICENSE file not found. Please add a LICENSE file."
        exit 1
    fi

    # 生成 Conda 配方 # brew install grayskull
    grayskull pypi bioat --output staged-recipes/recipes

    # 更新 meta.yaml 中的版本和 sha256
    META_YAML_PATH="staged-recipes/recipes/bioat/meta.yaml"
    if [ -f "$META_YAML_PATH" ]; then
        echo "Updating version and sha256 in meta.yaml..."
        # 适配 macOS 和 Linux 的 sed 语法
        sed -i '' "s/version: .*/version: \"$VERSION\"/" $META_YAML_PATH
        SHA256=$(curl -sL https://pypi.io/packages/source/b/bioat/bioat-$VERSION.tar.gz | sha256sum | awk '{ print $1 }')
        sed -i '' "s/sha256: .*/sha256: \"$SHA256\"/" $META_YAML_PATH
    else
        echo "Error: meta.yaml file not found. Please ensure Grayskull generated it correctly."
        exit 1
    fi
}

# 提交到 Conda Forge 或更新现有 PR
function submit_to_conda_forge {
    
    cd staged-recipes
    echo "现在位于 `pwd`"
    # 检查是否有 myfork 远程仓库存在于我的 fork 仓库中
    if ! git remote | grep -q myfork; then
        echo "对本地的staged-recipes仓库,添加我fork的远程仓库,即$FORKED_REPO_URL,并命名为 myfork 远程仓库..."
        git remote add myfork $FORKED_REPO_URL
    fi

    echo "向 conda-forge提交中..."

    echo "签出分支 add-$PACKAGE_NAME-$VERSION"
    echo 
    echo "=========== git ===========↓↓↓↓↓↓↓↓"
    git checkout add-$PACKAGE_NAME-$VERSION || git checkout -b add-$PACKAGE_NAME-$VERSION
    echo "=========== git ===========↑↑↑↑↑↑↑↑"
    echo 
    echo "复制配方文件到 staged-recipes/recipes/$PACKAGE_NAME 文件夹下"
    mkdir -p recipes/$PACKAGE_NAME
    
    echo "现在在位置 `pwd`下"
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

    cd ..
    echo "现在位于 `pwd`"
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
}

# 主函数，发布到 PyPI 和 Conda Forge
function main {
    # 检查是否安装了 poetry
    if ! command -v poetry &> /dev/null; then
        echo "Error: poetry is not installed. Please install poetry first (brew install poetry)."
        exit 1
    fi
    # 检查是否安装了 gh (GitHub CLI)
    if ! command -v gh &> /dev/null; then
        echo "Error: GitHub CLI (gh) is not installed. Please install gh first (brew install gh)."
        exit 1
    fi
    # 检查是否安装了 grayskull
    if ! command -v grayskull &> /dev/null; then
        echo "Error: grayskull is not installed. Please install grayskull first (brew install grayskull)."
        exit 1
    fi

    # 检查 PyPI 版本是否存在
    check_pypi_version
    if [ $? -eq 0 ]; then
        # 如果 PyPI 上没有相同版本，则发布到 PyPI
        publish_to_pypi
    fi
    echo "==============================>>>>>"
    echo "waiting for the PyPI release..."
    echo "==========>>>>>"
    echo "Waiting for 120 seconds."
    echo "==============================>>>>>"

    for ((i=120; i>0; i--)); do
        read -t 1 -n 1 input
        if [ $? -eq 0 ]; then
            echo "Skipped waiting!"
            break
        else
            echo -ne "Waiting: $i seconds remaining. <Press Enter to skip.> \033[0K\r"
        fi
    done

    # 生成 Conda 配方并提交到 Conda Forge
    create_conda_recipe
    # submit_to_conda_forge
    echo "now at:"
    echo `pwd`
    echo
    echo
    echo "==============================>>>>>"
    echo "Conda 配方文件已经生成, 参考下文自行提交新版本."
    echo "==============================>>>>>"
    echo "https://conda-forge.org/docs/maintainer/updating_pkgs/"

    # 检查并获取 bioat-feedstock
    if [ ! -d "bioat-feedstock" ]; then
        echo "~~~~~~~~~~~~~~>>>>>"
        echo "本地缺少 bioat-feedstock 仓库, 重新从 conda-forge 官方仓库[的我的fork仓库]克隆到本地..."
        echo "git clone git@github.com:hermanzhaozzzz/bioat-feedstock.git"
        git clone git@github.com:hermanzhaozzzz/bioat-feedstock.git
        echo "本地存在 bioat-feedstock 仓库, 尝试拉取官方仓库更新..."
        cd bioat-feedstock
        git checkout main
        git remote add upstream https://github.com/conda-forge/bioat-feedstock.git >/dev/null 2>&1
        git fetch upstream
        git rebase upstream/main
        cd ..
    else
        echo "~~~~~~~~~~~~~~>>>>>"
        echo "本地存在 bioat-feedstock 仓库, 尝试拉取官方仓库更新..."
        cd bioat-feedstock
        git checkout main
        git remote add upstream https://github.com/conda-forge/bioat-feedstock.git >/dev/null 2>&1
        git fetch upstream
        git rebase upstream/main
        cd ..
    fi
    echo "~~~~~~~~~~~~~~>>>>>"
    echo "cd bioat-feedstock"
    cd bioat-feedstock
    echo "try to delete old branch $VERSION"
    git branch -D $VERSION
    echo "create new branch $VERSION"
    echo "git checkout -b $VERSION"
    git checkout -b $VERSION
    echo "~~~~~~~~~~~~~~>>>>>"
    # 将 staged-recipes 中的最新 recipe 复制到 bioat-feedstock 的 recipe 目录
    echo "cp ../staged-recipes/recipes/bioat/meta.yaml recipe/meta.yaml"
    cp ../staged-recipes/recipes/bioat/meta.yaml recipe/meta.yaml
    echo "==============================>>>>>"
    echo "git add+commit+push the new meta.yaml to hermanzhaozzzz/bioat-feedstock branch $VERSION"
    echo "==============================>>>>>"
    # 更新内容后提交更改
    echo "git add ."
    git add .
    echo "git commit -m \"Update to $VERSION\""
    git commit -m "Update to $VERSION"
    echo "git push origin $VERSION"
    git push origin $VERSION -f
    echo "${VERSION} 已提交,切换到主分支"
    git checkout main
    echo "删除本地的 $VERSION 分支,因为已经提交到官方仓库"
    git branch -D $VERSION
    echo "~~~~~~~~~~~~~~>>>>>"
    cd ..
    echo "cd `pwd`"
    echo "==============================>>>>>"
    echo "bioat-feedstock 仓库已经更新, 请自行提交 PR 至官方仓库."
    echo "⭐️⭐️⭐️⭐️⭐️⭐️⭐️⭐️⭐️⭐️⭐️⭐️⭐️⭐️⭐️"
    echo "访问 https://github.com/hermanzhaozzzz/bioat-feedstock/pulls 创建PR请求"
    echo "[step1] 点击 New pull request"
    echo "[step2] [❎❎❎注意一般直接就是默认的❎❎❎]选择 conda-forge/bioat-feedstock 的main分支作为 base repository"
    echo "[step3] 选择 hermanzhaozzzz/bioat-feedstock 的$VERSION 分支作为 compare branch"
    echo "[step4] 点击 Create pull request"
    echo "[step5] 自动化测试通过后, 点击 Merge pull request"
    echo "[step6] 在PR下添加一个Comments, 内容 [⭐️ @conda-forge-admin, please rerender ⭐️] 通知 conda-forge 重新渲染"
    echo "[step7] 选择 hermanzhaozzzz/bioat-feedstock 的$VERSION 分支, 点击 Delete branch按钮删除分支, 完成发布."
    echo "⭐️"
    echo "按照Pull request successfully merged and closed提示点击Delete branch按钮删除 hermanzhaozzzz/bioat-feedstock/$VERSION 分支"
    echo "⭐️"
    echo "Conda 包已经发布, 可以半小时后运行 conda search -c conda-forge 自行检查是否成功更新版本."
    echo "⭐️⭐️⭐️⭐️⭐️⭐️⭐️⭐️⭐️⭐️⭐️⭐️⭐️⭐️⭐️"
}

main
