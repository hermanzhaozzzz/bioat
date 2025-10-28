#!/usr/bin/env bash
set -euo pipefail

########################################
#            简化固定配置              #
########################################
# 固定包名/维护者/feedstock：
PACKAGE_NAME="bioat"
MAINTAINERS=("hermanzhaozzzz")
FEEDSTOCK_GIT="git@github.com:hermanzhaozzzz/bioat-feedstock.git"
FEEDSTOCK_DIR="bioat-feedstock"

# 基本文件：
PYPROJECT_FILE="pyproject.toml"
LICENSE_FILE="LICENSE"

########################################
#          工具函数（跨平台）          #
########################################
SED_INPLACE() {
  if [[ "$OSTYPE" == darwin* ]]; then
    sed -i '' "$@"
  else
    sed -i "$@"
  fi
}

calc_sha256() {
  local url="$1"
  if command -v shasum >/dev/null 2>&1; then
    curl -sL "$url" | shasum -a 256 | awk '{print $1}'
  else
    curl -sL "$url" | sha256sum | awk '{print $1}'
  fi
}

rewrite_maintainers_block() {
  local meta_yaml="$1"

  # 删除占位符
  SED_INPLACE '/- *AddYourGitHubIdHere/d' "$meta_yaml"

  # 去掉旧的 recipe-maintainers 块
  awk '
    BEGIN{drop=0}
    {
      if ($0 ~ /^[[:space:]]*recipe-maintainers[[:space:]]*:/) { drop=1; next }
      if (drop==1) {
        if ($0 ~ /^[[:space:]]*-/) { next }
        if ($0 ~ /^[[:space:]]*[A-Za-z0-9_-]+[[:space:]]*:/ || $0 ~ /^$/) { drop=0 }
        if (drop==1) { next }
      }
      print
    }
  ' "$meta_yaml" > "${meta_yaml}.tmp" && mv "${meta_yaml}.tmp" "$meta_yaml"

  # 确保 extra: 存在
  grep -q '^extra:' "$meta_yaml" || printf "\nextra:\n" >> "$meta_yaml"

  # 追加干净块
  SED_INPLACE '/^extra:/a\
\
  recipe-maintainers:
' "$meta_yaml"

  for m in "${MAINTAINERS[@]}"; do
    SED_INPLACE "/^  recipe-maintainers:/a\\
\\
    - ${m}
" "$meta_yaml"
  done
}

########################################
#           依赖检查 & 版本取值         #
########################################
for bin in poetry grayskull git curl; do
  command -v "$bin" >/dev/null 2>&1 || { echo "Error: 需要安装 $bin"; exit 1; }
done

if [[ ! -f "$PYPROJECT_FILE" ]]; then
  echo "Error: 未找到 $PYPROJECT_FILE"
  exit 1
fi
if [[ ! -f "$LICENSE_FILE" ]]; then
  echo "Error: 未找到 $LICENSE_FILE"
  exit 1
fi

VERSION=$(grep -E '^version = "' "$PYPROJECT_FILE" | sed 's/version = "//;s/"//')
if [[ -z "${VERSION:-}" ]]; then
  echo "Error: 未能从 $PYPROJECT_FILE 解析出 version"
  exit 1
fi

########################################
#           步骤 1：发布到 PyPI         #
########################################
echo "检查 PyPI 上是否已有 ${PACKAGE_NAME}==${VERSION} ..."
code=$(curl -s -o /dev/null -w "%{http_code}" "https://pypi.org/pypi/${PACKAGE_NAME}/${VERSION}/json")
if [[ "$code" != "200" ]]; then
  echo "未发现该版本，开始发布至 PyPI ..."
  poetry build
  poetry publish
else
  echo "PyPI 已存在 ${PACKAGE_NAME}==${VERSION}，跳过发布。"
fi

# 轮询 PyPI 可见性（最多 5 分钟）
echo "等待 PyPI 展示 ${PACKAGE_NAME}==${VERSION} ..."
for i in {1..30}; do
  code=$(curl -s -o /dev/null -w "%{http_code}" "https://pypi.org/pypi/${PACKAGE_NAME}/${VERSION}/json")
  if [[ "$code" == "200" ]]; then
    echo "PyPI 已可见。"
    break
  fi
  echo "  未就绪（$i/30），10s 后再试..."
  sleep 10
done

########################################
#        步骤 2：用 grayskull 产配方     #
########################################
echo "用 grayskull 生成配方 ..."
/bin/rm -rf staged-recipes 2>/dev/null || true
git clone https://github.com/conda-forge/staged-recipes.git >/dev/null
grayskull pypi "${PACKAGE_NAME}" --output staged-recipes/recipes

META_YAML="staged-recipes/recipes/${PACKAGE_NAME}/meta.yaml"
[[ -f "$META_YAML" ]] || { echo "Error: 未生成 $META_YAML"; exit 1; }

# 修正 version / sha256（以源码 tarball 为准）
echo "修正 version / sha256 ..."
SED_INPLACE "s/^\\([[:space:]]*version:[[:space:]]*\\).*/\\1\"${VERSION}\"/" "$META_YAML"
TARBALL_URL="https://pypi.io/packages/source/${PACKAGE_NAME:0:1}/${PACKAGE_NAME}/${PACKAGE_NAME}-${VERSION}.tar.gz"
SHA256=$(calc_sha256 "$TARBALL_URL")
SED_INPLACE "s/^\\([[:space:]]*sha256:[[:space:]]*\\).*/\\1\"${SHA256}\"/" "$META_YAML"

# 修补维护者
rewrite_maintainers_block "$META_YAML"

########################################
#   步骤 3：更新本地 feedstock 并推分支  #
########################################
if [[ ! -d "$FEEDSTOCK_DIR" ]]; then
  echo "克隆你的 feedstock fork：$FEEDSTOCK_GIT"
  git clone "$FEEDSTOCK_GIT" "$FEEDSTOCK_DIR"
fi

echo "同步上游 conda-forge/${PACKAGE_NAME}-feedstock ..."
(
  cd "$FEEDSTOCK_DIR"
  git checkout main
  git remote add upstream "https://github.com/conda-forge/${PACKAGE_NAME}-feedstock.git" >/dev/null 2>&1 || true
  git fetch upstream
  # 保守起见使用 reset 确保一致
  git reset --hard upstream/main
)

# 新建以版本命名的工作分支（存在则先删）
(
  cd "$FEEDSTOCK_DIR"
  git branch -D "${VERSION}" >/dev/null 2>&1 || true
  git checkout -b "${VERSION}"
)

# 拷贝 meta.yaml
cp "$META_YAML" "${FEEDSTOCK_DIR}/recipe/meta.yaml"

# 提交并推送
(
  cd "$FEEDSTOCK_DIR"
  git add recipe/meta.yaml
  git commit -m "Update to ${VERSION}"
  git push origin "${VERSION}" -f
  git checkout main
  git branch -D "${VERSION}" >/dev/null 2>&1 || true
)

########################################
#            完成与提示                 #
########################################
echo "✅ 已推送分支：${PACKAGE_NAME}-feedstock:${VERSION}"
echo "下一步：到 GitHub 上创建 PR："
echo "  base：conda-forge/${PACKAGE_NAME}-feedstock 的 main"
echo "  compare：你的 fork ${PACKAGE_NAME}-feedstock 的 ${VERSION}"
echo "创建 PR 后在评论中 @conda-forge-admin, please rerender"
echo "CI 通过后合并即可。"