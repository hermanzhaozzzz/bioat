from importlib.metadata import PackageNotFoundError, metadata
from pathlib import Path

try:
    import tomllib
except ModuleNotFoundError:  # pragma: no cover - Python < 3.11 fallback
    tomllib = None

__PKG_NAME__ = "bioat"  # !don't change this line unless you know what you are doing


def _load_pyproject_fallback() -> dict[str, str]:
    pyproject = Path(__file__).resolve().parents[2] / "pyproject.toml"
    fallback = {
        "Author": "Huanan Herman Zhao",
        "Author-email": "hermanzhaozzzz@gmail.com",
        "Summary": (
            "bioat, a python package & command line toolkit for Bioinformatics "
            "and data science!"
        ),
        "License": "Apache-2.0",
        "Version": "0.0.0",
    }

    if tomllib is None or not pyproject.exists():
        return fallback

    with pyproject.open("rb") as handle:
        project = tomllib.load(handle).get("project", {})

    author = next(iter(project.get("authors", [])), {})
    license_info = project.get("license", {})
    fallback.update(
        {
            "Author": author.get("name", fallback["Author"]),
            "Author-email": author.get("email", fallback["Author-email"]),
            "Summary": project.get("description", fallback["Summary"]),
            "License": license_info.get("text", fallback["License"]),
            "Version": project.get("version", fallback["Version"]),
        }
    )
    return fallback


try:
    __META__ = metadata(__PKG_NAME__)
except PackageNotFoundError:
    __META__ = _load_pyproject_fallback()

__AUTHOR__ = __META__["Author"]
__AUTHOR_EMAIL__ = __META__["Author-email"]
__DESCRIPTION__ = __META__["Summary"]
__DOC_FORMAT__ = "google"
__DOC_PAGE__ = f"https://{__PKG_NAME__}.readthedocs.io/en/latest/"
__HOME_PAGE__ = f"https://github.com/hermanzhaozzzz/{__PKG_NAME__}"
__ISSUE_PAGE__ = f"{__HOME_PAGE__}/issues"
__LICENSE__ = __META__["License"]
__VERSION__ = __META__["Version"]


if __name__ == "__main__":
    print(f"__PKG_NAME__ = {__PKG_NAME__}")
    print(f"__AUTHOR__ = {__AUTHOR__}")
    print(f"__AUTHOR_EMAIL__ = {__AUTHOR_EMAIL__}")
    print(f"__DESCRIPTION__ = {__DESCRIPTION__}")
    print(f"__DOC_FORMAT__ = {__DOC_FORMAT__}")
    print(f"__DOC_PAGE__ = {__DOC_PAGE__}")
    print(f"__HOME_PAGE__ = {__HOME_PAGE__}")
    print(f"__ISSUE_PAGE__ = {__ISSUE_PAGE__}")
    print(f"__LICENSE__ = {__LICENSE__}")
    print(f"__VERSION__ = {__VERSION__}")
