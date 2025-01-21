from importlib.metadata import metadata

__PKG_NAME__ = "bioat"  # !don't change this line unless you know what you are doing
__META__ = metadata(__PKG_NAME__)
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
