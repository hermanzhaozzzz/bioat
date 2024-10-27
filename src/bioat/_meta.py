from importlib.metadata import metadata


def _get_bam_parser_backend():
    backend = ""

    try:
        import pysam

        backend = "pysam"
    except ImportError:
        pass

    if not backend:
        try:
            import bamnostic as pysam

            backend = "bamnostic"
        except ImportError:
            backend = ""
    return pysam, backend


__PKG_NAME__ = "bioat"  # !don't change this line unless you know what you are doing
__META__ = metadata(__PKG_NAME__)
__AUTHOR__ = __META__["Author"]
__AUTHOR_EMAIL__ = __META__["Author-email"]
pysam, __BAM_PARSER_BACKEND__ = _get_bam_parser_backend()
__DESCRIPTION__ = __META__["Summary"]
__DOC_FORMAT__ = "google"
__DOC_PAGE__ = f"https://{__PKG_NAME__}.readthedocs.io/en/latest/"
__HOME_PAGE__ = __META__["Home-page"]
__ISSUE_PAGE__ = f"{__HOME_PAGE__}/issues"
__LICENSE__ = __META__["License"]
__VERSION__ = __META__["Version"]


if __name__ == "__main__":
    print(f"__PKG_NAME__ = {__PKG_NAME__}")
    print(f"__AUTHOR__ = {__AUTHOR__}")
    print(f"__AUTHOR_EMAIL__ = {__AUTHOR_EMAIL__}")
    print(f"__BAM_PARSER_BACKEND__ = {__BAM_PARSER_BACKEND__}")
    print(f"__DESCRIPTION__ = {__DESCRIPTION__}")
    print(f"__DOC_FORMAT__ = {__DOC_FORMAT__}")
    print(f"__DOC_PAGE__ = {__DOC_PAGE__}")
    print(f"__HOME_PAGE__ = {__HOME_PAGE__}")
    print(f"__ISSUE_PAGE__ = {__ISSUE_PAGE__}")
    print(f"__LICENSE__ = {__LICENSE__}")
    print(f"__VERSION__ = {__VERSION__}")
