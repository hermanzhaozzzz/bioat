import os


HOME = os.path.expanduser("~")


def exists_in_PATH(x) -> bool:
    """Check is there a `x` command in the PATH var, return True or False"""
    return any(
        os.access(os.path.join(path, x), os.X_OK)
        for path in os.environ["PATH"].split(os.pathsep)
    )
