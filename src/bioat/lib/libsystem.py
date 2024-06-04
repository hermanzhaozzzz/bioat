import os
import sys

from bioat.logger import get_logger

__module_name__ = "bioat.lib.libsystem"


def check_cmd(x):
    """
    Check if a command is available in system PATH var.
    """
    logger = get_logger(
        level="debug", module_name=__module_name__, func_name="check_cmd"
    )
    logger.info("Checking command '%s'", x)
    return any(
        os.access(os.path.join(path, x), os.X_OK)
        for path in os.environ["PATH"].split(os.pathsep)
    )
