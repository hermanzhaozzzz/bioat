import os
from pathlib import Path

from bioat.exceptions import (
    BioatInvalidOptionError,
    BioatInvalidParameterError,
    BioatMissingDependencyError,
)
from bioat.logger import LoggerManager

lm = LoggerManager(mod_name="bioat.lib.libpath")

HOME = os.path.expanduser("~")


def check_cmd(x, log_level="WARNING") -> bool:
    """Check if a command is available in the system's PATH.

    Args:
        x (str): The command name to check.

    Returns:
        bool: True if the command is executable and found in PATH, False otherwise.
    """
    lm.set_names(func_name="check_cmd")
    lm.set_level(log_level)

    lm.logger.info("Checking command '%s'", x)
    result = any(
        os.access(os.path.join(path, x), os.X_OK)
        for path in os.environ["PATH"].split(os.pathsep)
    )
    if result:
        lm.logger.info("Command '%s' is available", x)
    else:
        lm.logger.warning("Command '%s' is not available", x)
    return result


def check_executable(
    x: str | None,
    name: str | None,
    log_level: str = "WARNING",
) -> None:
    if not x:
        if not name:
            msg = "Either x or name must be provided only one."
            raise BioatInvalidParameterError(
                msg,
            )
        if not check_cmd(name, log_level):
            msg = f"{name} not found in PATH"
            raise BioatMissingDependencyError(msg)
    elif not name:
        if not os.path.exists(x) or not os.access(x, os.X_OK):
            msg = f"{x} not found or not executable"
            raise BioatInvalidOptionError(msg)
    else:
        msg = "Either x or name must be provided only one."
        raise BioatInvalidParameterError(
            msg,
        )
