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


def is_path(string):
    return Path(string).exists()


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
    x: str | None, name: str | None, log_level: str = "WARNING"
) -> None:
    if not x:
        if not name:
            raise BioatInvalidParameterError(
                "Either x or name must be provided only one."
            )
        else:
            if not check_cmd(name, log_level):
                raise BioatMissingDependencyError(f"{name} not found in PATH")
    else:
        if not name:
            if not os.path.exists(x) or not os.access(x, os.X_OK):
                raise BioatInvalidOptionError(f"{x} not found or not executable")
        else:
            raise BioatInvalidParameterError(
                "Either x or name must be provided only one."
            )