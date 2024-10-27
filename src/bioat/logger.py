"""
Module for Logging.

This module provides a logging utility that allows for configurable
logging levels and formats. It uses the `logging` standard library
and enhances it with colorized output for better visibility in
command-line interfaces.

Usage:
    You can create a logger object by calling the `get_logger()`
    function, specifying the desired logging level, module name,
    and optionally a function name.

Examples:
    # Create a logger for the main module
    logger = get_logger("DEBUG")
    logger.debug("Debug message")

    # Create a logger for a specific function
    logger = get_logger("ERROR", func_name="my_function")
    logger.error("Error message")
"""

import logging
import logging.handlers
from logging import CRITICAL, DEBUG, ERROR, INFO, NOTSET, WARNING, Logger

import coloredlogs

from ._meta import __PKG_NAME__

LEVEL = dict(
    zip(
        ["CRITICAL", "ERROR", "WARNING", "INFO", "DEBUG", "NOTSET"],
        [CRITICAL, ERROR, WARNING, INFO, DEBUG, NOTSET],
    )
)


def get_logger(
    level="ERROR",
    module_name=__PKG_NAME__,
    func_name: None | str = None,
) -> Logger:
    """Return a logger object.

    This function creates and returns a logger object configured with the specified logging level,
    module name, and optionally, function name.

    Args:
        level (str): Logger level. Can be set to CRITICAL, ERROR, WARNING, INFO, DEBUG, NOTSET.
            Defaults to "ERROR".
        module_name (str): The name of the module for which the logger is created.
            Defaults to __PKG_NAME__.
        func_name (str, optional): The name of the function for which the logger is created.
            Defaults to None.

    Returns:
        Logger: A logger object configured with the provided parameters.
    """

    level = LEVEL[level.upper()]

    # define the name var in fmt: %(name)
    if func_name:
        logger = logging.getLogger(f"{module_name}.{func_name}")
    else:
        logger = logging.getLogger(f"{module_name}")
    # 设置日志格式
    fmt = "%(asctime)s - [%(name)s] - %(filename)s[line:%(lineno)4d] - %(levelname)+8s: %(message)s"
    # add color feature to logger obj (cli color)
    coloredlogs.DEFAULT_FIELD_STYLES = dict(
        asctime=dict(color="green"),
        name=dict(color="blue"),
        filename=dict(color="magenta"),
        lineno=dict(color="cyan"),
        levelname=dict(color="blue", bold=True),
    )
    coloredlogs.install(fmt=fmt, level=level, logger=logger, reconfigure=True)
    return logger


if __name__ == "__main__":
    logger = get_logger(
        "DEBUG",
        # module_name="bioat.logger",
        # func_name="test_function"
    )
    logger.debug("Python debug")
    logger.info("Python info")
    logger.warning("Python warning")
    logger.error("Python Error")
    logger.critical("Python critical")
