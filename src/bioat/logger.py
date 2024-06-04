"""Module for logging

author: Herman Huanan Zhao
email: hermanzhaozzzz@gmail.com
homepage: https://github.com/hermanzhaozzzz

Logging format and color

example 1: # TODO
    bioat list
        <in shell>:
            $ bioat list
        <in python consolo>:
            >>> from bioat.cli import Cli
            >>> bioat = Cli()
            >>> bioat.list()
            >>> print(bioat.list())

example 2:
    _example_
"""

import logging
import logging.handlers
import sys
from logging import CRITICAL, DEBUG, ERROR, INFO, NOTSET, WARNING, Logger

import coloredlogs

# !dont move this -> __PKG_NAME__, referenced by many functions
__PKG_NAME__ = "bioat"

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

    :param level: logger level, can be set to CRITICAL, ERROR, WARNING, INFO, DEBUG, NOTSET, defaults to "ERROR"
    :type level: str, optional
    :param module_name: _description_, defaults to __PKG_NAME__
    :type module_name: _type_, optional
    :param func_name: _description_, defaults to None
    :type func_name: _type_, optional
    :return: a logger object
    :rtype: Logger
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
