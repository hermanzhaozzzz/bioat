import logging
import os
import logging.handlers
import coloredlogs
from bioat import BioatParameterFormatError
from logging import (CRITICAL, ERROR, WARNING, INFO, DEBUG, NOTSET)

LEVEL = {
    k: v
    for k, v in zip(
        ['CRITICAL', 'ERROR', 'WARNING', 'INFO', 'DEBUG', 'NOTSET'],
        [CRITICAL, ERROR, WARNING, INFO, DEBUG, NOTSET]
    )
}


def get_logger(level='ERROR', module_name='bioat', func_name=''):
    level = LEVEL[level.upper()]

    # define the name var in fmt: %(name)
    if func_name:
        logger = logging.getLogger(f'{module_name}.{func_name}')
    else:
        logger = logging.getLogger(f'{module_name}')
    # 设置日志格式
    fmt = '%(asctime)s - [%(name)+34s] - %(filename)+14s[line:%(lineno)5d] - %(levelname)+8s: %(message)s'
    # add color feature to logger obj (cli color)
    coloredlogs.DEFAULT_FIELD_STYLES = dict(
        asctime=dict(color='green'),
        name=dict(color='blue'),
        filename=dict(color='magenta'),
        lineno=dict(color='cyan'),
        levelname=dict(color='blue', bold=True),
    )
    coloredlogs.install(fmt=fmt, level=level, logger=logger, reconfigure=True)
    return logger


if __name__ == '__main__':
    logger = get_logger('DEBUG', module_name='bioat.logger', func_name='test_function')
    logger.debug('Python debug')
    logger.info('Python info')
    logger.warning('Python warning')
    logger.error('Python Error')
    logger.critical('Python critical')
