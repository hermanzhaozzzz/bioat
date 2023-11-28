import logging
import os
import logging.handlers
import coloredlogs
from bioat import BioatParameterFormatError
from logging import (CRITICAL, ERROR, WARNING, INFO, DEBUG, NOTSET)


def get_logger(level='ERROR', module_name='bioat', func_name='', log_mode='s', log_path=None):
    level = level.upper()

    dt_level = {
        k: v
        for k, v in zip(
            ['CRITICAL', 'ERROR', 'WARNING', 'INFO', 'DEBUG', 'NOTSET'],
            [CRITICAL, ERROR, WARNING, INFO, DEBUG, NOTSET]
        )
    }

    if func_name:
        logger = logging.getLogger(f'{module_name}.{func_name}')
    else:
        logger = logging.getLogger(f'{module_name}')
    # 设置日志格式
    fmt = '%(asctime)s - [%(name)+34s] - %(filename)+14s[line:%(lineno)5d] - %(levelname)+8s: %(message)s'
    formatter = logging.Formatter(fmt)

    # 创建Handler, 输出到控制台-s 或者文件-f
    if log_mode == 's' and not log_path:
        logger_handler = logging.StreamHandler()
        logger_handler.setFormatter(formatter)
    elif log_mode == 'f' and log_path:
        # 创建日志文件
        if not os.path.isdir(log_path):
            os.mkdir(log_path)
        logger_handler = logging.FileHandler(log_path)
        logger_handler.setFormatter(formatter)
    else:
        raise BioatParameterFormatError

    logger.addHandler(logger_handler)
    logger.setLevel(level=dt_level[level])
    # 当日志输出到控制台时，会带有颜色
    coloredlogs.DEFAULT_FIELD_STYLES = dict(
        asctime=dict(color='green'),
        name=dict(color='blue'),
        filename=dict(color='magenta'),
        lineno=dict(color='cyan'),
        levelname=dict(color='blue', bold=True),
    )
    coloredlogs.install(fmt=fmt, level=level, logger=logger)
    #
    # if level == 'DEBUG':
    #     logger.info(f'set log level = {level}')

    return logger


if __name__ == '__main__':
    logger = get_logger('DEBUG', module_name='bioat.logger', func_name='test_function')
    logger.debug('Python debug')
    logger.info('Python info')
    logger.warning('Python warning')
    logger.error('Python Error')
    logger.critical('Python critical')
