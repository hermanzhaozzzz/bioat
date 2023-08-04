import logging
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

    logger.setLevel(level=dt_level[level])
    if log_mode == 's':
        logger_handler = logging.StreamHandler()
        logger_formatter = logging.Formatter(
            fmt='%(levelname)-5s @ %(asctime)s %(name)s: %(message)s',
            datefmt='%Y-%m-%d %H:%M:%S'
        )
    elif log_mode == 'f':
        logger_handler = logging.FileHandler(log_path)
        logger_formatter = logging.Formatter(
            fmt='%(levelname)-5s @ %(asctime)s %(name)s: %(message)s',
            datefmt='%Y-%m-%d %H:%M:%S'
        )
    else:
        raise KeyError

    logger_handler.setFormatter(logger_formatter)
    logger.addHandler(logger_handler)

    if level not in ('NOTSET', 'INFO'):
        logger.info(f'set log level = {level}')

    return logger


if __name__ == '__main__':
    #  ['CRITICAL', 'ERROR', 'WARNING', 'INFO', 'DEBUG', 'NOTSET'],
    # 一个比一个轻
    # logger = set_logging_level('CRITICAL')  # 只打印 CRITICAL
    # logger = set_logging_level('ERROR')  # DEBUG、INFO、WARNING 没打印
    # logger = set_logging_level('WARNING')  # DEBUG、INFO 没打印
    logger = get_logger('INFO')  # DEBUG 没打印
    # logger = set_logging_level('DEBUG')  # 全打印
    # logger = set_logging_level('NOTSET')  # 全打印
    # CRITICAL > ERROR > WARNING > INFO > DEBUG = NOTSET
    # 也就是说,如果设置为WARNING,就只打印WARNING及左边,的CRITICAL和ERROR
    logger.debug('Python debug')
    logger.info('Python info')
    logger.warning('Python warning')
    logger.error('Python Error')
    logger.critical('Python critical')
