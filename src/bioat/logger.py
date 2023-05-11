import sys
import logging
from logging import (CRITICAL, ERROR, WARNING, INFO, DEBUG, NOTSET)


def set_logging_level(level='ERROR'):
    level = level.upper()

    dt_level = {
        k: v
        for k, v in zip(
            ['CRITICAL', 'ERROR', 'WARNING', 'INFO', 'DEBUG', 'NOTSET'],
            [CRITICAL, ERROR, WARNING, INFO, DEBUG, NOTSET]
        )
    }
    logging.basicConfig(
        level=dt_level[level],
        format='%(levelname)-5s @ %(asctime)s: %(message)s ',
        datefmt='%Y-%m-%d %H:%M:%S',
        stream=sys.stderr,
        filemode="w"
    )
    # set logger
    lib_name = __name__
    function_name = sys._getframe().f_code.co_name
    logger = logging.getLogger(f'{lib_name}.{function_name} ==> ')

    if level not in ('NOTSET', 'INFO'):
        logger.info(f'set logging level = {level}')


if __name__ == '__main__':
    #  ['CRITICAL', 'ERROR', 'WARNING', 'INFO', 'DEBUG', 'NOTSET'],
    # 一个比一个轻
    set_logging_level('CRITICAL')  # 只打印 CRITICAL
    # set_logging_level('ERROR')  # DEBUG、INFO、WARNING 没打印
    # set_logging_level('WARNING')  # DEBUG、INFO 没打印
    # set_logging_level('INFO')  # DEBUG 没打印
    # set_logging_level('DEBUG')  # 全打印
    # set_logging_level('NOTSET')  # 全打印
    # CRITICAL > ERROR > WARNING > INFO > DEBUG = NOTSET
    # 也就是说,如果设置为WARNING,就只打印WARNING及左边,的CRITICAL和ERROR
    logging.debug('Python debug')
    logging.info('Python info')
    logging.warning('Python warning')
    logging.error('Python Error')
    logging.critical('Python critical')
