from __future__ import absolute_import
import pandas as pd
import logging
import sys
from bioat.logger import set_logging_level

__all__ = ['set_option']


def set_option(max_colwidth: int = 40, display_width: int = 120, display_max_columns: int = None,
               display_max_rows: int = 50, log_level: str = 'INFO', **kwargs):
    # set logger
    if log_level in ['CRITICAL', 'ERROR', 'WARNING', 'INFO', 'DEBUG', 'NOTSET']:
        set_logging_level(level=log_level)
    else:
        raise KeyError("log_level must be element in ['CRITICAL', 'ERROR', 'WARNING', 'INFO', 'DEBUG', 'NOTSET']")

    lib_name = __name__
    function_name = sys._getframe().f_code.co_name
    logger = logging.getLogger(f'{lib_name}.{function_name} ==> ')
    logger.info(f'set pandas.set_options: max_colwidth={max_colwidth}')
    pd.set_option("max_colwidth", max_colwidth)  # column最大宽度
    logger.info(f'set pandas.set_options: display.width={display_width}')
    pd.set_option("display.width", display_width)  # dataframe宽度
    logger.info(f'set pandas.set_options: display.max_columns={display_max_columns}')
    pd.set_option("display.max_columns", display_max_columns)  # column最大显示数
    logger.info(f'set pandas.set_options: display.max_rows={display_max_rows}')
    pd.set_option("display.max_rows", display_max_rows)  # row最大显示数
    return
