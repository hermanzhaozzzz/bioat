import pandas as pd
import sys
from bioat import get_logger

__all__ = ['set_option']
__module_name__ = 'bioat.lib.libpandas'


def set_option(max_colwidth: int = 40, display_width: int = 120, display_max_columns: int = None,
               display_max_rows: int = 50, log_level='INFO'):
    logger = get_logger(level=log_level, module_name=__module_name__, func_name=sys._getframe().f_code.co_name)

    logger.info(f'set pandas: max_colwidth={max_colwidth}')
    pd.set_option("max_colwidth", max_colwidth)  # column最大宽度
    logger.info(f'set pandas: display.width={display_width}')
    pd.set_option("display.width", display_width)  # dataframe宽度
    logger.info(f'set pandas: display.max_columns={display_max_columns}')
    pd.set_option("display.max_columns", display_max_columns)  # column最大显示数
    logger.info(f'set pandas: display.max_rows={display_max_rows}')
    pd.set_option("display.max_rows", display_max_rows)  # row最大显示数
    return None


if __name__ == '__main__':
    set_option()