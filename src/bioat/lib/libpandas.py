"""_summary_.

author: Herman Huanan Zhao
email: hermanzhaozzzz@gmail.com
homepage: https://github.com/hermanzhaozzzz

_description_

example 1:
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

import pandas as pd

from bioat.logger import LoggerManager

lm = LoggerManager(mod_name="bioat.lib.libpandas")

__all__ = ["set_option"]


def set_option(
    progres_bar: bool = True,
    max_colwidth: int = 40,
    display_width: int = 120,
    display_max_columns: int | None = None,
    display_max_rows: int = 50,
    log_level="INFO",
):
    lm.set_names(func_name="set_option")
    lm.set_level(log_level)

    # 通过判断是否在notebook环境中，来选择不同的进度条
    try:
        from IPython.display import display

        lm.logger.info("set pandas from a notebook environment")
        from tqdm.notebook import tqdm
    except ImportError:
        lm.logger.info("set pandas from a terminal environment")
        from tqdm import tqdm

    if progres_bar:
        lm.logger.info("set pandas: tqdm.pandas")
        tqdm.pandas()
    lm.logger.info(f"set pandas: max_colwidth={max_colwidth}")
    pd.set_option("max_colwidth", max_colwidth)  # column最大宽度
    lm.logger.info(f"set pandas: display.width={display_width}")
    pd.set_option("display.width", display_width)  # dataframe宽度
    lm.logger.info(f"set pandas: display.max_columns={display_max_columns}")
    pd.set_option("display.max_columns", display_max_columns)  # column最大显示数
    lm.logger.info(f"set pandas: display.max_rows={display_max_rows}")
    pd.set_option("display.max_rows", display_max_rows)  # row最大显示数


def move_column(
    df: pd.DataFrame,
    to_move: str,
    after: str | None = None,
    before: str | None = None,
    first: bool | None = None,
    last: bool | None = None,
    inplace=False,
) -> pd.DataFrame | None:
    """Move a specified column in a DataFrame to a new position.

    :param df: The DataFrame to operate on.
    :param to_move: The name of the column to move.
    :param after: (Optional) If provided, move the `to_move` column after the `after` column.
    :param before: (Optional) If provided, move the `to_move` column before the `before` column.
    :param first: (Optional) If True, move the column to the first position.
    :param last: (Optional) If True, move the column to the last position.
    :param inplace: (Optional) If True, modify the original DataFrame without returning a new one.
    :return: A new DataFrame with the moved column, or None if inplace=True.
    """
    # Validate inputs for to_move
    if to_move not in df.columns:
        msg = f"Column '{to_move}' not found in the DataFrame."
        raise ValueError(msg)
    # Validate inputs for first and last
    if first and last:
        msg = "Cannot specify both 'first' and 'last' as True. Please choose one."
        raise ValueError(
            msg
        )
    if first is not None and not isinstance(first, bool):
        msg = "'first' must be a boolean or None."
        raise ValueError(msg)
    if last is not None and not isinstance(last, bool):
        msg = "'last' must be a boolean or None."
        raise ValueError(msg)
    # Validate inputs for before and after
    if after is not None:
        if after not in df.columns:
            msg = f"Column '{after}' not found in the DataFrame. Cannot use it as a reference."
            raise ValueError(
                msg
            )
        if to_move == after:
            msg = f"to_move({to_move}) is not allowed to be the same with after({after})"
            raise ValueError(
                msg
            )
    if before is not None:
        if before not in df.columns:
            msg = f"Column '{before}' not found in the DataFrame. Cannot use it as a reference."
            raise ValueError(
                msg
            )
        if to_move == before:
            msg = f"to_move({to_move}) is not allowed to be the same with before({before})"
            raise ValueError(
                msg
            )
    if after is not None and before is not None:
        msg = "Cannot specify both 'after' and 'before'. Please choose one."
        raise ValueError(msg)

    current_columns = df.columns.tolist()
    if first:
        new_columns = [to_move] + [col for col in current_columns if col != to_move]
    elif last:
        new_columns = [col for col in current_columns if col != to_move] + [to_move]
    else:
        # Find the index of the column to move
        index_to_move = current_columns.index(to_move)
        # Remove the column from its current position
        current_columns.pop(index_to_move)
        # Determine the new position
        if after is not None:
            index_after = current_columns.index(after) + 1
            new_columns = (
                current_columns[:index_after]
                + [to_move]
                + current_columns[index_after:]
            )
        elif before is not None:
            index_before = current_columns.index(before)
            new_columns = (
                current_columns[:index_before]
                + [to_move]
                + current_columns[index_before:]
            )
        else:
            msg = "you must define one of the params(first/last/before/after) for operation"
            raise ValueError(
                msg
            )

    # Return the modified DataFrame or None if inplace=True
    if inplace:
        df[current_columns] = df[new_columns]
        return None
    return df[new_columns].copy()


if __name__ == "__main__":
    set_option()

    data = pd.DataFrame(
        {
            "a": [1, 2, 3, 4],
            "b": [2, 2, 3, 4],
            "c": [3, 2, 3, 4],
            "d": [4, 2, 3, 4],
        },
    )
    # print(move_column(data, to_move='d'))  # ValueError: you must define one of the params(first/last/before/after)
    # for operation
    # print(move_column(data, to_move='a', before='a'))  # to_move(a) is not allowed to be the same with before(a)
    # print(move_column(data, to_move='a', before='b'))
    # print(move_column(data, to_move='a', before='c'))
    # print(move_column(data, to_move='a', before='d'))
    # print(move_column(data, to_move='d', after='a'))
    # print(move_column(data, to_move='d', after='b'))
    # print(move_column(data, to_move='d', after='c'))
    # print(move_column(data, to_move='d', after='d'))  # to_move(d) is not allowed to be the same with after(d)
    # print(move_column(data, to_move='d', before='a'))
    # print(move_column(data, to_move='d', before='b'))
    # print(move_column(data, to_move='d', before='c'))
    # print(move_column(data, to_move='d', before='d'))  # to_move(d) is not allowed to be the same with before(d)
    # print(move_column(data, 'b', first=True))
    # print(move_column(data, to_move='b', last=True))
    print(f"data = \n{data}, \nid = {id(data)}")
    # data_new = move_column(data, to_move='a', last=True, inplace=False)
    data_new = move_column(data, to_move="a", last=True, inplace=True)
    print(f"data_new = \n{data_new}, \nid = {id(data_new)}")
    print(f"data = \n{data}, \nid = {id(data)}")
