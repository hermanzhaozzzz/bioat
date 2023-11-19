import numpy as np
import sys
import math
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from bioat import get_logger

__module_name__ = 'bioat.lib.libcolor'

def plot_colortable(colors, *, ncols=4, sort_colors=True, labels=None):

    cell_width = 212
    cell_height = 22
    swatch_width = 48
    margin = 12

    names = list(colors)

    if not labels:
        labels = names

    n = len(names)
    nrows = math.ceil(n / ncols)

    width = cell_width * 4 + 2 * margin
    height = cell_height * nrows + 2 * margin
    dpi = 72

    fig, ax = plt.subplots(figsize=(width / dpi, height / dpi), dpi=dpi)
    fig.subplots_adjust(margin / width, margin / height,
                        (width - margin) / width, (height - margin) / height)
    ax.set_xlim(0, cell_width * 4)
    ax.set_ylim(cell_height * (nrows - 0.5), -cell_height / 2.)
    ax.yaxis.set_visible(False)
    ax.xaxis.set_visible(False)
    ax.set_axis_off()

    for i, name in enumerate(names):
        row = i % nrows
        col = i // nrows
        y = row * cell_height

        swatch_start_x = cell_width * col
        text_pos_x = cell_width * col + swatch_width + 7

        ax.text(text_pos_x, y, labels[i], fontsize=14,
                horizontalalignment='left',
                verticalalignment='center')

        ax.add_patch(
            Rectangle(xy=(swatch_start_x, y - 9), width=swatch_width,
                      height=18, facecolor=name, edgecolor='0.7')
        )

    return fig

def convert_hex_to_rgb(hex_color: str) -> tuple:
    """Convert HEX color to RGB color.

    :param hex_color: str, like '#FFFFAA'
    :return: tuple, like (255, 255, 170)
    """
    value = hex_color.lstrip('#')
    value_len = len(value)
    return tuple(int(value[index: index + 2], 16) for index in range(0, value_len, 2))


def convert_rgb_to_hex(rgb_color: tuple) -> str:
    """Convert RGB color to HEX color.

    :param rgb_color: tuple, like (255, 255, 170)
    :return: str, like '#FFFFAA'
    """
    rgb_color = ('#%02x%02x%02x' % rgb_color).upper()
    return rgb_color


def map_color(value_vec, breaks, color_list):
    """
    INPUT:
        <value_vec>
            np.array or a list of values

        <breaks>
            A sorted value list, which can split all num into len(color_list) intervals.
            e.g. [0.01, 0.1, 0.5, 1] make all real num into 5 intervals, (-Inf,0.01], (0.01,0.1], (0.1, 0.5],  (0.5, 1], (1, +Inf]

        <color_list>
            A hex-format color list, which have to match with breaks

    RETURN
        <value_color_vec>
            A list map the value_vec with breaks
    """
    value_idx_list = []

    for value in value_vec:
        match_state = False
        for index, break_value in enumerate(breaks):
            if value <= break_value:
                value_idx_list.append(index)
                match_state = True
                break

        if not match_state:
            value_idx_list.append(index + 1)

    return tuple(color_list[col_idx] for col_idx in value_idx_list)


def make_color_list(low_color_RGB, high_color_RGB, length_out=20, return_fmt="HEX", log_level='DEBUG'):
    """
    INPUT
        <low_color_RGB> <high_color_RGB>
            Format like (210, 179, 150), tuple, list, or np.array

        <back_format>
            Hex OR RGB

    RETURN
        <color_list>
    """
    # set logger
    logger = get_logger(level=log_level, module_name=__module_name__, func_name=sys._getframe().f_code.co_name)

    return_fmt = return_fmt.upper()
    supported_fmt = ('HEX', 'RGB')

    if return_fmt not in supported_fmt:
        logger.critical(
            f'not supported color format: {return_fmt}\n'
            f'supported_fmt = {supported_fmt}'
        )
    low_color = np.array(low_color_RGB)
    high_color = np.array(high_color_RGB)

    color_list = []
    for index in range(0, length_out + 1):
        rgb_color = [abs(i) for i in low_color + (high_color - low_color) // length_out * index]
        if return_fmt == "HEX":
            color_list.append(convert_rgb_to_hex(tuple(rgb_color)))
        else:
            color_list.append(tuple(rgb_color))

    return color_list


if __name__ == '__main':
    # %%% test func plot_colortable
    colors = ['#64C1E8',
              '#80CED7',
              '#63C7B2',
              '#8E6C88',
              '#CA61C3',
              '#FF958C',
              '#883677']
    plot_colortable(colors, ncols=1, labels=[1, 2, 3, 4, 5, 6, 7])
    # plt.show()