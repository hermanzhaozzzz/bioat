import numpy as np
import sys
import logging

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


def make_color_list(low_color_RGB, high_color_RGB, length_out=20, return_fmt="HEX"):
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
    lib_name = __name__
    function_name = sys._getframe().f_code.co_name
    logger = logging.getLogger(f'{lib_name}.{function_name} ==> ')

    supported_fmt = ('HEX', 'RGB')
    return_fmt = return_fmt.upper()
    if return_fmt not in supported_fmt:
        logger.critical(
            f'not supported color format: {return_fmt}\n'
            f'supported_fmt = {supported_fmt}'
        )
    low_color = np.array(low_color_RGB)
    high_color = np.array(high_color_RGB)

    color_list = []
    for index in range(0, length_out + 1):
        rgb_color = low_color + (high_color - low_color) // length_out * index
        if return_fmt == "HEX":
            color_list.append(convert_rgb_to_hex(tuple(rgb_color)))
        else:
            color_list.append(tuple(rgb_color))

    return color_list
