"""bioat.lib.libplot

author: Herman Huanan Zhao
email: hermanzhaozzzz@gmail.com
homepage: https://github.com/hermanzhaozzzz

This module provides functions for plotting.

example 1:
    init_matplotlib
        <in python consolo>:
            >>> import matplotlib.pyplot as plt
            >>> from bioat.lib.libplot import init_matplotlib
            >>> init_matplotlib(log_level='info')
            >>> plt.plot([1, 2, 3], [4, 5, 6])
            >>> plt.show()

example 2:
    plot_colortable
        <in python consolo>:
            >>> from bioat.lib.libplot import plot_colortable
            >>> colors = ['#64C1E8', '#80CED7', '#63C7B2', '#8E6C88', '#CA61C3', '#FF958C', '#883677']
            >>> plot_colortable(colors, ncols=1, labels=[1, 2, 3, 4, 5, 6, 7])
            >>> plt.show()
"""

import math
import os
import shutil
import sys

import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.patches import Rectangle

from bioat.logger import get_logger

__all__ = ["init_matplotlib", "plot_colortable"]
__module_name__ = "bioat.lib.libplot"

DATAPATH = os.path.join(os.path.dirname(__file__), "libplot")


def _copy_fonts(log_level):
    logger = get_logger(
        level=log_level,
        module_name=__module_name__,
        func_name="_copy_fonts",
    )
    try:
        to_path = None
        for i in sys.path:
            if i.endswith("site-packages"):
                to_path = i
                break
        if not to_path:
            raise FileNotFoundError("site-packages not found in sys.path")
        to_path = [i for i in sys.path if i.endswith("site-packages")][0]
        to_path = os.path.join(to_path, "matplotlib/mpl-data/fonts/ttf")
        from_path = DATAPATH
        logger.debug(f"Copying fonts from {from_path} to {to_path}")
        fonts = [
            "Helvetica-Bold.ttf",
            "Helvetica-BoldOblique.ttf",
            "Helvetica-Light.ttf",
            "Helvetica-Oblique.ttf",
            "Helvetica.ttf",
        ]
        for font in fonts:
            shutil.copyfile(
                os.path.join(from_path, font),
                os.path.join(to_path, font),
            )
        logger.debug("Fonts copied successfully")
    except Exception as e:
        logger.error(f"Failed to copy fonts: {e}")


def init_matplotlib(log_level='INFO'):
    logger = get_logger(
        level=log_level,
        module_name=__module_name__,
        func_name="init_matplotlib",
    )
    logger.info('Initializing matplotlib')
    _copy_fonts(log_level)
    logger.debug('ref: https://matplotlib.org/stable/api/style_api.html')
    logger.info("set: plt.style.use('ggplot')  # use ggplot style")
    plt.style.use('ggplot')

    logger.info("set: plt.rcParams['font.family'] = 'Helvetica'")
    plt.rcParams['font.family'] = 'Helvetica'
    logger.info("set: plt.rcParams['font.sans-serif'] = ['Helvetica']")
    plt.rcParams['font.sans-serif'] = ['Helvetica']

    logger.debug('ref: https://matplotlib.org/stable/api/matplotlib_configuration_api.html#matplotlib.rcParams')
    logger.info("set: plt.rcParams['pdf.use14corefonts'] = True  # trigger core fonts for PDF backend")
    plt.rcParams['pdf.use14corefonts'] = True
    logger.info("set: plt.rcParams['ps.useafm'] = True  # trigger core fonts for PS backend")
    plt.rcParams['ps.useafm'] = True
    logger.info("set: plt.rcParams['svg.fonttype'] = 'none'  # change 'path' to 'none' (use font)for SVG backend")
    plt.rcParams['svg.fonttype'] = 'none'


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


if __name__ == '__main__':
    # test for init_matplotlib
    init_matplotlib(log_level='debug')
    # %%% test func plot_colortable
    colors = ['#64C1E8',
              '#80CED7',
              '#63C7B2',
              '#8E6C88',
              '#CA61C3',
              '#FF958C',
              '#883677']
    plot_colortable(colors, ncols=1, labels=[1, 2, 3, 4, 5, 6, 7])
    plt.show()
