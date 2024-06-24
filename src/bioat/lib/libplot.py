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
from glob import glob

import matplotlib.pyplot as plt
import numpy as np
from dna_features_viewer import GraphicFeature, GraphicRecord
from matplotlib.patches import Rectangle

from bioat.exceptions import BioatException, BioatRuntimeWarning
from bioat.lib.libpath import HOME
from bioat.logger import get_logger

__all__ = ["init_matplotlib", "plot_colortable"]
__module_name__ = "bioat.lib.libplot"

BIOAT_DEFAULT_FONTS_DATAPATH = os.path.join(os.path.dirname(__file__), "libplot")
BIOAT_DEFAULT_FONTS = [
    "Helvetica-Bold.ttf",
    "Helvetica-BoldOblique.ttf",
    "Helvetica-Light.ttf",
    "Helvetica-Oblique.ttf",
    "Helvetica.ttf",
]


def _copy_fonts(refresh=False, log_level="warning"):
    """copy fonts from bioat package to matplotlib package

    :param log_level: log level for logger
    :type log_level: string
    :raises FileNotFoundError: if not found which site-packages to copy fonts to
    """
    logger = get_logger(
        level=log_level,
        module_name=__module_name__,
        func_name="_copy_fonts",
    )
    if refresh:
        files_to_remove = glob(os.path.join(HOME, ".cache", "matplotlib", "fontlist*"))
        for file in files_to_remove:
            os.remove(file)
    try:
        to_path = None
        for i in sys.path:
            if i.endswith("site-packages"):
                to_path = i
                break
        if not to_path:
            logger.warning(BioatRuntimeWarning("site-packages not found in sys.path"))
            return
        to_path = [i for i in sys.path if i.endswith("site-packages")][0]
        to_path = os.path.join(
            to_path, "matplotlib", "mpl-data", "fonts", "ttf"
        )
        from_path = BIOAT_DEFAULT_FONTS_DATAPATH
        logger.debug(f"Copying fonts from {from_path} to {to_path}")
        for font in BIOAT_DEFAULT_FONTS:
            shutil.copyfile(
                os.path.join(from_path, font),
                os.path.join(to_path, font),
            )
        logger.debug("Fonts copied successfully")
    except Exception as e:
        logger.error(BioatException(f"Failed to copy fonts: {e}"))


def init_matplotlib(
    style="ggplot", font="Helvetica", refresh=False, log_level="INFO", **kwargs
):
    """easily set matplotlib style

    :param style: matplotlib style, defaults to 'ggplot'
    :type style: str, optional
    :param font: use what font in matplotlib, defaults to 'Helvetica'
    :type font: str, optional
    :param refresh: wether to remove matplotlib font cache and reset it, defaults to False
    :type refresh: bool, optional
    :param log_level: log level, defaults to 'INFO'
    :type log_level: str, optional
    :param set_backend_pdf: whether to use core fonts for the PDF backend, defaults to True
    :type set_backend_pdf: bool, optional
    :param set_backend_ps: whether to use core fonts for the PS backend, defaults to True
    :type set_backend_ps: bool, optional
    :param set_backend_svg: whether to use 'none' to replace 'path' (use font but not plot path for characters)for the SVG backend, defaults to 'none'
    :type set_backend_svg: str, optional
    :raises BioatException: if failed to copy fonts
    :raises BioatRuntimeWarning: if site-packages not found in sys.path
    """
    logger = get_logger(
        level=log_level,
        module_name=__module_name__,
        func_name="init_matplotlib",
    )
    logger.info('Initializing matplotlib')
    _copy_fonts(refresh=refresh, log_level=log_level)
    logger.info(
        f"set: plt.style.use('{style}')\n"
        "# set matplotlib style theme\n"
        "# ref: https://matplotlib.org/stable/api/style_api.html"
    )
    plt.style.use(style)

    logger.info(f"set: plt.rcParams['font.family'] = '{font}'")
    plt.rcParams["font.family"] = font
    logger.info(f"set: plt.rcParams['font.sans-serif'] = ['{font}']")
    plt.rcParams["font.sans-serif"] = [font]

    # set backends
    set_backend_pdf = kwargs.get("set_backend_pdf", True)
    set_backend_ps = kwargs.get("set_backend_ps", True)
    set_backend_svg = kwargs.get("set_backend_svg", "none")

    if set_backend_pdf:
        logger.info(
            f"set: plt.rcParams['pdf.use14corefonts'] = {set_backend_pdf} \n"
            "# whether to use core fonts for the PDF backend\n"
            "# ref: https://matplotlib.org/stable/api/matplotlib_configuration_api.html#matplotlib.rcParams)"
        )
        plt.rcParams["pdf.use14corefonts"] = set_backend_pdf
    if set_backend_ps:
        logger.info(
            f"set: plt.rcParams['ps.useafm'] = {set_backend_ps}\n"
            "# whether to use core fonts for the PS backend"
        )
        plt.rcParams["ps.useafm"] = set_backend_ps
    if set_backend_svg:
        logger.info(
            f"set: plt.rcParams['svg.fonttype'] = '{set_backend_svg}'\n"
            "# whether to use 'none' to replace 'path' (use font but not plot path for characters)for the SVG backend"
        )
        plt.rcParams["svg.fonttype"] = set_backend_svg
    logger.info("matplotlib initialized successfully")


def plot_colortable(colors, *, ncols=4):
    cell_width = 212
    cell_height = 22
    swatch_width = 48
    margin = 12

    names = list(colors)

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

        ax.text(
            text_pos_x,
            y,
            names[i],
            fontsize=14,
            horizontalalignment="left",
            verticalalignment="center",
        )

        ax.add_patch(
            Rectangle(xy=(swatch_start_x, y - 9), width=swatch_width,
                      height=18, facecolor=name, edgecolor='0.7')
        )
    plt.show()
    plt.close()


def set_figsize(figsize=(3.5, 2.5)):
    """Set the figure size for matplotlib.

    Defined in :numref:`sec_calculus`"""
    # use_svg_display()
    plt.rcParams["figure.figsize"] = figsize


def set_axes(axes, xlabel, ylabel, xlim, ylim, xscale, yscale, legend):
    """Set the axes for matplotlib.

    Defined in :numref:`sec_calculus`"""
    axes.set_xlabel(xlabel)
    axes.set_ylabel(ylabel)
    axes.set_xscale(xscale)
    axes.set_yscale(yscale)
    axes.set_xlim(xlim)
    axes.set_ylim(ylim)

    if legend:
        axes.legend(legend)
    axes.grid()


def plot(
    X,
    Y=None,
    xlabel=None,
    ylabel=None,
    legend=[],
    xlim=None,
    ylim=None,
    xscale="linear",
    yscale="linear",
    fmts=("-", "m--", "g-.", "r:"),
    figsize=(3.5, 2.5),
    axes=None,
):
    """Plot data points.

    Defined in :numref:`sec_calculus`"""

    def has_one_axis(X):  # True if X (tensor or list) has 1 axis
        return (
            hasattr(X, "ndim")
            and X.ndim == 1
            or isinstance(X, list)
            and not hasattr(X[0], "__len__")
        )

    if has_one_axis(X):
        X = [X]
    if Y is None:
        X, Y = [[]] * len(X), X
    elif has_one_axis(Y):
        Y = [Y]
    if len(X) != len(Y):
        X = X * len(Y)

    set_figsize(figsize)
    if axes is None:
        axes = plt.gca()
    axes.cla()
    for x, y, fmt in zip(X, Y, fmts):
        axes.plot(x, y, fmt) if len(x) else axes.plot(y, fmt)
    set_axes(axes, xlabel, ylabel, xlim, ylim, xscale, yscale, legend)
    # return axes


if __name__ == '__main__':
    # test for init_matplotlib
    init_matplotlib(log_level="warning")
    plot_colortable(
        ["#64C1E8", "#80CED7", "#63C7B2", "#8E6C88", "#CA61C3", "#FF958C", "#883677"],
        ncols=1,
    )
    def normal(x, mu, sigma):
        p = 1 / math.sqrt(2 * math.pi * sigma**2)
        return p * np.exp(-0.5 / sigma**2 * (x - mu) ** 2)

    x = np.arange(-7, 7, 0.01)
    params = [(0, 1), (0, 2), (3, 1)]
    plot(
        x,
        [normal(x, mu, sigma) for mu, sigma in params],
        xlabel="x",
        ylabel="p(x)",
        figsize=(4.5, 2.5),
        legend=[f"mean {mu}, std {sigma}" for mu, sigma in params],
    )
    plt.show()
    plt.close()
    # !test
    import pandas as pd

    df = pd.DataFrame(
        {
            "name": ["AAaAAAAAAAAAAAAAAAA", "B", "C"] + ["DR"] * 4 + ["spacer"] * 4,
            "type": ["CDS", "CDS", "CDS"] + ["DR"] * 4 + ["spacer"] * 4,
            "start": (
                [0, 100, 200] + list(range(300, 500, 66)) + list(range(336, 536, 66))
            ),
            "end": (
                [50, 130, 280] + list(range(336, 540, 66)) + list(range(366, 566, 66))
            ),
            "strand": [1, -1, 1] + [1] * 4 + [1] * 4,
            "color": (
                [
                    "#ffd700",
                    "#ffcccc",
                    "#cffccc",
                ]
                + ["black"] * 4
                + ["red"] * 4
            ),
        }
    )
    print(df)
    features = []
    for i, row in df.iterrows():
        label = row.get("name", "?")
        feature = row.get("type", "?")
        start = row.get("start", 0)
        end = row.get("end", 0)
        strand = row.get("strand", 0)
        color = row.get("color", "black")
        feature = GraphicFeature(
            start=start, end=end, strand=strand, color=color, label=label
        )
        features.append(feature)

    record = GraphicRecord(sequence_length=600, features=features)
    record.plot(
        figure_width=20,
        draw_line=True,
        with_ruler=True,
        ruler_color=None,
        plot_sequence=False,
        annotate_inline=True,
        max_label_length=50,
        max_line_length=30,
        level_offset=0,
        strand_in_label_threshold="default",
        elevate_outline_annotations="default",
        x_lim=None,
        figure_height=None,
        sequence_params=None,
    )
    plt.show()
