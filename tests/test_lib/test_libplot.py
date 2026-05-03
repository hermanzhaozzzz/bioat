import inspect

import matplotlib

matplotlib.use("Agg", force=True)

import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.colors import to_hex

from bioat.lib.libplot import BIOAT_MORANDI_PALETTE, init_matplotlib


def test_init_matplotlib_publication_defaults():
    init_matplotlib(log_level="ERROR")

    assert plt.rcParams["axes.unicode_minus"] is False
    assert plt.rcParams["figure.dpi"] == 300
    assert plt.rcParams["svg.fonttype"] == "path"
    assert plt.rcParams["pdf.use14corefonts"] is True
    assert plt.rcParams["ps.useafm"] is True
    assert plt.rcParams["axes.grid"] is False
    assert [
        to_hex(color).upper()
        for color in sns.color_palette()[: len(BIOAT_MORANDI_PALETTE)]
    ] == list(BIOAT_MORANDI_PALETTE)


def test_init_matplotlib_keeps_illustrator_backend_options_in_kwargs():
    params = inspect.signature(init_matplotlib).parameters

    assert "kwargs" in params
    assert "refresh" in params
    assert "style" not in params
    assert "axes_unicode_minus" not in params
    assert "set_backend_pdf" not in params
    assert "set_backend_ps" not in params
    assert "set_backend_svg" not in params


def test_init_matplotlib_keeps_advanced_kwargs_overridable():
    init_matplotlib(
        log_level="ERROR",
        axes_unicode_minus=True,
        set_backend_pdf=False,
        set_backend_ps=False,
        set_backend_svg="none",
    )

    assert plt.rcParams["axes.unicode_minus"] is True
    assert plt.rcParams["pdf.use14corefonts"] is False
    assert plt.rcParams["ps.useafm"] is False
    assert plt.rcParams["svg.fonttype"] == "none"


def test_init_matplotlib_respects_seaborn_parameters():
    init_matplotlib(log_level="ERROR")
    default_font_size = plt.rcParams["font.size"]

    init_matplotlib(sns_style="whitegrid", sns_font_scale=2.0, log_level="ERROR")

    assert plt.rcParams["axes.grid"] is True
    assert plt.rcParams["font.size"] > default_font_size
