import inspect

import matplotlib

matplotlib.use("Agg", force=True)

import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.colors import to_hex

from bioat.lib.libplot import BIOAT_MORANDI_PALETTE, init_matplotlib


def test_init_matplotlib_publication_defaults():
    plt.rcdefaults()
    init_matplotlib(log_level="ERROR")

    assert plt.rcParams["axes.unicode_minus"] is False
    assert plt.rcParams["figure.dpi"] == 300
    assert plt.rcParams["svg.fonttype"] == "none"
    assert plt.rcParams["pdf.use14corefonts"] is True
    assert plt.rcParams["ps.useafm"] is True
    assert plt.rcParams["mathtext.fontset"] == "dejavusans"
    assert plt.rcParams["mathtext.default"] == "it"
    assert plt.rcParams["font.family"] == ["Helvetica"]
    assert plt.rcParams["axes.grid"] is False
    assert [
        to_hex(color).upper()
        for color in sns.color_palette()[: len(BIOAT_MORANDI_PALETTE)]
    ] == list(BIOAT_MORANDI_PALETTE)


def test_init_matplotlib_keeps_backend_options_in_kwargs():
    params = inspect.signature(init_matplotlib).parameters

    assert "kwargs" in params
    assert "style" not in params
    assert "refresh" in params
    assert "sns_context" in params
    assert "sns_style" in params
    assert "sns_palette" in params
    assert "sns_font_scale" in params
    assert "figure_dpi" in params
    assert "axes_unicode_minus" not in params
    assert "pdf_use14corefonts" not in params
    assert "ps_useafm" not in params
    assert "split_simple_math_titles" not in params
    assert "set_backend_svg" not in params


def test_init_matplotlib_keeps_advanced_kwargs_overridable():
    plt.rcdefaults()
    init_matplotlib(
        log_level="ERROR",
        axes_unicode_minus=True,
        pdf_use14corefonts=False,
        ps_useafm=False,
        set_backend_svg="none",
    )

    assert plt.rcParams["axes.unicode_minus"] is True
    assert plt.rcParams["pdf.use14corefonts"] is False
    assert plt.rcParams["ps.useafm"] is False
    assert plt.rcParams["svg.fonttype"] == "none"


def test_init_matplotlib_respects_seaborn_parameters():
    plt.rcdefaults()
    init_matplotlib(sns_style="white", log_level="ERROR")

    assert plt.rcParams["axes.grid"] is False


def test_init_matplotlib_patches_simple_script_titles():
    plt.rcdefaults()
    init_matplotlib(log_level="ERROR")
    fig, ax = plt.subplots()

    title = ax.set_title("Protein$_{mini}$ title", fontsize=7)

    assert title.get_text() == ""
    assert hasattr(ax, "_bioat_simple_script_title_artist")
    plt.close(fig)


def test_init_matplotlib_patches_multiple_simple_script_titles():
    plt.rcdefaults()
    init_matplotlib(log_level="ERROR")
    fig, ax = plt.subplots()

    title = ax.set_title("CFTR$_{wt}$ 16HBE14o$^-$ cell", fontsize=7)

    assert title.get_text() == ""
    assert hasattr(ax, "_bioat_simple_script_title_artist")
    plt.close(fig)


def test_init_matplotlib_keeps_complex_mathtext_native():
    plt.rcdefaults()
    init_matplotlib(log_level="ERROR")
    fig, ax = plt.subplots()

    title = ax.set_title(r"$\alpha$ title", fontsize=7)

    assert title.get_text() == r"$\alpha$ title"
    assert not hasattr(ax, "_bioat_simple_script_title_artist")
    plt.close(fig)


def test_init_matplotlib_can_disable_simple_script_title_patch():
    plt.rcdefaults()
    init_matplotlib(split_simple_math_titles=False, log_level="ERROR")
    fig, ax = plt.subplots()

    title = ax.set_title("Protein$_{mini}$ title", fontsize=7)

    assert title.get_text() == "Protein$_{mini}$ title"
    assert not hasattr(ax, "_bioat_simple_script_title_artist")
    plt.close(fig)
