import io
import inspect
import logging
import xml.etree.ElementTree as ET

import matplotlib

matplotlib.use("Agg", force=True)

import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.colors import to_hex

from bioat.lib.libplot import BIOAT_MORANDI_PALETTE, init_matplotlib


def _svg_tspan_texts(svg_path):
    root = ET.parse(svg_path).getroot()
    return [
        element.text
        for element in root.iter()
        if element.tag.endswith("}tspan") or element.tag == "tspan"
    ]


def test_init_matplotlib_publication_defaults():
    init_matplotlib(log_level="ERROR")

    assert plt.rcParams["axes.unicode_minus"] is False
    assert plt.rcParams["figure.dpi"] == 300
    assert plt.rcParams["svg.fonttype"] == "none"
    assert plt.rcParams["mathtext.default"] == "regular"
    assert plt.rcParams["pdf.use14corefonts"] is False
    assert plt.rcParams["pdf.fonttype"] == 42
    assert plt.rcParams["ps.useafm"] is False
    assert plt.rcParams["ps.fonttype"] == 42
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
    assert "pdf_fonttype" not in params
    assert "ps_fonttype" not in params
    assert "mathtext_default" not in params
    assert "patch_savefig_for_illustrator" not in params


def test_init_matplotlib_keeps_advanced_kwargs_overridable():
    init_matplotlib(
        log_level="ERROR",
        axes_unicode_minus=True,
        mathtext_default="it",
        pdf_fonttype=3,
        ps_fonttype=3,
        set_backend_svg="path",
    )

    assert plt.rcParams["axes.unicode_minus"] is True
    assert plt.rcParams["mathtext.default"] == "it"
    assert plt.rcParams["pdf.fonttype"] == 3
    assert plt.rcParams["ps.fonttype"] == 3
    assert plt.rcParams["svg.fonttype"] == "path"


def test_init_matplotlib_respects_seaborn_parameters():
    init_matplotlib(log_level="ERROR")
    default_font_size = plt.rcParams["font.size"]

    init_matplotlib(sns_style="whitegrid", sns_font_scale=2.0, log_level="ERROR")

    assert plt.rcParams["axes.grid"] is True
    assert plt.rcParams["font.size"] > default_font_size


def test_init_matplotlib_patches_savefig_once_and_logs():
    logger = logging.getLogger("bioat.lib.libplot.init_matplotlib")
    stream = io.StringIO()
    handler = logging.StreamHandler(stream)
    logger.addHandler(handler)
    try:
        init_matplotlib(log_level="INFO")
        first_savefig = plt.savefig

        init_matplotlib(log_level="INFO")
        second_savefig = plt.savefig
    finally:
        logger.removeHandler(handler)

    assert first_savefig is second_savefig
    assert getattr(plt.savefig, "_bioat_illustrator_patch", False) is True
    assert "Illustrator-friendly text spans" in stream.getvalue()


def test_svg_mathtext_spans_are_merged_for_illustrator(tmp_path):
    init_matplotlib(log_level="ERROR")

    fig, ax = plt.subplots(figsize=(3, 1.2))
    ax.axis("off")
    labels = ["CFTR$_{wt}$", "CFTR$_{SP101}$", "CFTR$_{Y122X}$"]
    for index, label in enumerate(labels):
        ax.text(0.1, 0.8 - index * 0.25, label, fontsize=12)

    svg_path = tmp_path / "mathtext.svg"
    plt.savefig(svg_path)
    plt.close(fig)

    tspan_texts = _svg_tspan_texts(svg_path)
    assert tspan_texts.count("CFTR") == 3
    assert "wt" in tspan_texts
    assert "SP101" in tspan_texts
    assert "Y122X" in tspan_texts
    assert "w" not in tspan_texts
    assert "t" not in tspan_texts
