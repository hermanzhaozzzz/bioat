import pytest

from bioat.foldtools import FoldTools

from ._pytest_meta import DATA_PATH


@pytest.fixture
def foldtools_fixture(tmp_path):
    return {
        "fold_tools": FoldTools(),
        "ref_seq": DATA_PATH / "pdb/dCas9.fa",
        "ref_pdb": DATA_PATH / "pdb/dCas9.pdb",
        "cut_seq": DATA_PATH / "pdb/dCas9_rm_HNH.fa",
        "cut_pdb": DATA_PATH / "pdb/dCas9_rm_HNH.pdb",
        "output_fig": tmp_path / "output.html",
        "output_fasta": tmp_path / "test.fa",
    }


def test_show_ref_cut(foldtools_fixture):
    ft = foldtools_fixture
    ft["fold_tools"].show_ref_cut(
        ref_seq=ft["ref_seq"],
        ref_pdb=ft["ref_pdb"],
        cut_seq=ft["cut_seq"],
        cut_pdb=ft["cut_pdb"],
        output_fig=ft["output_fig"],
        ref_color="blue",
        cut_color="green",
        gap_color="red",
        ref_style="cartoon",
        cut_style="cartoon",
        gap_style="cartoon",
        log_level="DEBUG",
    )
    assert ft["output_fig"].exists()


def test_pdb2fasta(foldtools_fixture):
    ft = foldtools_fixture
    ft["fold_tools"].pdb2fasta(
        input_pdb=ft["ref_pdb"],
        output_fasta=ft["output_fasta"],
        log_level="INFO",
    )
    assert ft["output_fasta"].exists()
    content = ft["output_fasta"].read_text()
    assert content.startswith(">")
