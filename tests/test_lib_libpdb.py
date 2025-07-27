from functools import partial

import numpy as np
import pytest
from Bio import SeqIO
from Bio.PDB.mmcifio import MMCIFIO
from Bio.PDB.PDBParser import PDBParser
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from bioat.exceptions import (
    BioatFileFormatError,
    BioatFileNotFoundError,
    BioatInvalidParameterError,
)
from bioat.lib.libpdb import (
    _load_seq,
    _map_ref_colors,
    load_structure,
    pdb2fasta,
    show_ref_cut,
    structure2string,
)

from ._pytest_meta import DATA_PATH

ref_fasta_file = DATA_PATH / "pdb/EGFP_4EUL_ref.fa"
ref_pdb_file = DATA_PATH / "pdb/EGFP_4EUL_ref.pdb"
cut_fasta_file = DATA_PATH / "pdb/EGFP_4EUL_cut1.fa"
cut_pdb_file = DATA_PATH / "pdb/EGFP_4EUL_cut1.pdb"


def test_load_seq_from_file():
    seq, label = _load_seq(ref_fasta_file)
    assert isinstance(seq, Seq)
    assert isinstance(label, str)
    assert len(seq) > 0
    assert label != ""


def test_load_seq_from_seq_object():
    raw_seq = Seq("MKTG")
    seq, label = _load_seq(raw_seq, label="myseq")
    assert seq == raw_seq
    assert label == "myseq"


def test_load_seq_from_temp_file(tmp_path):
    record = SeqRecord(Seq("ACDEFGHIK"), id="temp_test")
    tmp_file = tmp_path / "temp.fa"
    SeqIO.write(record, tmp_file, "fasta")

    seq, label = _load_seq(tmp_file)
    assert str(seq) == "ACDEFGHIK"
    assert label == "temp_test"


def test_load_structure_from_file():
    structure, label = load_structure(ref_pdb_file, "test_label")
    assert hasattr(structure, "id")
    assert isinstance(label, str)
    assert any(structure.get_chains()), "Structure should have chains"


def test_load_structure_from_structure_object():
    parser = PDBParser(QUIET=True)
    raw_structure = parser.get_structure("test_label", ref_pdb_file)
    structure, label = load_structure(raw_structure)  # type: ignore
    assert structure is raw_structure
    assert label == "test_label"


def test_load_structure_invalid_input():
    with pytest.raises(BioatInvalidParameterError):
        load_structure(12345)  # type: ignore


def test_structure2string_valid():
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("test_id", ref_pdb_file)

    pdb_str = structure2string(structure)  # type: ignore

    assert isinstance(pdb_str, str)
    assert pdb_str.startswith("HEADER") or "ATOM" in pdb_str
    assert (
        "END" in pdb_str
        or pdb_str.strip().endswith("TER")
        or pdb_str.strip().endswith("END")
    )
    assert len(pdb_str) > 100  # 简单长度判断


def test_pdb2fasta_return_records():
    """Test .pdb input returns SeqRecord list."""
    records = pdb2fasta(ref_pdb_file, func_return=True, log_level="WARNING")
    assert isinstance(records, list)
    assert all(r.__class__.__name__ == "Seq" for r in records)


def test_pdb2fasta_write_file(tmp_path):
    """Test .pdb input writes fasta file."""
    out = tmp_path / "out.fa"
    result = pdb2fasta(ref_pdb_file, output_fasta=out, log_level="WARNING")
    assert out.exists()
    records = list(SeqIO.parse(out, "fasta"))
    assert len(records) > 0
    assert result is None


def test_pdb2fasta_func_return_and_file(tmp_path):
    """Test return + file write consistency."""
    out = tmp_path / "combined.fa"
    records = pdb2fasta(
        ref_pdb_file, output_fasta=out, func_return=True, log_level="WARNING"
    )
    assert out.exists()
    records2 = list(SeqIO.parse(out, "fasta"))
    assert len(records) == len(records2)  # type: ignore


def test_pdb2fasta_with_cif(tmp_path):
    """Test generated .cif input works."""
    cif_path = tmp_path / "test.cif"
    # Convert pdb -> cif via Biopython
    parser = PDBParser(QUIET=True)
    io = MMCIFIO()
    structure = parser.get_structure("4EUL", ref_pdb_file)
    io.set_structure(structure)
    io.save(str(cif_path))

    assert cif_path.exists()

    records = pdb2fasta(cif_path, func_return=True, log_level="WARNING")
    assert isinstance(records, list)
    assert all(r.__class__.__name__ == "Seq" for r in records)


def test_pdb2fasta_invalid_file_raises():
    with pytest.raises(BioatFileNotFoundError):
        pdb2fasta("not_exist.pdb")


def test_pdb2fasta_invalid_suffix_raises(tmp_path):
    bad_file = tmp_path / "bad.xyz"
    bad_file.write_text("HEADER FAKE")
    with pytest.raises(BioatFileFormatError):
        pdb2fasta(bad_file)


def test_pdb2fasta_invalid_object_type():
    with pytest.raises(BioatInvalidParameterError):
        pdb2fasta(12345)  # type: ignore


def test__map_ref_colors_normal():
    colors, names, idxes = _map_ref_colors(
        ref_map_colors=("#000000", "#FFFFFF"),
        ref_map_values={"A_10": 0.0, "G_11": 1.0},
        ref_map_value_random=False,
    )
    assert names == ["A", "G"]
    assert idxes == [10, 11]
    assert colors[0].lower() == "#000000"
    assert colors[1].lower() == "#ffffff"


def test__map_ref_colors_random(monkeypatch):
    # Patch RNG to deterministic

    monkeypatch.setattr(
        np.random, "default_rng", partial(np.random.default_rng, seed=42)
    )

    ref_map_colors = ("red", "blue")
    ref_map_values = {
        "A_1": 1.0,
        "V_2": 0.5,
        "G_3": 0.0,
    }

    colors, names, idxes = _map_ref_colors(
        ref_map_colors=ref_map_colors,
        ref_map_values=ref_map_values,
        ref_map_value_random=True,
    )

    assert names == ["A", "V", "G"]
    assert idxes == [1, 2, 3]
    assert len(colors) == 3


def test__map_ref_colors_invalid_color_tuple():
    with pytest.raises(BioatInvalidParameterError):
        _map_ref_colors(
            ref_map_colors=("#000000",),  # 长度错误
            ref_map_values={"A_10": 0.1},
            ref_map_value_random=False,
        )


def test__map_ref_colors_invalid_residue_name():
    with pytest.raises(BioatInvalidParameterError):
        _map_ref_colors(
            ref_map_colors=("#000000", "#FFFFFF"),
            ref_map_values={"XYZ_9": 0.5},  # 不合法的残基名
            ref_map_value_random=False,
        )


def test_show_ref_cut_only_ref(tmp_path):
    out_file = tmp_path / "out.html"
    show_ref_cut(
        ref_seq=ref_fasta_file,
        ref_pdb=ref_pdb_file,
        annotate=True,
        output_fig=out_file,
    )
    assert out_file.exists()
    assert out_file.stat().st_size > 1000


def test_show_ref_cut_ref_and_only_cut_seq(tmp_path):
    out_file = tmp_path / "out.html"
    show_ref_cut(
        ref_seq=ref_fasta_file,
        cut_seq=cut_fasta_file,
        ref_pdb=ref_pdb_file,
        output_fig=out_file,
    )
    assert out_file.exists()
    assert out_file.stat().st_size > 1000


def test_show_ref_cut_ref_and_cut(tmp_path):
    out_file = tmp_path / "out.html"
    show_ref_cut(
        ref_seq=ref_fasta_file,
        cut_seq=cut_fasta_file,
        ref_pdb=ref_pdb_file,
        cut_pdb=cut_pdb_file,
        cut_labels="cut1",
        output_fig=out_file,
    )
    assert out_file.exists()
    assert out_file.stat().st_size > 1000


def test_show_ref_cut_ref_and_cuts_list(tmp_path):
    out_file = tmp_path / "out.html"
    show_ref_cut(
        ref_seq=ref_fasta_file,
        cut_seq=[cut_fasta_file, cut_fasta_file],
        ref_pdb=ref_pdb_file,
        cut_pdb=[cut_pdb_file, cut_pdb_file],
        cut_labels=["cut1", "cut2"],
        output_fig=out_file,
    )
    assert out_file.exists()
    assert out_file.stat().st_size > 1000


def test_show_ref_cut_ref_and_cuts_comma_separated_list(tmp_path):
    out_file = tmp_path / "out.html"
    show_ref_cut(
        ref_seq=ref_fasta_file,
        cut_seq=f"{cut_fasta_file},{cut_fasta_file}",
        ref_pdb=ref_pdb_file,
        cut_pdb=f"{cut_pdb_file},{cut_pdb_file}",
        cut_labels="cut1,cut2",
        output_fig=out_file,
    )
    assert out_file.exists()
    assert out_file.stat().st_size > 1000
