"""TODO."""

import math
import os
import shutil
import subprocess
import sys
import tempfile
from io import StringIO
from pathlib import Path
from typing import cast

import numpy as np
import py3Dmol
from Bio import SeqIO
from Bio.PDB.MMCIFParser import MMCIFParser
from Bio.PDB.PDBIO import PDBIO
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.Structure import Structure as BiopythonStructure
from Bio.PDB.Superimposer import Superimposer
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from bioat.exceptions import (
    BioatFileFormatError,
    BioatFileNotFoundError,
    BioatInvalidParameterError,
)
from bioat.lib.libalignment import (
    get_aligned_seq,
    get_best_alignment,
    instantiate_pairwise_aligner,
)
from bioat.lib.libcolor import map_colors_between_two
from bioat.lib.libfastx import (
    AMINO_ACIDS_1CODE,
    AMINO_ACIDS_3CODE,
    AMINO_ACIDS_3CODE_EXTEND,
)
from bioat.logger import LoggerManager

__all__ = [
    "get_cut2ref_aln_info",
    "load_structure",
    "show_ref_cut",
    "structure2string",
]

lm = LoggerManager(mod_name="bioat.lib.libpdb")
lm.set_level("DEBUG")


def _load_seq(seq: str | Path | Seq, label=None):
    if isinstance(seq, str):
        seq = Path(seq)

    if (
        isinstance(seq, Path)
        and seq.is_file()
        and seq.suffix.lower() in [".fasta", ".fa"]
    ):
        seq = next(SeqIO.parse(seq, "fasta"))
        label = seq.id if label is None else label
        seq = seq.seq
    elif isinstance(seq, Seq):
        label = "no_name" if label is None else label
    else:
        msg = f"Invalid sequence format. seq: {seq}, type={type(seq)}, seq.is_file()={seq.is_file()}, seq.suffix.lower()={seq.suffix.lower()}"
        raise BioatInvalidParameterError(msg)
    return seq, label


def load_structure(
    structure: str | Path | BiopythonStructure, label: str | None = None
) -> tuple[BiopythonStructure, str]:
    """Load a PDB structure file or a BiopythonStructure object.

    Args:
        structure (str | Path | BiopythonStructure): Structure file path or BiopythonStructure object.
        label (str | None, optional): Structure label. If None, the label will be set to the file name. Defaults to None.

    Raises:
        BioatInvalidParameterError: If the input structure is not a valid file path or a BiopythonStructure object.

    Returns:
        tuple[BiopythonStructure, str]: A tuple of BiopythonStructure object and label.
    """
    if isinstance(structure, str):
        structure = Path(structure)

    if isinstance(structure, Path) and structure.is_file():
        parser = PDBParser(QUIET=True)
        structure = cast(BiopythonStructure, parser.get_structure(label, structure))
    if isinstance(structure, BiopythonStructure):
        label = cast(str, structure.id if label is None else label)
        return structure, label
    msg = f"Invalid structure format. structure: {structure}"
    raise BioatInvalidParameterError(msg)


def structure2string(structure: BiopythonStructure) -> str:
    """Convert BiopythonStructure object to string format.

    Args:
        structure (BiopythonStructure): Structure object to convert.

    Returns:
        str: PDB context.
    """
    pdb_io = PDBIO()
    pdb_io.set_structure(structure)
    # Use StringIO to save to memory string
    pdb_string = StringIO()
    pdb_io.save(pdb_string)
    # Return string
    return pdb_string.getvalue()


def pdb2fasta(
    input_pdb: str | Path | BiopythonStructure,
    output_fasta: str | Path | None = None,
    func_return: bool = False,
    log_level: str = "DEBUG",
) -> list[Seq] | None:
    """Converts a PDB file to a FASTA file.

    Details:
        1. **Proteins**:
            The protein sequence for each chain will be extracted as "Chain X Protein".
        2. **DNA and RNA**:
            Bases for DNA (A, T, G, C) will be saved as "Chain X DNA", and bases for RNA (A, U, G, C) will be saved as "Chain X RNA".
        3. **Other molecules**:
            Any unrecognized molecules (e.g., ions, modified molecules) will be labeled as [residue] and stored as "Chain X Other molecules".
        4. **Multi-chain complexes**:
            The program supports multi-chain structures in complexes, and the content of each chain will be recorded separately.

    Args:
        input_pdb (str or Path or BiopythonStructure):
            Path to the input PDB/CIF file or Biopython Structure.
        output_fasta (str or Path, optional):
            Output file path. If None, the output file will be named as the
            basename of the input file with a ".fa" extension. Defaults to None.
        func_return: (bool, optional)
            Whether to return a list of SeqRecord objects, useful when used as a function but not for command line. Defaults to False.
        log_level (str, optional):
            Logging level. Defaults to "WARNING".

    Returns:
        List of SeqRecord if func_return is True, otherwise None.
    """
    lm.set_names(func_name="pdb2fasta")
    lm.set_level(log_level)

    structure = None
    structure_id = "structure"

    # === Handle input file or structure object ===
    if isinstance(input_pdb, str | Path):
        input_pdb = Path(input_pdb)

        if not input_pdb.exists():
            msg = f"Input file not found: {input_pdb}"
            raise BioatFileNotFoundError(msg)

        suffix = input_pdb.suffix.lower()
        parser = (
            MMCIFParser(QUIET=True)
            if suffix.endswith(".cif") or suffix.endswith(".cif.gz")
            else PDBParser(QUIET=True)
            if suffix.endswith(".pdb") or suffix.endswith(".pdb.gz")
            else None
        )
        if parser is None:
            raise BioatFileFormatError("Only .pdb and .cif files are allowed.")

        structure_id = input_pdb.stem
        with input_pdb.open("rt") as f:
            structure = parser.get_structure(structure_id, f)

        if output_fasta is None:
            output_fasta = input_pdb.with_suffix(".fa")
            lm.logger.debug(f"Output is not defined, set output to {output_fasta}")
        else:
            output_fasta = Path(output_fasta)
            if output_fasta.suffix.lower() not in (".fa", ".fasta"):
                raise ValueError("Output file must have a '.fa' / '.fasta' suffix")

    elif isinstance(input_pdb, BiopythonStructure):
        structure = input_pdb
        structure_id = getattr(structure, "id", "structure")
    else:
        msg = f"Invalid input type. Expect Path or BiopythonStructure, got {type(input_pdb)}"
        raise BioatInvalidParameterError(msg)

    lm.logger.debug(f"Processing structure: {structure_id}")

    records = []

    if isinstance(structure, BiopythonStructure):
        for model in structure:
            for chain in model:
                lm.logger.debug(f"Processing chain {chain.id}")
                residues = list(chain)
                if not residues:
                    continue

                res_ids = [res.id[1] for res in residues]
                offset = min(res_ids)
                length = max(res_ids) - offset + 1

                protein = ["-"] * length
                dna = ["-"] * length
                rna = ["-"] * length
                other = ["-"] * length

                for res in residues:
                    idx = res.id[1] - offset
                    name = res.resname.strip().upper()

                    if name in AMINO_ACIDS_3CODE:
                        protein[idx] = AMINO_ACIDS_3CODE[name]
                    elif name in AMINO_ACIDS_3CODE_EXTEND:
                        protein[idx] = AMINO_ACIDS_3CODE_EXTEND[name]
                    elif name in ("DA", "DT", "DG", "DC"):
                        dna[idx] = name[1]
                    elif name in ("A", "U", "G", "C"):
                        rna[idx] = name
                    else:
                        other[idx] = f"[{name}]"

                def _append(seq, label, chain):
                    if any(c != "-" for c in seq):
                        records.append(
                            SeqRecord(
                                Seq("".join(seq)),
                                id=f"{structure_id}|Chain_{chain.id}|{label}",
                                description=f"Chain {chain.id} {label}",
                            )
                        )

                _append(protein, "Protein", chain)
                _append(dna, "DNA", chain)
                _append(rna, "RNA", chain)
                _append(other, "UNKNOWN", chain)
    else:
        raise BioatInvalidParameterError("Invalid structure format.")
    if output_fasta:
        SeqIO.write(records, output_fasta, "fasta")
        lm.logger.info(f"FASTA saved to: {output_fasta}")

    return [s.seq for s in records] if func_return else None


def _map_ref_colors(
    ref_map_colors: tuple,
    ref_map_values: dict,
    ref_map_value_random: bool,
) -> tuple:
    """Map reference residue colors based on values.

    This function maps residue-level scalar values (e.g., importance scores) to colors
    along a linear gradient defined by `ref_map_colors`. Each residue is identified by
    a string key in the format "X_123", where X is the 1-letter amino acid code and 123 is the residue index.

    Args:
        ref_map_colors (tuple): A tuple of two color values (e.g., hex or RGB), representing the color gradient endpoints.
                                For example: ("blue", "red") or ("#0000FF", "#FF0000").
        ref_map_values (dict): A dictionary mapping residue keys to float values.
                               Each key must be in the format "<residue>_<index>", such as "A_25".
                               The value represents a scalar quantity (e.g., attention weight) to be color-mapped.
        ref_map_value_random (bool): If True, the actual mapped values will be replaced with random numbers in [0, 1],
                                     ignoring the values in `ref_map_values`. Useful for debugging or visualization testing.

    Returns:
        tuple:
            - colors (List[str or tuple]): List of color values corresponding to each residue, generated by interpolation.
            - res_names (List[str]): List of 1-letter residue codes.
            - res_idxes (List[int]): List of residue indices (integers parsed from the keys).

    Raises:
        BioatInvalidParameterError:
            - If `ref_map_colors` is not a tuple of length 2.
            - If a residue name in `ref_map_values` is not a valid 1-letter amino acid code.

    Example:
        >>> ref_map_colors = ("blue", "red")
        >>> ref_map_values = {"A_25": 0.1, "V_26": 0.9}
        >>> colors, names, idxes = _map_ref_colors(ref_map_colors, ref_map_values, False)
        >>> print(colors)  # List of interpolated colors between blue and red
        >>> print(names)  # ['A', 'V']
        >>> print(idxes)  # [25, 26]
    """
    if len(ref_map_colors) != 2:
        raise BioatInvalidParameterError("ref_map_colors must be a tuple of length 2")

    color1, color2 = ref_map_colors

    values, res_names, res_idxes = [], [], []
    for k, v in ref_map_values.items():
        res_name, res_idx = k.split("_")
        if res_name not in AMINO_ACIDS_1CODE:
            msg = f"Invalid residue name: {res_name}"
            raise BioatInvalidParameterError(msg)
        res_names.append(res_name)
        res_idxes.append(int(res_idx))
        values.append(v)

    # 映射颜色
    rng = np.random.default_rng()
    mapped_values = rng.random(len(values)) if ref_map_value_random else values

    colors = map_colors_between_two(color1, color2, mapped_values)
    return colors, res_names, res_idxes


def show_ref_cut(
    ref_seq: str | Path | Seq,
    ref_pdb: str | Path | BiopythonStructure,
    cut_seq: list[str | Path | Seq] | str | Path | Seq | None = None,
    cut_pdb: list[str | Path | BiopythonStructure]
    | str
    | Path
    | BiopythonStructure
    | None = None,
    cut_labels: list[str] | str | None = None,
    ref_color: str = "red",
    ref_map_colors: tuple[str, str] | None = None,
    ref_map_values: dict | None = None,
    cut_color="lightgray",
    gap_color="purple",
    ref_style="cartoon",
    cut_style="cartoon",
    gap_style="cartoon",
    ref_map_value_random: bool = False,
    output_fig: str | Path | None = None,
    col: int = 4,
    scale: float = 1.0,
    annotate: bool = True,
    text_interval: int = 5,
    log_level="WARNING",
):
    """Visualizes the alignment of sequences and highlights changes in PDB structures using py3Dmol.

    Args:
        ref_seq (str or Seq): Amino acid sequence content for the ref protein.
        ref_pdb (str or BiopythonStructure): Path to the PDB file of the reference structure.
        cut_seq (str, Seq or None, optional): Amino acid sequence content for the cut protein.
        cut_pdb (str, BiopythonStructure or None, optional): Path to the PDB file of the cut structure.
        cut_labels (list[str] or str or None, optional): Label for the cut proteins. If None, the label will be set to "cut".
        ref_color (str, optional): Color for reference residues.
        ref_map_colors (tuple[str, str] or None, optional): ref_map_colors will be used as color bar from ref_map_colors[0] to ref_map_colors[1]. If None, do not apply color mapping. Defaults to None.
        ref_map_values (dict or None, optional): A dictionary of values for the ref color map, it will be normalized to the range of [0 - 1]. If None, all residues will be colored with the same color. e.g. ref_map_values = {'V_0': 0.4177215189873418, 'S_1': 0.8185654008438819, 'K_2': 0.9915611814345991, 'G_3': 0.42616033755274263, ...}
        cut_color (str, optional): Color for cut residues.
        gap_color (str, optional): Color for gaps or removed residues.
        ref_style (str, optional): "stick", "sphere", "cartoon", or "line"
        cut_style (str, optional): "stick", "sphere", "cartoon", or "line"
        gap_style (str, optional): "stick", "sphere", "cartoon", or "line"
        ref_map_value_random (bool, optional): If True, ref_map_values will be randomly generated. Defaults to False.
        output_fig (str or Path or None, optional): Output figure file path. If None, the figure will not be saved in html format. Defaults to None.
        col (int, optional): Number of columns for the visualization. Defaults to 3.
        scale (float, optional): Scale factor for the visualization. Defaults to 1.0.
        annotate (bool, optional): Whether to annotate the visualization with labels. Defaults to True.
        text_interval (int, optional): The interval between text annotations. Defaults to 5.
        log_level (str, optional): Log level. Defaults to "WARNING".
    """
    lm.set_names(func_name="show_ref_cut")
    lm.set_level(log_level)
    if isinstance(cut_seq, str) and "," in cut_seq:
        cut_seq = cast(list[str | Path | Seq], cut_seq.split(","))
    if isinstance(cut_pdb, str) and "," in cut_pdb:
        cut_pdb = cast(list[str | Path | BiopythonStructure], cut_pdb.split(","))
    ref_seq, ref_label = _load_seq(ref_seq, "ref")
    ref_pdb, ref_label = load_structure(ref_pdb, ref_label)

    aligner = instantiate_pairwise_aligner(
        scoring_match=5,
        penalty_mismatch=-100,
        penalty_gap_open=-1,
        penalty_gap_extension=0,
        penalty_query_left_gap_score=0,
        penalty_query_right_gap_score=0,
        log_level=log_level,
    )

    if cut_seq:  # not None
        cut_seqs = [cut_seq] if not isinstance(cut_seq, list) else cut_seq
        cut_pdbs = [cut_pdb] if not isinstance(cut_pdb, list) else cut_pdb

        ref_seqs_fix = []
        ref_pdbs_fix = []
        ref_labels = []
        cut_seqs_fix = []
        cut_pdbs_fix = []
        cut_labels_fix = []
        gap_indices = []
        rmsds = []
        alignment_dicts = []

        for idx, _cut_seq in enumerate(cut_seqs):
            ref_seqs_fix.append(ref_seq)
            ref_labels.append(ref_label)
            # 执行Sequence Alignment并得到gaps的位置
            _cut_seq, cut_label = _load_seq(_cut_seq, "cut")
            cut_seqs_fix.append(_cut_seq)

            if cut_labels:
                if isinstance(cut_labels, list):
                    cut_labels_fix.append(cut_labels[idx])
                elif isinstance(cut_labels, str):
                    cut_labels_fix.append(cut_labels)
            else:
                cut_labels_fix.append(cut_label)

            # 获得gaps的位置
            # Perform alignment
            alignment = get_best_alignment(
                seq_a=ref_seq,
                seq_b=_cut_seq,
                aligner=aligner,
                consider_strand=False,
            )
            alignment_dict = get_aligned_seq(alignment, reverse=False, letn_match=False)
            alignment_dicts.append(alignment_dict)
            gap_indices.append(  # 不指定 cut pdb 时必须在这里计算而不是get_cut2ref_aln_info中
                [i for i, n in enumerate(alignment_dict["target_seq"]) if n == "-"],
            )

            cut_pdb = cut_pdbs[idx]

            if cut_pdb:
                # 执行PDB Alignment并得对齐后的PDB文件和RMSD
                res = get_cut2ref_aln_info(
                    ref_pdb,
                    cut_pdb,
                    cal_rmsd=True,
                    cal_tmscore=False,
                    label1=ref_label,
                    label2=cut_label,
                    log_level=log_level,
                )
                ref_pdbs_fix.append(res[ref_label])
                cut_pdbs_fix.append(res[cut_label])
                rmsds.append(res["RMSD"])

                if not isinstance(res[cut_label], BiopythonStructure):
                    msg = "Invalid input format. cut_pdb must be (one or list of) BiopythonStructure"
                    raise BioatInvalidParameterError(msg)
            else:
                ref_pdbs_fix.append(ref_pdb)
                cut_pdbs_fix.append(None)
                rmsds.append(None)
    else:
        lm.logger.info(
            "Cut sequence is not provided, skip alignment and annotation for gaps, and just show ref visualization.",
        )
        ref_seqs_fix = [ref_seq]
        ref_pdbs_fix = [ref_pdb]
        ref_labels = [ref_label]
        cut_seqs_fix = [None]
        cut_pdbs_fix = [None]
        cut_labels_fix = [None]
        gap_indices = [None]
        rmsds = [None]
        alignment_dicts = [None]

    n = len(ref_seqs_fix)
    col = min(n, col)
    row = math.ceil(n / col)

    view = py3Dmol.view(
        width=int(500 * scale * col),
        height=int(500 * scale * row),
        viewergrid=(row, col) if n > 1 else None,
        linked=True,  # 是否同步网格间的分子运动
    )
    lm.logger.debug(f"viewergrid: ({row}, {col})")
    viewers = [(i // col, i % col) for i in range(n)]
    lm.logger.debug(f"viewers = {viewers}")

    for i, viewer in enumerate(viewers):
        lm.logger.debug(f"viewer = {viewer}")
        _add_one_submodel(
            view=view,
            ref_pdb=ref_pdbs_fix[i],
            ref_label=ref_labels[i],
            cut_pdb=cut_pdbs_fix[i],
            cut_label=cut_labels_fix[i],
            ref_color=ref_color,
            ref_map_colors=ref_map_colors,
            ref_map_values=ref_map_values,
            ref_map_value_random=ref_map_value_random,
            cut_color=cut_color,
            gap_color=gap_color,
            gap_indices=gap_indices[i],
            ref_style=ref_style,
            cut_style=cut_style,
            gap_style=gap_style,
            alignment_dict=alignment_dicts[i],
            rmsd=rmsds[i],
            annotate=annotate,
            text_interval=text_interval,
            viewer=viewer,
        )

    if output_fig:
        output_fig = Path(output_fig)
        if output_fig.suffix != ".html":
            lm.logger.warning(
                f"Output file must have a '.html' extension, but got {output_fig}"
            )
            sys.exit(1)
        else:
            with output_fig.open("w") as f:
                f.write(view._make_html())  # noqa: SLF001
            lm.logger.info(f"Visualization saved to: {output_fig}")
    return view.show()


def _add_one_submodel(
    view,
    ref_pdb,
    ref_label,
    cut_pdb,
    cut_label,
    ref_color,
    ref_map_colors,
    ref_map_values,
    ref_map_value_random,
    cut_color,
    gap_color,
    gap_indices,
    ref_style,
    cut_style,
    gap_style,
    alignment_dict,
    rmsd,
    annotate,
    text_interval,
    viewer,
):
    # for annotation
    label_xpos = -70
    label_ypos = 75
    label_zpos = -30

    if ref_map_colors is not None:
        colors, res_names, res_idxes = _map_ref_colors(
            ref_map_colors,
            ref_map_values,
            ref_map_value_random,
        )
    else:
        colors, res_names, res_idxes = None, None, None
    color1, color2 = None, None
    if ref_map_colors is not None:
        color1, color2 = ref_map_colors
    if ref_map_colors and ref_map_values:
        values = []
        res_names = []
        res_idxes = []
        for k, v in ref_map_values.items():
            res_name, res_idx = k.split("_")
            if res_name not in AMINO_ACIDS_1CODE:
                msg = f"Invalid residue name: {res_name}"
                raise BioatInvalidParameterError(msg)
            res_idx = int(res_idx)
            values.append(v)
            res_names.append(res_name)
            res_idxes.append(res_idx)
        if not ref_map_value_random:
            colors = map_colors_between_two(color1, color2, values)
        else:
            rng = np.random.default_rng()
            mapped_values = rng.random(len(values))
            colors = map_colors_between_two(color1, color2, mapped_values)

    view.setViewStyle({"style": "outline", "color": "black", "width": 0.1})

    view.addModel(
        structure2string(ref_pdb),
        "pdb",
        viewer=viewer,
    )
    view.setStyle(
        {"model": 0},
        {
            ref_style: {
                "arrows": True,
                "tubes": True,
                "style": "oval",
                "color": ref_color,
            },
        },
        viewer=viewer,
    )
    label_ypos -= text_interval
    view.addLabel(
        f"{ref_label}: {ref_color}",
        {
            "fontColor": ref_color,
            "fontSize": 14,
            "position": {"x": label_xpos, "y": label_ypos, "z": label_zpos},
        },
        viewer=viewer,
    )

    # 标记ref 值的颜色
    if (
        ref_map_colors is not None
        and colors is not None
        and res_names is not None
        and res_idxes is not None
    ):
        for i in range(len(colors)):
            idx = res_idxes[i]
            res = res_names[i]
            color = colors[i]
            view.addStyle(
                {"model": 0, "resi": idx + 1, "resn": AMINO_ACIDS_1CODE[res]},
                {ref_style: {"color": color, "radius": 0.5}},
                viewer=viewer,
            )
        if annotate:
            label_ypos -= text_interval
            view.addLabel(
                f"{ref_label} values: {ref_color} (0) to {ref_color} (1)",
                {
                    "fontColor": ref_color,
                    "fontSize": 14,
                    "position": {"x": label_xpos, "y": label_ypos, "z": label_zpos},
                },
                viewer=viewer,
            )
    # 标记gaps的颜色
    if gap_indices:
        for idx, res in enumerate(alignment_dict["reference_seq"]):
            # 标记 gap位置的颜色
            if idx in gap_indices:
                view.addStyle(
                    {"model": 0, "resi": idx + 1, "resn": AMINO_ACIDS_1CODE[res]},
                    {gap_style: {"color": gap_color, "radius": 0.5}},
                    viewer=viewer,
                )
        if annotate:
            label_ypos -= text_interval
            view.addLabel(
                f"gaps: {gap_color}",
                {
                    "fontColor": gap_color,
                    "fontSize": 14,
                    "position": {"x": label_xpos, "y": label_ypos, "z": label_zpos},
                },
                viewer=viewer,
            )
    if cut_pdb and rmsd:
        view.addModel(structure2string(cut_pdb), "pdb", viewer=viewer)
        view.setStyle(
            {"model": 1},
            {
                cut_style: {
                    "arrows": True,
                    "tubes": True,
                    "style": "oval",
                    "color": cut_color,
                },
            },
            viewer=viewer,
        )
        if annotate:
            label_ypos -= text_interval
            view.addLabel(
                f"{cut_label}: {cut_color}\tRMSD: {rmsd:.2f}",
                {
                    "fontColor": cut_color,
                    "fontSize": 14,
                    "position": {"x": label_xpos, "y": label_ypos, "z": label_zpos},
                },
                viewer=viewer,
            )

    # 聚焦显示
    # Finalize the view
    view.zoomTo(viewer=viewer)
    # view.render(viewer=viewer)


# 对齐函数
def get_cut2ref_aln_info(
    ref: str | Path | BiopythonStructure,
    cut: str | Path | BiopythonStructure,
    cal_rmsd=True,
    cal_tmscore=False,
    label1="ref",
    label2="cut",
    usalign_bin: str | Path = "usalign",
    log_level="WARNING",
) -> dict:
    """Align cutted pdb to ref pdb using the CA atoms.

    Aligns a truncated protein structure (cut) to its full-length reference structure (ref)
    using Ca atoms and Biopython's Superimposer.

    This function:
    - Extracts all Ca atoms from `ref` and `cut`
    - Removes atoms from `ref` at the indices listed in `gap_indices`
    - Aligns the remaining atoms from `cut` to the corresponding positions in `ref`
    - Modifies the `cut` structure in-place to match the aligned orientation
    - Returns both structures and the RMSD value of the alignment

    It assumes:
    - One-to-one correspondence between residues after gap removal
    - Structures are predicted by AlphaFold2 / ESMFold (no missing atoms)

    Args:
        ref (str or Bio.PDB.Structure.Structure): Reference structure path or loaded Structure.
        cut (str or Bio.PDB.Structure.Structure): Truncated structure path or loaded Structure.
        cal_rmsd (bool, optional): Whether to calculate RMSD. Default is True.
        cal_tmscore (bool, optional): Whether to calculate TM-score using USalign. Default is False.
        label1 (str, optional): Name for the reference structure. Default is "ref".
        label2 (str, optional): Name for the cut structure. Default is "cut".
        usalign_bin (str or Path, optional): Path to the USalign binary for TM-score calculation. Default is "usalign".
        log_level (str, optional): Logging level. Default is "WARNING".

    Returns:
        dict: {
                "{label1}": aln label1 structure,  # if cal_rmsd is True, unaltered label1 structure
                "{label2}}": fixed label2 structure,  # if cal_rmsd is True, fix label2 coords in-place
                "RMSD": 0.123  # if cal_rmsd is True, the RMSD value between label1 and label2
                f"{label1}_seq": ref_seq,  # if cal_rmsd is True, the sequence of label1 structure
                f"{label2}_seq": cut_seq,  # if cal_rmsd is True, the sequence of label2 structure
                "alignment_dict": alignment_dict,  # if cal_rmsd is True, the alignment dict of label1 and label2
                "gap_indices": gap_indices,  # if cal_rmsd is True, the indices of gaps in label1 structure
                "TM-score:mean": 0.623,  # if cal_tmscore is True, the mean TM-score value
                "TM-score:TM1": 0.456,  # if cal_tmscore is True, use label1 as ref <L_N> in calculation
                "TM-score:TM2": 0.789,  # if cal_tmscore is True, use label2 as ref <L_N> in calculation
                ...
            }
    """
    # 读取 PDB 文件
    ref, cut = (
        load_structure(ref, label1)[0],
        load_structure(cut, label2)[0],
    )
    aln_info = {}

    if cal_rmsd:
        lm.logger.debug(
            f"Calculating RMSD between {label1} and {label2} using CA atoms.",
        )
        ref_seq = cast(Seq, pdb2fasta(ref, func_return=True)[0])  # type: ignore
        cut_seq = cast(Seq, pdb2fasta(cut, func_return=True)[0])  # type: ignore
        lm.logger.debug(f"fetch seq info from structure: {ref.id} -> {ref_seq}")
        lm.logger.debug(f"fetch seq info from structure:  {cut.id} -> {cut_seq}")
        aligner = instantiate_pairwise_aligner(
            scoring_match=5,
            penalty_mismatch=-100,
            penalty_gap_open=-1,
            penalty_gap_extension=0,
            penalty_query_left_gap_score=0,
            penalty_query_right_gap_score=0,
            log_level=log_level,
        )

        # 获得gaps的位置
        # Perform alignment
        alignment = get_best_alignment(
            seq_a=ref_seq,
            seq_b=cut_seq,
            aligner=aligner,
            consider_strand=False,
        )
        alignment_dict = get_aligned_seq(alignment, reverse=False, letn_match=False)
        gap_indices = [
            i for i, n in enumerate(alignment_dict["target_seq"]) if n == "-"
        ]
        # 选择对齐的原子 (e.g., CA)
        atoms1 = [atom for atom in ref.get_atoms() if atom.get_name() == "CA"]
        atoms1 = [atoms1[i] for i in range(len(atoms1)) if i not in gap_indices]
        atoms2 = [atom for atom in cut.get_atoms() if atom.get_name() == "CA"]

        if len(atoms1) != len(atoms2):
            msg = (
                f"The number of CA atoms in {label1} and {label2} are not equal. "
                "Use ref and cut pdbs predicted by AlphaFold2/ESMFold or other tools "
                "but not pdb downloaded from <PDB database> to ensure the CA atoms are aligned."
            )
            raise BioatInvalidParameterError(msg)

        # 使用 Biopython 的 Superimposer 进行对齐
        super_imposer = Superimposer()
        super_imposer.set_atoms(atoms1, atoms2)
        super_imposer.apply(cut.get_atoms())  # 修改 pdb 的坐标
        rmsd = super_imposer.rms
        lm.logger.info(
            f"RMSD (strict, remove gaps in {label1}) between {label1} and {label2}: {rmsd:.3f}",
        )
        aln_info[label1] = ref
        aln_info[label2] = cut  # 坐标已对齐
        aln_info["RMSD"] = (
            rmsd  # it should be close to the result of pymol cmd.align("cut and name CA", "ref and name CA", cutoff=10.0, cycles=0)[0]
        )
        aln_info[f"{label1}_seq"] = str(ref_seq)
        aln_info[f"{label2}_seq"] = str(cut_seq)
        aln_info["alignment_dict"] = alignment_dict
        aln_info["gap_indices"] = gap_indices

    if cal_tmscore:
        lm.logger.debug(
            f"Calculating TM-score between {label1} and {label2} using CA atoms.",
        )
        # 判断 usalign_bin 是否可执行
        if (
            isinstance(usalign_bin, str)
            and "/" not in usalign_bin
            and shutil.which(usalign_bin) is None
        ) or (
            isinstance(usalign_bin, str | Path)
            and "/" in str(usalign_bin)
            and not (
                Path(usalign_bin).expanduser().resolve().is_file()
                and os.access(Path(usalign_bin).expanduser().resolve(), os.X_OK)
            )
        ):
            lm.logger.error(
                f"usalign binary not found in var 'PATH' or not executable: {usalign_bin}, please set the --usalign_bin <path/to/usalign> param or use conda install -c conda-forge usalign to install it!"
            )
            msg = f"usalign binary not found or not executable: {usalign_bin}"
            raise BioatFileNotFoundError(msg)

        def save_structure_to_tempfile(structure):
            io = PDBIO()
            io.set_structure(structure)
            with tempfile.NamedTemporaryFile(delete=False, suffix=".pdb") as tmp:
                io.save(tmp.name)
                return tmp.name

        ref_path = save_structure_to_tempfile(ref)
        cut_path = save_structure_to_tempfile(cut)
        # 使用 usalign 计算 TM-score
        cmd = [str(usalign_bin), ref_path, cut_path, "-outfmt", "2", "-mol", "prot"]
        result = subprocess.run(cmd, capture_output=True, text=True, check=False)  # noqa: S603

        if result.returncode != 0:
            msg = f"USalign 执行失败: \n{result.stderr}"
            raise RuntimeError(msg)

        lines = result.stdout.strip().splitlines()
        if len(lines) < 2:
            msg = "USalign 输出异常, 结果缺失。"
            raise ValueError(msg)

        fields = lines[1].split()
        tm1 = float(fields[2])
        tm2 = float(fields[3])
        aln_info["TM-score:mean"] = (tm1 + tm2) / 2.0
        aln_info["TM-score:TM1"] = tm1
        aln_info["TM-score:TM2"] = tm2
    return aln_info


if __name__ == "__main__":
    ref_seq = "data/pdb/EGFP_4EUL_ref.fa"
    cut_seq = "data/pdb/EGFP_4EUL_cut1.fa"
    ref_pdb = "data/pdb/EGFP_4EUL_ref.pdb"
    cut_pdb = "data/pdb/EGFP_4EUL_cut1.pdb"
    ref_map_values = {
        "M_0": 0.0,
        "V_1": 0.890295358649789,
        "S_2": 0.270042194092827,
        "K_3": 0.24050632911392406,
        "G_4": 0.8523206751054853,
        "E_5": 0.5991561181434599,
        "K_238": 0.5316455696202531,
    }
    show_ref_cut(
        ref_seq=ref_seq,
        ref_pdb=ref_pdb,
        cut_seq=cut_seq,
        cut_pdb=cut_pdb,
        ref_color="red",
        ref_map_colors=("black", "blue"),
        ref_map_values=ref_map_values,
        cut_color="green",
        gap_color="red",
        ref_style="cartoon",
        cut_style="cartoon",
        gap_style="cartoon",
        ref_map_value_random=False,
        output_fig=None,
        log_level="DEBUG",
    )
    aln_info = get_cut2ref_aln_info(
        ref=ref_pdb,
        cut=cut_pdb,
        cal_rmsd=True,
        cal_tmscore=True,
        label1="ref",
        label2="cut",
        usalign_bin="/Users/zhaohuanan/micromamba/envs/ProEvaluator/bin/usalign",
        log_level="DEBUG",
    )
    print(aln_info)
