import gzip
import math
import os
from io import StringIO
from typing import List

import Bio
import numpy as np
import py3Dmol
from Bio import SeqIO
from Bio.PDB import PDBIO, MMCIFParser, PDBParser, Superimposer
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
from bioat.lib.libpath import is_file
from bioat.logger import LoggerManager

__all__ = ["load_structure", "structure_to_string", "show_ref_cut"]

lm = LoggerManager(mod_name="bioat.lib.libpdb")
lm.set_level("DEBUG")

def _load_seq(seq, label=None, log_level="WARNING"):
    if isinstance(seq, str) and is_file(seq, log_level=log_level):
        seq = next(SeqIO.parse(seq, "fasta"))
        label = seq.id if label is None else label
        seq = seq.seq
    elif isinstance(seq, Seq):
        label = "no_name" if label is None else label
    else:
        raise BioatInvalidParameterError(f"Invalid sequence format. seq: {seq}")
    return seq, label


def load_structure(structure, label=None, log_level="WARNING"):
    if isinstance(structure, str) and is_file(structure, log_level=log_level):
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure(label, structure)
    elif isinstance(structure, Bio.PDB.Structure.Structure):
        label = structure.id if label is None else label
    else:
        raise BioatInvalidParameterError(
            f"Invalid structure format. structure: {structure}. "
            "file path or Bio.PDB.Structure.Structure is expected."
        )
    return structure, label


def structure_to_string(structure: Bio.PDB.Structure.Structure):
    """Convert Bio.PDB.Structure.Structure to string.

    Args:
        structure (Bio.PDB.Structure.Structure): Structure object to convert.

    Returns:
        str: PDB context。
    """
    pdb_io = PDBIO()
    pdb_io.set_structure(structure)

    # 使用 StringIO 保存到内存中的字符串
    pdb_string = StringIO()
    pdb_io.save(pdb_string)

    # 返回字符串
    return pdb_string.getvalue()


# 对齐函数
def _align_cut2ref(
    ref: str | Bio.PDB.Structure.Structure,
    cut: str | Bio.PDB.Structure.Structure,
    gap_indices: list[int],
    label1="ref",
    label2="cut",
    log_level="WARNING",
) -> dict:
    """Align cutted pdb to ref pdb using the CA atoms.
    Aligns a truncated protein structure (cut) to its full-length reference structure (ref)
    using Cα atoms and Biopython's Superimposer. Residues that are deleted in the cut structure
    must be specified using `gap_indices` (i.e., their positions in the reference structure).

    This function:
    - Extracts all Cα atoms from `ref` and `cut`
    - Removes atoms from `ref` at the indices listed in `gap_indices`
    - Aligns the remaining atoms from `cut` to the corresponding positions in `ref`
    - Modifies the `cut` structure in-place to match the aligned orientation
    - Returns both structures and the RMSD value of the alignment

    It assumes:
    - One-to-one correspondence between residues after gap removal
    - Structures are predicted by AlphaFold2 / ESMFold (no missing atoms)

    Example:
        Given:
            ref = "CRISPRCRISPR"
            cut = "CRISRCIPR"
            gap_indices = [4, 7, 9]

        Then:
            aln_ref = CRISPRCRISPR
            aln_cut = CRIS-RC-I-PR
                          *  * *
                      012345678901
                        (aligned by removing ref[4], ref[7], ref[9])
        Finally, the function will return:
            {
                "ref": aln_ref structure,  # unaltered ref structure
                "cut": fixed_cut structure,  # fix cut coords in-place
                "RMSD": 0.123  # example RMSD value
            }

    Args:
        ref (str or Bio.PDB.Structure.Structure): Reference structure path or loaded Structure.
        cut (str or Bio.PDB.Structure.Structure): Truncated structure path or loaded Structure.
        gap_indices (list[int]): Indices (0-based) in ref of residues that are deleted in cut.
        label1 (str, optional): Name for the reference structure. Default is "ref".
        label2 (str, optional): Name for the cut structure. Default is "cut".
        log_level (str, optional): Logging level. Default is "WARNING".

    Returns:
        dict: {
            label1: Bio.PDB.Structure.Structure (unaltered),
            label2: Bio.PDB.Structure.Structure (aligned in-place),
            "RMSD": float
        }
    """
    # 读取 PDB 文件
    ref, cut = (
        load_structure(ref, label1, log_level)[0],
        load_structure(cut, label2, log_level)[0],
    )

    # 选择对齐的原子 (e.g., CA)
    atoms1 = [atom for atom in ref.get_atoms() if atom.get_name() == "CA"]
    atoms1 = [atoms1[i] for i in range(len(atoms1)) if i not in gap_indices]
    atoms2 = [atom for atom in cut.get_atoms() if atom.get_name() == "CA"]

    assert len(atoms1) == len(atoms2), (
        "The number of CA atoms in pdb1 and pdb2 are not equal. "
        "Use ref and cut pdbs predicted by AlphaFold2/ESMFold or other tools "
        "but not pdb downloaded from <PDB database> to ensure the CA atoms are aligned."
    )
    # 使用 Biopython 的 Superimposer 进行对齐
    super_imposer = Superimposer()
    super_imposer.set_atoms(atoms1, atoms2)
    super_imposer.apply(cut.get_atoms())  # 修改 pdb 的坐标
    rmsd = super_imposer.rms
    lm.logger.info(
        f"RMSD (strict, remove gaps in {label1}) between {label1} and {label2}: {rmsd:.3f}"
    )

    return {
        label1: ref,
        label2: cut,  # 坐标已对齐
        "RMSD": rmsd,  # it should be close to the result of pymol cmd.align("cut and name CA", "ref and name CA", cutoff=10.0, cycles=0)[0]
    }

def _map_ref_colors(ref_map_colors, ref_map_values, ref_map_value_random):
    # fix color mapping
    assert (ref_map_colors is None and ref_map_values is None) or (
        isinstance(ref_map_colors, tuple)
        and isinstance(ref_map_values, dict)
        and len(ref_map_colors) == 2
    ), "ref_color_base and ref_map_values must be both None or both not None"
    color1, color2 = None, None
    if ref_map_colors is not None:
        color1, color2 = ref_map_colors

    if ref_map_colors and ref_map_values:
        values = []
        res_names = []
        res_idxes = []
        for k, v in ref_map_values.items():
            res_name, res_idx = k.split("_")
            assert res_name in AMINO_ACIDS_1CODE, f"Invalid residue name: {res_name}"
            res_idx = int(res_idx)
            values.append(v)
            res_names.append(res_name)
            res_idxes.append(res_idx)
        if not ref_map_value_random:
            colors = map_colors_between_two(color1, color2, values)
        else:
            colors = map_colors_between_two(color1, color2, np.random.rand(len(values)))
        return colors, res_names, res_idxes
    else:
        return None, None, None


def show_ref_cut(
    ref_seq: str | Seq,
    ref_pdb: str | Bio.PDB.Structure.Structure,
    cut_seq: List[str | Seq] | str | Seq | None = None,
    cut_pdb: List[str | Bio.PDB.Structure.Structure]
    | str
    | Bio.PDB.Structure.Structure
    | None = None,
    cut_labels: List[str] | None = None,
    ref_color: str = "red",
    ref_map_colors: tuple[str] | None = None,
    ref_map_values: dict | None = None,
    cut_color="lightgray",
    gap_color="purple",
    ref_style="cartoon",
    cut_style="cartoon",
    gap_style="cartoon",
    ref_map_value_random: bool = False,
    output_fig: str | None = None,
    col: int = 4,
    scale: float = 1.0,
    annotate: bool = True,
    text_interval: int = 5,
    log_level="WARNING",
):
    """
    Visualizes the alignment of sequences and highlights changes in PDB structures using py3Dmol.

    Args:
        ref_seq (str or Seq): Amino acid sequence content for the ref protein.
        ref_pdb (str or Bio.PDB.Structure.Structure): Path to the PDB file of the reference structure.
        cut_seq (str, Seq or None, optional): Amino acid sequence content for the cut protein.
        cut_pdb (str, Bio.PDB.Structure.Structure or None, optional): Path to the PDB file of the cut structure.
        cut_labels (str or None, optional): Label for the cut proteins. If None, the label will be set to "cut".
        ref_color (str, optional): Color for reference residues.
        ref_map_colors (tuple[str, str] or None, optional): ref_map_colors will be used as color bar from ref_map_colors[0] to ref_map_colors[1]. If None, do not apply color mapping. Defaults to None.
        ref_map_values (dict or None, optional): A dictionary of values for the ref color map, it will be normalized to the range of [0 - 1]. If None, all residues will be colored with the same color. e.g. ref_map_values = {'V_0': 0.4177215189873418, 'S_1': 0.8185654008438819, 'K_2': 0.9915611814345991, 'G_3': 0.42616033755274263, ...}
        cut_color (str, optional): Color for cut residues.
        gap_color (str, optional): Color for gaps or removed residues.
        ref_style (str, optional): "stick", "sphere", "cartoon", or "line"
        cut_style (str, optional): "stick", "sphere", "cartoon", or "line"
        gap_style (str, optional): "stick", "sphere", "cartoon", or "line"
        ref_map_value_random (bool, optional): If True, ref_map_values will be randomly generated. Defaults to False.
        output_fig (str or None, optional): Output figure file path. If None, the figure will not be saved in html format. Defaults to None.
        col (int, optional): Number of columns for the visualization. Defaults to 3.
        scale (float, optional): Scale factor for the visualization. Defaults to 1.0.
        annotate (bool, optional): Whether to annotate the visualization with labels. Defaults to True.
        text_interval (int, optional): The interval between text annotations. Defaults to 5.
        log_level (str, optional): Log level. Defaults to "WARNING".
    """
    lm.set_names(func_name="show_ref_cut")
    lm.set_level(log_level)
    lm.logger.debug(
        f"""\
Params:
-------
ref_seq: {ref_seq}
ref_pdb: {ref_pdb}
cut_seq: {cut_seq}
cut_pdb: {cut_pdb}
ref_color: {ref_color}
ref_map_colors: {ref_map_colors}
ref_map_values: {ref_map_values}
cut_color: {cut_color}
gap_color: {gap_color}
ref_style: {ref_style}
cut_style: {cut_style}
gap_style: {gap_style}
ref_map_value_random: {ref_map_value_random},
output_fig : {output_fig},
col: {col},
scale: {scale},
annotate: {annotate},
text_interval: {text_interval},
log_level: {log_level}"""
    )

    ref_seq, ref_label = _load_seq(ref_seq, "ref", log_level=log_level)
    ref_pdb, ref_label = load_structure(ref_pdb, ref_label, log_level=log_level)
    assert isinstance(ref_pdb, Bio.PDB.Structure.Structure), (
        "Invalid input format. ref_pdb must be Bio.PDB.Structure"
    )

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
        if not isinstance(cut_seq, list):
            cut_seqs = [cut_seq]
        else:
            cut_seqs = cut_seq
        if not isinstance(cut_pdb, list):
            cut_pdbs = [cut_pdb]
        else:
            cut_pdbs = cut_pdb

        ref_seqs_fix = []
        ref_pdbs_fix = []
        ref_labels = []
        cut_seqs_fix = []
        cut_pdbs_fix = []
        cut_labels_fix = []
        gap_indices = []
        rmsds = []
        alignment_dicts = []

        for idx, cut_seq in enumerate(cut_seqs):
            ref_seqs_fix.append(ref_seq)
            ref_labels.append(ref_label)
            # 执行Sequence Alignment并得到gaps的位置
            cut_seq, cut_label = _load_seq(cut_seq, "cut")
            cut_seqs_fix.append(cut_seq)

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
                seq_a=ref_seq, seq_b=cut_seq, aligner=aligner, consider_strand=False
            )
            alignment_dict = get_aligned_seq(alignment, reverse=False, letn_match=False)
            alignment_dicts.append(alignment_dict)
            gap_indices.append(
                [i for i, n in enumerate(alignment_dict["target_seq"]) if n == "-"]
            )

            cut_pdb = cut_pdbs[idx]

            if cut_pdb:
                # 执行PDB Alignment并得对齐后的PDB文件和RMSD
                res = _align_cut2ref(
                    ref_pdb, cut_pdb, gap_indices[idx], ref_label, cut_label, log_level
                )
                ref_pdbs_fix.append(res[ref_label])
                cut_pdbs_fix.append(res[cut_label])
                rmsds.append(res["RMSD"])

                assert isinstance(res[cut_label], Bio.PDB.Structure.Structure), (
                    "Invalid input format. cut_pdb must be (one or list of) Bio.PDB.Structure"
                )
            else:
                ref_pdbs_fix.append(ref_pdb)
                cut_pdbs_fix.append(None)
                rmsds.append(None)
    else:
        lm.logger.info(
            "Cut sequence is not provided, skip alignment and annotation for gaps, and just show ref visualization."
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
            ref_seq=ref_seqs_fix[i],
            ref_pdb=ref_pdbs_fix[i],
            ref_label=ref_labels[i],
            cut_seq=cut_seqs_fix[i],
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
        assert output_fig.endswith(".html"), "Output file must have a '.html' extension"
        with open(output_fig, "w") as f:
            f.write(view._make_html())
    return view.show()


def _add_one_submodel(
    view,
    ref_seq,
    ref_pdb,
    ref_label,
    cut_seq,
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

    colors, res_names, res_idxes = _map_ref_colors(
        ref_map_colors, ref_map_values, ref_map_value_random
    )
    color1, color2 = None, None
    if ref_map_colors is not None:
        color1, color2 = ref_map_colors
    if ref_map_colors and ref_map_values:
        values = []
        res_names = []
        res_idxes = []
        for k, v in ref_map_values.items():
            res_name, res_idx = k.split("_")
            assert res_name in AMINO_ACIDS_1CODE, f"Invalid residue name: {res_name}"
            res_idx = int(res_idx)
            values.append(v)
            res_names.append(res_name)
            res_idxes.append(res_idx)
        if not ref_map_value_random:
            colors = map_colors_between_two(color1, color2, values)
        else:
            colors = map_colors_between_two(color1, color2, np.random.rand(len(values)))

    view.setViewStyle({"style": "outline", "color": "black", "width": 0.1})

    view.addModel(
        structure_to_string(ref_pdb),
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
    if ref_map_colors is not None:
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
        view.addModel(structure_to_string(cut_pdb), "pdb", viewer=viewer)
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


def pdb2fasta(pdb_file, output_fasta, log_level="DEBUG"):
    """Converts a PDB / CIF file to a FASTA file.

        This function processes the provided PDB file and extracts protein, DNA,
        RNA sequences, and other molecules appropriately to create a FASTA file.
    e   1. Proteins:The protein sequence for each chain will be extracted as Chain X Protein.
        2. DNA and RNA: Bases for DNA (A, T, G, C) will be saved as Chain X DNA and bases for RNA (A, U, G, C) will be saved as Chain X RNA.
        3. Other molecules: Any unrecognized molecules (e.g., ions, modified molecules) will be labeled as [residue] and stored as Chain X Other molecules.
        4. Multi-chain complexes: The program supports multi-chain structures in complexes, and the content of each chain will be recorded separately.

        Args:
            pdb_file (str): input file path.
            output_fasta (str | None, optional): output file path. If None, the output file will be named as the basename of the input file with a ".fa" extension. Defaults to None.
            log_level (str, optional): log level. Defaults to "WARNING".
    """
    lm.set_names(func_name="pdb2fasta")
    lm.set_level(log_level)
    lm.logger.debug(
        f"""\
Params:
-------
input: {input}
output_fasta: {output_fasta}
log_level: {log_level}"""
    )

    # 设置输出文件名
    if output_fasta is None:
        if pdb_file.endswith(".gz"):
            idx = ".".join(os.path.basename(pdb_file).split(".")[:-2])
        else:
            idx = ".".join(os.path.basename(pdb_file).split(".")[:-1])
        output_fasta = f"{idx}.fa" if idx else "output.fa"
        lm.logger.debug(f"Output is not defined, set output to {output_fasta}")
    else:
        assert output_fasta.endswith(".fa"), "Output file must have a '.fa' extension"
        idx = ".".join(os.path.basename(output_fasta).split(".fa")[:-1])

    # 检查输入文件是否存在
    if not os.path.exists(pdb_file):
        lm.logger.error(f"Input file not found: {pdb_file}")
        raise BioatFileNotFoundError(f"Input file not found: {pdb_file}")

    # 检测文件格式并选择解析器
    if pdb_file.endswith(".cif") or pdb_file.endswith(".cif.gz"):
        parser = MMCIFParser(QUIET=True)
    elif pdb_file.endswith(".pdb") or pdb_file.endswith(".pdb.gz"):
        parser = PDBParser(QUIET=True)
    else:
        lm.logger.error("Only .pdb and .cif files are allowed.")
        raise BioatFileFormatError("Only .pdb and .cif files are allowed.")

    # 打开文件，支持 gzip 格式
    input_handler = (
        gzip.open(pdb_file, "rt") if pdb_file.endswith(".gz") else open(pdb_file, "rt")
    )
    structure = parser.get_structure(idx, input_handler)

    records = []

    # 遍历结构中的所有模型和链，提取序列
    for model in structure:
        for chain in model:
            lm.logger.debug(f"processing chain {chain.id}")
            # 找出当前链中所有氨基酸的编号范围
            residue_ids = [residue.id[1] for residue in chain]

            if not residue_ids:
                continue

            min_res_id, max_res_id = min(residue_ids), max(residue_ids)
            # 分别初始化序列列表，使用 "-" 占位符填充缺失位置
            protein_seq = ["-"] * (max_res_id - min_res_id + 1)
            dna_seq = ["-"] * (max_res_id - min_res_id + 1)
            rna_seq = ["-"] * (max_res_id - min_res_id + 1)
            other_seq = ["-"] * (max_res_id - min_res_id + 1)

            for residue in chain:
                residue_name = residue.resname.strip().upper()  # 标准化残基名称
                residue_id = residue.id[1] - min_res_id  # 计算相对于最小编号的偏移

                if residue_name in AMINO_ACIDS_3CODE:
                    # 蛋白质序列，填充氨基酸到相应位置
                    protein_seq[residue_id] = AMINO_ACIDS_3CODE[residue_name]
                elif residue_name in AMINO_ACIDS_3CODE_EXTEND:
                    # 蛋白质序列，填充氨基酸到相应位置
                    protein_seq[residue_id] = AMINO_ACIDS_3CODE_EXTEND[residue_name]
                elif residue_name in ("DA", "DT", "DG", "DC"):
                    # DNA 序列，填充 A, T, G, C
                    dna_seq[residue_id] = residue_name[1]
                elif residue_name in ("A", "U", "G", "C"):
                    # RNA 碱基
                    rna_seq[residue_id] = residue_name
                else:
                    # 其他分子类型，用 [] 标记
                    other_seq[residue_id] = f"[{residue_name}]"

            # 根据链类型创建 SeqRecord，并添加到 records 列表
            if any(res != "-" for res in protein_seq):
                records.append(
                    SeqRecord(
                        Seq("".join(protein_seq)),
                        id=f"{structure.id}|Chain_{chain.id}|Protein",
                        description=f"Chain {chain.id} Protein",
                    )
                )
            if any(res != "-" for res in dna_seq):
                records.append(
                    SeqRecord(
                        Seq("".join(dna_seq)),
                        id=f"{structure.id}|Chain_{chain.id}|DNA",
                        description=f"Chain {chain.id} DNA",
                    )
                )
            if any(res != "-" for res in rna_seq):
                records.append(
                    SeqRecord(
                        Seq("".join(rna_seq)),
                        id=f"{structure.id}|Chain_{chain.id}|RNA",
                        description=f"Chain {chain.id} RNA",
                    )
                )
            if any(res != "-" for res in other_seq):
                records.append(
                    SeqRecord(
                        Seq("".join(other_seq)),
                        id=f"{structure.id}|Chain_{chain.id}|UNKNOWN",
                        description=f"Chain {chain.id} Other molecules",
                    )
                )

    # 将所有链的序列保存为 FASTA 文件
    SeqIO.write(records, output_fasta, "fasta")
    input_handler.close()
    lm.logger.info(f"FASTA saved to: {output_fasta}")


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