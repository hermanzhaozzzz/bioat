import gzip
import os
from io import StringIO

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
from bioat.lib.libpath import is_path
from bioat.logger import LoggerManager

lm = LoggerManager(mod_name="bioat.lib.libpdb")
lm.set_level("DEBUG")


def structure_to_string(structure: Bio.PDB.Structure.Structure):
    """
    将 Bio.PDB.Structure.Structure 对象转换为 PDB 格式的字符串。

    参数:
        structure (Bio.PDB.Structure.Structure): 待转换的结构对象。

    返回:
        str: PDB 格式的字符串。
    """
    pdb_io = PDBIO()
    pdb_io.set_structure(structure)

    # 使用 StringIO 保存到内存中的字符串
    pdb_string = StringIO()
    pdb_io.save(pdb_string)

    # 返回字符串
    return pdb_string.getvalue()


# 对齐函数
def align_cut2ref(
    ref: str | Bio.PDB.Structure.Structure,
    cut: str | Bio.PDB.Structure.Structure,
    gap_indices: list[int],
    label1="ref",
    label2="cut",
    log_level="WARNING",
):
    """Align pdb to ref using the CA atoms.

    Args:
        ref (str or Bio.PDB.Structure.Structure): Path to the PDB file of the ref structure.
        pdb (str or Bio.PDB.Structure.Structure): Path to the PDB file of the pdb structure.
        gap_indices (list[int]): List of indices of gaps in the cut structure.
        label1 (str, optional): ref label. Defaults to "structure1".
        label2 (str, optional): pdb label. Defaults to "structure2".
        log_level (str, optional): Log level. Defaults to "WARNING".

    Returns:
        Bio.PDB.Structure.Structure: The aligned pdb structure.
    """
    lm.set_names(func_name="align_pdb2ref")
    lm.set_level(log_level)

    # 读取 PDB 文件
    parser = PDBParser(QUIET=True)

    ref = parser.get_structure(label1, ref) if isinstance(ref, str) else ref
    cut = parser.get_structure(label2, cut) if isinstance(cut, str) else cut

    assert (
        isinstance(ref, Bio.PDB.Structure.Structure)
        and isinstance(cut, Bio.PDB.Structure.Structure)
    ), "Invalid input format. ref and cut must be both str or Bio.PDB.Structure.Structure"

    # 选择对齐的原子 (e.g., CA)
    atoms1 = [atom for atom in ref.get_atoms() if atom.get_name() == "CA"]
    atoms1 = [atoms1[i] for i in range(len(atoms1)) if i not in gap_indices]
    atoms2 = [atom for atom in cut.get_atoms() if atom.get_name() == "CA"]

    assert len(atoms1) == len(
        atoms2
    ), "The number of CA atoms in ref and cut are not equal."
    # 使用 Biopython 的 Superimposer 进行对齐
    super_imposer = Superimposer()
    super_imposer.set_atoms(atoms1, atoms2)
    super_imposer.apply(cut.get_atoms())  # 修改 pdb 的坐标
    rmsd = super_imposer.rms

    return {
        label1: ref,
        label2: cut,
        "RMSD": rmsd,
    }


def show_ref_cut(
    ref_seq: str | Seq,
    ref_pdb: str | Bio.PDB.Structure.Structure,
    cut_seq: str | Seq | None = None,
    cut_pdb: str | Bio.PDB.Structure.Structure | None = None,
    ref_color: str = "blue",
    ref_color_base: str | None = None,
    ref_value_dict: dict | None = None,
    cut_color="green",
    gap_color="red",
    ref_style="cartoon",
    cut_style="cartoon",
    gap_style="cartoon",
    ref_value_random: bool = False,
    output_fig: str | None = None,
    log_level="WARNING",
):
    """
    Visualizes the alignment of sequences and highlights changes in PDB structures using py3Dmol.

    Args:
        ref_seq (str or Seq): Amino acid sequence content for the ref protein.
        ref_pdb (str or Bio.PDB.Structure.Structure): Path to the PDB file of the reference structure.
        cut_seq (str, Seq or None, optional): Amino acid sequence content for the cut protein.
        cut_pdb (str, Bio.PDB.Structure.Structure or None, optional): Path to the PDB file of the cut structure.
        ref_color (str, optional): Color for reference residues.
        ref_color_base (str or None, optional): ref_color_base will be used as base color, and ref_color will be used as target color. If None, do not apply color mapping. Defaults to None.
        ref_value_dict (dict or None, optional): A dictionary of values for the ref color map, it will be normalized to the range of [0 - 1]. If None, all residues will be colored with the same color. e.g. ref_value_dict = {'V_0': 0.4177215189873418, 'S_1': 0.8185654008438819, 'K_2': 0.9915611814345991, 'G_3': 0.42616033755274263, ...}
        cut_color (str, optional): Color for cut residues.
        gap_color (str, optional): Color for gaps or removed residues.
        ref_style (str, optional): "stick", "sphere", "cartoon", or "line"
        cut_style (str, optional): "stick", "sphere", "cartoon", or "line"
        gap_style (str, optional): "stick", "sphere", "cartoon", or "line"
        ref_value_random (bool, optional): If True, ref_value_dict will be randomly generated. Defaults to False.
        output_fig (str or None, optional): Output figure file path. If None, the figure will not be saved in html format. Defaults to None.
        log_level (str, optional): Log level. Defaults to "WARNING".
    """
    lm.set_names(func_name="show_ref_cut")
    lm.set_level(log_level)

    # fix color mapping
    assert (ref_color_base is None and ref_value_dict is None) or (
        isinstance(ref_color_base, str) and isinstance(ref_value_dict, dict)
    ), "ref_color_base and ref_value_dict must be both None or both not None"

    if ref_value_dict:
        values = []
        res_names = []
        res_idxes = []
        for k, v in ref_value_dict.items():
            res_name, res_idx = k.split("_")
            assert res_name in AMINO_ACIDS_1CODE, f"Invalid residue name: {res_name}"
            res_idx = int(res_idx)
            values.append(v)
            res_names.append(res_name)
            res_idxes.append(res_idx)
        if not ref_value_random:
            colors = map_colors_between_two(ref_color_base, ref_color, values)
        else:
            colors = map_colors_between_two(
                ref_color_base, ref_color, np.random.rand(len(values))
            )

    # load ref_seq
    def load_seq(seq, label=None):
        if isinstance(seq, str) and is_path(seq):
            seq = next(SeqIO.parse(seq, "fasta"))
            label = seq.id if label is None else label
            seq = seq.seq
        elif isinstance(seq, Seq):
            label = "no_name" if label is None else label
        else:
            raise BioatInvalidParameterError(f"Invalid sequence format. seq: {seq}")
        return seq, label

    ref_seq, ref_label = load_seq(ref_seq, "ref")

    # load ref_pdb
    parser = PDBParser(QUIET=True)
    ref_pdb = (
        ref_pdb
        if isinstance(ref_pdb, Bio.PDB.Structure.Structure)
        else parser.get_structure(ref_label, ref_pdb)
    )

    gap_indices, rmsd = [], None

    if cut_seq:
        # 执行Sequence Alignment并得到gaps的位置
        cut_seq, cut_label = load_seq(cut_seq, "cut")
        # Perform alignment
        aligner = instantiate_pairwise_aligner(
            scoring_match=5,
            penalty_mismatch=-1,
            penalty_gap_open=-1,
            penalty_gap_extension=0,
            penalty_query_left_gap_score=0,
            penalty_query_right_gap_score=0,
            log_level="WARNING",
        )
        alignment = get_best_alignment(
            seq_a=ref_seq, seq_b=cut_seq, aligner=aligner, consider_strand=False
        )
        alignment_dict = get_aligned_seq(alignment, reverse=False, letn_match=False)
        gap_indices = [
            i for i, n in enumerate(alignment_dict["target_seq"]) if n == "-"
        ]

        if cut_pdb:
            # 执行PDB Alignment并得对齐后的PDB文件和RMSD
            res = align_cut2ref(ref_pdb, cut_pdb, gap_indices, ref_label, cut_label)
            ref_pdb = res[ref_label]
            cut_pdb = res[cut_label]
            rmsd = res["RMSD"]
    else:
        lm.logger.info(
            "Cut sequence is not provided, skip alignment and annotation for gaps, and just show ref visualization."
        )
        cut_seq, cut_label, cut_pdb = None, None, None

    assert isinstance(
        ref_pdb, Bio.PDB.Structure.Structure
    ), "Invalid input format. ref_pdb must be Bio.PDB.Structure"

    if cut_pdb:
        assert isinstance(
            cut_pdb, Bio.PDB.Structure.Structure
        ), "Invalid input format. cut_pdb must be Bio.PDB.Structure"

    # show ref_pdb
    view = py3Dmol.view(width=1000, height=1000)
    label_ypos = -15
    view.addModel(structure_to_string(ref_pdb), "pdb")
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
    )
    label_ypos -= 5
    view.addLabel(
        f"[{ref_label}]: {ref_color}",
        {
            "fontColor": ref_color,
            "fontSize": 14,
            "position": {"x": -30, "y": label_ypos, "z": -30},
        },
    )
    # 标记ref 值的颜色
    if ref_color_base is not None:
        for i in range(len(colors)):
            idx = res_idxes[i]
            res = res_names[i]
            color = colors[i]
            view.addStyle(
                {"model": 0, "resi": idx + 1, "resn": AMINO_ACIDS_1CODE[res]},
                {ref_style: {"color": color, "radius": 0.5}},
            )
        label_ypos -= 5
        view.addLabel(
            f"[{ref_label}] values: {ref_color_base} (0) to {ref_color} (1)",
            {
                "fontColor": ref_color,
                "fontSize": 14,
                "position": {"x": -30, "y": label_ypos, "z": -30},
            },
        )
    # 标记gaps的颜色
    if gap_indices:
        for idx, res in enumerate(alignment_dict["reference_seq"]):
            # 标记 gap位置的颜色
            if idx in gap_indices:
                view.addStyle(
                    {"model": 0, "resi": idx + 1, "resn": AMINO_ACIDS_1CODE[res]},
                    {gap_style: {"color": gap_color, "radius": 0.5}},
                )
        label_ypos -= 5
        view.addLabel(
            f"gaps: {gap_color}",
            {
                "fontColor": gap_color,
                "fontSize": 14,
                "position": {"x": -30, "y": label_ypos, "z": -30},
            },
        )
    if cut_pdb and rmsd:
        view.addModel(structure_to_string(cut_pdb), "pdb")
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
        )
        label_ypos -= 5
        view.addLabel(
            f"{cut_label}: {cut_color}",
            {
                "fontColor": cut_color,
                "fontSize": 14,
                "position": {"x": -30, "y": label_ypos, "z": -30},
            },
        )
        label_ypos -= 5
        view.addLabel(
            f"RMSD: {rmsd:.2f}",
            {
                "fontColor": "white",
                "fontSize": 14,
                "position": {"x": -30, "y": label_ypos, "z": -30},
            },
        )

    # 聚焦显示
    # Finalize the view
    view.zoomTo()

    if output_fig:
        assert output_fig.endswith(".html"), "Output file must have a '.html' extension"

        # 保存为 SVG 文件
        with open(output_fig, "w") as f:
            f.write(view._make_html())
        return view.show()


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
        raise BioatFileNotFoundError(f"Input file not found: {pdb_file}")

    # 检测文件格式并选择解析器
    if pdb_file.endswith(".cif") or pdb_file.endswith(".cif.gz"):
        parser = MMCIFParser(QUIET=True)
    elif pdb_file.endswith(".pdb") or pdb_file.endswith(".pdb.gz"):
        parser = PDBParser(QUIET=True)
    else:
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
    # ref_seq
    pass
