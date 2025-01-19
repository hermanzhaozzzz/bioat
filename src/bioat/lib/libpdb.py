import gzip
import os

import py3Dmol
from Bio import SeqIO
from Bio.PDB import MMCIFParser, PDBParser
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from bioat.exceptions import BioatFileFormatError, BioatFileNotFoundError
from bioat.lib.libalignment import (
    get_aligned_seq,
    get_best_alignment,
    instantiate_pairwise_aligner,
)
from bioat.lib.libfastx import AMINO_ACIDS, AMINO_ACIDS_EXTEND
from bioat.logger import LoggerManager

lm = LoggerManager(mod_name="bioat.lib.libpdb")


def show_ref_cut(
    ref_seq: str,
    cut_seq: str,
    ref_pdb: str,
    cut_pdb: str | None = None,
    ref_color="blue",
    cut_color="green",
    gap_color="red",
    log_level="WARNING",
):
    """
    Visualizes the alignment of sequences and highlights changes in PDB structures using py3Dmol.

    Args:
        ref_seq (str): Amino acid sequence content for the ref protein.
        cut_seq (str): Amino acid sequence content for the cut protein
        ref_pdb (str): Path to the PDB file of the reference structure.
        cut_pdb (str): Path to the PDB file of the cut structure.
        gap_color (str): Color for gaps or removed residues.
        ref_color (str): Color for reference residues.
        cut_color (str): Color for cut residues.
    """
    lm.set_names(func_name="show_ref_cut")
    lm.set_level(log_level)

    def highlight_gaps_py3dmol(viewer, pdb_file, gaps, color, label=""):
        """
        Highlights specific residue ranges in a PDB file with a given color.

        Args:
            viewer: py3Dmol.View instance.
            pdb_file (str): Path to the PDB file.
            gaps (list of tuples): List of (start, end) residue ranges to highlight.
            color (str): Color for highlighting.
            label (str): Label for the highlighted residues.
        """
        with open(pdb_file, "r") as pdb:
            pdb_data = pdb.read()
        viewer.addModel(pdb_data, "pdb")

        for start, end in gaps:
            viewer.setStyle(
                {"resi": list(range(start, end))},
                {"stick": {"color": color}},
            )
            if label:
                viewer.addLabel(
                    label,
                    {"fontColor": color, "backgroundColor": "white"},
                    {"resi": start},
                )

    # Perform alignment
    aligner = instantiate_pairwise_aligner(
        scoring_match=5,
        penalty_mismatch=-1,
        penalty_gap_open=-1,
        penalty_gap_extension=0,
        penalty_query_left_gap_score=0,
        penalty_query_right_gap_score=0,
    )
    alignment = get_best_alignment(
        seq_a=ref_seq, seq_b=cut_seq, aligner=aligner, consider_strand=False
    )
    alignment_dict = get_aligned_seq(alignment, reverse=False, letn_match=False)

    # Create mapping of aligned indices and gaps
    ref_indices = alignment.aligned[0]
    cut_indices = alignment.aligned[1]
    gap_ranges = [
        (ref_indices[i][1], cut_indices[i][0])
        for i in range(len(ref_indices))
        if cut_indices[i][0] > ref_indices[i][1]
    ]

    # Visualize PDB structures using py3Dmol
    viewer = py3Dmol.view(width=400, height=300)
    viewer.setStyle({"cartoon": {"color": "spectrum"}})

    # Add reference and cut structures
    highlight_gaps_py3dmol(viewer, ref_pdb, [], ref_color, label="Reference")
    highlight_gaps_py3dmol(viewer, cut_pdb, [], cut_color, label="Cut")
    # Highlight gaps
    highlight_gaps_py3dmol(
        viewer, ref_pdb, gap_ranges, gap_color, label="Deleted Residues"
    )

    # Finalize the view
    viewer.zoomTo()
    return viewer.show()


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

                if residue_name in AMINO_ACIDS:
                    # 蛋白质序列，填充氨基酸到相应位置
                    protein_seq[residue_id] = AMINO_ACIDS[residue_name]
                elif residue_name in AMINO_ACIDS_EXTEND:
                    # 蛋白质序列，填充氨基酸到相应位置
                    protein_seq[residue_id] = AMINO_ACIDS_EXTEND[residue_name]
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