import gzip
import os

from Bio import SeqIO
from Bio.PDB import MMCIFParser, PDBParser
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from bioat.exceptions import BioatFileFormatError, BioatFileNotFoundError
from bioat.lib.libfastx import AMINO_ACIDS, AMINO_ACIDS_EXTEND
from bioat.logger import get_logger

__module_name__ = "bioat.lib.libpdb"


def pdb2fasta(pdb_file, output_fasta, log_level="DEBUG"):
    """Converts a PDB / CIF file to a FASTA file.

    This function processes the provided PDB file and extracts protein, DNA,
    RNA sequences, and other molecules appropriately to create a FASTA file.
    1. Proteins:
        The protein sequence for each chain will be extracted as Chain X Protein.
    2. DNA and RNA:
        Bases for DNA (A, T, G, C) will be saved as Chain X DNA and bases for RNA (A, U, G, C) will be saved as Chain X RNA.
    3. Other molecules:
        Any unrecognized molecules (e.g., ions, modified molecules) will be labeled as [residue] and stored as Chain X Other molecules.
    4. Multi-chain complexes:
        The program supports multi-chain structures in complexes, and the content of each chain will be recorded separately.

    Args:
        pdb_file (str): input file path.
        output_fasta (str | None, optional): output file path. If None, the output file will be named as the basename of the input file with a ".fa" extension. Defaults to None.
        log_level (str, optional): log level. Defaults to "WARNING".
    """
    logger = get_logger(
        level=log_level, module_name=__module_name__, func_name="pdb2fasta"
    )
    # 设置输出文件名
    if output_fasta is None:
        if pdb_file.endswith(".gz"):
            idx = ".".join(os.path.basename(pdb_file).split(".")[:-2])
        else:
            idx = ".".join(os.path.basename(pdb_file).split(".")[:-1])
        output_fasta = f"{idx}.fa" if idx else "output.fa"
        logger.debug(f"Output is not defined, set output to {output_fasta}")
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
            logger.debug(f"processing chain {chain.id}")
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
    logger.info(f"FASTA saved to: {output_fasta}")
