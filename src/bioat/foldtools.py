import Bio
from Bio.Seq import Seq

from bioat.lib.libpdb import get_cut2ref_aln_info, pdb2fasta, show_ref_cut
from bioat.logger import LoggerManager

lm = LoggerManager(mod_name="bioat.foldtools")


class FoldTools:
    """Folding toolbox."""

    def __init__(self):
        pass

    def show_ref_cut(
        self,
        ref_seq: str | Seq,
        ref_pdb: str | Bio.PDB.Structure.Structure,
        cut_seq: str | Seq | None = None,
        cut_pdb: str | Bio.PDB.Structure.Structure | None = None,
        ref_color: str = "blue",
        ref_map_colors: tuple[str] | None = None,
        ref_map_values: dict | None = None,
        cut_color="green",
        gap_color="red",
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
            ref_color (str, optional): Color for reference residues.
            ref_map_colors (tuple[str, str] or None, optional): ref_map_colors will be used as color bar from ref_map_colors[0] to ref_map_colors[1]. If None, do not apply color mapping. Defaults to None.
            ref_map_values (dict or None, optional): A dictionary of values for the ref color map, it will be normalized to the range of [0 - 1]. If None, all residues will be colored with the same color. e.g. ref_value_dict = {'V_0': 0.4177215189873418, 'S_1': 0.8185654008438819, 'K_2': 0.9915611814345991, 'G_3': 0.42616033755274263, ...}
            cut_color (str, optional): Color for cut residues.
            gap_color (str, optional): Color for gaps or removed residues.
            ref_style (str, optional): "stick", "sphere", "cartoon", or "line"
            cut_style (str, optional): "stick", "sphere", "cartoon", or "line"
            gap_style (str, optional): "stick", "sphere", "cartoon", or "line"
            ref_map_value_random (bool, optional): If True, ref_value_dict will be randomly generated. Defaults to False.
            output_fig (str or None, optional): Output figure file path. If None, the figure will not be saved in html format. Defaults to None.
            col (int, optional): Number of columns for the visualization. Defaults to 3.
            scale (float, optional): Scale factor for the visualization. Defaults to 1.0.
            annotate (bool, optional): Whether to annotate the visualization with labels. Defaults to True.
            text_interval (int, optional): The interval between text annotations. Defaults to 5.
            log_level (str, optional): Log level. Defaults to "WARNING".
        """
        show_ref_cut(
            ref_seq=ref_seq,
            cut_seq=cut_seq,
            ref_pdb=ref_pdb,
            cut_pdb=cut_pdb,
            ref_color=ref_color,
            ref_map_colors=ref_map_colors,
            ref_map_values=ref_map_values,
            cut_color=cut_color,
            gap_color=gap_color,
            ref_style=ref_style,
            cut_style=cut_style,
            gap_style=gap_style,
            ref_map_value_random=ref_map_value_random,
            output_fig=output_fig,
            col=col,
            scale=scale,
            annotate=annotate,
            text_interval=text_interval,
            log_level=log_level,
        )

    def pdb2fasta(self, input: str, output: str | None = None, log_level="WARNING"):
        """Converts a PDB file to a FASTA file.

        This function processes the provided PDB file and extracts protein, DNA,
        RNA sequences, and other molecules appropriately to create a FASTA file.

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
            input (str):
                Input file path.
            output (str, optional):
                Output file path. If None, the output file will be named as the
                basename of the input file with a ".fa" extension. Defaults to None.
            log_level (str, optional):
                Logging level. Defaults to "WARNING".

        Returns:
            None
        """
        pdb2fasta(pdb=input, output_fasta=output, log_level=log_level)

    def get_cut2ref_aln_info(
        self,
        ref: str | Bio.PDB.Structure.Structure,
        cut: str | Bio.PDB.Structure.Structure,
        cal_rmsd=True,
        cal_tmscore=False,
        label1="ref",
        label2="cut",
        usalign_bin: str = "usalign",
        log_level="WARNING",
    ):
        """Align cutted pdb to ref pdb using the CA atoms.
        Aligns a truncated protein structure (cut) to its full-length reference structure (ref)
        using Cα atoms and Biopython's Superimposer.

        This function:
        - Extracts all Cα atoms from `ref` and `cut`
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
            usalign_bin (str, optional): Path to the USalign binary for TM-score calculation. Default is "usalign".
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
        aln_info = get_cut2ref_aln_info(
            ref=ref,
            cut=cut,
            cal_rmsd=cal_rmsd,
            cal_tmscore=cal_tmscore,
            label1=label1,
            label2=label2,
            usalign_bin=usalign_bin,
            log_level=log_level,
        )
        for key, value in aln_info.items():
            print(f"{key}: {value}")