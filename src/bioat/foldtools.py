import inspect

import Bio
from Bio.Seq import Seq

from bioat.lib.libpdb import pdb2fasta, show_ref_cut
from bioat.logger import LoggerManager

lm = LoggerManager(mod_name="bioat.foldtools")


class FoldTools:
    """Folding toolbox."""

    lm.set_names(cls_name="FoldTools")

    def __init__(self):
        pass

    def show_ref_cut(
        self,
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
        **kwargs,
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
        print(f"Calling show_ref_cut with {kwargs}")
        print(f"Actual show_ref_cut module: {inspect.getmodule(show_ref_cut)}")
        lm.set_names(func_name="show_ref_cut")
        lm.set_level(log_level)
        show_ref_cut(
            ref_seq=ref_seq,
            cut_seq=cut_seq,
            ref_pdb=ref_pdb,
            cut_pdb=cut_pdb,
            ref_color=ref_color,
            ref_color_base=ref_color_base,
            ref_value_dict=ref_value_dict,
            cut_color=cut_color,
            gap_color=gap_color,
            ref_style=ref_style,
            cut_style=cut_style,
            gap_style=gap_style,
            ref_value_random=ref_value_random,
            output_fig=output_fig,
            log_level=log_level,
        )
        lm.logger.debug(
            f"""\
Params:
-------
ref_seq: {ref_seq}
ref_pdb: {ref_pdb}
cut_seq: {cut_seq}
cut_pdb: {cut_pdb}
ref_color: {ref_color}
ref_color_base: {ref_color_base}
cut_color: {cut_color}
gap_color: {gap_color}
ref_style: {ref_style}
cut_style: {cut_style}
gap_style: {gap_style}
ref_value_random: {ref_value_random},
output_fig : {output_fig},
log_level: {log_level}"""
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
        lm.set_names(func_name="pdb2fasta")
        lm.set_level(log_level)
        lm.logger.debug(
            f"""\
Params:
-------
input: {input}
output: {output}
log_level: {log_level}"""
        )
        pdb2fasta(pdb_file=input, output_fasta=output, log_level=log_level)
