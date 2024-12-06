from bioat.lib.libpdb import pdb2fasta
from bioat.logger import LoggerManager

lm = LoggerManager(mod_name="bioat.foldtools")


class FoldTools:
    """Folding toolbox."""

    lm.set_names(cls_name="FoldTools")

    def __init__(self):
        pass

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
        lm.logger.debug(f"""\
Params:
-------
input: {input}
output: {output}
log_level: {log_level}""")
        pdb2fasta(pdb_file=input, output_fasta=output, log_level=log_level)
