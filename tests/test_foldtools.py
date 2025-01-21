import os
import unittest

from bioat.foldtools import FoldTools


class TestFoldTools(unittest.TestCase):
    def setUp(self):
        self.fold_tools = FoldTools()
        self.ref_seq = "data/pdb/dCas9.fa"
        self.ref_pdb = "data/pdb/dCas9.pdb"
        self.cut_seq = "data/pdb/dCas9_rm_HNH.fa"
        self.cut_pdb = "data/pdb/dCas9_rm_HNH.pdb"
        self.output_fig = "/tmp/output.html"
        self.output_fasta = "/tmp/test.fa"

    def test_show_ref_cut(self):
        # Call the method directly
        self.fold_tools.show_ref_cut(
            ref_seq=self.ref_seq,
            ref_pdb=self.ref_pdb,
            cut_seq=self.cut_seq,
            cut_pdb=self.cut_pdb,
            output_fig=self.output_fig,
            ref_color="blue",
            cut_color="green",
            gap_color="red",
            ref_style="cartoon",
            cut_style="cartoon",
            gap_style="cartoon",
            log_level="DEBUG",
        )
        # Assert the output file is created
        self.assertTrue(os.path.exists(self.output_fig))

    def test_pdb2fasta(self):
        # Call the method directly
        self.fold_tools.pdb2fasta(
            input=self.ref_pdb, output=self.output_fasta, log_level="INFO"
        )
        # Assert the output file is created
        self.assertTrue(os.path.exists(self.output_fasta))

        # Additional validation (optional): Check the contents of the output FASTA file
        with open(self.output_fasta, "r") as fasta_file:
            content = fasta_file.read()
            self.assertTrue(content.startswith(">"))  # Ensure it contains FASTA format

    def tearDown(self):
        if os.path.exists(self.output_fig):
            os.remove(self.output_fig)
        if os.path.exists(self.output_fasta):
            os.remove(self.output_fasta)


if __name__ == "__main__":
    unittest.main()
