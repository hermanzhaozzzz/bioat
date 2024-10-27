# import argparse
# import sys

# from pyfasta import Fasta


# def get_opt():
#     """Check and parsing the opts"""
#     parser = argparse.ArgumentParser(
#         prog="seqFasta",
#         description="seqFasta: A program to get fastx file from BED file ",
#         usage="%(prog)s.py [options] -g genome -b bedfile",
#     )
#     parser.add_argument(
#         "-g",
#         "--genome",
#         nargs="?",
#         type=str,
#         default=sys.stdin,
#         help="The path of fastx genome file, the suffix must be '.fa' or '.fastx'. [Required]",
#         required=True,
#     )
#     parser.add_argument(
#         "-b",
#         "--bedfile",
#         nargs="?",
#         type=argparse.FileType("r"),
#         default=sys.stdin,
#         help="BED6 file. detail, http://genome.ucsc.edu/FAQ/FAQformat.html#format1 . [Required]",
#         required=True,
#     )
#     parser.add_argument(
#         "-n",
#         "--seqname",
#         nargs="?",
#         type=str,
#         default=None,
#         help="The prefix of fastx seq name, default is None, means to use the origin region_id of the bed",
#         required=False,
#     )
#     parser.add_argument(
#         "-f",
#         "--flank",
#         nargs="?",
#         type=int,
#         default=0,
#         help="The number of the upstream and downstream slip bases, default is 0",
#         required=False,
#     )
#     parser.add_argument(
#         "-s",
#         "--strand",
#         type=str,
#         help="Whether to fetch from strand information in bed, True/False.",
#         required=True,
#     )
#     parser.add_argument(
#         "-o",
#         "--outfile",
#         nargs="?",
#         type=argparse.FileType("w"),
#         default=sys.stdout,
#         help="Output file name for storing the results, default is stdout.",
#         required=False,
#     )
#     args = parser.parse_args()
#     return args


# def run(args):
#     genome = Fasta(args.genome)
#     bed = list(filter(lambda x: x.strip(), args.bedfile.readlines()))
#     bed_list = list(map(lambda x: x.strip().split(), bed))

#     result = map(
#         lambda i: ">{0}({1})\n{2}".format(
#             bed_list[i][3] if args.seqname is None else i + 1,
#             bed_list[i][5] if args.strand == "True" else "?",
#             genome.sequence(
#                 {
#                     "chr": bed_list[i][0],
#                     "start": int(bed_list[i][1]) - args.flank,
#                     "stop": int(bed_list[i][2]) + args.flank,
#                     "strand": bed_list[i][5] if args.strand == "True" else "+",
#                 }
#             ).upper(),
#         ),
#         range(len(bed_list)),
#     )

#     if args.outfile:
#         args.outfile.write("\n".join(result))
#     else:
#         print("".join(result))


# if __name__ == "__main__":
#     args = get_opt()
#     run(args)
