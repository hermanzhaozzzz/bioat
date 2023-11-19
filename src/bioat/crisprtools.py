class CrisprTools:
    """CRISPR mining toolbox."""
    def crispr_repeats_finder(self):
        pass
    def cas12_finder(self):
        pass
    def cas13_finder(self):
        pass



    # https://github.com/yszhou2016/Cas13/blob/master/0.Cas-Finder/1.Cas13-Finder.pl

    # import re
    #
    # file = sys.argv[1]  # Get filename from command line arguments
    # file2 = re.sub(r'.fasta$', '', file)  # Remove .fasta extension from file2
    #
    # with open(file, 'r') as fasta_file:
    #     with open(file2 + '.fa', 'w') as fa_file:
    #         sequence_flag = False
    #         for line in fasta_file:
    #             line = line.rstrip()
    #             if line.startswith('>'):
    #                 if sequence_flag:
    #                     fa_file.write('\n')
    #                 fa_file.write(line + '\n')
    #                 sequence_flag = True
    #             else:
    #                 fa_file.write(line)
    #         fa_file.write('\n')
    #
    # with open(file2 + '.fa', 'r') as fa_file:
    #     with open(file2 + '.RxxxxH.fa', 'w') as output_file:
    #         for line in fa_file:
    #             line2 = fa_file.readline().rstrip()
    #             length = len(line2)
    #             sequence = list(line2)
    #             RNH_flag = False
    #             RNQ_flag = False
    #             start = 0
    #             end = 0
    #             for i in range(length):
    #                 if sequence[i] == 'R' and sequence[i + 1] == 'H' and sequence[i + 5] == 'H':
    #                     if not RNH_flag:
    #                         start = i
    #                         RNH_flag = True
    #                     end = i
    #                 elif sequence[i] == 'R' and sequence[i + 1] == 'N' and sequence[i + 5] == 'H':
    #                     RNQ_flag = True
    #                 elif sequence[i] == 'R' and sequence[i + 1] == 'Q' and sequence[i + 5] == 'H':
    #                     RNQ_flag = True
    #             if RNH_flag and RNQ_flag and start * 2 < length and end * 2 > length:
    #                 output_file.write('>' + line + '_L' + str(length) + '\n' + line2 + '\n')