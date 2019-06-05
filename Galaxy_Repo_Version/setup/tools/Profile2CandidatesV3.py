#!/usr/bin/env python

__author__ = 'Ralf Hauenschild, Lukas Schmidt'
__maintainer__ = 'Lukas Schmidt'
__email__ = 'sclukas@students.uni-mainz.de'
__version__ = '3.0'


# This programme filters out possible modified bases depending on given input values for mismatch, arrest, position and
# coverage.

import sys

########################################################################################################################


def sub_calculate_misms(my_line, misms):
    """
    Calculates the mismatch values according to the base of interest.
    :param my_line: Line from the profile
    :param misms: Overall mismatch percentage.
    :return: Mismatch values for each possible type of mismatch
    """
    if misms == 0.0:
        gmism = tmism = cmism = amism = str(0.0)
    elif origbase == 'A':
        amism = str(0.0)
        gmism = str(float(int(my_line[7]) + int(my_line[12])) / misms)
        tmism = str(float(int(my_line[8]) + int(my_line[13])) / misms)
        cmism = str(float(int(my_line[9]) + int(my_line[14])) / misms)
    elif origbase == 'C':
        amism = str(float(int(my_line[6]) + int(my_line[11])) / misms)
        gmism = str(float(int(my_line[7]) + int(my_line[12])) / misms)
        tmism = str(float(int(my_line[8]) + int(my_line[13])) / misms)
        cmism = str(0.0)
    elif origbase == 'G':
        amism = str(float(int(my_line[6]) + int(my_line[11])) / misms)
        gmism = str(0.0)
        tmism = str(float(int(my_line[8]) + int(my_line[13])) / misms)
        cmism = str(float(int(my_line[9]) + int(my_line[14])) / misms)
    else:
        amism = str(float(int(my_line[6]) + int(my_line[11])) / misms)
        gmism = str(float(int(my_line[7]) + int(my_line[12])) / misms)
        tmism = str(0.0)
        cmism = str(float(int(my_line[9]) + int(my_line[14])) / misms)
    return amism, cmism, gmism, tmism

########################################################################################################################


def calculate_misms(my_line):
    """
    Calculates the mismatch percentage
    :param my_line: List containing all information on one line of the profile
    :return: Mismatch values (float) for each base.
    """
    if origbase == 'A':
        misms = float(my_line[7]) + float(my_line[12]) + float(my_line[8]) + float(my_line[13]) + float(my_line[9]) + \
                float(my_line[14])
    elif origbase == 'C':
        misms = float(my_line[6]) + float(my_line[11]) + float(my_line[7]) + float(my_line[12]) + float(my_line[8]) + \
                float(my_line[13])
    elif origbase == 'G':
        misms = float(my_line[6]) + float(my_line[11]) + float(my_line[8]) + float(my_line[13]) + float(my_line[9]) + \
                float(my_line[14])
    elif origbase == 'T':
        misms = float(my_line[6]) + float(my_line[11]) + float(my_line[7]) + float(my_line[12]) + float(my_line[9]) + \
                float(my_line[14])
    else:
        print('No/invalid base entered.')
        sys.exit()

    return sub_calculate_misms(my_line, misms)

########################################################################################################################

# infile format: ref_seg pos ref_base cov pre_base rel_mism A G T C N a g t c n jump_rates arrest_rate

min_cov = float(sys.argv[4])
min_cov3p = float(sys.argv[5])
minrelmism = float(sys.argv[6])
minarrestrate = float(sys.argv[7])
origbase = sys.argv[8]
posmin = int(sys.argv[9])
posmax = int(sys.argv[10])
# excludepos
'''
print("min_cov", min_cov)
print("min_cov3p", min_cov3p)
print("minrelmism", minrelmism)
print("minarrestrate", minarrestrate)
print("origbase", origbase)
print("posmin", posmin)
print("posmax", posmax)
'''

# Creates a dictionary containing the reference segment and the positions of the specified base.
m1a_dict = {}
mod_file = open(sys.argv[11], 'r')  # Reference Genome

line1 = (mod_file.readline())[0:]
while len(line1) > 0:
    # Find reference segment and add it to the dictionary as a key
    # Find position of first ' ', because reference segments containing spaces are not recognized by mpileup-function
    # example: >1 dna:chromosome chromosome:GRCh38:1:151008391:151035713:1 This reference segment will be saved as '1'
    #  curr_seg = str(line1[1:line1.find(' ')])  # Starts at index 1 because '>' is at index 0 in fasta file
    if '>' in line1 and line1 not in m1a_dict:
        curr_seg = str(line1[1:(line1.find(' ') if line1.__contains__(' ') else len(line1))])
        m1a_dict[curr_seg] = {}  # Add to dictionary
    # Go to next line without '>' (a line containing the sequence or part of the sequence)
    while '>' in line1:
        line1 = (mod_file.readline())
    temp_line = ''
    while '>' not in line1:  # Go through each line as long as it is not a new reference segment
        temp_line += line1.strip()  # Line containing the complete sequence
        line1 = mod_file.readline().strip()
        if not line1:  # If the next line is empty, stop. Applies only when file ends.
            break
    # Iterate over all bases. If the base is the specified base, save the position in the dictionary.
    for i in range(len(temp_line)):
        if temp_line[i] == origbase:
            (m1a_dict[curr_seg])[i + 1] = True
mod_file.close()

########################################################################################################################

infile = open(sys.argv[1], 'r')
outfile = open(sys.argv[2], 'w')
outfile2 = open(sys.argv[3], 'w')
header = 'ref_seg\tpos\tref_base\tcov\tpre_base\tarrest_rate\tmism_rate\tamism\tcmism\tgmism\ttmism' \
         '\tsingle_jump_rate_direct\tsingle_jump_rate_delayed\tdouble_jump_rate\n'

# For candidates
outfile.write(header)
# For non-candidates
outfile2.write(header)

# Iterate through all lines in the input file ('profile'-format, contains the positions in the sequence). For every
# position, the algorithm checks whether it fulfills the given requirements (mismatch-rate, arrest-rate, etc.). If the
# position satisfies the requirements, then it is written into an output-file containing all candidates for a given RNA-
# modification. If requirements are not met, the position is written into an output-file containing all bases that are
# not modified.
lines = infile.readlines()
for i in range(1, len(lines)):
    if i + 1 < len(lines):
        line1 = lines[i].rstrip().split('\t')
        line2 = lines[i + 1].rstrip().split('\t')
    ref_seg = line1[0]
    pos = line1[1]

    '''
    print(line1)
    print("Coverage ", int(line1[3]) >= min_cov and int(line2[3]) >= min_cov3p)
    print("Mism ", float(line1[4]) >= minrelmism)
    print("Arrest ", float(line1[-1]) >= minarrestrate)
    print("Base ", line1[2] == origbase)
    print("Position ", posmin <= int(line1[1]) <= posmax)
    print("Refseq, Pos in dict ", (ref_seg in m1a_dict and int(pos) in m1a_dict[ref_seg]))
    '''
    if (int(line1[3]) >= min_cov and int(line2[3]) >= min_cov3p and float(line1[5]) >= minrelmism and float(line1[-1])
            >= minarrestrate and line1[2] == origbase and (posmin <= int(line1[1]) <= posmax) and
            (ref_seg in m1a_dict and int(pos) in m1a_dict[ref_seg])):

        cov = int(line1[3])
        misms = calculate_misms(line1)
        outfile.write(line1[0] + '\t' + line1[1] + '\t' + line1[2] + '\t' + line1[3] + '\t' + line1[4] + '\t' +
                      line1[-1] + '\t' + line1[5] + '\t' + misms[0] + '\t' + misms[1] + '\t' + misms[2] + '\t' +
                      misms[3] + '\t' + line1[-4] + '\t' + line1[-3] + '\t' + line1[-2] + '\n')

    elif int(line1[3]) >= min_cov and line1[2] == origbase and (posmin <= int(line1[1]) <= posmax):
        cov = int(line1[3])
        misms = calculate_misms(line1)
        outfile2.write(line1[0] + '\t' + line1[1] + '\t' + line1[2] + '\t' + line1[3] + '\t' + line1[4] + '\t' +
                      line1[-1] + '\t' + line1[5] + '\t' + misms[0] + '\t' + misms[1] + '\t' + misms[2] + '\t' +
                      misms[3] + '\t' + line1[-4] + '\t' + line1[-3] + '\t' + line1[-2] + '\n')

infile.close()
outfile.close()
outfile2.close()

# The End

'''
/home/akhelm/Lukas/Ergebnisse/HeLa/2510_S10_L001_merged_001.fastq.pritr_vs_fusion.sam.sorted.bam.pileup.overhang-trimmed.pileup.profileV3
/home/akhelm/Lukas/Ergebnisse/HeLa/2510_S10_L001_merged_001.fastq.pritr_vs_fusion.sam.sorted.bam.pileup.overhang-trimmed.pileup.profileV3.A_candidatesV3
/home/akhelm/Lukas/Ergebnisse/HeLa/2510_S10_L001_merged_001.fastq.pritr_vs_fusion.sam.sorted.bam.pileup.overhang-trimmed.pileup.profileV3.nonA_candidatesV3
50
50
0.1
0.2
A
0
10000
/home/akhelm/Lukas/Ergebnisse/HeLa/Reference/Fasta_Seq_tRNA_genes_Homo_sapiens.fasta
'''

