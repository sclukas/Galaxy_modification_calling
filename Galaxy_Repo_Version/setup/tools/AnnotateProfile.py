#!/usr/bin/python

__author__ = 'Lukas Schmidt'
__maintainer__ = 'Lukas Schmidt'
__email__ = 'sclukas@students.uni-mainz.de'
__version__ = '1.0'

# Steps in this order: pileup2profile -> FilterByBase -> AnnotateProfile -> PredictModifications
# The algorithm prepares an input in profile format for prediction using data of a random forest.

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

# Infile format: ref_seg pos ref_base cov pre_base rel_mism A G T C N a g t c n single_jump_rate_direct single_jumprate_delayed
#  double_jump_rate arrest_rate
origbase = sys.argv[3]
with open(sys.argv[1], 'r') as infile, open(sys.argv[2], 'w') as outfile:
    modname = 'x'
    outfile.write('arrest_rate\tmism_rate\tpre_base\tamism\tcmism\tgmism\ttmism\tjump_rate_total\tmod_type\n')
    for line in infile:
        split_list = line.strip().split('\t')
        print(split_list)
        if split_list[0] == 'ref_seg':  # First line of input-file contains the header, which is not written to output
            continue
        else:
            misms = calculate_misms(split_list)
            outfile.write(
                split_list[-1] + '\t' + split_list[5] + '\t' + split_list[4] + '\t' + misms[0] + '\t' + misms[1] + '\t'
                + misms[2] + '\t' + misms[3] + '\t' + str(float(split_list[-4]) + float(split_list[-3]) +
                float(split_list[-2])) + '\t' + "'" + modname + "'" + '\n')

# The End

