#!/usr/bin/python

__author__ = 'Lukas Schmidt'
__maintainer__ = 'Lukas Schmidt'
__email__ = 'sclukas@students.uni-mainz.de'
__version__ = '2.0'

# This algorithm prepares a file of given modified RNA bases or non-modified bases for machine learning (random forest).
# Since the input file in profile-format contains data which cannot be used during machine learning, the irrelevant
# features have to be filtered out. This script removes the name, position and reference base information from the input
# and also adds the modification name (given as argument). The 'annotated' bases will then be written into the output.

import sys

# Infile format: ref_seg pos refbase cov pre_base relmism arrestrate gmism tmism cmism single_jump_rate_direct
# single_jump_rate_delayed double_jump_rate
with open(sys.argv[1], 'r') as infile, open(sys.argv[2], 'w') as outfile:
    modname = sys.argv[3]  # 'm1A'/'non-m1A'
    outfile.write('arrest_rate\tmism_rate\tpre_base\tamism\tcmism\tgmism\ttmism\tjump_rate_total\tmod_type\n')
    for line in infile:
        split_list = line.strip().split('\t')
        if split_list[0] == 'ref_seg':  # First line of input-file contains the header, which is not written to output
            continue
        else:
            outfile.write(
                split_list[5] + '\t' + split_list[6] + '\t' + split_list[4] + '\t' + split_list[7] + '\t' +
		        split_list[8] + '\t' + split_list[9] + '\t' + split_list[10] + '\t' + str(float(split_list[11])
		        + float(split_list[12]) + float(split_list[13])) + '\t' + "'" + modname + "'" + '\n')

# The End

