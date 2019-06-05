#!/usr/bin/python

__author__ = 'Lukas Schmidt'
__maintainer__ = 'Lukas Schmidt'
__email__ = 'sclukas@students.uni-mainz.de'
__version__ = '1.0'

# Steps in this order: pileup2profile -> FilterByBase -> AnnotateProfile -> PredictModifications -> Get_Predictions
# This algorithm filters out positions by modification type from an annotated profile.

import sys

# infile format: ref_seg pos ref_base cov pre_base arrest_rate mism_rate amism cmism gmism tmism jump_rate_total
#  predicted_mod_type
with open(sys.argv[1], 'r') as fin:
    for line in fin:
        split_line = line.strip().split()
        if 'non' in split_line[-1]:  # Write header into file
            mod_type = split_line[-1].replace('non', '')
        else:
            mod_type = split_line[-1]

with open(sys.argv[1], 'r') as fin, open(sys.argv[2], 'w') as fout_mod, open(sys.argv[3], 'w') as fout_non_mod:
    for line in fin:
        split_line = line.strip().split()
        if split_line[0] == 'ref_seg':
            fout_mod.write(line)
            fout_non_mod.write(line)
            continue

        # Write line if it matches the base and the minimum coverage requirement is fulfilled
        if split_line[-1] == mod_type:
            fout_mod.write(line)
        else:
            fout_non_mod.write(line)

# The End
