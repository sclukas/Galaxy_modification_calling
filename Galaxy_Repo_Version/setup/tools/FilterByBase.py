#!/usr/bin/python

__author__ = 'Lukas Schmidt'
__maintainer__ = 'Lukas Schmidt'
__email__ = 'sclukas@students.uni-mainz.de'
__version__ = '1.0'

# Steps in this order: Pileup2Profile -> FilterByBase -> AnnotateProfile -> PredictModifications
# This algorithm filters out positions from a profile file based on a given reference base.

import sys

# infile format:  ref_seg pos ref_base cov pre_base rel_mism A G T C N a g t c n arrest_rate
with open(sys.argv[1], 'r') as fin, open(sys.argv[2], 'w') as fout:
    refbase = sys.argv[3]
    for line in fin:
        split_line = line.strip().split()
        if split_line[0] == 'ref_seg':  # Write header into file
            fout.write(line)
	# Write line if it matches the base and the minimum coverage requirement is fulfilled        
	if split_line[2] == refbase and int(split_line[3]) >= int(sys.argv[4]):
            fout.write(line)

# The End

