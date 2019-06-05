#!/usr/bin/python

__author__ = 'Lukas Schmidt'
__maintainer__ = 'Lukas Schmidt'
__email__ = 'sclukas@students.uni-mainz.de'
__version__ = '1.0'

# This programme writes a required positional table for the visualization of the mismatch pattern.
# The algorithm requires as input the name of a reference-sequence found in the profile file.

import sys

# input: Profile-file
ref = sys.argv[3]
with open(sys.argv[1], 'r') as fin, open(sys.argv[2], 'w') as fout:
    ref_seq = ''
    skip_header = False
    fout.write('position' + '\t' + 'refbase' + '\t' + 'coverage' + '\t' + 'mismatch_rate' + '\t' + 'a_mism' +
               '\t' + 'g_mism' + '\t' + 't_mism' + '\t' + 'c_mism' + '\t' + 'arrest_rate' + '\n')
    for line in fin:
        line1 = line.strip().split()
        if line1[0] == ref.replace('__tt__', '|'):
            fout.write(line1[1] + '\t' + line1[2] + '\t' + line1[3] + '\t' + line1[4] + '\t' + str(int(line1[5]) +
                       int(line1[10])) + '\t' + str(int(line1[6]) + int(line1[11])) + '\t' + str(int(line1[7]) +
                       int(line1[12])) + '\t' + str(int(line1[8]) + int(line1[13])) + '\t' + line1[-1] + '\n')

# The End
