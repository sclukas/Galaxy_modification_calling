#!/usr/bin/python

__author__ = 'Lukas Schmidt'
__maintainer__ = 'Lukas Schmidt'
__email__ = 'sclukas@students.uni-mainz.de'
__version__ = '1.0'

# This algorithm calculates the relative change in the arrest and mismatch rates for two independent samples for the
# same reference and writes instances with a high level of arrest/mismatch increase or decrease into an output file.
# The first profile given contains the normal run, the second profile contains the result of the demethylated run.

import sys

base_of_interest = sys.argv[4]

relative_change_threshold = sys.argv[5]
if float(relative_change_threshold) > 1.0:
    relative_change_threshold = float(relative_change_threshold) / 100

absolut_change_threshold = sys.argv[6]
if float(absolut_change_threshold) > 1.0:
    absolut_change_threshold = float(absolut_change_threshold) / 100
minimum_coverage = sys.argv[7]

with open(sys.argv[1], 'r') as fin_1, \
     open(sys.argv[2], 'r') as fin_2, \
     open(sys.argv[3], 'w') as fout:
    # Relevant variables for relative and absolute changes
    rel_mismatch_change, absolute_mismatch_change, rel_arrest_change, absolute_arrest_change = 0.0, 0.0, 0.0, 0.0

    # Simultaneously read in both profiles
    line_file_1 = fin_1.readlines()
    line_file_2 = fin_2.readlines()
    i, j = 1, 1
    while i < len(line_file_1) and j < len(line_file_2):
        processed, assigned = False, False

        # In rare cases, the files are offset by multiple positions. This needs to be taken into account when
        # processing the information. In order to not erroneously calculate relative and absolute changes for
        # unaligning positions, the positions from both files have to match. This is tested after every iteration of
        # the loop. The boolean 'processed' ensures that nothing is falsely processed and that the loop is stopped
        # after the relevant information has been calculated.

        while not processed:
            if not assigned:  # If lines have no yet been assigned --> assign
                line1 = line_file_1[i].strip().split('\t')
                line2 = line_file_2[j].strip().split('\t')
            if line1[0] + line1[1] == line2[0] + line2[1]:  # In case both lines contain the same sequence-position
                mism_orig = float(line1[5])
                mism_new = float(line2[5])
                arrest_orig = float(line1[-1])
                arrest_new = float(line2[-1])
                try:
                    rel_mismatch_change = (float(mism_new) - float(mism_orig)) / float(mism_orig)
                    rel_arrest_change = (float(arrest_new) - float(arrest_orig)) / float(arrest_orig)
                    absolute_mismatch_change = abs(float(mism_new) - float(mism_orig))
                    absolute_arrest_change = abs(float(arrest_new) - float(arrest_orig))
                except:  # Avoid division by zero
                    pass

                if rel_mismatch_change <= float(relative_change_threshold) * -1 and absolute_mismatch_change >= float(absolut_change_threshold) and line1[2] == base_of_interest \
                        and int(line2[3]) >= int(minimum_coverage) and int(line1[3]) >= int(minimum_coverage):
                    fout.write('\t'.join(line2) + '\n')
                i += 1
                j += 1
                processed = True  # Set to True so that loop ends

            # These two if-clauses make sure that the lines of the original profile and the demethylation profile
            # contain the same position within the sequence
            elif int(line1[1]) > int(line2[1]):
                while int(line1[1]) != int(line2[1]):
                    j += 1
                    try:
                        line2 = line_file_2[j].strip().split()
                    except:
                        break
                assigned = True  # Change to True since both lines now contain the same position of the sequence
            else:
                while int(line1[1]) != int(line2[1]):
                    i += 1
                    try:
                        line1 = line_file_1[i].strip().split('\t')
                    except:
                        break
                assigned = True


# The End

