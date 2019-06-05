#!/usr/bin/python

__author__ = 'Lukas Schmidt'
__maintainer__ = 'Lukas Schmidt'
__email__ = 'sclukas@students.uni-mainz.de'
__version__ = '1.0'

# This programme allows you to filter out known modifications from profile files
# This requires an input file containing known modifications
# Format (Tab-separated):
# Reference-sequence    Position
# Example:
# tdbR00000193|Saccharomyces_cerevisiae|4932|Lys|TTT	58

import os
import sys


modification = {'m1A': 'm1A', 'm2,2G': 'm2,2G', 'm1G': 'm1G', 'm1G_m2,2G': 'm1G_m2,2G'}


########################################################################################################################

def find_pre_base(line, index):
    """
    :param line:
    :param index:
    :return:
    """
    pre_base = ''
    if index + 1 < len(lines):
        pre_base = lines[index + 1].rstrip().split('\t')[2]
    return pre_base

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

# File containing known modifications
#mod_file = open(sys.argv[1], 'r')
#mod_name = sys.argv[5]
# Base of interest
origbase = sys.argv[5]
    #'/home/akhelm/Lukas/Ergebnisse/all_tRNAs/Known_m1As/Known_m1A_sites_yeast'
#mod_file = '/home/akhelm/Lukas/Ergebnisse/all_tRNAs/Raw/HeLa_RT03_m1As'
#mod_file = '/home/akhelm/Lukas/Ergebnisse/Mouse_tRNA/Filtered_Candidates/Mouse_likely_m1As_RT03'
#mod_file = '/home/akhelm/Lukas/Ergebnisse/all_tRNAs/Known_m2,2Gs/Known_m2,2Gs_yeast'
#mod_file = '/home/akhelm/Lukas/Ergebnisse/all_tRNAs/Known_m1Gs/Known_m1Gs_sites_yeast'
#mod_file = '/home/akhelm/Lukas/Ergebnisse/E_coli/Profiles/m1G_E_coli'
#mod_file = '/home/akhelm/Lukas/Ergebnisse/all_tRNAs/Known_m2,2Gs/Known_m1G_m2,2G_yeast'

#infile = open(mod_file, 'r')

# Get positions of the confirmed modified base
header = ''
counter = 0
modified_positions = []
with open(sys.argv[1], 'r') as infile:
    for line in infile:
        if counter == 0:
            header = line
            counter += 1
            continue
        else:
            counter += 1
            modified_positions.append(line.strip().split()[0] + line.strip().split()[1])


profile = open(sys.argv[2], 'r')

# if x.endswith('_prevBasestats'):
# print(path + '/' + x[0:6] + '_m1A_positions')

outfile = open(sys.argv[3], 'w')
outfile2 = open(sys.argv[4], 'w')

header = 'ref_seg' + '\t' + 'pos' + '\t' + 'ref_base' + '\t' + 'cov' + '\t' + 'pre_base' + '\t' + 'arrest_rate' \
         + '\t' + 'mism_rate' + '\t' + 'amism' + '\t' + 'cmism' + '\t' + 'gmism' + '\t' + 'tmism' + '\t' + \
         'single_jump_rate_direct' + '\t' + 'single_jump_rate_delayed' + '\t' + 'double_jump_rate' + '\n'
outfile.write(header)
outfile2.write(header)

lines = profile.readlines()
for i in range(1, len(lines)):
    if i + 1 < len(lines):
        line1 = lines[i].rstrip().split('\t')
        line2 = lines[i + 1].rstrip().split('\t')

    # print(line1)
    split_line = line1

    if line1[0] + line1[1] in modified_positions and int(split_line[3]) >= counter:
        misms = calculate_misms(line1)
        pre_base = find_pre_base(line1, i)
        outfile.write(
            line1[0] + '\t' + line1[1] + '\t' + line1[2] + '\t' + line1[3] + '\t' + line1[4] + '\t' +
            line1[-1] + '\t' + line1[5] + '\t' + misms[0] + '\t' + misms[1] + '\t' + misms[2] + '\t' + misms[3]
            + '\t' + line1[-4] + '\t' + line1[-3] + '\t' + line1[-2] + '\n')

    elif line1[2] == origbase and int(line1[3]) >= counter:
        misms = calculate_misms(line1)
        pre_base = find_pre_base(line1, i)
        outfile2.write(
            line1[0] + '\t' + line1[1] + '\t' + line1[2] + '\t' + line1[3] + '\t' + line1[4] + '\t' +
            line1[-1] + '\t' + line1[5] + '\t' + misms[0] + '\t' + misms[1] + '\t' + misms[2] + '\t' + misms[3]
            + '\t' + line1[-4] + '\t' + line1[-3] + '\t' + line1[-2] + '\n')

outfile.close()
outfile2.close()

# The End


#/home/akhelm/Lukas/Ergebnisse/all_tRNAs/Known_m1As/Known_m1A_sites_yeast
#/home/akhelm/Lukas/Ergebnisse/S_cerevisiae/Tailing_G/BAMs_Profiles_Pileups/Profiles/All/MH1503.profile
#/home/akhelm/Lukas/Ergebnisse/S_cerevisiae/Tailing_G/BAMs_Profiles_Pileups/Profiles/All/MH1503.profile_m1A_test
#/home/akhelm/Lukas/Ergebnisse/S_cerevisiae/Tailing_G/BAMs_Profiles_Pileups/Profiles/All/MH1503.profile_non_m1A_test
#m1A
#A
