#!/usr/bin/python

__author__ = 'Ralf Hauenschild, Lukas Schmidt'
__maintainer__ = 'Lukas Schmidt'
__email__ = 'sclukas@students.uni-mainz.de'
__version__ = '5.0'

# This programme extracts numerical values for each base from a pileup-file and writes them into an output-file in
# profile format. The values extracted are: mismatch, arrest, c-mismatch, g-mismatch, t-mismatch, multiple jump-rates
# and the preBase.

import sys
import copy


########################################################################################################################

def check_if_intron(pileup_line):
    """
    Function checks whether a given line from the pileup file only contains sequence skips (e.g. '<' and '>')
    :param pileup_line: Given line from the pileup file
    :return: Boolean value 'is_intron'. True if pileup_line only consists of skips, else false.
    """
    is_intron = True
    pileup_line = pileup_line.split()
    if len(pileup_line) > 4:  # Check if pileup is empty
        for symbol in pileup_line[4]:
            if symbol not in ['<', '>']:
                return False
    else:
        return False
    return is_intron


########################################################################################################################

def parse_pile(pile, splitted_line, curr_line):
    """
    Counts the instances for each base given the pileup data. Also, counts the number of reference skips in the pile,
    since '<' and '>' are included in the coverage and have to be removed.
    :param pile: Line containing the pileup data.
    :param curr_line: Boolean determining whether the line that is processed is the current line or the following line.
    :return: skip_counter: Integer containing the number of reference skips.
    """
    skip_counter = 0
    i = 0
    while i < len(pile):
        if pile[i] not in '+-^$<>':
            if curr_line:
                count_dict[pile[i]] += 1  # Counts all aligned bases within the line
            else:
                count_dict_2[pile[i]] += 1
        else:
            if pile[i] in '+-':  # Indel
                number = pile[i + 1]
                j = i + 1
                while j + 1 < len(pile) and pile[j + 1] in '1234567890':
                    number += pile[j + 1]
                    j += 1
                if pile[i] == '-' and int(number) <= 3 and curr_line:
                    (jump_dict['0'])[int(number)] += 1  # Count jumps caused by indel.
                i += 1 + int(number)
            elif pile[i] == '^':  # $-case requires no action
                i += 1
            elif pile[i] in ['<', '>']:
                skip_counter += 1
        i += 1  # Move along the pile

    if curr_line:
        count_dict[splitted_line[2]] += count_dict['.'] + count_dict[',']
    else:
        count_dict_2[splitted_line[2]] += count_dict_2['.'] + count_dict_2[',']
    '''
    try:
        count_dict[splitted_line[2]] += count_dict['.'] + count_dict[',']
    except:
        count_dict[splitted_line[2]] = count_dict['.'] + count_dict[',']
        print('Unexpected symbol in split_list[2]: ', splitted_line[2])
    '''
    return skip_counter


########################################################################################################################

def initialize_first_line(fin):
    """
    Returns the first line from the pileup file.
    :param fin: Path to the pileup file.
    :return: First line from the pileup file.
    """
    with open(fin, 'r') as infile:
        line1 = ''
        for line in infile:
            line1 = line
            break
    return line1


########################################################################################################################

def write_output_pre_base(splitted_line, coverage, pre_base, rel_mism, single_jump_rate_direct, single_jump_rate_delayed
                          , double_jump_rate, arrest):
    """
    Writes the given values to the output-file.
    """
    outfile.write(splitted_line[0] + '\t' + splitted_line[1] + '\t' + splitted_line[2] + '\t' + str(coverage) + '\t' +
                  str(pre_base) + '\t'
                  + rel_mism + '\t' + str(count_dict['A']) + '\t' + str(count_dict['G']) + '\t' + str(count_dict['T']) +
                  '\t' + str(count_dict['C']) + '\t' + str(count_dict['N']) + '\t' + str(count_dict['a']) + '\t' +
                  str(count_dict['g']) + '\t' + str(count_dict['t']) + '\t' + str(count_dict['c']) + '\t' +
                  str(count_dict['n']) + '\t' + str(single_jump_rate_direct) + '\t' + str(single_jump_rate_delayed) +
                  '\t' + str(double_jump_rate) + '\t' + str(arrest) + '\n')


########################################################################################################################

def write_output(splitted_line, coverage, rel_mism, single_jump_rate_direct, single_jump_rate_delayed, double_jump_rate,
                 arrest):
    """
    Writes the given values to the output-file.
    """
    outfile.write(splitted_line[0] + '\t' + splitted_line[1] + '\t' + splitted_line[2] + '\t' + str(coverage) + '\t'
                  + rel_mism + '\t' + str(count_dict['A']) + '\t' + str(count_dict['G']) + '\t' + str(count_dict['T']) +
                  '\t' + str(count_dict['C']) + '\t' + str(count_dict['N']) + '\t' + str(count_dict['a']) + '\t' +
                  str(count_dict['g']) + '\t' + str(count_dict['t']) + '\t' + str(count_dict['c']) + '\t' +
                  str(count_dict['n']) + '\t' + str(single_jump_rate_direct) + '\t' + str(single_jump_rate_delayed) +
                  '\t' + str(double_jump_rate) + '\t' + str(arrest) + '\n')


########################################################################################################################

def calculate_features(splitted_line_1, splitted_line_2, pair):
    """
    Calculates mismatch-rate, arrest-rate and jump-rate for the line contained in the first argument (splitted_line_1)
    and calls a function to write the relevant values into the output-file.
    :param splitted_line_1: List containing the 'above'-line from the pileup-file.
    :param splitted_line_2: List containing the 'below'-line from the pileup-file. Arrest/jump-rates for the 'above'-
    line can only be calculated, if the 'below'-line is given.
    :param pair: Boolean value. True if given lines are successive, false if they are not.
    """
    pile = splitted_line_1[4]
    skips = parse_pile(pile, splitted_line_1, True)
    if pair:
        skips_2 = parse_pile(splitted_line_2[4], splitted_line_2, False)
        cov2 = float(splitted_line_2[3]) - float(skips_2) #- count_dict_2['*']
    cov = float(splitted_line_1[3]) - float(skips)  # Subtract the number of skips from the overall coverage.

    rel_mism = '0.0'
    single_jump_rate_direct = '0.0'
    single_jump_rate_delayed = '0.0'
    double_jump_rate = '0.0'
    arrest = 0.0

    if cov > 0.0:
        try:
            rel_mism = str((cov - count_dict['.'] - count_dict[','] - count_dict['*']) / (cov - count_dict['*']))
        except:
            rel_mism = str(0.0)

        single_jump_rate_direct = (jump_dict['-1'])[1] / cov
        single_jump_rate_delayed = (jump_dict['-2'])[1] / cov
        double_jump_rate = (jump_dict['-2'])[2] / cov
    cov = cov - count_dict['*']

    # Calculate the arrest-rate
    if float(splitted_line_2[3]) > 0.0 and pair:  # If the given lines are successive --> calculate
        if cov2 == 0:
            arrest = str(0.0)
        else:
            arrest = str(((splitted_line_2[4]).count('^')) / cov2)  #float(splitted_line_2[3]))
    elif not pair:  # Else --> arrest = 0
        arrest = str(0.0)

    if pre_base_feature:
        pre_base = splitted_line_2[2]

        # does not count indels
        # Output format: ref_seg  pos  ref_base  cov pre_base rel_mism  A  G  T  C  N  a  g  t  c  n  arrest_rate
        write_output_pre_base(splitted_line_1, int(cov), pre_base, rel_mism, single_jump_rate_direct,
                              single_jump_rate_delayed, double_jump_rate, arrest)

    else:
        write_output(splitted_line_1, int(cov), rel_mism, single_jump_rate_direct, single_jump_rate_delayed,
                     double_jump_rate, arrest)

    # Move the jump-dicts by one position
    jump_dict['-3'] = copy.deepcopy(jump_dict['-2'])
    jump_dict['-2'] = copy.deepcopy(jump_dict['-1'])
    jump_dict['-1'] = copy.deepcopy(jump_dict['0'])
    jump_dict['0'] = {1: 0, 2: 0, 3: 0}

    # Reset the counter-dictionaries
    for char_x in char_string:
        count_dict[char_x] = 0
        count_dict_2[char_x] = 0


########################################################################################################################

# MAIN #

count_dict = {}
count_dict_2 = {}
char_string = 'agtcnAGTCNXRKS.,*_'

# Initialize count-dicts with zeroes
for char_x in char_string:
    count_dict[char_x] = 0
    count_dict_2[char_x] = 0

pile = ''
ref_seg = ''
pos = ''
cov = 0.0
split_list1 = []
split_list2 = []

# Dictionary containing information on the jumps over multiple positions.
jump_dict = {'0': {1: 0, 2: 0, 3: 0}, '-1': {1: 0, 2: 0, 3: 0}, '-2': {1: 0, 2: 0, 3: 0}, '-3': {1: 0, 2: 0, 3: 0}}

# Initialize the first two lines
line1 = initialize_first_line(sys.argv[1])
line2 = None

post_intron_line = ''

#s = sys.argv[2]
#pre_base_feature = True if s == 'True' else False

pre_base_feature = True

with open(sys.argv[1], 'r') as infile, open(sys.argv[2], 'w') as outfile:
    if pre_base_feature:
        outfile.write('ref_seg\tpos\tref_base\tcov\tpre_base\tmism_rate\tA\tG\tT\tC\tN\ta\tg\tt\tc\tn'
                      '\tsingle_jump_rate_direct\tsingle_jumprate_delayed\tdouble_jump_rate\tarrest_rate\n')
    else:
        outfile.write('ref_seg\tpos\tref_base\tcov\tmism_rate\tA\tG\tT\tC\tN\ta\tg\tt\tc\tn\tsingle_jump_rate_direct'
                      '\tsingle_jumprate_delayed\tdouble_jump_rate\tarrest_rate\n')
    pre_intron_line = ''
    temp_bool = True
    for line in infile:
        line2 = line
        if temp_bool:  # First line is already initialized -> skip
            temp_bool = False
            continue
        elif check_if_intron(line):  # Is the current line part of an intron?
            pre_intron_line = line1  # Save the information of the line preceding the intron.
            for line in infile:
                if check_if_intron(line):  # Check if the current line only contains '<' and '>'. If yes -> intron.
                    continue
                line1 = pre_intron_line  # If intron ends, reinitialize line1.
                post_intron_line = line  # Save the first line following an intron
                calculate_features(line1.split('\t'), line2.split('\t'), True)
                break  # No intron anymore -> end loop
        else:
            if pre_intron_line != '':  # Special case where the first line after the intron is written into the output
                line1 = post_intron_line
                line2 = line
                calculate_features(line1.split('\t'), line2.split('\t'), True)
                line1 = line2
                # Reset variables
                pre_intron_line = ''
                post_intron_line = ''
                continue
            split_list1 = line1.split('\t')
            split_list2 = line2.split('\t')
            ref_seg = split_list1[0]
            next_ref_seg = split_list2[0]

            # Check if current line still lies within the same reference segment. Additionally, check whether the two
            # positions are next to each other in the reference. If that is not the case, check if the two positions are
            # separated by an intron
            if (ref_seg == next_ref_seg and (int(split_list1[1]) == int(split_list2[1]) - 1) or (
                    pre_intron_line.split('\t')
                    [0] == next_ref_seg)):
                calculate_features(split_list1, split_list2, True)
            else:
                calculate_features(split_list1, split_list2, False)

            pre_intron_line = ''
            line1 = line2

    # Add last line to the profile
    if line2 is not None:
        calculate_features(line2.split('\t'), line2.split('\t'), True)

# The End

