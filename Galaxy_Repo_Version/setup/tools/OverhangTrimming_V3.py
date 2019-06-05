#!/usr/bin/python

__author__ = 'Ralf Hauenschild, Lukas Schmidt'
__maintainer__ = 'Lukas Schmidt'
__email__ = 'sclukas@students.uni-mainz.de'
__version__ = '3.0'

# This programme trims ligation overhang relics from files in pileup-format. It works in two steps:
# First step: The programme calculates the statistical distribution of elongated overhangs.
# Example: Case 1
# 5' - GTAACCTGC - 3' Read
# 5' - ATAACCTGC - 3' Reference
# The 5' G in the read is very likely to be an overhang relic, since it is located at the end of the read and does not
# match with the base in the reference.
#
# Case 2
# GGGTAGCCT Read
# GACTAGCCT Reference
# Here, the 5' G in the read matches with the base in the reference. However, since both Gs next to it do not match with
# the reference, there is a high likelihood that all three Gs stem from ligation overhang and can be trimmed. The
# algorithm checks the probability of encountering these kind of overhang relics.
# Typically, the tailing efficiency of G determined by this tool must be much lower than the other ones, because much
# more reads generally start with G than with other bases. This lowers the resulting percentage for G. As this statistic
# is only applied on GG spots for the second G, we expect a much higher fraction of true Gs on the first G of GG,
# because of course some reads simply start untailed at this G und should not be trimmed.
#
# Second step: The programme cuts the overhangs from each line in the pileup file and writes the trimmed outputs to an
# output file. The trimming is dependent on the probabilities calculated in the first step.


# ATTENTION! This program requires a minimum read length of 5 to work properly

# IMPORTANT: Python offers multiple functions to read input data. Most notably readline() and readlines(). However,both
# functions read the whole input file into the working memory. Since files in pileup-format can be very large for an
# alignment to a whole genome or transcriptome, these functions can be very slow and can cause the computer to crash or
# slow it down so much that working becomes impossible. Therefore, in this algorithm, the input is read line by line. On
# one hand, this makes it more difficult to access individual lines by index or to manipulate loops. On the other hand,
# it makes the algorithm more resource efficient and much faster for large input files.

import sys
sys.path.append('/usr/local/lib/python2.7/dist-packages')
import copy
import numpy
from scipy.stats import sem
import datetime

a = datetime.datetime.now()


########################################################################################################################

def determine_overhang_prolongation_probabilities(fin):
    """
    Determines statistical data of positions which are likely overhangs.
    :param fin: Pileup file.
    :return: Array with statistical values (standard deviation, etc.).
    """
    base_dependent_prolongation_stats = {'A': [], 'G': [], 'T': [], 'C': [], 'N': []}

    potential_trimming_indices = []
    immediate_trimming_indices = []
    ending_indices = []
    starting_indices = []
    last_ref_seg = ''
    last_ref_base = ''
    with open(fin, 'r') as infile:
        for line in infile:
            split_list = line.split('\t')

            ref_seg = split_list[0]
            ref_base = split_list[2].capitalize()
            pileup = split_list[4]
            qual = split_list[5]

            if ref_seg == tailing_base_up:
                g_series += 1  # Value containing the number of successive Gs
            else:
                g_series = 0

            if ref_seg != last_ref_seg:
                potential_trimming_indices = []
                immediate_trimming_indices = []
                ending_indices = []
                starting_indices = []

            seq_words = []
            qual_words = []
            i = 0
            qual_index = 0
            word_index = 0
            while i < len(pileup):
                word = ''
                qual_word = ''
                if pileup[i] == '^':  # Read starts. This symbol is always followed by an ASCII-symbol for the quality
                    # score and the first symbol for the alignment
                    word = pileup[i:i + 3]  # Add '^', the quality score and the first symbol to the word
                    qual_word += qual[qual_index]
                    qual_index += 1
                    i += 3
                    starting_indices.append(word_index)  # Save the index for the word.
                    if i < len(pileup):
                        if pileup[i] in ['+', '-']:  # Indel
                            word += pileup[i]
                            i += 1
                            number = 0
                            while pileup[i] in ['0', '1', '2', '3', '4', '5', '6', '7', '8', '9']:  # Length of indel
                                number *= 10
                                number += int(pileup[i])
                                i += 1
                            indel = pileup[i:i + number]
                            word += str(number) + indel
                            i += number
                        if i < len(pileup):
                            if pileup[i] == '$':  # One-base read ends. Save index of position where read ends.
                                word += pileup[i]
                                i += 1
                                ending_indices.append(word_index)

                elif pileup[i] in ['*', ',', '.', 'A', 'G', 'T', 'C', 'N', 'a', 'g', 't', 'c', 'n', '<', '>']:
                    word += pileup[i]
                    qual_word += qual[qual_index]
                    qual_index += 1
                    i += 1
                    if i < len(pileup):
                        if pileup[i] == '$':  # Read ends. Save index of position where read ends.
                            word += pileup[i]
                            i += 1
                            ending_indices.append(word_index)
                        elif pileup[i] in ['+', '-']:  # Insertion / Deletion
                            word += pileup[i]
                            i += 1
                            number = 0
                            while pileup[i] in ['0', '1', '2', '3', '4', '5', '6', '7', '8', '9']:
                                number *= 10
                                number += int(pileup[i])
                                i += 1
                            indel = pileup[i:i + number]
                            word += str(number) + indel
                            i += number
                else:
                    i += 1

                seq_words.append(word)
                qual_words.append(qual_word)
                word_index += 1

            curr_prolongationcounts = {'A': 0.0000001, 'G': 0.0000001, 'T': 0.0000001, 'C': 0.0000001, 'N': 0.0000001}
            curr_startstats = {'A': 0.0000001, 'G': 0.0000001, 'T': 0.0000001, 'C': 0.0000001, 'N': 0.0000001}

            for i in starting_indices:
                word = seq_words[i]
                if ref_base != tailing_base_up:
                    if len(word) == 3:
                        # reference must be non-G to increase overhang evidence
                        if word[2] in [tailing_base_up, tailing_base_low]:
                            immediate_trimming_indices.append(i)
                if word[2] in ['.', ',']:
                    curr_startstats[ref_base] += 1
                else:
                    curr_startstats[(word[2]).upper()] += 1
                # This case is not defined in this statistics script due to assumed minimum read length
                if word.endswith('$'):
                    print('Word ends with $: ', word)

            pti_deletor = []
            for pti in potential_trimming_indices:  # Overhang trimming continuation
                if len(seq_words[pti]) == 1:
                    if seq_words[pti] not in [tailing_base_up, tailing_base_low] and not (
                            seq_words[pti] in ['.', ','] and ref_base == tailing_base_up):  # read base is non-G
                        try:
                            curr_prolongationcounts[(seq_words[pti]).upper()] += 1
                        except:
                            curr_prolongationcounts[ref_base] += 1
                        pti_deletor.append(pti)

                    elif ref_base == tailing_base_up and seq_words[pti] in ['.', ',']:
                        curr_prolongationcounts[tailing_base_up] += 1
                        pti_deletor.append(pti)
                    else:
                        pass  # G takes part in the statistics only if ref_base == G. Remember this, when summing up
                        # rb-spec tailing efficiencies.

            for pti in pti_deletor:
                potential_trimming_indices.remove(pti)

            for rb in bases:      # ['A', 'T', 'C', 'N']:
                for i in range(0, int(curr_startstats[rb] + curr_prolongationcounts[rb])):
                    (base_dependent_prolongation_stats[rb]).append(
                     float(curr_prolongationcounts[rb]) / (float(curr_startstats[rb] + curr_prolongationcounts[rb])))
            if ref_base == tailing_base_up and last_ref_base != tailing_base_up:
                for i in range(0, int(curr_startstats[tailing_base_up] + curr_prolongationcounts[tailing_base_up])):
                    (base_dependent_prolongation_stats[tailing_base_up]).append(
                     float(curr_prolongationcounts[tailing_base_up]) / (float(curr_startstats[tailing_base_up] + curr_prolongationcounts[tailing_base_up])))

            # Here, the latest starting reads are trimmed. The indices of interest were collected before trimming old
            # potentially interesting indices. Therefore, the potentially gut-shot pile requires correction of top-level
            # indices of starting reads.
            for pti in immediate_trimming_indices:
                potential_trimming_indices.append(pti)  # Only relevant from next round on

            for i in range(0, len(potential_trimming_indices)):
                # Adjustment of potential_trimming_indices to upstream ending reads
                subtractor = 0
                for j in range(0, len(ending_indices)):
                    if potential_trimming_indices[i] > ending_indices[j]:
                        subtractor += 1
                potential_trimming_indices[i] -= subtractor

            immediate_trimming_indices = []
            ending_indices = []  # Old upstream ending reads are irrelevant now
            starting_indices = []

            last_ref_seg = ref_seg
            last_ref_base = ref_base

    stats = []
    for rb in ['A', 'G', 'T', 'C']:
        print(rb, 'mean', numpy.mean(base_dependent_prolongation_stats[rb]))
        print(rb, 'stderr', sem(base_dependent_prolongation_stats[rb]))
        print(rb, 'sd', numpy.std(base_dependent_prolongation_stats[rb]))
        stats.append((rb, numpy.mean(base_dependent_prolongation_stats[rb])))
    return stats


########################################################################################################################

def check_if_intron(pileup_line):
    """
    Function checks whether a given line from the pileup file only contains sequence skips (e.g. '<' and '>')
    :param pileup_line: Given line from the pileup file
    :return: Boolean value 'is_intron'. True if pileup_line only consists of skips, else false.
    """
    '''
    is_intron = True
    pileup_line = pileup_line.split()
    if len(pileup_line) > 4:  # Check if pileup is empty
        for char_x in pileup_line[4]:
            if char_x not in ['<', '>']:
                is_intron = False
                return is_intron
    else:
        False
    return is_intron
    '''

    is_intron = True
    pileup_line = pileup_line.split()
    if len(pileup_line) > 4:  # Check if pileup is empty
        for char_x in pileup_line[4]:
            if char_x not in ['<', '>']:
                is_intron = False
                if not is_intron:
                    return is_intron
    else:
        is_intron = False
    return is_intron


########################################################################################################################

def reinitialize_packages(p1, p2, p3, p4, index, line):
    """
    This function is called when a new exon starts. It re-initializes the furthest forward of the four packages (p4)
    for the read frame at the new positions within the pileup and moves up the other packages by one position.
    :param p1-p4: Packages (dictionaries) for the read frame.
    :param line: Line from pileup file
    :return: All four packages (dictionaries) and the updated index containing the position of 'p4'.
    """
    p1 = copy.deepcopy(p2)
    p2 = copy.deepcopy(p3)
    p3 = copy.deepcopy(p4)

    p4 = parse_pileup(line, p3['ref_seg'], index)
    return p1, p2, p3, p4


########################################################################################################################

def parse_pileup(line, lastrefseq, curr_index):
    """
    This function creates a package given a line from the pileup file.
    :param line: Single line from pileup file.
    :param lastrefseq:
    :param curr_index: Index of the line within pileup file
    :return: Dictionary (Package) containing all relevant information on this line from the pileup.
    """
    split_list = line.split('\t')

    # Split different pileup-sections into lists
    ref_seg = split_list[0]
    pos = split_list[1]
    ref_base = split_list[2].capitalize()
    cov = int(split_list[3])
    pileup = split_list[4]
    qual = split_list[5]

    # Initialize lists for each category
    seq_words = []
    qual_words = []
    ending_indices = []
    starting_indices = []

    # Counts the number of times a read starts at the current position specifically for the three non-G bases
    starting_counts = {'A': 0, 'T': 0, 'C': 0, 'G': 0, 'N': 0}
    i = 0
    qual_index = 0
    len_pileup = len(pileup)
    word_index = 0
    while i < len_pileup:
        word = ''  # Sequence
        qual_word = ''  # Quality
        # In case the read starts ('^')
        if pileup[i] == '^':
            word = pileup[i:i + 3]
            qual_word += qual[qual_index]
            qual_index += 1
            i += 3
            starting_indices.append(word_index)
            starting_counts[(fit_base[ref_base])[word[2]]] += 1
            if i < len_pileup:
                if pileup[i] in ['+', '-']:  # Insertion (+) or deletion (-)
                    word += pileup[i]
                    i += 1
                    number = 0
                    # Check length of insertion/deletion. A 2bp insertion can look like this: +2AG
                    while pileup[i] in ['0', '1', '2', '3', '4', '5', '6', '7', '8', '9']:
                        number *= 10  # Insertion/deletion no longer in list
                        number += int(pileup[i])
                        i += 1
                    indel = pileup[i:i + number]  # Get insertion/deletion string
                    word += str(number) + indel
                    i += number  # Always increment i to get to next position
                if i < len(pileup):
                    if pileup[i] == '$':  # One-base read ends
                        word += pileup[i]
                        i += 1
        elif pileup[i] in ['*', ',', '.', 'A', 'G', 'T', 'C', 'N', 'a', 'g', 't', 'c', 'n', '>', '<']:
            word += pileup[i]
            qual_word += qual[qual_index]
            qual_index += 1
            i += 1
            if i < len_pileup:
                if pileup[i] == '$':  # Read ends
                    word += pileup[i]
                    i += 1
                    ending_indices.append(word_index)
                elif pileup[i] in ['+', '-']:  # Indel
                    word += pileup[i]
                    i += 1
                    number = 0
                    while pileup[i] in ['0', '1', '2', '3', '4', '5', '6', '7', '8', '9']:
                        number *= 10
                        number += int(pileup[i])
                        i += 1
                    indel = pileup[i:i + number]
                    word += str(number) + indel
                    i += number

        # Add relevant information to each category
        seq_words.append(word)
        word_index += 1
        qual_words.append(qual_word)

    return {'last_ref_seg': last_ref_seg, 'ref_seg': ref_seg, 'pos': pos, 'ref_base': ref_base, 'cov': cov, 'pileup':
            pileup, 'qual': qual, 'seq_words': seq_words, 'qual_words': qual_words, 'ending_indices': ending_indices,
            'starting_indices': starting_indices, 'starting_counts': starting_counts, 'curr_position': curr_index}


########################################################################################################################

def overhang_continuation_evidence(p1, p2, p3, p4, level):
    """
    Checks whether the G overhang is comprised of multiple bases????? I don't know
    :param p1, p2, p3, p4: Dictionaries containing information on the pileup file.
    :param level: wtf is level?
    :return: Boolean value. True if the G in the read is likely part of the overhang. False if it is unlikely.
    """
    subtractor = 0
    for index in p1['ending_indices']:
        if level > index:
            subtractor += 1
    level -= subtractor
    if (p2['seq_words'])[level] in [tailing_base_up, tailing_base_low]:  # G in read, but not in reference
        return True
    elif p2['ref_base'] == tailing_base_up and (p2['seq_words'])[level] in ['.', ',']:  # G in read, G in reference
        subtractor = 0
        for index in p2['ending_indices']:
            if level > index:
                subtractor += 1
        level -= subtractor
        if (p3['seq_words'])[level] in [tailing_base_up, tailing_base_low]:  # G in read, other in reference
            return True
        elif p3['ref_base'] == tailing_base_up and (p3['seq_words'])[level] in ['.', ',']:  # G in read, G in reference
            subtractor = 0
            for index in p3['ending_indices']:
                if level > index:
                    subtractor += 1
            level -= subtractor
            if (p4['seq_words'])[level] in [tailing_base_up, tailing_base_low]:  # G in read, other in reference
                return True
            else:
                return False
        else:
            return False
    else:
        return False


########################################################################################################################

def build_output(package):
    """
    Prepares the trimmed pileup and writes it into the output file.
    :param package: A dictionary containing all relevant information for the pileup format including the reference
    sequence, position within the reference, reference base, coverage, pile and quality scores.
    """
    for pti in potential_trimming_indices:
        (package['seq_words'])[pti] = '^!' + (package['seq_words'])[pti]
        elevator_dict[pti] = True  # Move starting read to current top of pile
        try:
            elevator_list.remove(pti)  # priority order of starting points is corrected
        except:
           pass
        elevator_list.append(pti)

    cov = package['cov']
    for pti in pti_deletor:
        potential_trimming_indices.remove(pti)
    for pti in immediate_trimming_indices:
        deletion_indices[pti] = True  # Only pseudo deletion
        cov -= 1

        potential_trimming_indices.append(pti)  # Only relevant from next round on

    # Adjustment of potential_trimming_indices to upstream ending reads
    for i in range(0, len(potential_trimming_indices)):
        subtractor = 0
        for j in range(0, len(package_1['ending_indices'])):
            if potential_trimming_indices[i] > package_1['ending_indices'][j]:
                subtractor += 1
        # Use subtractor to avoid in-place-overwriting inconsistencies (faster than sort alternative)
        potential_trimming_indices[i] -= subtractor

    new_pile = ''
    new_qual = ''
    # Create trimmed outputs
    for i in range(0, len(package['seq_words'])):
        if i not in elevator_dict:
            new_pile += (package['seq_words'])[i]
            new_qual += (package['qual_words'])[i]
    # print('Len seq_words: ', len(package['seq_words']))
    # print('Seq_words ', package['seq_words'])
    # print('Len qual_words: ', len(package['qual_words']))
    # print('Qual_words ', package['qual_words'])
    # print('Elevator_list: ', len(elevator_list))
    # print(elevator_list)
    for index in elevator_list:
        new_pile += (package['seq_words'])[index]
        new_qual += (package['qual_words'])[index]

    outfile.write(package['ref_seg'] + '\t' + package['pos'] + '\t' + package['ref_base'] + '\t' + str(package['cov']) +
                  '\t' + new_pile + '\t' + new_qual + '\n')


########################################################################################################################

def get_dropped_elevators(elevator_list, elevator_dict, ending_indices):
    drop_dict = {}
    dropped_elevator_dict = {}
    dropped_elevator_list = []
    for i in range(0, len(elevator_list)):
        drop_dict[elevator_list[i]] = 0
        for end in ending_indices:
            if elevator_list[i] > end:  # ending happens somewhere below to-be-elevated line
                drop_dict[elevator_list[i]] += 1
        if drop_dict[elevator_list[i]] > 0:
            dropped_elevator_dict[elevator_list[i] - drop_dict[elevator_list[i]]] = elevator_dict[elevator_list[i]]
            dropped_elevator_list.append(elevator_list[i] - drop_dict[elevator_list[i]])
        else:
            dropped_elevator_dict[elevator_list[i]] = elevator_dict[elevator_list[i]]
            dropped_elevator_list.append(elevator_list[i])
    return dropped_elevator_dict, dropped_elevator_list


########################################################################################################################

def get_active_elevators(elevator_list, elevator_dict, ending_indices):
    """
    Removes inactive elevators and returns an updated list containing all elevators.
    :param elevator_list:
    :param elevator_dict:
    :param ending_indices: List containing indices of ending reads ('$'-symbol in pileup format).
    :return:
    """
    for elevator in range(0, len(ending_indices)):
        try:
            elevator_list.remove(ending_indices[elevator])
            del elevator_dict[ending_indices[elevator]]  # inactivate elevator, if present
        except:
            pass
    return elevator_dict, elevator_list


########################################################################################################################

def initialize_packages(file_path):
    """
    Build the four initial packages (dictionaries) that serve as a read frame for the pileup format. Package 1 contains
    the first line of the pileup file, package 4 contains the fourth.
    :param file_path: Path to the input file.
    :return: packages 1-4
    """
    with open(file_path, 'r') as fin:
        temp_counter = 0
        list_lines = []
        for line in fin:
            if temp_counter < 4:
                list_lines.append(line)
                temp_counter += 1
            else:
                break
    # Initialize first four packages starting with the first line
    package_1 = parse_pileup(list_lines[0], last_ref_seg, 0)
    package_2 = parse_pileup(list_lines[1], package_1['ref_seg'], 1)
    package_3 = parse_pileup(list_lines[2], package_2['ref_seg'], 2)
    package_4 = parse_pileup(list_lines[3], package_3['ref_seg'], 3)

    return package_1, package_2, package_3, package_4


########################################################################################################################

# MAIN #
tailing_base_up = sys.argv[2]  # Base to be trimmed
tailing_base_low = tailing_base_up.lower()

if tailing_base_up == 'A':
    bases = ['C', 'G', 'T', 'N']
elif tailing_base_up == 'C':
    bases = ['A', 'G', 'T', 'N']
elif tailing_base_up == 'G':
    bases = ['A', 'C', 'T', 'N']
elif tailing_base_up == 'T':
    bases = ['A', 'C', 'G', 'N']

print(tailing_base_up, tailing_base_low)

# Get overhang statistics
statistics = determine_overhang_prolongation_probabilities(sys.argv[1])

a_e_o = (statistics[0])[1]
g_e_o = (statistics[1])[1]
t_e_o = (statistics[2])[1]
c_e_o = (statistics[3])[1]
n_e_o = 0.1  # Dummy value

# a_e_o = 0.5
# g_e_o = 0.9
# t_e_o = 0.6
# c_e_o = 0.7

# Dictionary containing the overhang statistics
expected_overhang = {
    'A': {'.': a_e_o, ',': a_e_o, 'G': g_e_o, 'g': g_e_o, 'T': t_e_o, 't': t_e_o, 'C': c_e_o, 'c': c_e_o, 'N': n_e_o,
          'n': n_e_o},
    'G': {'A': a_e_o, 'a': a_e_o, '.': g_e_o, ',': g_e_o, 'T': t_e_o, 't': t_e_o, 'C': c_e_o, 'c': c_e_o, 'N': n_e_o,
          'n': n_e_o},
    'T': {'A': a_e_o, 'a': a_e_o, 'G': g_e_o, 'g': g_e_o, '.': t_e_o, ',': t_e_o, 'C': c_e_o, 'c': c_e_o, 'N': n_e_o,
          'n': n_e_o},
    'C': {'A': a_e_o, 'a': a_e_o, 'G': g_e_o, 'g': g_e_o, 'T': t_e_o, 't': t_e_o, '.': c_e_o, ',': c_e_o, 'N': n_e_o,
          'n': n_e_o},
    'N': {'A': a_e_o, 'a': a_e_o, 'G': g_e_o, 'g': g_e_o, 'T': t_e_o, 't': t_e_o, 'C': c_e_o, 'c': c_e_o, '.': n_e_o,
          ',': n_e_o, 'N': n_e_o, 'n': n_e_o}}

# Dictionary containing all possible symbols and their meaning that can occur at the given reference base
fit_base = {'A': {'.': 'A', ',': 'A', 'G': 'G', 'g': 'G', 'T': 'T', 't': 'T', 'C': 'C', 'c': 'C', 'n': 'N', 'N': 'N'},
            'G': {'A': 'A', 'a': 'A', '.': 'G', ',': 'G', 'T': 'T', 't': 'T', 'C': 'C', 'c': 'C', 'n': 'N', 'N': 'N'},
            'T': {'A': 'A', 'a': 'A', 'G': 'G', 'g': 'G', '.': 'T', ',': 'T', 'C': 'C', 'c': 'C', 'n': 'N', 'N': 'N'},
            'C': {'A': 'A', 'a': 'A', 'G': 'G', 'g': 'G', 'T': 'T', 't': 'T', '.': 'C', ',': 'C', 'n': 'N', 'N': 'N'},
            'N': {'A': 'A', 'a': 'A', 'G': 'G', 'g': 'G', 'T': 'T', 't': 'T', 'C': 'C', 'c': 'C', ',': 'N', '.': 'N',
                  'N': 'N', 'n': 'N'}}

line = None
index_intron_start = 0   # Variable to store starting index of an intron.
last_ref_seg = ''
package_1, package_2, package_3, package_4 = initialize_packages(sys.argv[1])
line_counter = 2
immediate_trimming_indices = []
# with open(sys.argv[1], 'r') as infile, open(sys.argv[1] + '.overhang-trimmed', 'w') as outfile:
with open(sys.argv[1], 'r') as infile, open(sys.argv[3], 'w') as outfile:
    # First cycle of new reference
    auxiliary_line = ''
    intron_lines = []
    potential_trimming_indices = []
    immediate_trimming_indices = []
    elevator_dict = {}
    elevator_list = []
    deletion_indices = {}

    # Iterate over lines of pileup file.
    for my_index, line in enumerate(infile):
        if my_index > 3:  # First 4 lines can be skipped, since packages have already been initialized.
            line_counter += 1
            cov = package_1['cov']
            len_ending_indices = len(package_1['ending_indices'])

            # Examine read starts and activate elevators for pileup definition maintenance
            # ToDo: Warum ist hier kein G mit dabei????
            applied_one_base_overhang = {'A': 0.0, 'T': 0.0, 'C': 0.0, 'N': 0.0, 'G': 0.0}

            for i in package_1['starting_indices']:
                word = (package_1['seq_words'])[i]
                elevator_dict[i] = True  # Move starting read to current top of pile
                elevator_list.append(i)
                if word[2] in [tailing_base_up, tailing_base_low]:  # Refbase can't be G, otherwise word[2] is . or ,
                    immediate_trimming_indices.append(i)
                elif package_1['ref_base'] == tailing_base_up and word[2] in ['.', ',']:
                    if overhang_continuation_evidence(package_1, package_2, package_3, package_4, i):
                        immediate_trimming_indices.append(i)
                    else:
                        # print(package_2['seq_words'])
                        # print(len_ending_indices)
                        # print(fit_base[package_2['ref_base']])
                        # print(i - len_ending_indices)
                        # print(package_2)
                        if ((package_2['seq_words'])[i - len_ending_indices])[0] not in ['<', '>']:
                            responsible_base = (fit_base[package_2['ref_base']])[((package_2['seq_words'])[i - len_ending_indices])[0]]
                        else:
                            continue

                        #print('Package 2 seq_words: ', package_2['seq_words'])
                        #print('Package 2 ref_base: ', expected_overhang[package_2['ref_base']])
                        #print(i - len_ending_indices)
                        #print(((expected_overhang[package_2['ref_base']])[((package_2['seq_words'])[i - len_ending_indices])[0]]))
                        tailing_efficiency = ((expected_overhang[package_2['ref_base']])[((package_2['seq_words'])[i - len_ending_indices])[0]])

                        # rule of 3: the quotient calculates 1 % of all bases by dividing the untailed bases by the
                        # expected fraction.
                        if responsible_base != tailing_base_up and applied_one_base_overhang[responsible_base] < (
                                (package_2['starting_counts'])[responsible_base] / (100.0 * (
                                1.0 - tailing_efficiency))) * (100 * tailing_efficiency):
                            immediate_trimming_indices.append(i)
                            applied_one_base_overhang[responsible_base] += 1.0

            # Overhang trimming continuation
            applied_one_base_overhang = {'A': 0.0, 'T': 0.0, 'C': 0.0, 'N': 0.0, 'G': 0.0}
            pti_deletor = []
            for pti in potential_trimming_indices:
                if (package_1['seq_words'])[pti] in [tailing_base_up, tailing_base_low]:  # Active overhang
                    deletion_indices[pti] = True
                    cov -= 1  # Coverage adjustment
                elif package_1['ref_base'] == tailing_base_up and (package_1['seq_words'])[pti] in ['.', ',']:
                    if overhang_continuation_evidence(package_1, package_2, package_3, package_4, pti):
                        deletion_indices[pti] = True
                        cov -= 1  # Coverage adjustment
                    else:
                        if ((package_2['seq_words'])[pti - len_ending_indices])[0] not in ['<', '>']:
                            responsible_base = (fit_base[package_2['ref_base']])[
                                ((package_2['seq_words'])[pti - len_ending_indices])[0]]
                        else:
                            continue
                        tailing_efficiency = ((expected_overhang[package_2['ref_base']])[((package_2['seq_words'])
                        [pti - len_ending_indices])[0]])

                        if responsible_base != tailing_base_up and applied_one_base_overhang[responsible_base] < (
                                (package_2['starting_counts'])[responsible_base] / (100.0 * (
                                1.0 - tailing_efficiency))) * (100 * tailing_efficiency):  # Apply overhang statistics
                            deletion_indices[pti] = True
                            cov -= 1
                            applied_one_base_overhang[responsible_base] += 1.0

                        else:
                            (package_1['seq_words'])[pti] = '^!' + (package_1['seq_words'])[pti]
                            # This case should also be prepared for
                            """
                            if (seq_words[pti])[-1] == "$":
                            """
                            elevator_dict[pti] = True  # Move starting read to current top of pile
                            try:
                                elevator_list.remove(pti)  # priority order of starting points is corrected
                            except:
                                pass
                            elevator_list.append(pti)
                            pti_deletor.append(pti)
                else:
                    (package_1['seq_words'])[pti] = '^!' + (package_1['seq_words'])[pti]
                    # This case should also be prepared for
                    """
                    if (seq_words[pti])[-1] == "$":
                    """
                    elevator_dict[pti] = True  # Move starting read to current top of pile
                    try:
                        elevator_list.remove(pti)  # priority order of starting points is corrected
                    except:
                        pass
                    elevator_list.append(pti)
                    pti_deletor.append(pti)
            for pti in pti_deletor:
                potential_trimming_indices.remove(pti)

            # Here, the latest starting reads are trimmed. The indices of interest were collected before trimming old
            # potentially interesting indices. Therefore, the potentially gut-shot pile requires correction of top-level
            # indices of starting reads.
            for pti in immediate_trimming_indices:
                deletion_indices[pti] = True  # Only pseudo deletion
                cov -= 1
                potential_trimming_indices.append(pti)  # Only relevant from next round on

            # Adjustment of potential_trimming_indices to upstream ending reads
            for i in range(0, len(potential_trimming_indices)):
                subtractor = 0
                for j in range(0, len(package_1['ending_indices'])):
                    if potential_trimming_indices[i] > package_1['ending_indices'][j]:
                        subtractor += 1
                # Use subtractor to avoid in-place-overwriting inconsistencies (faster than sort alternative)
                potential_trimming_indices[i] -= subtractor

            # Build trimmed output strings
            new_pile = ''
            new_qual = ''
            for i in range(0, len(package_1['seq_words'])):
                if i not in elevator_dict and i not in deletion_indices:
                    new_pile += (package_1['seq_words'])[i]
                    new_qual += (package_1['qual_words'])[i]
            for index in elevator_list:
                if index not in deletion_indices:
                    new_pile += (package_1['seq_words'])[index]
                    new_qual += (package_1['qual_words'])[index]

            # Write to trimmed pileup file
            outfile.write(package_1['ref_seg'] + '\t' + package_1['pos'] + '\t' + package_1['ref_base'] + '\t' +
                          str(cov) + '\t' + new_pile + '\t' + new_qual + '\n')

            # Ended reads don't need elevators anymore
            elevator_dict, elevator_list = get_active_elevators(elevator_list, elevator_dict,
                                                                package_1['ending_indices'])

            if elevator_dict:
                # Adjustment of elevators according to ending reads
                elevator_dict, elevator_list = get_dropped_elevators(elevator_list, elevator_dict,
                                                                     package_1['ending_indices'])

            # Reset temporary variables
            immediate_trimming_indices = []
            deletion_indices = {}
            applied_one_base_overhang = 0.0

            # Shift frame 1 refbase forward
            if check_if_intron(line):  # Check if line is intron (e.g. only contains '<' and '>')
                intron_lines = []
                index_intron_start = line_counter
                # Save all lines containing only skips in list. These lines cannot be written into the output file
                # immediately, because in each loop only package_1 is written into the output. Here, only package_4
                # makes the 'jump' beyond the intron, while the other packages are still at positions before the intron.
                # Accordingly, the loop has to be repeated until package_1 also makes the 'jump'. Only then the lines
                # containing the intron are written to the output file.
                intron_lines.append(line)
                line_counter += 1
                for line in infile:
                    if not check_if_intron(line):
                        auxiliary_line = line  # This line is usually jumped over -> temporary variable
                        line_counter += 1
                        break
                    else:
                        intron_lines.append(line)
                    line_counter += 1
                # Reinitialize package_4 at the new position after the intron
                package_1, package_2, package_3, package_4 = reinitialize_packages(
                    package_1, package_2, package_3, package_4, line_counter, auxiliary_line)
            else:
                package_1, package_2, package_3, package_4 = reinitialize_packages(
                    package_1, package_2, package_3, package_4, line_counter, line)

            # Write intron lines to the output
            if package_1['curr_position'] >= index_intron_start and len(intron_lines) > 0:
                for i in range(len(intron_lines)):
                    outfile.write(intron_lines[i])
                intron_lines = []
                index_intron_start = 0

    if line is not None:  # When function has reached the end of the file write remaining packages into output.
        build_output(package_1)
        # Ended reads don't need elevators anymore
        elevator_dict, elevator_list = get_active_elevators(elevator_list, elevator_dict, package_1['ending_indices'])
        # Adjustment of elevators according to ending reads
        elevator_dict, elevator_list = get_dropped_elevators(elevator_list, elevator_dict, package_1['ending_indices'])
        immediate_trimming_indices = []
        deletion_indices = {}
        applied_one_base_overhang = 0.0
        build_output(package_2)
        elevator_dict, elevator_list = get_active_elevators(elevator_list, elevator_dict, package_2['ending_indices'])
        elevator_dict, elevator_list = get_dropped_elevators(elevator_list, elevator_dict, package_2['ending_indices'])
        immediate_trimming_indices = []
        deletion_indices = {}
        applied_one_base_overhang = 0.0
        build_output(package_3)
        elevator_dict, elevator_list = get_active_elevators(elevator_list, elevator_dict, package_3['ending_indices'])
        elevator_dict, elevator_list = get_dropped_elevators(elevator_list, elevator_dict, package_3['ending_indices'])
        immediate_trimming_indices = []
        deletion_indices = {}
        applied_one_base_overhang = 0.0
        build_output(package_4)

b = datetime.datetime.now()
print('Overall time: ', b-a)

# The End

