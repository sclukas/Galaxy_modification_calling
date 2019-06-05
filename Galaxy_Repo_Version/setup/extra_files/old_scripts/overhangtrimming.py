#!/usr/bin/python

import os
import sys
sys.path.append('/usr/local/lib/python2.7/dist-packages')
import copy
import subprocess
import numpy
from scipy.stats import sem

# ATTENTION! This program requires a minimum read length of 5 to work properly


# Typically, the tailing efficiency of G determined by this tool must be much lower than the other ones, because much more reads generally start with G than with other bases.
# This lowers the resulting percentage for G. As this statistic is only applied on GG spots for the second G, we expect a much higher fraction of true Gs on the first G of GG,
# because of course some reads simply start untailed at this G und should not be trimmed.


# Determination of overhang prolongation probability
def determine_overhang_prolongation_probabilities(fin):

    curr_startstats = {'A': 0.0000001, 'G': 0.0000001, 'T': 0.0000001, 'C': 0.0000001, 'N': 0.0000001}
    curr_prolongationcounts = {'A': 0.0000001, 'G': 0.0000001, 'T': 0.0000001, 'C': 0.0000001, 'N': 0.0000001}
    basedependent_prolongation_stats = {'A': [], 'G': [], 'T': [], 'C': [], 'N': []}

    infile = open(fin, "r")

    potential_trimming_indices = []
    immediate_trimming_indices = []
    elevatordict = {}
    elevatorlist = []
    ending_indices = []
    starting_indices = []

    lastrefseg = ""
    lastrefbase = ""
    line = infile.readline()
    while len(line) > 0:
        splitlist = line.split("\t")

        refseg = splitlist[0]
        pos = splitlist[1]
        refbase = splitlist[2]
        cov = int(splitlist[3])
        pileup = splitlist[4]
        qual = splitlist[5]

        if refseg == 'G':
            g_series += 1
        else:
            g_series = 0

        if refseg != lastrefseg:
            print(refseg)
            potential_trimming_indices = []
            immediate_trimming_indices = []
            ending_indices = []
            starting_indices = []

        seqwords = []
        qualwords = []
        i = 0
        qualindex = 0
        wordindex = 0
        while i < len(pileup):
            word = ""
            qualword = ""
            if pileup[i] == "^":  # Read starts
                word = pileup[i:i + 3]
                qualword += qual[qualindex]
                qualindex += 1
                i += 3
                starting_indices.append(wordindex)
                if i < len(pileup):

                    if pileup[i] in ["+", "-"]:
                        word += pileup[i]
                        i += 1
                        number = 0
                        while pileup[i] in ["0", "1", "2", "3", "4", "5", "6", "7", "8", "9"]:
                            number *= 10
                            number += int(pileup[i])
                            i += 1
                        indel = pileup[i:i + number]
                        word += str(number) + indel
                        i += number
                    if i < len(pileup):
                        if pileup[i] == "$":  # One-base read ends
                            word += pileup[i]
                            i += 1
                            ending_indices.append(wordindex)

            elif pileup[i] in ["*", ",", ".", "A", "G", "T", "C", "N", "a", "g", "t", "c", "n"]:

                word += pileup[i]
                qualword += qual[qualindex]
                qualindex += 1
                i += 1
                if i < len(pileup):
                    if pileup[i] == "$":  # Read ends
                        word += pileup[i]
                        i += 1
                        ending_indices.append(wordindex)
                    elif pileup[i] in ["+", "-"]:  # Indel
                        word += pileup[i]
                        i += 1
                        number = 0
                        while pileup[i] in ["0", "1", "2", "3", "4", "5", "6", "7", "8", "9"]:
                            number *= 10
                            number += int(pileup[i])
                            i += 1
                        indel = pileup[i:i + number]
                        word += str(number) + indel
                        i += number

            else:
                "Alarm!!! This should never happen."

            seqwords.append(word)
            qualwords.append(qualword)
            wordindex += 1

        curr_prolongationcounts = {'A': 0.0000001, 'G': 0.0000001, 'T': 0.0000001, 'C': 0.0000001, 'N': 0.0000001}
        curr_startstats = {'A': 0.0000001, 'G': 0.0000001, 'T': 0.0000001, 'C': 0.0000001, 'N': 0.0000001}

        for i in starting_indices:
            word = seqwords[i]
            if refbase != 'G':
                if len(word) == 3:
                    if word[2] in ["G", "g"]:  # reference must be non-G to increase overhang evidence
                        immediate_trimming_indices.append(i)
            if word[2] in ['.', ',']:
                curr_startstats[refbase] += 1
            else:
                curr_startstats[(word[2]).upper()] += 1

            if word.endswith("$"):  # This case is not defined in this statistics script due to assumed minimum read length
                print(word)

        ptideletor = []
        for pti in potential_trimming_indices:  # Overhang trimming continuation

            if len(seqwords[pti]) == 1:
                if seqwords[pti] not in ['G', 'g'] and not (
                        seqwords[pti] in ['.', ','] and refbase == 'G'):  # read base is non-G

                    try:
                        curr_prolongationcounts[(seqwords[pti]).upper()] += 1
                    except:
                        curr_prolongationcounts[refbase] += 1
                    ptideletor.append(pti)

                elif refbase == "G" and seqwords[pti] in ['.', ',']:
                    curr_prolongationcounts['G'] += 1
                    ptideletor.append(pti)
                else:
                    pass  # G takes part in the statistics only if refbase == G. Remember this, when summing up rb-spec tailing efficiencies.
            # print refbase
            # else :
            # ptideletor.append(pti)

        for pti in ptideletor:
            potential_trimming_indices.remove(pti)

        for rb in ['A', 'T', 'C', 'N']:  # curr_prolongationcounts:
            # if int(curr_startstats[rb] + curr_prolongationcounts[rb]) >= 20:

            for i in range(0, int(curr_startstats[rb] + curr_prolongationcounts[rb])):
                (basedependent_prolongation_stats[rb]).append(
                    float(curr_prolongationcounts[rb]) / (float(curr_startstats[rb] + curr_prolongationcounts[rb])))
        if refbase == 'G' and lastrefbase != 'G':
            # if int(curr_startstats['G'] + curr_prolongationcounts['G']) >= 1:
            for i in range(0, int(curr_startstats['G'] + curr_prolongationcounts['G'])):
                (basedependent_prolongation_stats['G']).append(
                    float(curr_prolongationcounts['G']) / (float(curr_startstats['G'] + curr_prolongationcounts['G'])))
            # print line
            # print "pos", pos, "curr_startstats['G']", curr_startstats['G'], "curr_prolongationcounts[rb]", curr_prolongationcounts['G'], "ratio", float(curr_prolongationcounts['G']) /  (  float(curr_startstats['G'] + curr_prolongationcounts['G'])   )
            # print "pos", pos, "curr_startstats['A']", curr_startstats['A'], "curr_prolongationcounts[rb]", curr_prolongationcounts['A'], "ratio", float(curr_prolongationcounts['A']) /  (  float(curr_startstats['A'] + curr_prolongationcounts['A'])   )
            # print "pos", pos, "curr_startstats['T']", curr_startstats['T'], "curr_prolongationcounts[rb]", curr_prolongationcounts['T'], "ratio", float(curr_prolongationcounts['T']) /  (  float(curr_startstats['T'] + curr_prolongationcounts['T'])   )

        # print line
        # print "pos", pos, "curr_startstats['G']", curr_startstats['G'], "curr_prolongationcounts[rb]", curr_prolongationcounts['G'], "ratio", float(curr_prolongationcounts['G']) /  (  float(curr_startstats['G'] + curr_prolongationcounts['G'])   )
        # print "pos", pos, "curr_startstats['A']", curr_startstats['A'], "curr_prolongationcounts[rb]", curr_prolongationcounts['A'], "ratio", float(curr_prolongationcounts['A']) /  (  float(curr_startstats['A'] + curr_prolongationcounts['A'])   )
        # print "pos", pos, "curr_startstats['T']", curr_startstats['T'], "curr_prolongationcounts[rb]", curr_prolongationcounts['T'], "ratio", float(curr_prolongationcounts['T']) /  (  float(curr_startstats['T'] + curr_prolongationcounts['T'])   )

        # if lnongseries == 1: # only first base in each non-G-series is respected

        # Here, the latest starting reads are trimmed. The indices of interest were collected before trimming old potentially interesting indices.
        # Therefore, the potentially gut-shot pile requires correction of top-level indices of starting reads.
        for pti in immediate_trimming_indices:
            potential_trimming_indices.append(pti)  # Only relevant from next round on

        for i in range(0, len(
                potential_trimming_indices)):  # Adjustment of potential_trimming_indices to upstream ending reads
            subtractor = 0
            for j in range(0, len(ending_indices)):
                if potential_trimming_indices[i] > ending_indices[j]:
                    subtractor += 1
            potential_trimming_indices[i] -= subtractor

        immediate_trimming_indices = []
        ending_indices = []  # Old upstream ending reads are irrelevant now
        starting_indices = []

        lastrefseg = refseg
        lastrefbase = refbase

        line = infile.readline()

    infile.close()

    # print basedependent_prolongation_stats

    stats = []
    #outfile = open(sys.argv[1] + "_tailing_stats", "w")
    for rb in ['A', 'G', 'T', 'C']:
        print(rb, "mean", numpy.mean(basedependent_prolongation_stats[rb]))
        print(rb, "stderr", sem(basedependent_prolongation_stats[rb]))
        print(rb, "sd", numpy.std(basedependent_prolongation_stats[rb]))

        stats.append((rb, numpy.mean(basedependent_prolongation_stats[rb])))
        #outfile.write(rb + "\tmean\t" + str(numpy.mean(basedependent_prolongation_stats[rb])) + "\n")
        #outfile.write(rb + "\tstderr\t" + str(sem(basedependent_prolongation_stats[rb])) + "\n")
        #outfile.write(rb + "\tsd\t" + str(numpy.std(basedependent_prolongation_stats[rb])) + "\n")
    #outfile.close()
    return stats
    #outfile.close()
########################################################################################################################

"""
p = subprocess.Popen("python /home/ralf/Desktop/Code_PhD/DetermineOverhangProlongationProbabilities.py" + " " + sys.argv[1], stdout=subprocess.PIPE)
result = p.communicate()[0]

print "#" + result + "#"
asfasf

stats = result.split(" ")

expected_overhang['A':float(stats[1]), 'G':float(stats[3]), 'T':float(stats[5]), 'C':float(stats[7])]
"""


def stripLastNPathElements(path, n):  # Remove n rightmost Elements of a path (w.r.t. "/")

    return (path.rsplit("/", n))[0]


#os.system("python " + stripLastNPathElements(sys.argv[0], 1) + "/DetermineOverhangProlongationProbabilities.py " + sys.argv[
     #   1])




statistics = determine_overhang_prolongation_probabilities(sys.argv[1])

#statsfile = open('/home/akhelm/Downloads/1501_Ludi2_10000000.pileup_tailing_stats', 'r')

a_e_o = (statistics[0])[1]
g_e_o = (statistics[1])[1]
t_e_o = (statistics[2])[1]
c_e_o = (statistics[3])[1]

'''
print(a_e_o)
a_e_o = float((((statsfile.readline())[:-1]).split("\t"))[-1])
print(a_e_o)
statsfile.readline()
statsfile.readline()
g_e_o = float((((statsfile.readline())[:-1]).split("\t"))[-1])
print(g_e_o)
statsfile.readline()
statsfile.readline()
t_e_o = float((((statsfile.readline())[:-1]).split("\t"))[-1])
print(t_e_o)
statsfile.readline()
statsfile.readline()
c_e_o = float((((statsfile.readline())[:-1]).split("\t"))[-1])
print(c_e_o)
statsfile.close()
'''
#print("a_e_o", a_e_o, "g_e_o", g_e_o, "t_e_o", t_e_o, "c_e_o", c_e_o)

"""
a_e_o = 0.5
g_e_o = 0.9
t_e_o = 0.6
c_e_o = 0.7

"""

n_e_o = 0.6

expected_overhang = {
    'A': {'.': a_e_o, ',': a_e_o, 'G': g_e_o, 'g': g_e_o, 'T': t_e_o, 't': t_e_o, 'C': c_e_o, 'c': c_e_o, 'N': n_e_o,
          'n': n_e_o},
    'G': {'A': a_e_o, 'a': a_e_o, '.': g_e_o, ',': g_e_o, 'T': t_e_o, 't': t_e_o, 'C': c_e_o, 'c': c_e_o, 'N': n_e_o,
          'n': n_e_o},
    'T': {'A': a_e_o, 'a': a_e_o, 'G': g_e_o, 'g': g_e_o, '.': t_e_o, ',': t_e_o, 'C': c_e_o, 'c': c_e_o, 'N': n_e_o,
          'n': n_e_o},
    'C': {'A': a_e_o, 'a': a_e_o, 'G': g_e_o, 'g': g_e_o, 'T': t_e_o, 't': t_e_o, '.': c_e_o, ',': c_e_o, 'N': n_e_o,
          'n': n_e_o}}

fitbase = {'A': {'.': 'A', ',': 'A', 'G': 'G', 'g': 'G', 'T': 'T', 't': 'T', 'C': 'C', 'c': 'C', 'n': 'N', 'N': 'N'},
           'G': {'A': 'A', 'a': 'A', '.': 'G', ',': 'G', 'T': 'T', 't': 'T', 'C': 'C', 'c': 'C', 'n': 'N', 'N': 'N'},
           'T': {'A': 'A', 'a': 'A', 'G': 'G', 'g': 'G', '.': 'T', ',': 'T', 'C': 'C', 'c': 'C', 'n': 'N', 'N': 'N'},
           'C': {'A': 'A', 'a': 'A', 'G': 'G', 'g': 'G', 'T': 'T', 't': 'T', '.': 'C', ',': 'C', 'n': 'N', 'N': 'N'},
           'N': {'A': 'A', 'a': 'A', 'G': 'G', 'g': 'G', 'T': 'T', 't': 'T', 'C': 'C', 'c': 'C', ',': 'N', '.': 'N'}}


def parsePileup(line, lastrefseq):
    splitlist = line.split("\t")

    refseg = splitlist[0]
    pos = splitlist[1]
    refbase = splitlist[2]
    cov = int(splitlist[3])
    pileup = splitlist[4]
    qual = splitlist[5]

    seqwords = []
    qualwords = []
    ending_indices = []
    starting_indices = []
    starting_counts = {'A': 0, 'T': 0, 'C': 0, 'G': 0, 'N': 0}  # Counts the number of times a read starts at the current position specifically for the three non-G bases
    i = 0
    qualindex = 0
    len_pileup = len(pileup)
    wordindex = 0
    while i < len_pileup:
        word = ""  # Sequence
        qualword = ""  # Quality
        if pileup[i] == "^":  # Read starts
            word = pileup[i:i + 3]
            qualword += qual[qualindex]
            qualindex += 1
            i += 3
            starting_indices.append(wordindex)
            starting_counts[(fitbase[refbase])[word[2]]] += 1
            if i < len_pileup:
                if pileup[i] in ["+", "-"]:
                    word += pileup[i]
                    i += 1
                    number = 0
                    while pileup[i] in ["0", "1", "2", "3", "4", "5", "6", "7", "8", "9"]:
                        number *= 10
                        number += int(pileup[i])
                        i += 1
                    indel = pileup[i:i + number]
                    word += str(number) + indel
                    i += number
                if i < len(pileup):
                    if pileup[i] == "$":  # One-base read ends
                        word += pileup[i]
                        i += 1
        elif pileup[i] in ["*", ",", ".", "A", "G", "T", "C", "N", "a", "g", "t", "c", "n"]:
            word += pileup[i]
            qualword += qual[qualindex]
            qualindex += 1
            i += 1
            if i < len_pileup:
                if pileup[i] == "$":  # Read ends
                    word += pileup[i]
                    i += 1
                    ending_indices.append(wordindex)
                elif pileup[i] in ["+", "-"]:  # Indel
                    word += pileup[i]
                    i += 1
                    number = 0
                    while pileup[i] in ["0", "1", "2", "3", "4", "5", "6", "7", "8", "9"]:
                        number *= 10
                        number += int(pileup[i])
                        i += 1
                    indel = pileup[i:i + number]
                    word += str(number) + indel
                    i += number

        else:
            print("Alarm!!! This should never happen. Invalid word start!")
        seqwords.append(word)
        qualwords.append(qualword)
        wordindex += 1

    return {"lastrefseg": lastrefseg, "refseg": refseg, "pos": pos, "refbase": refbase, "cov": cov, "pileup": pileup,
            "qual": qual, "seqwords": seqwords, "qualwords": qualwords, "ending_indices": ending_indices,
            "starting_indices": starting_indices, "starting_counts": starting_counts}


def overhangContinuationEvidence(package_1, package_2, package_3, package_4, level):
    substractor = 0
    for index in package_1['ending_indices']:
        if level > index:
            substractor += 1
    level -= substractor
    if (package_2['seqwords'])[level] in ['G', 'g']:  # G in read, other in reference
        return True
    elif package_2['refbase'] == 'G' and (package_2['seqwords'])[level] in ['.', ',']:  # G in read, G in reference
        substractor = 0
        for index in package_2['ending_indices']:
            if level > index:
                substractor += 1
        level -= substractor
        if (package_3['seqwords'])[level] in ['G', 'g']:  # G in read, other in reference
            return True
        elif package_3['refbase'] == 'G' and (package_3['seqwords'])[level] in ['.', ',']:  # G in read, G in reference
            substractor = 0
            for index in package_3['ending_indices']:
                if level > index:
                    substractor += 1
            level -= substractor
            if (package_4['seqwords'])[level] in ['G', 'g']:  # G in read, other in reference
                return True
            else:
                return False
        else:
            return False
    else:
        return False


def getDroppedElevators(elevatorlist, elevatordict, ending_indices):
    dropdict = {}
    droppedelevatordict = {}
    droppedelevatorlist = []
    for i in range(0, len(elevatorlist)):

        dropdict[elevatorlist[i]] = 0
        for end in ending_indices:
            if elevatorlist[i] > end:  # ending happens somewhere below to-be-elevated line
                dropdict[elevatorlist[i]] += 1

        if dropdict[elevatorlist[i]] > 0:
            droppedelevatordict[elevatorlist[i] - dropdict[elevatorlist[i]]] = elevatordict[elevatorlist[i]]
            droppedelevatorlist.append(elevatorlist[i] - dropdict[elevatorlist[i]])
        else:
            droppedelevatordict[elevatorlist[i]] = elevatordict[elevatorlist[i]]
            droppedelevatorlist.append(elevatorlist[i])

    return droppedelevatordict, droppedelevatorlist


def getActiveElevators(elevatorlist, elevatordict, ending_indices):
    for i in range(0, len(ending_indices)):
        try:
            elevatorlist.remove(ending_indices[i])
            del elevatordict[ending_indices[i]]  # inactivate elevator, if present
        except:
            pass

    return elevatordict, elevatorlist


infile = open(sys.argv[1], "r")
outfile = open(sys.argv[2], "w")

# First cycleof new reference
potential_trimming_indices = []
immediate_trimming_indices = []
elevatordict = {}
elevatorlist = []
deletion_indices = {}

lastrefseg = ""
line = infile.readline()
line2 = infile.readline()
line3 = infile.readline()
line4 = infile.readline()
package_1 = parsePileup(line, lastrefseg)
package_2 = parsePileup(line2, package_1['refseg'])
package_3 = parsePileup(line3, package_2['refseg'])
package_4 = parsePileup(line4, package_3['refseg'])

# Iterate through lines of pileup file
while len(line) > 0:

    cov = package_1['cov']

    len_ending_indices = len(package_1['ending_indices'])

    # Examine read starts and activate elevators for pileup definition maintenance
    # one_base_overhang_candidates = {'A':0.00000001, 'T':0.00000001, 'C':0.00000001}
    applied_onebase_overhang = {'A': 0.0, 'T': 0.0, 'C': 0.0, 'N': 0.0}
    for i in package_1['starting_indices']:
        word = (package_1['seqwords'])[i]
        elevatordict[i] = True  # Move starting read to current top of pile
        elevatorlist.append(i)
        if word[2] in ["G", "g"]:  # Refbase can't be G, otherwise word[2] is . or ,
            immediate_trimming_indices.append(i)
        elif package_1['refbase'] == "G" and word[2] in [".", ","]:
            if overhangContinuationEvidence(package_1, package_2, package_3, package_4, i):
                immediate_trimming_indices.append(i)
            else:
                responsible_base = (fitbase[package_2['refbase']])[((package_2['seqwords'])[i - len_ending_indices])[0]]

                # one_base_overhang_candidates[responsible_base] += 1.0 # consider incrementing for non-G responsibles only
                # if applied_onebase_overhang / one_base_overhang_candidates < (expected_overhang[package_2['refbase']])[((package_2['seqwords'])[i-len_ending_indices])[0]]:    # Apply overhang statistics
                tailing_efficiency = ((expected_overhang[package_2['refbase']])[((package_2['seqwords'])[i - len_ending_indices])[0]])

                # if responsible_base != "G" and applied_onebase_overhang[responsible_base] / one_base_overhang_candidates[responsible_base] < (tailing_efficiency  * (package_2['starting_counts'])[responsible_base] ) / ((1.0-tailing_efficiency)*one_base_overhang_candidates[responsible_base]):    # Apply overhang statistics
                if responsible_base != "G" and applied_onebase_overhang[responsible_base] < (
                        (package_2['starting_counts'])[responsible_base] / (100.0 * (1.0 - tailing_efficiency))) * (
                        100 * tailing_efficiency):  # rule of 3: the quotient calculates 1 % of all bases by dividing the untailed bases by the expected fraction. the product then
                    immediate_trimming_indices.append(i)
                    applied_onebase_overhang[responsible_base] += 1.0
                    # print "te", tailing_efficiency, "applied", applied_onebase_overhang[responsible_base], "starting", ((package_2['starting_counts'])[responsible_base], "target", ((package_2['starting_counts'])[responsible_base] / (100.0*(1.0-tailing_efficiency))) *  (100*tailing_efficiency))

    # Overhang trimming continuation
    # one_base_overhang_candidates = {'A':0.00000001, 'T':0.00000001, 'C':0.00000001}
    applied_onebase_overhang = {'A': 0.0, 'T': 0.0, 'C': 0.0, 'N': 0.0}
    ptideletor = []
    for pti in potential_trimming_indices:
        if (package_1['seqwords'])[pti] in ["G", "g"]:  # Active overhang
            deletion_indices[pti] = True
            cov -= 1  # Coverage adjustment
        elif package_1['refbase'] == "G" and (package_1['seqwords'])[pti] in [".", ","]:
            if overhangContinuationEvidence(package_1, package_2, package_3, package_4, pti):
                deletion_indices[pti] = True
                cov -= 1  # Coverage adjustment
            else:
                responsible_base = (fitbase[package_2['refbase']])[
                    ((package_2['seqwords'])[pti - len_ending_indices])[0]]

                # one_base_overhang_candidates[responsible_base] += 1.0
                # if applied_onebase_overhang / one_base_overhang_candidates < (expected_overhang[package_2['refbase']])[((package_2['seqwords'])[pti-len_ending_indices])[0]]:    # Apply overhang statistics
                tailing_efficiency = (
                (expected_overhang[package_2['refbase']])[((package_2['seqwords'])[pti - len_ending_indices])[0]])

                # if responsible_base != "G" and applied_onebase_overhang[responsible_base] / one_base_overhang_candidates[responsible_base] < (tailing_efficiency  * (package_2['starting_counts'])[responsible_base] ) / ((1.0-tailing_efficiency)*one_base_overhang_candidates[responsible_base]):    # Apply overhang statistics
                if responsible_base != "G" and applied_onebase_overhang[responsible_base] < (
                        (package_2['starting_counts'])[responsible_base] / (100.0 * (1.0 - tailing_efficiency))) * (
                        100 * tailing_efficiency):  # Apply overhang statistics
                    deletion_indices[pti] = True
                    cov -= 1
                    applied_onebase_overhang[responsible_base] += 1.0
                    # print "te", tailing_efficiency, "applied", applied_onebase_overhang[responsible_base], "starting", ((package_2['starting_counts'])[responsible_base], "target", ((package_2['starting_counts'])[responsible_base] / (100.0*(1.0-tailing_efficiency))) *  (100*tailing_efficiency))

                else:
                    (package_1['seqwords'])[pti] = "^!" + (package_1['seqwords'])[pti]
                    # This case should also be prepared for
                    """
                    if (seqwords[pti])[-1] == "$":
                        dfdfhdfh
                    """
                    elevatordict[pti] = True  # Move starting read to current top of pile
                    try:
                        elevatorlist.remove(pti)  # priority order of starting points is corrected
                    except:
                        pass
                    elevatorlist.append(pti)
                    ptideletor.append(pti)
        else:
            (package_1['seqwords'])[pti] = "^!" + (package_1['seqwords'])[pti]
            # This case should also be prepared for
            """
            if (seqwords[pti])[-1] == "$":
                dfdfhdfh
            """
            elevatordict[pti] = True  # Move starting read to current top of pile
            try:
                elevatorlist.remove(pti)  # priority order of starting points is corrected
            except:
                pass
            elevatorlist.append(pti)
            ptideletor.append(pti)
    for pti in ptideletor:
        potential_trimming_indices.remove(pti)

    # Here, the latest starting reads are trimmed. The indices of interest were collected before trimming old potentially interesting indices.
    # Therefore, the potentially gut-shot pile requires correction of top-level indices of starting reads.
    for pti in immediate_trimming_indices:
        deletion_indices[pti] = True  # only pseudo deletion
        cov -= 1

        potential_trimming_indices.append(pti)  # Only relevant from next round on

    # Adjustment of potential_trimming_indices to upstream ending reads
    for i in range(0, len(potential_trimming_indices)):
        subtractor = 0
        for j in range(0, len(package_1['ending_indices'])):
            if potential_trimming_indices[i] > package_1['ending_indices'][j]:
                subtractor += 1
        potential_trimming_indices[
            i] -= subtractor  # use substractor to avoid in-place-overwriting inconsistencies (faster than sort alternative)

    # Build trimmed output strings
    newpile = ""
    newqual = ""
    for i in range(0, len(package_1['seqwords'])):
        if i not in elevatordict and i not in deletion_indices:
            newpile += (package_1['seqwords'])[i]
            newqual += (package_1['qualwords'])[i]
    for index in elevatorlist:
        if index not in deletion_indices:
            newpile += (package_1['seqwords'])[index]
            newqual += (package_1['qualwords'])[index]

    # Write to trimmed pileup file
    outfile.write(package_1['refseg'] + "\t" + package_1['pos'] + "\t" + package_1['refbase'] + "\t" + str(
        cov) + "\t" + newpile + "\t" + newqual + "\n")

    # Ended reads don't need elevators anymore
    elevatordict, elevatorlist = getActiveElevators(elevatorlist, elevatordict, package_1['ending_indices'])

    # Adjustment of elevators according to ending reads
    elevatordict, elevatorlist = getDroppedElevators(elevatorlist, elevatordict, package_1['ending_indices'])

    # Reset temporary variables
    immediate_trimming_indices = []
    deletion_indices = {}
    applied_onebase_overhang = 0.0

    # Shift frame 1 refbase forward
    line = line2
    package_1 = copy.deepcopy(package_2)
    line2 = line3
    package_2 = copy.deepcopy(package_3)
    line3 = line4
    package_3 = copy.deepcopy(package_4)
    line4 = infile.readline()
    if len(line4) > 0:
        package_4 = parsePileup(line4, package_3['refseg'])

    # New reference starts soon
    if len(line4) == 0 or package_4['refseg'] != package_1['refseg']:

        print("Done:", package_1['refseg'])

        ####################
        # n-2

        # Terminate all running trimmings
        for pti in potential_trimming_indices:

            (package_1['seqwords'])[pti] = "^!" + (package_1['seqwords'])[pti]

            elevatordict[pti] = True  # Move starting read to current top of pile
            try:
                elevatorlist.remove(pti)  # priority order of starting points is corrected
            except:
                pass
            elevatorlist.append(pti)

        # Build trimmed output strings
        newpile = ""
        newqual = ""
        for i in range(0, len(package_1['seqwords'])):
            if i not in elevatordict:
                newpile += (package_1['seqwords'])[i]
                newqual += (package_1['qualwords'])[i]
        for index in elevatorlist:
            newpile += (package_1['seqwords'])[index]
            newqual += (package_1['qualwords'])[index]

        outfile.write(package_1['refseg'] + "\t" + package_1['pos'] + "\t" + package_1['refbase'] + "\t" + str(
            cov) + "\t" + newpile + "\t" + newqual + "\n")

        # Ended reads don't need elevators anymore
        elevatordict, elevatorlist = getActiveElevators(elevatorlist, elevatordict, package_1['ending_indices'])
        # Adjustment of elevators according to ending reads
        elevatordict, elevatorlist = getDroppedElevators(elevatorlist, elevatordict, package_1['ending_indices'])

        ####################
        # n-1

        # Build trimmed output strings
        newpile = ""
        newqual = ""
        for i in range(0, len(package_2['seqwords'])):
            if i not in elevatordict:
                newpile += (package_2['seqwords'])[i]
                newqual += (package_2['qualwords'])[i]
        for index in elevatorlist:
            newpile += (package_2['seqwords'])[index]
            newqual += (package_2['qualwords'])[index]

        outfile.write(package_2['refseg'] + "\t" + package_2['pos'] + "\t" + package_2['refbase'] + "\t" + str(
            package_2['cov']) + "\t" + newpile + "\t" + newqual + "\n")

        # Ended reads don't need elevators anymore
        elevatordict, elevatorlist = getActiveElevators(elevatorlist, elevatordict, package_2['ending_indices'])
        # Adjustment of elevators according to ending reads
        elevatordict, elevatorlist = getDroppedElevators(elevatorlist, elevatordict, package_2['ending_indices'])

        ####################
        # n

        # Build trimmed output strings
        newpile = ""
        newqual = ""
        for i in range(0, len(package_3['seqwords'])):
            if i not in elevatordict:
                newpile += (package_3['seqwords'])[i]
                newqual += (package_3['qualwords'])[i]
        for index in elevatorlist:
            newpile += (package_3['seqwords'])[index]
            newqual += (package_3['qualwords'])[index]

        outfile.write(package_3['refseg'] + "\t" + package_3['pos'] + "\t" + package_3['refbase'] + "\t" + str(
            package_3['cov']) + "\t" + newpile + "\t" + newqual + "\n")

        # Relaunch on new reference
        # First cycle of new reference
        potential_trimming_indices = []
        immediate_trimming_indices = []
        elevatordict = {}
        elevatorlist = []
        deletion_indices = {}

        lastrefseg = ""
        line = line4
        line2 = infile.readline()
        line3 = infile.readline()
        line4 = infile.readline()
        if len(line) > 0:
            package_1 = parsePileup(line, lastrefseg)
        if len(line2) > 0:
            package_2 = parsePileup(line2, package_1['refseg'])
        if len(line3) > 0:
            package_3 = parsePileup(line3, package_2['refseg'])
        if len(line4) > 0:
            package_4 = parsePileup(line4, package_3['refseg'])

infile.close()
outfile.close()
