#!/usr/bin/python

import os
import sys
import copy

# This V2 version of Pileup2Profile.py is minimized and optimized for publication


print "### Pileup2Profile..."


countdict = {}
charstring = "agtcnAGTCNXRKS.,*_"

for charx in charstring:
    countdict[charx] = 0


pile = ""
refseg = ""
pos = ""
cov = 0.0
splitlist1 = []
splitlist2 = []

jumpdict = {'0':{1:0, 2:0, 3:0}, '-1':{1:0, 2:0, 3:0}, '-2':{1:0, 2:0, 3:0}, '-3':{1:0, 2:0, 3:0}}

def parsePile():
    i = 0
    
    while (i < len(pile)):
        if pile[i] not in "+-^$":
            countdict[pile[i]] += 1
        else:
            if pile[i] in "+-":
                number = pile[i+1]
                j = i+1
                while(j+1 < len(pile) and pile[j+1] in "1234567890"):
                    number += pile[j+1]
                    j += 1
		if pile[i] == "-" and int(number) <= 3:
                    (jumpdict['0'])[int(number)] += 1
                i += 1 + int(number)
            elif pile[i] == "^": # $-case requires no action
                i += 1
        i += 1
    
    try:
        countdict[splitlist1[2]] += countdict["."] + countdict[","]
    except:
        countdict[splitlist1[2]] = countdict["."] + countdict[","]
        print "Unexpected symbol in splitlist[2]:", splitlist[2]


#with open(sys.argv[1], 'r') as infile, open(sys.argv[2], 'w') as outfile:
infile = open(sys.argv[1], "r")
outfile = open(sys.argv[2], "w")
outfile.write("ref_seg\tpos\tref_base\tcov\trel_mism\tA\tG\tT\tC\tN\ta\tg\tt\tc\tn\tsingle_jump_rate_direct\tsingle_jumprate_delayed\tdouble_jump_rate\tarrest_rate\n")



# reference_segment                                            position    reference_base      coverage    bases     quals
# Tb427.03.4393_brucei_Lister_427_serine_genomic_forward       65          G                   3           ,$,$,$    AB1

line1 = infile.readline()
splitlist1 = line1.split("\t")
line2 = infile.readline()
splitlist2 = line2.split("\t")

refseg = splitlist1[0]
nextrefseg = splitlist2[0]

while len(line1) > 0 and len(line2) > 0:
    
    while len(line1) > 0 and len(line2) > 0 and refseg == nextrefseg and int(splitlist1[1]) == int(splitlist2[1])-1:
        
        
        pile = splitlist1[4]
        refseg = splitlist1[0]
        pos = splitlist1[1]
        cov = float(splitlist1[3])
        
        parsePile()
        
        relmism = "0.0"
        single_jump_rate_direct = "0.0"
        single_jump_rate_delayed = "0.0"
        double_jump_rate = "0.0"
        if cov > 0.0:
            relmism = str((cov-countdict["."]-countdict[","]-(jumpdict['-1'])[1]-(jumpdict['-2'])[2])/cov)
            single_jump_rate_direct = (jumpdict['-1'])[1]/cov
            single_jump_rate_delayed = (jumpdict['-2'])[1]/cov
            double_jump_rate = (jumpdict['-2'])[2]/cov
        arrest = "0.0"
        if float(splitlist2[3]) > 0.0:
            arrest = str(((splitlist2[4]).count("^"))/float(splitlist2[3]))
        
        #                                                                                    does not count indels
        #             ref_seg\t               pos\t                    ref_base\t            cov\t                rel_mism\t                                                        A\t                                G\t                        T\t                        C\t                            N\t                            a\t                        g\t                        t\t                            c\t                            n\t            arrest_rate
        outfile.write(splitlist1[0] + "\t" + splitlist1[1] + "\t" + splitlist1[2] + "\t" + splitlist1[3] + "\t" + relmism + "\t" + str(countdict["A"]) + "\t" + str(countdict["G"]) + "\t" + str(countdict["T"]) + "\t" + str(countdict["C"]) + "\t" + str(countdict["N"]) + "\t" + str(countdict["a"]) + "\t" + str(countdict["g"]) + "\t" + str(countdict["t"]) + "\t" + str(countdict["c"]) + "\t" + str(countdict["n"]) + "\t" + str(single_jump_rate_direct) + "\t" + str(single_jump_rate_delayed) + "\t" + str(double_jump_rate) + "\t" + arrest + "\n")
        
        jumpdict['-3'] = copy.deepcopy(jumpdict['-2'])
        jumpdict['-2'] = copy.deepcopy(jumpdict['-1'])
	jumpdict['-1'] = copy.deepcopy(jumpdict['0'])
	jumpdict['0'] = {1:0, 2:0, 3:0}

        for charx in charstring:
            countdict[charx] = 0
        splitlist1 = splitlist2
        
        nextrefseg = splitlist2[0]
        refseg = splitlist1[0]
        line1 = line2
        line2 = infile.readline()
        splitlist2 = line2.split("\t")
    
    pile = splitlist1[4]
    refseg = splitlist1[0]
    pos = splitlist1[1]
    cov = float(splitlist1[3])
    
    parsePile()
    
    relmism = "0.0"
    single_jump_rate_direct = "0.0"
    single_jump_rate_delayed = "0.0"
    double_jump_rate = "0.0"
    if cov > 0.0:
        relmism = str((cov-countdict["."]-countdict[","]-(jumpdict['-1'])[1]-(jumpdict['-2'])[2])/cov)
        single_jump_rate_direct = (jumpdict['-1'])[1]/cov
        single_jump_rate_delayed = (jumpdict['-2'])[1]/cov
        double_jump_rate = (jumpdict['-2'])[2]/cov
    # last line of old contig
    #                                                                                    does not count indels
    #             ref_seg\t               pos\t                    ref_base\t            cov\t                rel_mism\t                                                        A\t                                G\t                        T\t                        C\t                            N\t                            a\t                        g\t                        t\t                            c\t                            n\t            arrest_rate
    
    
    outfile.write(splitlist1[0] + "\t" + splitlist1[1] + "\t" + splitlist1[2] + "\t" + splitlist1[3] + "\t" + relmism + "\t" + str(countdict["A"]) + "\t" + str(countdict["G"]) + "\t" + str(countdict["T"]) + "\t" + str(countdict["C"]) + "\t" + str(countdict["N"]) + "\t" + str(countdict["a"]) + "\t" + str(countdict["g"]) + "\t" + str(countdict["t"]) + "\t" + str(countdict["c"]) + "\t" + str(countdict["n"]) + "\t" + str(single_jump_rate_direct) + "\t" + str(single_jump_rate_delayed) + "\t" + str(double_jump_rate) + "\t" + "0.0" + "\n")
    
    jumpdict['-3'] = copy.deepcopy(jumpdict['-2'])
    jumpdict['-2'] = copy.deepcopy(jumpdict['-1'])
    jumpdict['-1'] = copy.deepcopy(jumpdict['0'])
    jumpdict['0'] = {1:0, 2:0, 3:0}
    
    for charx in charstring:
        countdict[charx] = 0
    
    splitlist1 = splitlist2
    
    nextrefseg = splitlist2[0]
    refseg = splitlist1[0]
    line1 = line2
    line2 = infile.readline()
    splitlist2 = line2.split("\t")
    

infile.close()
outfile.close()


