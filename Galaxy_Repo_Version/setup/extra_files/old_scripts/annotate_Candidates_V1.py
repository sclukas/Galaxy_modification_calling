
import os
import sys
from random import shuffle

# sys.argv[1] /home/akhelm/Lukas/Bioinformatics_Ralf/C_elegans/run2100/testLukas_newReference/hypoxia/profiles_renamed
modname = sys.argv[3]  # 'm1A' / 'nonm1A'


with open(sys.argv[1]) as fin, open(sys.argv[2], 'w') as fout:
	line_list = []
        
        fout.write("arrest_rate\tmism_rate\tgmism\ttmism\tcmism\tjump_rate_total\tmod_type\n")

        line = fin.readline()
        if "mism" in line:
            line = fin.readline()
        while len(line) > 0:
            split_list = (line[:-1]).split("\t")
            if len(split_list) >= 11:
# ref_seg	pos	refbase	cov	relmism	arrestrate	gmism	tmism	cmism	single_jump_rate_direct	single_jump_rate_delayed	double_jump_rate
                line_list.append(
                    split_list[5] + "\t" + split_list[4] + "\t" + split_list[6] + "\t" + split_list[7] + "\t" +
                    split_list[8] + "\t" + str(float(split_list[9]) + float(split_list[10]) + float(split_list[11]))
                    + "\t" + "'" + modname + "'" + "\n")
            line = fin.readline()

        shuffle(line_list)
        i = 0
        while i < len(line_list):
            fout.write(line_list[i])
            i += 1

# The End

