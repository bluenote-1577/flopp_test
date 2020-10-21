import numpy as np
import pickle
from collections import defaultdict
import scipy.stats
import itertools
import sys

frags_file_short = sys.argv[1]
corr_frags_file = sys.argv[2]
num_variants = sys.argv[3]

l = 2
epsilon = 0.01
ploidy = 3
pval = 0.05

#frags_file_short= "/home/jshaw/summer_project_phasing/pot_dos_files/illumina_10x_frags.p"
#corr_frags_file = "qual_illumina_10x_frags.txt"

f_l = open(frags_file_short,'rb')
R_l = pickle.load(f_l)
ex = pickle.load(f_l)
R_ld = pickle.load(f_l)
R_haplo_list = pickle.load(f_l)
Q_ld = pickle.load(f_l)



key_to_blockinfo = dict()
counter = 0

for i,q_d in enumerate(Q_ld):
    pos_to_del = set()
    for pos in q_d:
        if ord(q_d[pos]) < 50:
            pos_to_del.add(pos)

    for pos in pos_to_del:
        q_d.pop(pos)
        R_ld[i].pop(pos)


for r_d in R_ld:
    blocks = []
    prev_block_index = -1
    current_block = []
    for pos in r_d:
        if pos > prev_block_index+1 and prev_block_index != -1:
            blocks.append(current_block)
            current_block = [pos]
            prev_block_index = pos
        else:
            current_block.append(pos)
            prev_block_index = pos

    if len(current_block) > 0:
        blocks.append(current_block)
    if len(blocks) == 0:
        counter+=1
        continue
    if len(blocks) < 2 and len(blocks[0]) == 1:
        counter+=1
        continue
    else:
        key_to_blockinfo[counter] = blocks

    counter+=1


with open(corr_frags_file,'w') as f:
    for key in key_to_blockinfo:
        blocks = key_to_blockinfo[key]
        blockseq = R_ld[key]
        s3 = str(len(blocks)) + "\t" + str(key) + "hap" + str(R_haplo_list[key])
        qual_str = "\t"
        for block in blocks:
            s3+= "\t" + str(block[0]) + "\t"
            for index in block:
             #   print(blockseq[index])
                s3 += str(blockseq[index])
        for i in range(len(blockseq)):
            qual_str += 'I'
        f.write(s3+qual_str)
        f.write('\n')

