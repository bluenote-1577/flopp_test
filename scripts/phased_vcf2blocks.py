import pickle
import numpy as np
from collections import defaultdict
import sys
import itertools

ploidy = int(sys.argv[2])
input= sys.argv[1]
output= input + ".txt"
total_ham_rate = 0
total_swer = 0

blocks = []
res = open(output,'w')
inp = open(input,'r')
current_tetra_block = list()
for i in range(ploidy):
    current_tetra_block.append([])
tag = None

for line in inp:
    if line[0] == '#':
        continue
    else:
        genops = line.split()[9]
        haplotag = genops.split(':')
        new_tag = haplotag[-1]
        haplotypes = haplotag[0]
        if '/' in haplotypes:
            for hap in current_tetra_block:
                hap.append(-1)
        else:
            if new_tag != tag and tag != None:
                blocks.append(current_tetra_block)
                current_tetra_block = list()
                for i in range(ploidy):
                    current_tetra_block.append([])
            tag = new_tag
            haplotypes = haplotag[0]

            haplotypes = haplotypes.split('|')
            for i in range(ploidy):
                current_tetra_block[i].append(int(haplotypes[i]))
        

blocks.append(current_tetra_block)
total_len = 0
num_blocks = 0
hist = []
for block in blocks:
    total_len += len(block[0])
    #if total_len > 1280:
    #    break
    num_blocks +=1
    hist.append(len(block[0]))

print(num_blocks)
print(hist)
blocks = np.array(blocks)

N50_inter = 0
N50 = 0
truth_coords = 0
sorted_lengths = sorted(hist,reverse = True)
for l in sorted_lengths:
    N50_inter += l
    if N50_inter > 9513/2:
        N50 = l
        break

print(N50,'N50!!')

res.write("***WHP\n")
start_pos = 1
for block in blocks:
    res.write("***\n")
    for i in range(len(block[0])):
        res.write(str(start_pos) + '\t')
        for hap in block:
            res.write(str(hap[i]) + '\t')
        res.write('\n')
        start_pos += 1
    res.write("***\n")


