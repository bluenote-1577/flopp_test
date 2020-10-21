import numpy as np
import sys
from glob import glob

##TODO CHANGE FILE NAME FROM TEST

pswitch = 0.01
hap_files_dir = sys.argv[1]
variant_haplos_file = glob(hap_files_dir+"/*varianthaplos.txt")[0]
test_str = ''

hap_file = open(variant_haplos_file,'r')
header = next(hap_file)
hap_list = []
for line in hap_file:
    splitted = line.split()
    pos_all = splitted[2]
    ref_all = splitted[3]
    alt_all = splitted[4]
    hap_list.append((pos_all,ref_all,alt_all))


current_pos = 0
switch_pos1 = []

while current_pos < len(hap_list):
    z = np.random.geometric(p=pswitch)
    current_pos += z
    if current_pos < len(hap_list):
        switch_pos1.append(z)

current_pos = 0
switch_pos2 = []
while current_pos < len(hap_list):
    z = np.random.geometric(p=pswitch)
    current_pos += z
    if current_pos < len(hap_list):
        switch_pos2.append(z)

h1 = []
h2 = []
h3 = []
h4 = []
H =[h1,h2,h3,h4]
for (k,switch) in enumerate(switch_pos1):
    for i in range(switch):
        if k%2 == 0:
            h1.append(0)
            h2.append(1)
        else:
            h1.append(1)
            h2.append(0)

for k in range(len(hap_list)-len(h1)):
    h1.append(0)
    h2.append(1)

for (k,switch) in enumerate(switch_pos2):
    for i in range(switch):
        if k%2 == 0:
            h3.append(0)
            h4.append(1)
        else:
            h3.append(1)
            h4.append(0)

for k in range(len(hap_list)-len(h3)):
    h3.append(0)
    h4.append(1)

hap_file = open(variant_haplos_file+test_str,'w')

hap_file.write(header)
for i in range(len(hap_list)):
    hap_file.write(str(i+1)+'\tchr01\t')
    hap_file.write(hap_list[i][0] + '\t' + hap_list[i][1] + '\t' + hap_list[i][2] + '\t')
    for hap in H:
        hap_file.write(str(hap[i])+'\t')
    hap_file.write('\n')

hap_fastas = glob(hap_files_dir+"/*.fa")
hap_fastas = sorted(hap_fastas)
print(hap_fastas)

for (i,hap_fasta) in enumerate(hap_fastas):
    f = open(hap_fasta,'r')
    head = next(f)
    seq = next(f)
    new_seq = ""
    f.close()
    hap_encode = H[i]
    list_of_change_pos = set()
    change_pos_allles = dict()
    for j in range(len(hap_encode)):
        zero_ind_genome_pos = int(hap_list[j][0]) - 1
        list_of_change_pos.add(zero_ind_genome_pos)
        if hap_encode[j] == 1:
            change_pos_allles[zero_ind_genome_pos] = hap_list[j][2]
        else:
            change_pos_allles[zero_ind_genome_pos] = hap_list[j][1]

    for (j,char) in enumerate(seq):
        if j not in list_of_change_pos:
            new_seq += char
        else:
            new_seq += change_pos_allles[j]
    f = open(hap_fasta+test_str,'w')
    f.write(head)
    f.write(new_seq)

