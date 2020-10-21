import numpy as np
import pickle
import sys
haplo_string = sys.argv[1]
out_haplo_string = sys.argv[2]
ploidy = int(sys.argv[3])

print("THER ENEEDS TO BE A FILLER LINE ON TOP!! ASSUME START AT LINE 2")
with open(haplo_string,'r') as file:
    next(file)
    H = []
    for i in range(ploidy):
        H.append([])
    for row in file:
        sp = row.split()
        #print(sp)
        for i in range(5,5+ploidy):
            #print(sp[i],row)
            if sp[i] != '0' and sp[i] != '1' and sp[i] != '2' and sp[i] != '3':
                H[i-5].append(-1)
            else:
                H[i-5].append(int(sp[i]))

    #for i in range(length_of_haplotype - len(H[0])):
    #    for j in range(ploidy):
    #        H[j].append(-1)
    #print(H)

#print(np.array(H)[:,:30])
pickle.dump(H,open(out_haplo_string,'wb'))


