import sys

ranbow_file = sys.argv[1]
output = sys.argv[2]

f = open(ranbow_file,'r')
next(f)
H = []
for line in f:
    hap = []
    hap_str = line.split()[3]
    Q_str = line.split()[6]
    print(Q_str)
    start = Q_str[3]
    for i in range(int(start) - 1):
        hap.append('-')
    for char in hap_str:
        hap.append(char)
    H.append(hap)

ploidy = len(H)
print(ploidy)
out = open(output,'w')
out.write("****\n")
for i in range(len(H[1])):
    out.write(str(i+1)+'\t')
    for k in range(ploidy):
        print(i,k)
        out.write(str(H[k][i]) + '\t')

    out.write('\n')
out.write("******")
