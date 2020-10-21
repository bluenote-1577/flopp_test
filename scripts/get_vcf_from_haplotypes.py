from collections import defaultdict
import sys

files = sys.argv[3:]
ref = sys.argv[1]
outfile = sys.argv[2]

#files = ['./hap_files//pot_dos.fa_hap1.fa','./hap_files//pot_dos.fa_hap2.fa','./hap_files//pot_dos.fa_hap3.fa','./hap_files//pot_dos.fa_hap4.fa']

ref_file = open(ref,'r')
ref_name = next(ref_file)[1:].split()[0]
ref_string = next(ref_file)

strings = []

ploidy = len(files)

for file in files:
    f = open(file,'r')
    next(f)
    s = next(f)
    strings.append(s)

prev = 1
hist = []
num = 0
variants = []
for i in range(len(strings[0])):
    pos = i+1
    list_chars = []
    chars = defaultdict(int)
    for j in range(ploidy):
        chars[strings[j][i]] += 1
        list_chars.append(strings[j][i])

    if len(chars) > 1:
#        print(chars,'VAR',pos)
#        print(pos - prev)
        hist.append(pos-prev)
        prev = pos
        num+=1
        ref = ref_string[i]
        genotype = []
        for key in chars:
            if key != ref:
                alt = key
                ac = chars[alt]
        for char in list_chars:
            if char == alt:
                genotype.append(1)
            else:
                genotype.append(0)


        variants.append({'pos' : pos,'ref' : ref, 'alt' : alt, 'AC' : ac, 'gt' : genotype})

print(sum(hist)/float(num))
print(variants[0])

out = open(outfile,'w')
out.write("##fileformat=VCFv4.1\n")
out.write('##FILTER=<ID=PASS,Description="All filters passed">\n')
out.write('##contig=<ID=%s>\n' %ref_name)
outstr = '##INFO=<ID=AC,Number=A,Type=Integer,Description="Total number of alternate alleles in called genotypes">\n##INFO=<ID=AN,Number=1,Type=Integer,Description="Total number of alleles in called genotypes">\n##INFO=<ID=CIGAR,Number=A,Type=String,Description="The extended CIGAR representation of each alternate allele, with the exception that = is replaced by M to ease VCF parsing.  Note that INDEL alleles do not have the first matched base (which is provided by default, per the spec) referred to by the CIGAR.">\n##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
out.write(outstr)
out.write('#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	NAMETBD\n')


for variant in variants:
    out_string = '%s\t' %(ref_name)
    out_string += str(variant['pos']) + '\t'
    out_string += '.'+ '\t'
    out_string += variant['ref']+ '\t'
    out_string += variant['alt']+ '\t'
    out_string += str(100)+ '\t'
    out_string += 'PASS'+ '\t'
    out_string += 'AC=' + str(variant['AC']) +';'+'AN=4;' + 'CIGAR=1X'+'\t'
    out_string += 'GT\t'
    for g in variant['gt']:
        out_string += str(g) + '|'
    out_string = out_string[:-1]
    out_string += '\n'
    out.write(out_string)

