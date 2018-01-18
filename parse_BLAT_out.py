import sys

input=open(sys.argv[1],"r")  # blat output
input2=open(sys.argv[2],"r") # novel peptide gene cor file

output=open(sys.argv[3],"w")

# parse blat output and collect peptides mapped to multiple loci
pep_dic={}
pep_map={}
pep_multi_loci={}
pep_unique_loci={}

n=0
for line in input:
    row=line.strip().split("\t")
    n+=1
    if row[0].isdigit():
        pep=row[-2].replace(",","")
        pep_dic[pep]=1
        match=int(row[0])
        strand=row[8][1]
        qsize=int(row[10])
        qstart=int(row[11])
        qend=int(row[12])
        chr=row[13]
        chr_start=row[15]
        chr_end=row[16]
        if match==qsize and qend-qstart==qsize:
            if pep not in pep_map:
                pep_map[pep]=[chr+"_"+chr_start+"_"+chr_end+"_"+strand]
            else:
                pep_map[pep].append(chr+"_"+chr_start+"_"+chr_end+"_"+strand)

print ("the number of peptides input for blat:"),n
print ("blat qsize=match peptides"),len(pep_map)
print ("blat return peptides"),len(pep_dic)

for pep in pep_map:
    if len(pep_map[pep])>1:
        pep_multi_loci[pep]=pep_map[pep]
    elif len(pep_map[pep])==1:
        pep_unique_loci[pep]=pep_map[pep]

print ("blat multi map peptides"),len(pep_multi_loci)
print ("blat unique map peptides"),len(pep_unique_loci)

output.write(input2.readline())
for line in input2:
    row=line.strip().split("\t")
    pep=row[0]

    if pep in pep_multi_loci:
        print "multimapped peptide",line
    elif pep in pep_unique_loci:
        output.write("\t".join(row)+"\n")
    else:
        output.write("\t".join(row)+"\n")]

input.close()
input2.close()
output.close()
