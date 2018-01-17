import sys
import re

input1=open(sys.argv[1],"r") # novel peptide list

input2=open(sys.argv[2],"r") # var search novel PSM table
input3=open(sys.argv[3],"r") # 6RF search novel psm table

output=open(sys.argv[4],"w") # Final novel peptide PSM

input1.readline()
pep_dic={}
for line in input1:
    row=line.strip().split("\t")
    pep=row[0]
    pep_dic[pep]=1

print(len(pep_dic))

output.write(input2.readline())

input3.readline()
spectra_dic={}
for line in input2:
    row=line.strip().split("\t")
    pep=re.sub("[\W\d]","",row[11].strip())
    specID=row[3]+row[5] # spectra file + scan number
    if pep in pep_dic:
        if specID not in spectra_dic:
            spectra_dic[specID]=1
            output.write(line)

for line in input3:
    row=line.strip().split("\t")
    pep=re.sub("[\W\d]","",row[11].strip())
    specID=row[3]+row[5] # spectra file + scan number
    if pep in pep_dic:
        if specID not in spectra_dic:
            spectra_dic[specID]=1
            output.write(line)


print(len(spectra_dic))
input1.close()
input2.close()
input3.close()
output.close()
