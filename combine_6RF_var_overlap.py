import sys
import re
import numpy as np

input1=open(sys.argv[1],"r") # Galaxy varDB search peptide table
input2=open(sys.argv[2],"r") # Galaxy 6RF search peptide table
output=open(sys.argv[3],"w")

output.write(input1.readline())
pepdic={}
for line in input1:
    row=line.strip().split("\t")
    pep=re.sub("[\W\d]","",row[0].strip())
    acc=row[1]
    row[0]=pep
    if acc[:3]=="lnc":
        if pep not in pepdic:
            pepdic[pep]=1
            output.write("\t".join(row)+"\n")
    elif acc[:6]=="PGOHUM":
        if pep not in pepdic:
            pepdic[pep]=1
            output.write("\t".join(row)+"\n")

input2.readline()
for line in input2:
    row=line.strip().split("\t")
    pep=re.sub("[\W\d]","",row[0].strip())
    acc=row[1]
    row[0]=pep
    if acc.count("chr")>1:
        print(pep,acc)
        continue
    if pep not in pepdic:
        pepdic[pep]=1
        output.write("\t".join(row)+"\n")


print("total number of unique novel peptides:",len(pepdic))
input1.close()
input2.close()
output.close()




