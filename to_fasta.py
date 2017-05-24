import sys

input=open(sys.argv[1],'r') # Novel peptide list
output=open(sys.argv[2],'w')
input.readline()

dic={}

for line in input:
    row=line.strip().split("\t")
    pep=row[0]

    output.write(">%s\n%s\n" % (pep,pep))

input.close()
output.close()
