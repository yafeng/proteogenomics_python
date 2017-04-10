import sys
from collections import OrderedDict

class Peptide(object):
    def __init__(self,chr=None,strand=None,start=0,end=0,psm=0,content=None):
        self.start=start
        self.end=end
        self.strand=strand
        self.chr=chr
        self.content=content

input1=open(sys.argv[1],"r")
output=open(sys.argv[2],"w")

peplist=[]
output.write("loci_tag\t"+input1.readline())

for line in input1:
    row=line.strip().split("\t")
    chr=row[2][3:].replace("X","23").replace("Y","24")
    peplist.append(Peptide(chr=chr,start=int(row[3]),end=int(row[4]),strand=row[5],content=line))


## sort by chr and start cor
peplist.sort(key=lambda x:map(int,(x.chr,x.start)))

print "total peptide number",len(peplist)

loci_group=OrderedDict()
loci_tag=0

for pep in peplist:
    found_neighbor="NO"
    for loci in loci_group:
        if pep.chr==loci_group[loci][0].chr:
            if abs(pep.start-loci_group[loci][0].end)<10000:
                loci_group[loci].append(pep)
                found_neighbor="YES"
            elif abs(pep.end-loci_group[loci][0].start)<10000:
                loci_group[loci].append(pep)
                found_neighbor="YES"

    if found_neighbor=="NO":
        loci_tag+=1
        loci_group[loci_tag]=[pep]

print "total coding loci",len(loci_group)
pep_count=[]
loci_list=[]

for loci in loci_group:
    peptide_groups=loci_group[loci]
    pep_count.append(len(peptide_groups))
    for pep in peptide_groups:
        loci_list.append(loci)
        output.write(str(loci)+"\t"+pep.content)

print 'loci Peptide count'
for count in sorted(list(set(pep_count))):
    print count,pep_count.count(count)


print 'the following Loci supported by more than 2 peptides'

dic={}
for loci in loci_list:
    if loci_list.count(loci)>=2:
        if loci not in dic:
            print loci,loci_list.count(loci)
            dic[loci]=1


input1.close()
output.close()








