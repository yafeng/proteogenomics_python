import sys
import re
import os
import getopt

def get_psm(infile):
    dic={}
    
    for line in infile:
        row=line.strip().split("\t")
        pep=re.sub("[\W\d]","",row[index].strip())#strip modifications
        if pep not in dic:
            dic[pep]=[line]
        else:
            dic[pep].append(line)
    return dic


################  Comand-line arguments ################
peptide_column="Peptide" #default name to look for peptide column in input file

if len(sys.argv[1:])<=1:  ### Indicates that there are insufficient number of command-line arguments
    print "Warning! wrong command, please read the mannual in Readme.txt."
    print "Example: python parse_blastp_output.py --input_psm PSM_filename --input_table novpep_table.txt --output output_filename"
else:
    options, remainder = getopt.getopt(sys.argv[1:],'', ['input_psm=','peptide_column=','input_pep=','output='])
    for opt, arg in options:
        if opt == '--input_psm': input_file=arg
        elif opt == '--peptide_column': peptide_column=arg
        elif opt == '--input_pep': table_file=arg
        elif opt == '--output': output_file=arg
        else:
            print "Warning! Command-line argument: %s not recognized. Exiting..." % opt; sys.exit()

input1=open(input_file,"r") # Novel peptide psm table with TMT intensity

header1=input1.readline().strip().split("\t")
print header1
if peptide_column not in header1:
    print "%s not in header" % peptide_column;sys.exit()
else:
    index=header1.index(peptide_column)

dic1=get_psm(input1)
input2=open(table_file,"r") # 1 mismatch novel peptide list
output=open(output_file,"w") # 1 mismatch novel peptide PSM
newheader=header1+["sub_type","sub_pos"]
output.write("\t".join(newheader)+"\n")

input2_header=input2.readline().split("\t")
index1=input2_header.index("category")
index2=input2_header.index("sub_type")
index3=input2_header.index("sub_pos")

for line in input2:
    row=line.strip().split("\t")
    pep=row[0]
    category=row[index1]
    sub_type=row[index2]
    sub_pos=row[index3]
    if "1 mismatch" in category:
        if sub_type=="N>D": # check if it is N>D deamination
            print line
            continue
        else:
            if pep in dic1:
                for entry in dic1[pep]:
                    s=entry.strip().split("\t")
                    s+=[sub_type,sub_pos]
                    output.write("\t".join(s)+"\n")

input1.close()
input2.close()
output.close()
