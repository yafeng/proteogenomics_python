import sys
import os
import getopt
import re


################  Comand-line arguments ################


if len(sys.argv[1:])<=1:  ### Indicates that there are insufficient number of command-line arguments
    print("Warning! wrong command, please read the manual in Readme.txt.")
    print("Example: python to_fasta.py --input input_file.txt --pep_seq_column 6 --output output.fasta")
    sys.exit()
else:
    options, remainder = getopt.getopt(sys.argv[1:],'', ['input=','pep_seq_column=','output='])
    for opt, arg in options:
        if opt == '--input': input_file=arg
        elif opt == '--pep_seq_column': pep_seq_column=int(arg)
        elif opt == '--output': out_file=arg
        else:
            print("Warning! Command-line argument: %s not recognized. Exiting...") % opt; sys.exit()

input=open(input_file,'r') # tab demilited text file including peptide seqeunces
output=open(out_file,'w')
input.readline()

dic={}


for line in input:
    row=line.strip().split("\t")
    pep=re.sub('[^A-Z]','',row[pep_seq_column-1])
    if pep not in dic:
        output.write(">%s\n%s\n" % (pep,pep))
        dic[pep]=1

input.close()
output.close()
