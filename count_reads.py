import sys
import os
import getopt
import numpy as np
import re
from pysam import *

class PEPTIDE(object):
    def __init__(self,ID=None,seq=None,type=None,chr=None,strand=None,start=0,end=0,splice_start=0,splice_end=0):
        self.ID=ID
        self.seq=seq
        self.type=type
        self.start=start  #chromosome start coordinate
        self.end=end  #chromosome end coordinate
        self.strand=strand
        self.splice_start=splice_start
        self.splice_end=splice_end
        self.chr=chr

################  Comand-line arguments ################


if len(sys.argv[1:])<=1:  ### Indicates that there are insufficient number of command-line arguments
    print "Warning! wrong command, please read the mannual in Readme.txt."
    print "Example: python scanBam.py --input_gff novelpep.gff3 --bam_files bam_files_list.txt --output novelpep_readcount.txt"
else:
    options, remainder = getopt.getopt(sys.argv[1:],'', ['gff_input=',
                                                         'bam_files=',
                                                         'output='])
     for opt, arg in options:
        if opt == '--gff_input': gff_file=arg
        elif opt == '--bam_files': bam_files=arg
        elif opt == '--output': out_file=arg
        else:
            print "Warning! Command-line argument: %s not recognized. Exiting..." % opt; sys.exit()


### pair-end reads need to be properly paired
### reads map to single location
### reads without mismatch
### reads mapped to the same strand of the peptide
### splice juncton peptides are only counted supported when found with corresponding
### reads overlap with peptide with minimum 1 nucleotide

input1=open(gff_file,"r")
input2=open(bam_files,"r")

pep_dic={} # store peptide object with sequence as key
for line in input1:
    if line[0]!="#":
        row=line.strip().split("\t")
        seq=row[8].replace("ID=","")
        if seq not in pep_dic:
            pep_dic[seq]=PEPTIDE(ID=seq,seq=seq,chr=row[0],start=row[3],end=row[4],strand=row[6],type="continous")
        else:
            pep_dic[seq].type="spliced"
            if pep_dic[seq].start==row[3]:
                pep_dic[seq].splice_start=row[4]
            elif pep_dic[seq].end==row[4]:
                pep_dic[seq].splice_start=row[3]


print len(pep_dic)

for line in input2:
    # add code to process bam files
    
    for seq in pep_dic:
        peptide=pep_dic[seq]
        if peptide.type="continous":

        elif peptide.type="spliced":




input1.close()
input2.close()








                                                         

                                                         
