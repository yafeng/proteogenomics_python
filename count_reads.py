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

## Function definitions.

# Determine if an alignment has few enough mismatches to pass.
def mismatches_ok(aln, max_mismatches=1):
    try:
        nm = aln.get_tag('nM') 
    except:
        try:
            nm = aln.get_tag('NM') 
        except:
            return(-1) # Could not find a tag for number of mismatches
    return (nm <= max_mismatches)

# Is pairing OK? Single-end passes automatically; paired-end passes if properly paired.
def pairing_ok(aln):
    if not aln.is_paired: 
        return True
    elif aln.is_proper_pair:
        return True
    else:
        return False
    
# Is the read multi-mapping? Fail if so
def multimapping_ok(aln, max_loci=1):
    try:
        if aln.get_tag('NH') > max_loci: 
            return False
        else:
            return True
    except:
        try:
            if aln.get_tag('XT') == 'R': # XT:A:U means unique, XT:A:R means repeat, apparently 
                return False 
            else: 
                return True
        except:
            return(-1) # Could not find a tag for number of loci


def find_eligible_alns(region, bamfile, max_mismatches=1):
    good_alns = []
    try:
        iter = bamfile.fetch(region[0], region[1], region[2])
    except:
        sys.exit("Region " + region[0] + ' ' + str(region[1]) + ' ' + str(region[2]) + '\nBAM file' + bamfile + '\nMake sure that you have an indexed BAM file!')
    for x in iter:
        if mismatches_ok(x) and pairing_ok(x) and multimapping_ok(x):        
            good_alns.append(x)
    return(good_alns)




################  Comand-line arguments ################


if len(sys.argv[1:])<=1:  ### Indicates that there are insufficient number of command-line arguments
    print("Warning! wrong command, please read the manual in Readme.txt.")
    print("Example: python scanBam.py --gff_input novelpep.gff3 --bam_files bam_files_list.txt --output novelpep_readcount.txt")
    sys.exit()
else:
    options, remainder = getopt.getopt(sys.argv[1:],'', ['gff_input=',
                                                         'bam_files=',
                                                         'output='])
    for opt, arg in options:
        if opt == '--gff_input': gff_file=arg
        elif opt == '--bam_files': bam_files=arg
        elif opt == '--output': out_file=arg
        else:
            print("Warning! Command-line argument: %s not recognized. Exiting...") % opt; sys.exit()


### pair-end reads need to be properly paired
### reads map to single location
### reads without mismatch
### reads mapped to the same strand of the peptide
### splice juncton peptides are only counted supported when found with corresponding
### reads overlap with peptide with minimum 1 nucleotide

input1=open(gff_file,"r")
input2=open(bam_files,"r")
output=open(out_file,"w")

pep_dic={} # store peptide object with sequence as key
for line in input1:
    if line[0]!="#":
        row=line.strip().split("\t")
        seq=row[8].replace("ID=","")
        if seq not in pep_dic:
            pep_dic[seq]=PEPTIDE(ID=seq,seq=seq,chr=row[0],start=row[3],end=row[4],strand=row[6],type="continuous")
        else:
            pep_dic[seq].type="spliced"
            if pep_dic[seq].start==row[3]:
                pep_dic[seq].splice_start=row[4]
            elif pep_dic[seq].end==row[4]:
                pep_dic[seq].splice_start=row[3]


print(len(pep_dic))

aln_table = {}

for line in input2:
    # add code to process bam files
    bam = line.strip()
    sys.stderr.write(bam + '\n')
    aln_count = {} # A dictionary that will collect eligible alignment counts by peptide.
    bamfile = AlignmentFile(bam,"rb")
    for seq in pep_dic:
        peptide=pep_dic[seq]
        region = [peptide.chr, int(peptide.start), int(peptide.end)]
        if peptide.type=="continuous":
            eligible = find_eligible_alns(region, bamfile)
        elif peptide.type=="spliced":
            pass
        aln_count[peptide.seq]=len(eligible)
    aln_table[bam] = aln_count

input1.close()
input2.close()

bam_files = sorted(list(aln_table.keys()))

# Write output file header.
output.write('sequence\t')
for i in range(0,len(bam_files)):
        output.write(bam_files[i].split('.')[0].split('/')[-1])
        if (i < (len(bam_files)-1)): output.write('\t')
        else: output.write('\n')



# And the counts.
for pep in sorted(aln_table[bam_files[0]].keys()):
    output.write(pep + '\t')
    for i in range(0,len(bam_files)): 
        output.write(str(aln_table[bam_files[i]][pep]))
        if (i < (len(bam_files)-1)): output.write('\t')
    output.write('\n')








                                                         

                                                         
