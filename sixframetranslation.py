import sys
import os
import getopt
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC

def OUTPUT_PEPTIDE(a,b,c,d,pep,r): # only output peptides in a certain length range defined by r
    if len(pep) in r:
        output.write("%s_%s_%s_%s\t%s\n" % (a,b,c,d,pep))

def translate_trypsin(seq,strand):
    peptide=""
    start=1
    chr_len=len(seq)
    end=chr_len
    for i in range(0,chr_len,3):
        codon=seq[i:i+3]
        if len(codon)!=3:
            continue
        if chr=="MT":
            aa=codon.translate(table=mito_trans_table) # The Mitochondrial Code
        elif chr=="M":
            aa=codon.translate(table=mito_trans_table) # The Mitochondrial Code
        else:
            aa=codon.translate(table=nuclear_trans_table) # The standard Code
        peptide+=aa
        if (str(aa) in cut_site):
            if strand=="-":
                start=end-len(peptide)*3+1
                OUTPUT_PEPTIDE(chr,start,end,strand,peptide,length_range)
                peptide=""
                end=start-1
            
            else:
                end = start + len(peptide)*3-1
                OUTPUT_PEPTIDE(chr,start,end,strand,peptide,length_range)
                peptide=""
                start=end+1

    if len(peptide)>1: #peptides at the terminal without cuting site
        if strand=="-":
            start=end-len(peptide)*3+1
            OUTPUT_PEPTIDE(chr,start,end,strand,peptide,length_range)
            
        else:
            end = start + len(peptide)*3-1
            OUTPUT_PEPTIDE(chr,start,end,strand,peptide,length_range)



################  Comand-line arguments ################


if len(sys.argv[1:])<=1:  ### Indicates that there are insufficient number of command-line arguments
    print "Warning! wrong command, please read the mannual in Readme.txt."
    print "Example: python sixframetranslation.py --input genome.fasta --output genome.6FT.txt --nuclear_trans_table 1 --mito_trans_table 2 --min_length 8 --max_length 30"
else:
    options, remainder = getopt.getopt(sys.argv[1:],'', ['input=',
                                                         'output=',
                                                         'nuclear_trans_table=',
                                                         'mito_trans_table=',
                                                         'min_length=',
                                                         'max_length='])
    for opt, arg in options:
        if opt == '--input': input_file=arg
        elif opt == '--output':output_file=arg
        elif opt == '--nuclear_trans_table': nuclear_trans_table=int(arg)
        elif opt == '--mito_trans_table': mito_trans_table=int(arg)
        elif opt == '--min_length':min_len=int(arg)
        elif opt == '--max_length':max_len=int(arg)
        else:
            print "Warning! Command-line argument: %s not recognized. Exiting..." % opt; sys.exit()

genome_dict = SeqIO.index(input_file, "fasta") #whole genome sequence
output=open(output_file,'w') #tab demiliated.txt

cut_site={"K":1,"R":1,"*":1}
length_range=range(min_len,max_len+1)

for record in genome_dict:
    coding_dna=genome_dict[record].seq
    
    f1_seq=coding_dna
    r1_seq=coding_dna.reverse_complement()
    
    translate_trypsin(f1_seq,"+")
    translate_trypsin(r1_seq,"-")

    f2_seq=f1_seq[1::]
    r2_seq=r1_seq[1::]
    
    translate_trypsin(f2_seq,"+")
    translate_trypsin(r2_seq,"-")


    f3_seq=f1_seq[2::]
    r3_seq=r1_seq[2::]
    
    translate_trypsin(f3_seq,"+")
    translate_trypsin(r3_seq,"-")

output.close()
