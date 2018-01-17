import sys
import os
import getopt
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC

def OUTPUT_PEPTIDE(a,b,c,d,pep,r): # only output peptides in a certain length range defined by r
    if len(pep) in r:
        output.write("%s_%s_%s_%s\t%s\n" % (a,b,c,d,str(pep).replace("*","")))

def translate_trypsin(seq,strand,start,end):
    peptide=""
    for i in range(0,len(seq),3):
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
    print("Warning! wrong command, please read the mannual in Readme.txt.")
    print("Example: python sixframetranslation.py --input genome.fasta --output genome.6FT.txt --nuclear_trans_table 1 --mito_trans_table 2 --min_length 8 --max_length 30")
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
            print("Warning! Command-line argument: %s not recognized. Exiting..." % opt); sys.exit()

genome_dict = SeqIO.index(input_file, "fasta") #whole genome sequence
output=open(output_file,'w') #tab demiliated.txt

cut_site={"K":1,"R":1,"*":1}
length_range=list(range(min_len,max_len+1))

for record in genome_dict:
    coding_dna=genome_dict[record].seq
    chr=record
    chr_len=len(coding_dna)
    
    f1_seq=coding_dna
    r1_seq=coding_dna.reverse_complement()
    
    N_count_5prime=0  # count NNN from 5 prime end
    N_count_3prime=0  # count NNN from 3 prime end
    for base in f1_seq:
        if str(base)=="N":
            N_count_5prime+=1
        else:
            break
    start=N_count_5prime+1

    for base in r1_seq:
        if str(base)=="N":
            N_count_3prime+=1
        else:
            break

    end=chr_len-N_count_3prime

    ORF1_seq=f1_seq[N_count_5prime:chr_len-N_count_3prime] #remove NNNN from DNA sequence
    ORF4_seq=r1_seq[N_count_3prime:chr_len-N_count_5prime] #remove NNNN from DNA sequence
    translate_trypsin(ORF1_seq,"+",start,end)
    translate_trypsin(ORF4_seq,"-",start,end)


    ORF2_seq=f1_seq[N_count_5prime+1:chr_len-N_count_3prime]
    ORF5_seq=r1_seq[N_count_3prime+1:chr_len-N_count_5prime]
    start+=1
    end=end-1

    translate_trypsin(ORF2_seq,"+",start,end)
    translate_trypsin(ORF5_seq,"-",start,end)


    ORF3_seq=f1_seq[N_count_5prime+2:chr_len-N_count_3prime]
    ORF6_seq=r1_seq[N_count_3prime+2:chr_len-N_count_5prime]
    start+=1
    end=end-1
    
    translate_trypsin(ORF3_seq,"+",start,end)
    translate_trypsin(ORF6_seq,"-",start,end)

output.close()
