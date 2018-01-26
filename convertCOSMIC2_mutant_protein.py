'''
    this script is make a fasta file of mutant protein sequence from COSMIC database.
    you need to download CosmicMutantExport.tsv and All_COSMIC_Genes.fasta from http://cancer.sanger.ac.uk/cosmic/download
'''

import sys
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
import re

class SNP(object):
    def __init__(self,gene=None,mRNA=None,dna_mut=None,aa_mut=None,type=None):
        self.gene=gene
        self.mRNA=mRNA
        self.aa_mut=aa_mut
        self.type=type
        self.dna_mut=dna_mut

COSMIC_CDS_DB=SeqIO.index(sys.argv[1],'fasta') # All_COSMIC_Genes.fasta

cosmic_input=open(sys.argv[2],'r') # CosmicMutantExport.tsv
cosmic_input.readline()

output=open(sys.argv[3],'w')

mutation_dic={}

print len(COSMIC_CDS_DB)
line_counter=1
for line in cosmic_input:
    if line_counter%10000==0:
        print line_counter
    line_counter+=1
    row=line.strip().split("\t")
    if "coding silent" in row[15]:
        continue;
    if "?" not in row[13]:
        snp=SNP(gene=row[0],mRNA=row[1],dna_mut=row[13],aa_mut=row[14],type=row[15])
        header="COSMIC_%s_%s_%s_%s_%s" % (snp.mRNA,snp.gene,snp.dna_mut,snp.aa_mut,snp.type.replace(" ",""))
        if header not in mutation_dic: #skip duplicates
            mutation_dic[header]=1
            try:
                seq=COSMIC_CDS_DB[snp.gene].seq
            except KeyError:
                print snp.mRNA
                continue;

            positions=re.findall(r'\d+',snp.dna_mut)
            if ">" in snp.dna_mut and len(positions)==1: # Substitution
                mut_dna=snp.dna_mut[-1]
                index=int(positions[0])-1
                seq_mut=seq[:index]+mut_dna+seq[index+1:]
            elif "ins" in snp.dna_mut:
                index=snp.dna_mut.index("ins")
                insert_dna=snp.dna_mut[index+3:]
                if insert_dna.isalpha():
                    ins_index1=int(positions[0])
                    seq_mut=seq[:ins_index1]+insert_dna+seq[ins_index1:]

            elif "del" in snp.dna_mut:
                index=snp.dna_mut.index("del")
                if len(positions)>1:
                    del_index1=int(positions[0])-1
                    del_index2=int(positions[1])
                    seq_mut=seq[:del_index1]+seq[del_index2:]
                else:
                    del_index1=int(positions[0])-1
                    seq_mut=seq[:del_index1]+seq[del_index1+1:]

            mut_pro_seq=seq_mut.translate(to_stop=True)

            entry=">%s\n%s\n" % (header,mut_pro_seq)
            output.write(entry)
    else:
        if "?" not in row[14]:
            snp=SNP(gene=row[0],mRNA=row[1],dna_mut=row[13],aa_mut=row[14],type=row[15])
            header="COSMIC_%s_%s_%s_%s_%s" % (snp.mRNA,snp.gene,snp.dna_mut,snp.aa_mut,snp.type.replace(" ",""))
            if header not in mutation_dic: #skip duplicates
                mutation_dic[header]=1
                try:
                    seq=COSMIC_CDS_DB[snp.gene].seq
                except KeyError:
                    print snp.mRNA
                    continue;

                positions=re.findall(r'\d+',snp.aa_mut)

                protein_seq=str(seq.translate(to_stop=True))

                if "Missense" in snp.type:
                    mut_aa=snp.aa_mut[-1]
                    index=int(positions[0])-1
                    mut_pro_seq=protein_seq[:index]+mut_aa+protein_seq[index+1:]
                elif "Nonsense" in snp.type:
                    index=int(positions[0])-1
                    mut_pro_seq=protein_seq[:index]
                elif "Insertion - In frame" in snp.type:
                    index=snp.aa_mut.index("ins")
                    insert_aa=snp.aa_mut[index+3:]
                    if insert_aa.isalpha():
                        ins_index1=int(positions[0])
                        mut_pro_seq=protein_seq[:ins_index1]+insert_aa+protein_seq[ins_index1:]
                elif "Deletion - In frame" in snp.type:
                    try:
                        index=snp.aa_mut.index("del")
                    except ValueError:
                        print snp.gene,snp.mRNA,snp.dna_mut,snp.aa_mut,snp.type
                        continue;
                    if len(positions)>1:
                        del_index1=int(positions[0])-1
                        del_index2=int(positions[1])
                        mut_pro_seq=protein_seq[:del_index1]+protein_seq[del_index2:]
                    else:
                        del_index1=int(positions[0])-1
                        mut_pro_seq=protein_seq[:del_index1]+protein_seq[del_index1+1:]
                elif "Complex" in snp.type and "frameshift" not in snp.type:
                    try:
                        index=snp.aa_mut.index(">")
                    except ValueError:
                        print snp.gene,snp.mRNA,snp.dna_mut,snp.aa_mut,snp.type
                        continue;
                    mut_aa=snp.aa_mut[index+1:]
                    if "deletion" in snp.type:
                        try:
                            del_index1=int(positions[0])-1
                            del_index2=int(positions[1])
                            mut_pro_seq=protein_seq[:del_index1]+mut_aa+protein_seq[del_index2:]
                        except IndexError:
                            print snp.gene,snp.mRNA,snp.dna_mut,snp.aa_mut,snp.type
                            continue;
                    elif "insertion" in snp.type:
                        try:
                            ins_index1=int(positions[0])-1
                        except IndexError:
                            print snp.gene,snp.mRNA,snp.dna_mut,snp.aa_mut,snp.type
                            continue;
                        mut_pro_seq=protein_seq[:ins_index1]+mut_aa+protein_seq[ins_index1+1:]
                    elif "compound substitution" in snp.type:
                        if "*" not in mut_aa:
                            try:
                                del_index1=int(positions[0])-1
                                del_index2=int(positions[1])
                                mut_pro_seq=protein_seq[:del_index1]+mut_aa+protein_seq[del_index2:]
                            except IndexError:
                                print snp.gene,snp.mRNA,snp.dna_mut,snp.aa_mut,snp.type
                                continue;
                        else:
                            try:
                                del_index1=int(positions[0])-1
                                del_index2=int(positions[1])
                                mut_pro_seq=protein_seq[:del_index1]+mut_aa.replace("*","")
                            except IndexError:
                                print snp.gene,snp.mRNA,snp.dna_mut,snp.aa_mut,snp.type
                                continue;

                entry=">%s\n%s\n" % (header,mut_pro_seq)
                output.write(entry)

print ("COSMIC contains in total",len(mutation_dic),"non redundant mutations")
cosmic_input.close()
output.close()
