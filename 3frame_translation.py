import sys
from Bio import Seq,SeqIO
from Bio.Alphabet import IUPAC

handle=open(sys.argv[1],'r')
output=open(sys.argv[2],'w')

for record in SeqIO.parse(handle,'fasta'):
    seq=record.seq
    ORF1=seq.translate()
    ORF2=seq[1::].translate()
    ORF3=seq[2::].translate()
    
    if record.id=="":
        print(record.description, record.seq)
        continue
    output.write("%s\n%s\n" % ('>'+record.id+'_ORF1',ORF1))
    output.write("%s\n%s\n" % ('>'+record.id+'_ORF2',ORF2))
    output.write("%s\n%s\n" % ('>'+record.id+'_ORF3',ORF3))
            
handle.close()

output.close()
