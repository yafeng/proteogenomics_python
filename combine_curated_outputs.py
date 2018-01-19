import sys
import pandas as pd
import os

if len(sys.argv)<2:
    sys.exit("USAGE: python " + sys.argv[0] + "<file1> <file2> ... <fileN>  > output.txt")

pepTable = pd.read_csv(sys.argv[1],sep="\t")

dfmerge = pepTable
for idx in range(2,len(sys.argv)):
    if os.path.exists(sys.argv[idx]):
        df = pd.read_csv(sys.argv[idx],sep="\t")
        dfmerge = pd.merge(left=dfmerge,right=df, how='left', left_on='Peptide', right_on='Peptide')
    else:
        print ("file %s doesn't exist" % sys.argv[idx])

dfmerge.to_csv("novel_peptides.combined.results.txt",sep="\t",index=F)
