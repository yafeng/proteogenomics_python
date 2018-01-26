'''
    the script is modified from Mikael Hussius @ SciLifeLab, https://github.com/hussius/gff-phylocsf-human
    
    download the following bigwig files first
    # wget https://data.broadinstitute.org/compbio1/PhyloCSFtracks/hg19/latest/PhyloCSF+0.bw
    # wget https://data.broadinstitute.org/compbio1/PhyloCSFtracks/hg19/latest/PhyloCSF+1.bw
    # wget https://data.broadinstitute.org/compbio1/PhyloCSFtracks/hg19/latest/PhyloCSF+2.bw
    # wget https://data.broadinstitute.org/compbio1/PhyloCSFtracks/hg19/latest/PhyloCSF-0.bw
    # wget https://data.broadinstitute.org/compbio1/PhyloCSFtracks/hg19/latest/PhyloCSF-1.bw
    # wget https://data.broadinstitute.org/compbio1/PhyloCSFtracks/hg19/latest/PhyloCSF-2.bw
'''

import sys
import os
import pyBigWig as pw

def predict_coding(vec):
    coding = "no"
    for v in vec:
        if not v: continue
        if v > 0: coding = "yes"
    return(coding)

if len(sys.argv)<4:
    sys.exit("USAGE: python " + sys.argv[0] + "<GFF file> <BigWig file path> <output file>")

infile = sys.argv[1]
bw_file_path = sys.argv[2]
outfile = sys.argv[3]

regs = []
chrom={}
starts={}
ends={}
peptide={}

for line in open(infile):
    if not line.startswith("chr"):
        continue
    fields = line.strip().split()
    (chr, start, end, pept) = (fields[0], fields[3], fields[4], fields[8])
    if not pept.startswith("Parent="): continue
    name = chr+":"+start+"-"+end
    chrom[name]=chr
    starts[name]=int(start)
    ends[name]=int(end)
    peptide[name]=pept.split("=")[1]
    regs.append(name)

scores = {}

rpathbase = os.path.join(bw_file_path,"PhyloCSF")

for rf in ["+0","+1","+2","-0","-1","-2"]:
    sys.stderr.write("Searching PhyloCSF reading frame " + rf + "\n")
    rpath = rpathbase + rf + ".bw"
    bw = pw.open(rpath)
    frame_score = {}
    count = 0
    for r in regs:
        count += 1
        if(count % 50 ==0): sys.stderr.write('\tProcessed ' + str(count) + " peptides out of " + str(len(regs)) + "\n")
        sys.stderr.flush()
        try:
            score = bw.stats(chrom[r], starts[r], ends[r])[0]
        except RuntimeError:
            pass
        frame_score[r] = score
        scores[rf] = frame_score


output = open(outfile,"w")
output.write("\t".join(["Peptide","chromosome","start","end","PhyloCSF+0.score","PhyloCSF+1.score","PhyloCSF+2.score","PhyloCSF-0.score","PhyloCSF-1.score","PhyloCSF-2.score","PhyloCSF_prediction"])+"\n")
for r in regs:
    scoreList = [scores["+0"][r], scores["+1"][r], scores["+2"][r], scores["-0"][r], scores["-1"][r], scores["-2"][r]]
    coords = [chrom[r], str(starts[r]), str(ends[r])]
    row = [peptide[r]]+coords+ ['NA' if x is None else str(x) for x in scoreList] + [predict_coding(scoreList)]
    output.write('\t'.join(row) + '\n')

