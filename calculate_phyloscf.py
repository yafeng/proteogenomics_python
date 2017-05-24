'''
    the script is written by Mikael Hussius @ SciLifeLab, https://github.com/hussius/gff-phylocsf-human
'''


import pyBigWig as pw
import sys

def predict_coding(vec):
    coding = "no"
    for v in vec:
        if not v: continue
        if v > 0: coding = "yes"
    return(coding)

if len(sys.argv)<2:
    sys.exit("USAGE: python " + sys.argv[0] + " <GFF file>")

regs = []
chrom={}
starts={}
ends={}
peptide={}

for line in open(sys.argv[1]):
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

rpathbase = "https://data.broadinstitute.org/compbio1/PhyloCSFtracks/hg19/latest/PhyloCSF"

scores = {}

for rf in ["+0","+1","+2","-0","-1","-2"]:
    sys.stderr.write("Reading frame " + rf + "\n")
    rpath = rpathbase + rf + ".bw"
    bw = pw.open(rpath)
    frame_score = {}
    count = 0
    for r in regs:
        count += 1
        if(count % 50 ==0): sys.stderr.write('\tProcessed ' + str(count) + " peptides out of " + str(len(regs)) + "\n")
        sys.stderr.flush()
        score = bw.stats(chrom[r], starts[r], ends[r])[0]
        frame_score[r] = score
        scores[rf] = frame_score

print("\t".join(["chromosome","start","end","peptide_seq","+0","+1","+2","-0","-1","-2","overall_prediction"]))
for r in regs:
    scoreList = [scores["+0"][r], scores["+1"][r], scores["+2"][r], scores["-0"][r], scores["-1"][r], scores["-2"][r]]
    coords = [chrom[r], str(starts[r]), str(ends[r])]
    print("\t".join(coords) + "\t" + peptide[r] + "\t" + "\t".join(map(str,scoreList)) + "\t" + predict_coding(scoreList))

