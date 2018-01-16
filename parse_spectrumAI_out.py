import sys
import re
import getopt


if len(sys.argv[1:])<=1:  ### Indicates that there are insufficient number of command-line arguments
    print "Warning! wrong command, please read the mannual in Readme.txt."
    print "Example: python parse_spectrumAI_out.py --spectrumAI_out specAI_file --input input_psms.txt --output output_filename"
else:
    options, remainder = getopt.getopt(sys.argv[1:],'', ['spectrumAI_out=',
                                                         'input=',
                                                         'output='])
    for opt, arg in options:
        if opt == '--spectrumAI_out': input1_file=arg
        elif opt == '--input':input2_file=arg
        elif opt == '--output': output_file=arg
        else:
            print "Warning! Command-line argument: %s not recognized. Exiting..." % opt; sys.exit()

input1=open(input1_file,"r")  ## SpectrumAI output
input2=open(input2_file,"r")  ## novelpep table, peptide sequence is in first column
output=open(output_file,"w")

header1=input1.readline().split("\t")
index1=header1.index("Peptide")
index2=header1.index("flanking_ions_support")


pep_yes={}  # found with b,y ion support

for line in input1:
    row=line.strip().split("\t")
    try:
        if row[index2]=="YES":
            pep=re.sub("[\W\d]","",row[index1].strip())
            pep_yes[pep]=1
    except IndexError:
        print line

header2=input2.readline().split("\t")
index3=header2.index("category")

output.write("\t".join(header2))

n=0
for line in input2:
    row=line.strip().split("\t")
    category=row[index3]
    pep=row[0]
    if "1 mismatch" in category:
        n+=1
        if pep in pep_yes:
            output.write(line)
    else:
        output.write(line)

input1.close()
input2.close()
output.close()

print "%d out of %d single substitution novel peptides passed SpectrumAI curation" % (len(pep_yes),n)



