import sys

if len(sys.argv[1:])<=1:  ### Indicates that there are insufficient number of command-line arguments
    print "Warning! wrong command, please read the mannual in Readme.txt."
    print "Example: python prepare_annovar_input.py --input example_vardb_6rf_novpep.hg19cor.txt --output example_novpep_avinput.txt"
else:
    options, remainder = getopt.getopt(sys.argv[1:],'', ['input=',
                                                         'output='])
    for opt, arg in options:
        if opt == '--input': input_file=arg
        elif opt == '--output': output_file=arg
        else:
            print "Warning! Command-line argument: %s not recognized. Exiting..." % opt; sys.exit()

input=open(input_file,"r")
output=open(output_file,"w")

header=input.readline().split("\t")

index1=header.index("chr")
index2=header.index("start")
index3=header.index("end")
index4=header.index("strand")

for line in input:
    row=line.strip().split("\t")
    try:
        pep=row[0]
        chr=row[index1]
        start=row[index2]
        end=row[index3]
        strand=row[index4]
    
        output.write("%s\t%s\t%s\t%s\t%s\tComments:Seq=%s\n" % (chr,start,start,"A","-",pep))
        output.write("%s\t%s\t%s\t%s\t%s\tComments:Seq=%s\n" % (chr,end,end,"A","-",pep))
    except IndexError:
        print "IndexError",line

input.close()
output.close()



