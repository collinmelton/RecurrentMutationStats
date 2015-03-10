from optparse import OptionParser


def getPositions(inputFile):
    f = open(inputFile, 'r')
    f.readline()
    f.readline()
    header=f.readline().strip().split("\t")
    line=f.readline()
    results={}
    while line !="":
        vals=line.strip().split("\t")
        if len(vals)>2:
            chrom, pos=vals[0], int(vals[1])
            if chrom not in results:
                results[chrom]=[]
            results[chrom].append(pos)
        line=f.readline()
    f.close()
    return results

def writeWigs(positions, outputFilePrefix, bychrom=True):
    chromstouse=['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', 'X', 'Y', 'MT']
    if bychrom:
        for chrom in chromstouse:
            if chrom in positions:
                print "writing", chrom
                positions[chrom].sort()
                f=open(outputFilePrefix+"_"+chrom+".wig", "w")
                for pos in positions[chrom]:
                    if f.tell()!=0: f.write("\n")
                    f.write("fixedStep chrom="+chrom+" start="+str(pos)+" step=1\n1")
                f.close()
    else:
        f=open(outputFilePrefix, "w")
        for chrom in chromstouse:
            if chrom in positions:
                print "writing", chrom
                positions[chrom].sort()    
                for pos in positions[chrom]:
                    if f.tell()!=0: f.write("\n")
                    f.write("fixedStep chrom="+chrom+" start="+str(pos)+" step=1\n1")
        f.close()

def run(inputFile, outputFilePrefix, bychrom=True):
    positions=getPositions(inputFile)
    writeWigs(positions, outputFilePrefix, bychrom)

# this functions gets the command line options for running the program
def getOptions():
    parser = OptionParser()
    parser.add_option("--I", dest = "inputFile", help = "input file in merged format",
                      metavar = "FILE", type = "string", default = "./TestData/01eef340-598c-4205-a990-cec190ac2ca5.unannotated.varscan.mutect.merge")
    parser.add_option("--O", dest="outputFilePrefix", help = "output file is a tab delimited file of key and count",
                      metavar = "FILE", type = "string", default = "./TestData/01eef340-598c-4205-a990-cec190ac2ca5.unannotated.varscan.mutect.merge.wig")
    parser.add_option("--BC", dest="bychrom", help = "whether to treat outputfileprefix as prefix then append chrom or as filename for all chroms",
                      metavar = "STRING", type = "string", default = "F")
    (options, args) = parser.parse_args()
    return options

# this is the main function
def runMain():
    # get options and defaults
    print "starting run"
    options = getOptions()
    print "about to run"
    print options.inputFile, options.outputFilePrefix
    run(options.inputFile, options.outputFilePrefix, bychrom=(options.bychrom=="T"))
    
runMain()
    