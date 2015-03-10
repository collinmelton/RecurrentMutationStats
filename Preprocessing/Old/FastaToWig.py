from optparse import OptionParser

def processWigHeader(line):
    chrom=line.split()[2].split(":")[2]
    start=line.split()[2].split(":")[3]
    step=line.split()[2].split(":")[5]
    result="fixedStep chrom="+chrom+" start="+start+" step="+step
    return(chrom, result)

def run(input, output):
    print "running with input, output:", (input, output)
    f=open(input, "r")
    g=open(output+"0.wig", "w")
    line=f.readline().strip()
    while line!="":
        if ">" in line:
            g.close()
            chrom, header=processWigHeader(line)
            print "writing file: ", output+chrom+".wig"
            g=open(output+chrom+".wig", "w")
            g.write(header)
        else:
            #if g.tell()!=: g.write("\n")
            for char in line:
                if char=="C" or char=="G":
                    g.write("\n2")
                elif char=="T" or char=="A":
                    g.write("\n1")
                else:
                    g.write("\n0")
        line=f.readline().strip()
    f.close()
    g.close()

def getOptions():
    parser = OptionParser()
    parser.add_option("--I", dest="input", help = "path plus input file prefix",
                      metavar = "FILE", type = "string", default = "/srv/gs1/projects/snyder/collinmelton/bundle/2.3/b37/d5/hs37d5.fa")
    parser.add_option("--O", dest="output", help="path and name of output file",
                      metavar="FILE", type="string", default="/srv/gs1/projects/snyder/collinmelton/bundle/2.3/b37/d5/hs37d5.CG_TA.")
    (options, args) = parser.parse_args()
    return options

# this is the main function
def runMain():
    options = getOptions()
    run(options.input, options.output)    

runMain()