import SourceWig, OverlapWig
from optparse import OptionParser

def run(coverageFile, overlapWigFile, resultFile):
    print coverageFile
    print overlapWigFile
    cov=SourceWig.SourceWig(coverageFile)
    overlap=OverlapWig.OverlapWig(overlapWigFile)
    
    p=cov.walk()
    result={}
    i=0
    while p!=None:
        i+=1
        val, chrom, pos=p
        if i%100000==0: print chrom, pos 
        annotation=overlap.getWigResults(chrom, pos)
        if annotation==None:
            annotation="NA_NA_NA"
        #print chrom, pos, annotation
        key=str(val)+"_"+annotation
        if "\x00" not in key:
            if key not in result:
                result[key]=0
            result[key]+=1
        else:
            print "weirdo key:", chrom, pos, val, key, str(val), str(key)
        p=cov.walk()
        #if i==100000: break
    f=open(resultFile, "w")
    f.write("key\tcovered\tbpType\treplicationTiming\ttranscriptType\tcount")
    for key in result:
        f.write("\n"+key+"\t"+"\t".join(key.split("_"))+"\t"+str(result[key]))
    f.close()
    

# this functions gets the command line options for running the program
def getOptions():
    parser = OptionParser()
    parser.add_option("--I", dest = "inputFile", help = "input file in wig format",
                      metavar = "FILE", type = "string", default = "./TestData/coverage.wig.corrected.22.wig")
    parser.add_option("--O", dest="outputFile", help = "output file is a tab delimited file of key and count",
                      metavar = "FILE", type = "string", default = "./TestData/Chrom22Covariates.tsv")
    parser.add_option("--V", dest = "overlapWigFile", help = "path to overlap wig file"+
                      "concurrently", metavar = "FILE", default = "./TestData/mergedBPRepTimingTranscript.22.wig", type = "string")
    (options, args) = parser.parse_args()
    return options

# this is the main function
def runMain():
    # get options and defaults
    print "starting run"
    options = getOptions()
    print "about to run"
    print options.inputFile, options.overlapWigFile, options.outputFile
    run(options.inputFile, options.overlapWigFile, options.outputFile)
    
runMain()
    
