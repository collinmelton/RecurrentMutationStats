import os, SourceWig, OverlapWig
import operator, os, Queue, time, threading, shutil
from optparse import OptionParser

class ThreadSplitWig(threading.Thread):
    def __init__(self, wigFile, outputPrefix, runEvent):
        threading.Thread.__init__(self)
        self.wigFile=wigFile
        self.chrom=""
        self.outputPrefix=outputPrefix
        self.runEvent=runEvent

    # this method runs through the pipeline for each job dictionary in the queue
    def run(self):
        print "start"
        sourceWig=SourceWig.SourceWig(self.wigFile)
        i=0
        next=sourceWig.walk()
        oldc, oldp = "0", 0
        f=open(self.outputPrefix+"."+self.chrom+".wig", "w")
        while next!=None:
            r, c, p = next
            if c!=oldc:
                self.chrom=c
                f.close()
                f=open(self.outputPrefix+"."+self.chrom+".wig", "w")
            if c!=oldc:
                f.write("fixedStep chrom="+c+" start="+str(p)+" step="+str(sourceWig.currentHeader.step))
            elif p!=oldp+sourceWig.currentHeader.step:
                f.write("\nfixedStep chrom="+c+" start="+str(p)+" step="+str(sourceWig.currentHeader.step))
            towrite="\n"+r
            f.write(towrite)
            oldc, oldp = c, p
            i+=1
            if i%100000==0: print i, self.chrom
            next=sourceWig.walk()
        f.close()
        self.runEvent.set()     

# this functions gets the command line options for running the program
def getOptions():
    parser = OptionParser()
    parser.add_option("--I", dest = "inputFile", help = "input file path",
                      metavar = "FILE", type = "string", default = "/Users/cmelton/Documents/Lab/SnyderLab/LocalMutationRate/WigFiles/summary.wig")#0a2coverage.wig.corrected.txt")
    parser.add_option("--O", dest="outputPrefix", help = "output file is a tab delimited file of key and count",
                      metavar = "FILE", type = "string", default = "/Users/cmelton/Documents/Lab/SnyderLab/LocalMutationRate/WigFiles/summary")#0a2coverage.wig.corrected")
    (options, args) = parser.parse_args()
    return options

# this is the main function
def runMain():
    # get options and defaults
    print "starting run"
    options = getOptions()
    # run merge
    runEvent=threading.Event()
    t = ThreadSplitWig(options.inputFile, options.outputPrefix, runEvent)
    t.setDaemon(True)
    t.start()
    time.sleep(1)
    runEvent.wait()
    
    
#     files=["3prime_overlapping_ncrna_exonInfo",
#             "antisense_exonInfo",
#             "IG_C_gene_exonInfo",
#             "IG_D_gene_exonInfo",
#             "IG_J_gene_exonInfo",
#             "IG_V_gene_exonInfo",
#             "lincRNA_exonInfo",
#             "miRNA_exonInfo",
#             "misc_RNA_exonInfo",
#             "Mt_rRNA_exonInfo",
#             "Mt_tRNA_exonInfo",
#             "polymorphic_pseudogene_exonInfo",
#             "processed_transcript_exonInfo",
#             "protein_coding_exonInfo",
#             "pseudogene_exonInfo",
#             "rRNA_exonInfo",
#             "sense_intronic_exonInfo",
#             "sense_overlapping_exonInfo",
#             "snoRNA_exonInfo",
#             "snRNA_exonInfo",
#             "TR_C_gene_exonInfo",
#             "TR_D_gene_exonInfo",
#             "TR_J_gene_exonInfo",
#             "TR_V_gene_exonInfo"]
#     directory="/Users/cmelton/Documents/Lab/SnyderLab/LocalMutationRate/WigFiles/"
#     for f in files:
#         # run merge
#         runEvent=threading.Event()
#         t = ThreadSplitWig(directory+f+".wig", directory+f, runEvent)
#         #t = ThreadSplitWig(options.inputFile, options.outputPrefix, runEvent)
#         t.setDaemon(True)
#         t.start()
#         time.sleep(1)
#         runEvent.wait()
    
runMain()