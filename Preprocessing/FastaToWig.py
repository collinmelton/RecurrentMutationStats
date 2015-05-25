from optparse import OptionParser

def processWigHeader(line):
    #chrom="_".join(line.split()[2].split(":")[0:3])
    chrom=line.split()[2].split(":")[2]
    start=line.split()[2].split(":")[3]
    stop=line.split()[2].split(":")[4]
    step=line.split()[2].split(":")[5]
    result="fixedStep chrom="+chrom+" start="+start+" step="+step
    return(chrom, result)

def run(input, output, tripletsVersusSinglets=True):
    print "running with input, output:", (input, output)
    f=open(input, "r")
    g=open(output+"0.wig", "w")
    line=f.readline().strip()
    tripletDict=dict([('GGT', 0), ('ACC', 0), ('CAT', 1), ('ATG', 1), ('CTT', 2), ('AAG', 2), ('TTT', 3), ('AAA', 3), ('GAT', 4), ('ATC', 4), ('GTT', 5), ('AAC', 5), ('TAT', 6), ('ATA', 6), ('CCT', 7), ('AGG', 7), ('GAG', 8), ('CTC', 8), ('GCT', 9), ('AGC', 9), ('TGT', 10), ('ACA', 10), ('TCT', 11), ('AGA', 11), ('ATT', 12), ('AAT', 12), ('CTA', 13), ('TAG', 13), ('AGT', 14), ('ACT', 14), ('GTG', 15), ('CAC', 15), ('CGT', 16), ('ACG', 16), ('CAA', 17), ('TTG', 17), ('CCA', 18), ('TGG', 18), ('CGG', 19), ('CCG', 19), ('GGG', 20), ('CCC', 20), ('CGA', 21), ('TCG', 21), ('CTG', 22), ('CAG', 22), ('GCG', 23), ('CGC', 23), ('GGA', 24), ('TCC', 24), ('TTA', 25), ('TAA', 25), ('GTC', 26), ('GAC', 26), ('GAA', 27), ('TTC', 27), ('TGA', 28), ('TCA', 28), ('GCA', 29), ('TGC', 29), ('GTA', 30), ('TAC', 30), ('GGC', 31), ('GCC', 31)])
    first="N"
    second="N"
    third="N"
    while line!="":
        if ">" in line:
            if tripletsVersusSinglets:
                g.write("\n-1")
            g.close()
            chrom, header=processWigHeader(line)
            print "writing file: ", output+chrom+".wig"
            g=open(output+chrom+".wig", "w")
            g.write(header)
            first="N"
            second="N"
            third="N"
            i=0
        else:
            if tripletsVersusSinglets:
                for char in line:
                    first=second
                    second=third
                    third=char
                    triplet=first+second+third
                    if i>0:
                        if triplet not in tripletDict: 
                            g.write("\n-1")#+"\t"+triplet)
                        else:
                            g.write("\n"+str(tripletDict[triplet]))#+"\t"+triplet)
                    i+=1
            else:
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
                      metavar = "FILE", type = "string", default="/Users/cmelton/Documents/GenomeReference/human_g1k_v37.fasta") #human_g1k_v37.fasta")#default = "/srv/gs1/projects/snyder/collinmelton/bundle/2.3/b37/d5/hs37d5.fa")#
    parser.add_option("--O", dest="output", help="path and name of output file",
                      metavar="FILE", type="string", default="./humang1kv37.singlets")#"/srv/gsfs0/projects/snyder/collinmelton/RecurrentMutationDecectionAlgorithm/WholeGenomeWigs/hs37d5.triplets.")#"/srv/gs1/projects/snyder/collinmelton/bundle/2.3/b37/d5/hs37d5.triplets.") #
    parser.add_option("--T", dest="triplets", help="True if tripets",
                      metavar="STRING", type="string", default="F") 
    (options, args) = parser.parse_args()
    return options

# this is the main function
def runMain():
    options = getOptions()
    run(options.input, options.output, options.triplets=="T")    

runMain()