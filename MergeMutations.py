'''
Created on Mar 5, 2015

@author: cmelton
'''

# Imports
from optparse import OptionParser

def GetMutationsFromFile(inputFile, pid, merged):
    f=open(inputFile, 'r')
    line=f.readline()
    while line[0]=="#":
        line=f.readline()
    header=line.strip().split("\t")
    line=f.readline()
    while line!="":
        vals=line.strip().split("\t")
        chrom, pos, covered, bpType, replicationTiming, transcriptType=vals[0:6]
        chrom_pos="chr"+chrom+"_"+pos
        if chrom_pos not in merged:
            merged[chrom_pos]={"pids":set(), "covariates":[]}
        merged[chrom_pos]["pids"].add(pid)
        if len(merged[chrom_pos]["pids"])==1:
            merged[chrom_pos]["covariates"].append([chrom, pos, covered, bpType, replicationTiming, transcriptType])
        line=f.readline()
    f.close()

def getMergedMutations(unparsedInputFiles):
    inputFiles=dict(map(lambda x:(x.split("|")[0], x.split("|")[1]), unparsedInputFiles))
    result={}
    for pid in inputFiles:
        GetMutationsFromFile(inputFiles[pid], pid, result)
    return inputFiles.keys(), result

def writeOutput(pids, mergedMutations, outputFile, splitByChrom):
    outputFiles={}
    chroms=set(map(lambda chrom_pos: chrom_pos.split("_")[0].replace("chr", ""), mergedMutations.keys()))
    if splitByChrom!="T":
        f=open(outputFile, 'w')
        f.write("\t".join(["chrom_pos", "chrom", "pos", "covered", "bpType", "replicationTiming", "transcriptType", "count", "patient_ids"]))
    for chrom in chroms:
        if splitByChrom=="T":
            outputFiles[chrom]=open(outputFile+"_"+chrom+".tsv", 'w')
            outputFiles[chrom].write("\t".join(["chrom_pos", "chrom", "pos", "covered", "bpType", "replicationTiming", "transcriptType", "count", "patient_ids"]))
        else:
            outputFiles[chrom]=f    
    for chrom_pos in mergedMutations:
        chrom=chrom_pos.split("_")[0].replace("chr", "")
        if len(mergedMutations[chrom_pos])>1:
            #chrom_pos.split("_")[0].replace("chr", ""), chrom_pos.split("_")[1]]
            for cov in mergedMutations[chrom_pos]["covariates"]:
                outputFiles[chrom].write("\n"+"\t".join([chrom_pos]+cov))
                outputFiles[chrom].write("\t"+str(len(mergedMutations[chrom_pos]["pids"]))+"\t"+"|".join(list(mergedMutations[chrom_pos]["pids"])))
    if splitByChrom=='T':
        for f in outputFiles.values(): f.close()
    else: f.close()

def getOptions():
    parser = OptionParser()
    parser.add_option("--I", dest = "inputs", help = "",
                      metavar = "STRING", type = "string", default = "01eef340-598c-4205-a990-cec190ac2ca5|./Temp/Chrom22MutCovariates.tsv||2_01eef340-598c-4205-a990-cec190ac2ca5|./Temp/Chrom22MutCovariates.tsv")
    parser.add_option("--O", dest = "output", help = "",
                      metavar = "STRING", type = "string", default = "./Temp/MergedMutations.tsv")
    parser.add_option("--S", dest = "splitByChrom", help = "",
                      metavar = "STRING", type = "string", default = "T")
    (options, args) = parser.parse_args()
    return options

if __name__ == '__main__':
    # get options
    options = getOptions()
    pids, mergedMutations=getMergedMutations(options.inputs.split("||"))
    writeOutput(pids, mergedMutations, options.output, options.splitByChrom)