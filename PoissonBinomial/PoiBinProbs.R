
# RScript PoiBinProbs.R '/Users/cmelton/Documents/Aptana Studio 3 Workspace/RegStatsCodeForPaper/PoissonBinomial/TestData/Chrom22MutCovariates_Prob.tsv' '/Users/cmelton/Documents/Aptana Studio 3 Workspace/RegStatsCodeForPaper/PoissonBinomial/TestData/WholeGenomeWindows_1a.1b.1d.1f.2a.2b.2c.3a.3b.4.5_chr22_1.txt' '/Users/cmelton/Documents/Aptana Studio 3 Workspace/RegStatsCodeForPaper/PoissonBinomial/TestData/PoiBinProbs.csv' '1'

library("poibin")

args<-commandArgs(TRUE)
probFile<-args[1] #"~/Documents/Lab/SnyderLab/LocalMutationRate/WindowProbs/RegionProbs_chr22_1.csv"
countsFile<-args[2] #"/Users/cmelton/Documents/Aptana Studio 3 Workspace/GenomeScanner/output/Mar5/OneBP/RawOutput/WholeGenomeWindows_1a.1b.1d.1f.2a.2b.2c.3a.3b.4.5_chr22_1.txt"
outputFile<-args[3]
size<-as.numeric(args[4])
print(args[4])
print(size)

# read in counts and probabilities for a chromosome
probs<-read.csv(probFile, header=T, sep=",", stringsAsFactors=F)
probs$regionID<-sapply(probs$regionID, function(prob) paste(strsplit(x=prob, split="_", fixed=T)[[1]][1],as.character(as.numeric(strsplit(x=prob, split="_", fixed=T)[[1]][2])-size+1), sep="_"))

#get rid of redudndant cols
realNames<-sapply(colnames(probs), function(prob) strsplit(x=prob, split=".[x,y]", fixed=F)[[1]])
probs<-probs[,!duplicated(realNames)]
colnames(probs)<-sapply(colnames(probs), function(prob) strsplit(x=prob, split=".[x,y]", fixed=F)[[1]])
pids<-colnames(probs)[2:length(probs)]
counts<-read.delim(countsFile, header=T, stringsAsFactors=F)
counts$regionID<-paste(counts$chrom, counts$pos, sep="_")

# counts and probs merged
m<-merge(probs, counts[,c("regionID", "number_of_tumors")], by="regionID", all.x=T, all.y=F)
m<-m[!is.na(m$number_of_tumors),]

# function to comput poisson binomial prob k>=K for a row of merged data

getProb<-function(data, i, pids) {
  probs<-as.numeric(data[i,pids])
  cts<-as.numeric(data[i,"number_of_tumors"])
   return(1-ppoibin(cts-1, probs))
}

# compute probabilities
probs<-sapply(1:length(m$regionID), function(i) getProb(m, i, pids))

# generate results data frame
result<-data.frame(regionID=m$regionID, prob=probs, number_of_tumors=m$number_of_tumors)

# write output to file
write.csv(result, file=outputFile, row.names=F, quote=F)
