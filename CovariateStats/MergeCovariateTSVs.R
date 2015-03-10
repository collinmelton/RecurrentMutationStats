

# RScript MergeCovariateTSVs.R '/Users/cmelton/Documents/Aptana Studio 3 Workspace/RegStatsCodeForPaper/CovariateStats/TestData/test_rand.tsv' '/Users/cmelton/Documents/Aptana Studio 3 Workspace/RegStatsCodeForPaper/CovariateStats/TestData/test_real.tsv' '/Users/cmelton/Documents/Aptana Studio 3 Workspace/RegStatsCodeForPaper/CovariateStats/TestData/RealRandMutsForOneSample/'

args<-commandArgs(TRUE)
randOutput<-args[1]
realOutput<-args[2]
inputFolder<-args[3]

allFiles<-list.files(path=inputFolder, pattern="*.covariates.tsv", full.names=T)
rands<-allFiles[grep(pattern="random_mutations", x=allFiles)]
reals<-allFiles[grep(pattern="unannotated.varscan.mutect", x=allFiles)]

# read in, merge, and write rand data
randData<-do.call(rbind, lapply(rands, function(x) read.delim(x)))
write.table(randData, file=randOutput, sep="\t", quote=F, row.names=F, col.names=T)

# read in, merge, and write real data
realData<-do.call(rbind, lapply(reals, function(x) read.delim(x)))
write.table(realData, file=realOutput, sep="\t", quote=F, row.names=F, col.names=T)
