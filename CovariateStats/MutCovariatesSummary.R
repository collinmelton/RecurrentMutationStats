library(reshape)
library(reshape2)

args<-commandArgs(TRUE)
input<-args[1] #"/Users/cmelton/Documents/AptanaStudio3Workspace/RecurrentMutationWrapper/Temp/Chrom22MutCovariates.tsv"
output<-args[2] #"/Users/cmelton/Documents/AptanaStudio3Workspace/RecurrentMutationWrapper/Temp/Chrom22MutCovariates.tsv.merged"

data<-read.delim(input, header = T, sep = "\t")
data$count<-1
data<-cast(melt(data, id.vars=c("covered", "bpType", "replicationTiming", "transcriptType"), measure.vars=c("count")), covered+bpType+replicationTiming+transcriptType~variable, fun.aggregate = sum)
write.table(data, file = output, sep = "\t", quote = F, row.names = F, col.names = T)
