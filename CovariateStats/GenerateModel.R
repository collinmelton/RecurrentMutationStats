library(reshape, verbose = F, warn.conflicts = F)
library(reshape2, verbose = F, warn.conflicts = F)
library(dplyr, verbose = F, warn.conflicts = F)

args<-commandArgs(TRUE)
modeldata<-args[1] #"/Users/cmelton/Documents/AptanaStudio3Workspace/RecurrentMutationWrapper/Temp/ModelData|/Users/cmelton/Documents/AptanaStudio3Workspace/RecurrentMutationWrapper/Temp/ModelData"#args[1] #"/Users/cmelton/Documents/AptanaStudio3Workspace/RecurrentMutationWrapper/Temp/ModelData|/Users/cmelton/Documents/AptanaStudio3Workspace/RecurrentMutationWrapper/Temp/ModelData"# 
output<-args[2] #"/Users/cmelton/Documents/AptanaStudio3Workspace/RecurrentMutationWrapper/Temp/Model"#args[2] # 
data<-rbind_all(lapply(strsplit(modeldata, split = "|", fixed = T)[[1]], function(filename) read.delim(filename, header=T)))
data$pid<-as.character(data$pid)
print(unique(data$pid))
data<-na.omit(data)
fittedmodel<-glm(cbind(mutant, normal) ~ pid+refPair*pid+refPair+replicationTiming+pid*replicationTiming+replicationTiming+transcriptType, family = binomial, data=data, x=F, y=F)
save(fittedmodel, file = output)
