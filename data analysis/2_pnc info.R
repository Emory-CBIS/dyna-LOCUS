# extract other information about the subjects
load(file = "/Users/scarlett/Dropbox/dFC/analysis/data/pnc\ data/PNC515.RData")
pnc.info = read.delim("/Users/scarlett/Dropbox/dFC/analysis/data/pnc\ data/PNC_Phenotypes_imaging_subset.txt", header = TRUE, sep = "\t", dec = ".")

# select interested columns
pnc.info = pnc.info[,c(2,4,5)]
# select subjects in subj.list 
pnc.info = pnc.info[pnc.info$SUBJID %in% subj.list$V1,] # after selection, 10 people lack info(doesn't matter)
# final data pnc information for use
pnc.info = merge(data.frame(SUBJID = subj.list$V1), pnc.info, all=T)

# store the information
setwd("/Users/scarlett/Dropbox/dFC/analysis/data/generated\ data")
# pnc.info
write.csv(pnc.info, file = "pnc\ info.csv")
