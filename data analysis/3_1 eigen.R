##########################
# calculate eigen values
##########################

# library
library(parallel)

# load file
load(file = "/projects/guo_lab/neurostat/project/dFC/data/raw_data/Ycombined15TR.Rdata")

# calculate eigen values
Y = sweep(Y,2,apply(Y,2,mean),"-") 
eigenA = eigen(Y%*%t(Y), T)

# save eigenA
save(eigenA, file = "/projects/guo_lab/neurostat/project/dFC/data/raw_data/eigenA/eigenA_15TR.RData")
