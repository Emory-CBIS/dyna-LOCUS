################################################
# calculate whiten matrix and dewhitening matrix
################################################

# load file
load(file = "/projects/guo_lab/neurostat/project/dFC/data/raw_data/Ycombined15TR.Rdata")

# calculate whiten matrix
Y = sweep(Y,2,apply(Y,2,mean),"-") 
N = dim(Y)[1]
qs = seq(5, 40, 5)

# eigenA = eigen(Y%*%t(Y), T)
load("/projects/guo_lab/neurostat/project/dFC/data/generated_data/eigenA/eigenA_15TR.RData")

# import environment parameter
args = commandArgs(trailingOnly = TRUE)
iter = as.numeric(args[1])

# get the value of X
q = qs[iter]

# calculate whitenmat
whitenmat = diag((eigenA$values[1:q] - mean(eigenA$values[(q+1):N]))^(-0.5)) %*% t(eigenA$vectors[,1:q])
# save whitenmat
save(whitenmat, file = paste0("/projects/guo_lab/neurostat/project/dFC/data/generated_data/whitenmat_q/whitenmat15TR_q", q, ".RData"))

# calculate dewhitening mat
# H_star = eigenA$vectors[,1:q] %*% diag((eigenA$values[1:q]-mean(eigenA$values[(q+1):N]))^(0.5))
# save(H_star, file = paste0("/projects/guo_lab/neurostat/project/LOCUS/dynamic/data/generated data/dewhitenmat_q/15TR/dewhitenmat15TR_q", q, ".RData"))
