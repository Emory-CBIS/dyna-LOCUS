#####################
# whitenmap bootstrap
#####################

# load data
load(file = "/projects/guo_lab/neurostat/project/dFC/data/raw_data/slidingwindow15TR.Rdata")

# import environment parameter
args = commandArgs(trailingOnly = TRUE)
iter = as.numeric(args[1])

# set seed
set.seed(iter)
# sample the subjects with replacement
index = sample(1:514, replace = TRUE)
Yraw.subj.new = Yraw.subj[index]

# combine all people's 91*p matrices together for decomposition [dimention: (91*514)*p]
Y = do.call(rbind, Yraw.subj.new)
Y = sweep(Y,2,apply(Y,2,mean),"-") 
N = dim(Y)[1]

# calculate eigen value
eigenA = eigen(Y%*%t(Y), T)

for(q in seq(5, 40, 5)){
  # calculate H (whiten mat)
  H = diag((eigenA$values[1:q]-mean(eigenA$values[(q+1):N]))^(-0.5)) %*% t(eigenA$vectors[,1:q])
  # calculate dewhitening matrix
  H_star = eigenA$vectors[,1:q] %*% diag((eigenA$values[1:q]-mean(eigenA$values[(q+1):N]))^(0.5))
  
  # store the result
  setwd(paste0("/projects/guo_lab/neurostat/project/dFC/data/generated_data/whitenmat_bootstrap/q", q))
  save(H, file = paste0("H_q", q, "_seed", iter,".Rdata"))
  save(H_star, file = paste0("H_star_q", q, "_seed", iter,".Rdata"))
}