# --------------------------
# decomposition using dlocus
# --------------------------

# library
library(dlocus)

# load data
load(file = "/projects/guo_lab/neurostat/project/dFC/data/raw_data/slidingwindow15TR.Rdata")
# load bootstrap samples
names = list.files("/projects/guo_lab/neurostat/project/dFC/data/generated_data/whitenmat_bootstrap", full.names = FALSE)

# import environment parameter
args = commandArgs(trailingOnly = TRUE)
iter = as.numeric(args[1])

# parameter grid
param_grid = expand.grid(
  # q
  q = seq(5, 40, 5),
  # index
  index = 1:100
)

# parameters
q = param_grid[iter,]$q
index = param_grid[iter,]$index

# load bootstrap samples
samples = list.files(paste0("/projects/guo_lab/neurostat/project/dFC/data/generated_data/whitenmat_bootstrap/q", q), full.names = FALSE)
seed = as.numeric(gsub(".Rdata", "", gsub(paste0("H_q", q, "_seed"), "", samples)))
seed = sort(seed[!is.na(seed)])[index]

# set seed
set.seed(seed)
# sample the subjects with replacement
index = sample(1:514, replace = TRUE)
Yraw.subj.new = Yraw.subj[index]

# combine all people's 91*p matrices together for decomposition [dimention: (91*514)*p]
Y = do.call(rbind, Yraw.subj.new)
Y = sweep(Y,2,apply(Y,2,mean),"-") 

# dlocus
result = dlocus(Y, q = q, V = 264, n_subject = 514,
                preprocess = T, penalty_S = "SCAD", phi = 2, penalty_A =TRUE, nu = exp(-2),
                maxIteration = 20000, approximation = T,
                espli1 = 0.01, espli2 = 0.01, rho = 0.95, silent = F)

# record the result
setwd("/projects/guo_lab/neurostat/project/dFC/result/decomposition_bootstrap")
save(result, file = paste0("result_q", q, "_seed", seed,".Rdata"))



