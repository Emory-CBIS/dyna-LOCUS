############################
# decomposition using dlocus
############################

# library dlocus package
library(dlocus)

# load data
load(file = "/projects/guo_lab/neurostat/project/dFC/data/raw_data/Ycombined15TR.Rdata")

# decomposition
q = 30
result = dlocus(Y, q = q, V = 264, n_subject = 514,
                preprocess = T, penalty_S = "SCAD", phi = 2, penalty_A =TRUE, nu = exp(-2),
                maxIteration = 20000, approximation = T,
                espli1 = 0.01, espli2 = 0.01, rho = 0.95, silent = F)
 
# save decomposition results
save(result, file = paste0("/projects/guo_lab/neurostat/project/dFC/result/decomposition/result15TR_q", q, ".RData"))
