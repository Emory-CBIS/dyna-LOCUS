# library
library(fastICA)
library(dlocus)
library(LOCUS)

# simulation parameters
parameter_grid = expand.grid(
  # type
  type = 1,
  # n_subj
  n_subj = c(20, 50),
  # sd
  sd = seq(0, 3, 0.5),
  # n_sims
  n_sims = 1:100
)

# import environment parameter
args = commandArgs(trailingOnly = TRUE)
iter = as.numeric(args[1])
# get the local parameters
type = parameter_grid$type[iter]
n_subj = parameter_grid$n_subj[iter]
sd = parameter_grid$sd[iter]
n_sims = parameter_grid$n_sims[iter]

# connICA and dlocus decomposition
# load the connectivity matrix
setwd(paste0("/projects/guo_lab/neurostat/project/dFC/simulation/sim_Y/type", type))
load(file = paste0("Y_n", n_subj, "_sd", sd, "_", n_sims, ".RData"))

# hyperparameters for dlocus
phi = 1.5; rho = 0.95; nu = 0.01
# decomposition via locus based on selected parameters
result = dlocus(Y, q = 8, V = 50, n_subject = n_subj,
                preprocess = T, penalty = "SCAD", phi = phi, penalty_A = TRUE, nu = nu,
                maxIteration = 10000, approximation = T,
                espli1 = 0.01, espli2 = 0.01, rho = rho, silent = F)

# decomposition via connICA
Y_demean = sweep(Y, 2, apply(Y, 2, mean), "-")
result_conn = fastICA(Y_demean, n.comp = 8, maxit = 10000)

# decomposition via locus
result_locus = LOCUS(Y, q = 8, V = 50, MaxIteration = 10000, penalty="SCAD", phi = phi,
                     approximation = TRUE, preprocess = TRUE, 
                     espli1 = 0.01, espli2 = 0.01, rho = rho, silent = FALSE)

# save decomposition results
setwd(paste0("/projects/guo_lab/neurostat/project/dFC/simulation/decomp/type", type))
save(result, file = paste0("dlocus_n", n_subj, "_sd", sd, "_", n_sims, ".RData"))
save(result_conn, file = paste0("connICA_n", n_subj, "_sd", sd, "_", n_sims, ".RData"))
save(result_locus, file =paste0("locus_n", n_subj, "_sd", sd, "_", n_sims, ".RData") )