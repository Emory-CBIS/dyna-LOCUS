# ---------------------
# generate Y on cluster
# ---------------------

# parameter grid
parms = expand.grid(type = 1,
                    n_subj = c(20, 50),
                    sd = seq(0, 3, 0.5),
                    n_sims = 1:100)

# import environment parameter
args = commandArgs(trailingOnly = TRUE)
iter = as.numeric(args[1])

# get parameters
type = parms[iter,]$type
n_subj = parms[iter,]$n_subj
sd = parms[iter,]$sd
n_sims = parms[iter,]$n_sims

# load latent sources
setwd("/projects/guo_lab/neurostat/project/dFC/simulation/sim_data")
load(paste0("S_type", type, ".RData"))
# load loading matrix
setwd("/projects/guo_lab/neurostat/project/dFC/simulation/sim_data")
load(paste0("Asim_smooth_n", n_subj, ".RData"))
# generate datasets
Ytruth = Asim %*% Struth

# set seed
set.seed(n_sims)
# add noise based on sd
Y = Ytruth + rnorm(prod(dim(Ytruth)), mean = 0, sd = sd)
setwd(paste0("/projects/guo_lab/neurostat/project/dFC/simulation/sim_Y/type", type))
save(Y, file = paste0("Y_n", n_subj, "_sd", sd, "_", n_sims, ".RData"))
