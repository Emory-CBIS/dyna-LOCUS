# ---------------------------------------------------
# merge decomposition result for dictionary learning
# ---------------------------------------------------

# simulation parameters
parms = expand.grid(
  # type
  type = 1,
  # n_subj
  n_subj = c(20, 50),
  # sd
  sd = seq(0, 3, 0.5)
)

# import environment parameter
args = commandArgs(trailingOnly = TRUE)
iter = as.numeric(args[1])

# get parameters
type = parms[iter,]$type
n_subj = parms[iter,]$n_subj
sd = parms[iter,]$sd

# load true S and true A
setwd("/projects/guo_lab/neurostat/project/dFC/simulation/sim_data")
# Struth
load(paste0("S_type", type, ".RData"))
# Asim
load(paste0("Asim_smooth_n", n_subj, ".RData"))

# set working directory
setwd(paste0("/projects/guo_lab/neurostat/project/dFC/simulation/decomp/dl/type", type))

# empty vector to store results
dl_cor_A = c(); dl_cor_S = c()
# merge dl results
for(n_sims in 1:100){
  if(sd %in% 0:3){
    # load dl decompostion results
    dl_S_filename = paste0("S_n", n_subj, "_sd", sd, ".0_", n_sims, ".csv")
    dl_A_filename = paste0("A_n", n_subj, "_sd", sd, ".0_", n_sims, ".csv")
  }else{
    # load dl decompostion results
    dl_S_filename = paste0("S_n", n_subj, "_sd", sd, "_", n_sims, ".csv")
    dl_A_filename = paste0("A_n", n_subj, "_sd", sd, "_", n_sims, ".csv")
  }
  if(file.exists(dl_S_filename)){
    S = read.csv(dl_S_filename, header = F)
    # record average correlation
    dl_cor_S = c(dl_cor_S, mean(apply(abs(cor(t(Struth), t(S))), 1, max)))
  }
  if(file.exists(dl_A_filename)){
    A = read.csv(dl_A_filename, header = F)
    # record average correlation
    dl_cor_A = c(dl_cor_A, mean(apply(abs(cor(Asim, A)), 1, max)))
  }
}
# save result
setwd(paste0("/projects/guo_lab/neurostat/project/dFC/simulation/summary/dl/type", type))
save(dl_cor_A, file = paste0("dl_cor_A_n", n_subj, "_sd", sd, ".RData"))
save(dl_cor_S, file = paste0("dl_cor_S_n", n_subj, "_sd", sd, ".RData"))

