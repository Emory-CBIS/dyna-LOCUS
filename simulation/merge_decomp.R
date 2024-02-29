# ---------------------------
# merge decomposition result
# ---------------------------

# simulation parameters
parms = expand.grid(
  # type
  type = 1,
  # n_subj
  n_subj = c(20, 50),
  # sd
  sd = seq(0, 3, 0.5),
  # method 
  method = c("dlocus", "connICA", "locus")
)

# import environment parameter
args = commandArgs(trailingOnly = TRUE)
iter = as.numeric(args[1])

# get parameters
type = parms[iter,]$type
n_subj = parms[iter,]$n_subj
sd = parms[iter,]$sd
method = parms[iter,]$method

# load true S and true A
setwd("/projects/guo_lab/neurostat/project/dFC/simulation/sim_data")
# Struth
load(paste0("S_type", type, ".RData"))
# Asim
load(paste0("Asim_smooth_n", n_subj, ".RData"))

# set working directory
setwd(paste0("/projects/guo_lab/neurostat/project/dFC/simulation/decomp/type", type))

if(method == "dlocus"){
  # empty vector to store results
  dlocus_cor_A = c(); dlocus_cor_S = c()
  # merge dlocus results
  for(n_sims in 1:100){
    # load dlocus decompostion results
    dlocus_filename = paste0("dlocus_n", n_subj, "_sd", sd, "_", n_sims, ".RData")
    if(file.exists(dlocus_filename)){
      load(dlocus_filename)
    }
    # extract loading and latent source
    A = result$A; S = result$S;
    # record average correlation
    dlocus_cor_A = c(dlocus_cor_A, mean(apply(abs(cor(Asim, A)), 1, max)))
    dlocus_cor_S = c(dlocus_cor_S, mean(apply(abs(cor(t(Struth), t(S))), 1, max)))
  }
  # save result
  setwd(paste0("/projects/guo_lab/neurostat/project/dFC/simulation/summary/type", type))
  save(dlocus_cor_A, file = paste0("dlocus_cor_A_n", n_subj, "_sd", sd, ".RData"))
  save(dlocus_cor_S, file = paste0("dlocus_cor_S_n", n_subj, "_sd", sd, ".RData"))
}else if(method == "connICA"){
  # empty vector to store results
  connICA_cor_A = c(); connICA_cor_S = c()
  # merge dlocus results
  for(n_sims in 1:100){
    # load dlocus decompostion results
    conn_filename = paste0("connICA_n", n_subj, "_sd", sd, "_", n_sims, ".RData")
    if(file.exists(conn_filename)){
      load(conn_filename)
    }
    # extract loading and latent source
    A = result_conn$S; S = result_conn$A;
    # record average correlation
    connICA_cor_A = c(connICA_cor_A, mean(apply(abs(cor(Asim, A)), 1, max)))
    connICA_cor_S = c(connICA_cor_S, mean(apply(abs(cor(t(Struth), t(S))), 1, max)))
  }
  # save result
  setwd(paste0("/projects/guo_lab/neurostat/project/dFC/simulation/summary/type", type))
  save(connICA_cor_A, file = paste0("connICA_cor_A_n", n_subj, "_sd", sd, ".RData"))
  save(connICA_cor_S, file = paste0("connICA_cor_S_n", n_subj, "_sd", sd, ".RData"))
}else{
  # empty vector to store results
  locus_cor_A = c(); locus_cor_S = c()
  # merge dlocus results
  for(n_sims in 1:100){
    # load dlocus decompostion results
    locus_filename = paste0("locus_n", n_subj, "_sd", sd, "_", n_sims, ".RData")
    if(file.exists(locus_filename)){
      load(locus_filename)
    }
    # extract loading and latent source
    A = result_locus$A; S = result_locus$S;
    # record average correlation 
    locus_cor_A = c(locus_cor_A, mean(apply(abs(cor(Asim, A)), 1, max)))
    locus_cor_S = c(locus_cor_S, mean(apply(abs(cor(t(Struth), t(S))), 1, max)))
  }
  # save result
  setwd(paste0("/projects/guo_lab/neurostat/project/dFC/simulation/summary/type", type))
  save(locus_cor_A, file = paste0("locus_cor_A_n", n_subj, "_sd", sd, ".RData"))
  save(locus_cor_S, file = paste0("locus_cor_S_n", n_subj, "_sd", sd, ".RData"))
}

