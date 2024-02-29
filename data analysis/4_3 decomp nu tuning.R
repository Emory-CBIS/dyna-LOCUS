# --------------------------
# decomposition using dlocus
# --------------------------

# library dlocus package
library(dlocus)

# load file
load(file = "/projects/guo_lab/neurostat/project/dFC/data/raw_data/Ycombined15TR.Rdata")

# simulation parameters
parameter_grid = expand.grid(
  # nu
  log_nu = c(1, seq(-7, 0, by = 0.5))
)

# import environment parameter
args = commandArgs(trailingOnly = TRUE)
iter=as.numeric(args[1])
# get the local parameters
log_nu = parameter_grid$log_nu[iter]
nu = exp(log_nu)

# decomposition
if(log_nu != 1){
  result = dlocus(Y, q = 30, V = 264, n_subject = 514,
                  preprocess = T, penalty_S = "SCAD", phi = 2, penalty_A = TRUE, nu = nu,
                  maxIteration = 20000, approximation = T,
                  espli1 = 0.01, espli2 = 0.01, rho = 0.95, silent = F)
}else{
  result = dlocus(Y, q = 30, V = 264, n_subject = 514,
                  preprocess = T, penalty_S = "SCAD", phi = 2, penalty_A = FALSE, nu = nu,
                  maxIteration = 20000, approximation = T,
                  espli1 = 0.01, espli2 = 0.01, rho = 0.95, silent = F)
}


# save decomposition results
setwd("/projects/guo_lab/neurostat/project/dFC/result/nu_tuning")
save(result, file = paste0("nu_", log_nu, ".RData"))

# calculate the log reconstruction err
Y_demean = sweep(Y, 2, apply(Y, 2, mean), "-") 
A = result$A
S = result$S
err = log(sum((Y_demean - A %*% S)^2))

# save the error
save(err, file = paste0("err_", log_nu, ".RData"))


# ---------------------------------
# merge decomposition results
# ---------------------------------

# load the results
setwd("/projects/guo_lab/neurostat/project/dFC/result/nu_tuning")
# simulation parameters
log_nu = c(1, seq(-7, 0, by = 0.5))

# create a new vector
err_comb = numeric(length(log_nu))

for(i in 1:length(log_nu)){
  # save the error
  load(paste0("err_", log_nu[i], ".RData"))
  
  err_comb[i] = err
}

# check the result
df = cbind(log_nu, err_comb)

