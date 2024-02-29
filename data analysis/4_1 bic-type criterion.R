##############################################
## BIC-type criterion calculation
##############################################

# library
library(LOCUS)

# load file
load(file = "/projects/guo_lab/neurostat/project/dFC/data/raw_data/Ycombined15TR.Rdata")

# simulation parameters
parameter_grid = expand.grid(
  # phi
  phi = seq(0.5, 4, 0.5),
  # rho
  rho = seq(0.7, 0.95, 0.05)
)

# import environment parameter
args = commandArgs(trailingOnly = TRUE)
iter=as.numeric(args[1])

# define a norm function
norm_vec = function(x) sum(x^2)
phi = parameter_grid$phi[iter]
rho = parameter_grid$rho[iter]

# use LOCUS for decomposition
result = LOCUS(Y, q = 30, V = 264, MaxIteration = 100000, phi = phi,
               preprocess = TRUE, rho = rho, espli1 = 0.01, espli2 = 0.01, approximation = TRUE)
# save result
setwd("/projects/guo_lab/neurostat/project/dFC/result/bic")
save(result, file = paste0("result_phi", phi, "_rho", rho, ".RData"))

# get S
p = dim(Y)[2]
theta = result$theta
S = result$S

# calculate preprocessed Y
Ynew = sweep(Y, 2, apply(Y, 2, mean), "-")
rm(Y)
load(file = "/projects/guo_lab/neurostat/project/dFC/data/generated_data/whitenmat_q/whitenmat15TR_q30.RData")
Ynew = whitenmat%*%Ynew 
Ynew = Ynew/sd(Ynew)*5

# calculate A
A = Ynew %*% t(S) %*% solve(S %*% t(S))

# calculate BIC
n = dim(Ynew)[1]
p = dim(Ynew)[2]
sigma = sqrt(1/(n*p)*sum(apply(Ynew - A%*%S, MARGIN = 1, norm_vec)))

# calculate loglikelihood
loglike = 0
for (i in 1:n){
  mean = as.vector(A[i,] %*% S)
  loglike = loglike - 2 * sum(log(dnorm(Ynew[i,], mean, sigma)))
}

# calculate logN*sum(||S||0)
L11 = log(n)*sum(abs(S)>1e-1)

# calculate bic
bic1 = loglike + L11

BIC = c(phi, bic1, loglike, L11)

# store the result
setwd("/projects/guo_lab/neurostat/project/dFC/result/bic")
save(BIC, file = paste0("bic_phi", phi, "_rho", rho, ".RData"))

