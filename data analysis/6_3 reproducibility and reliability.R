##############################################
# reproducibility dynamic 
##############################################

# qs
qs = seq(5, 40, 5)

# import environment parameter
args = commandArgs(trailingOnly = TRUE)
iter = as.numeric(args[1]) 
# get the filename
q = qs[iter]

# load the original results
load(file = paste0("/projects/guo_lab/neurostat/project/dFC/result/decomposition/result15TR_q", q, ".RData"))
# extract the latent sources
S_A = result$S
# seed
names = list.files("/projects/guo_lab/neurostat/project/dFC/result/decomposition_bootstrap", full.names = FALSE)
names = names[grepl(paste0("result_q", q), names, fixed=TRUE)]
seeds = as.numeric(gsub(".Rdata", "", gsub(paste0("result_q", q, "_seed"), "", names)))

# record matched latent source correlation
matching = matrix(0, nrow = q, ncol = length(seeds))
# record index of the matched latent source
indexing = matrix(0, nrow = q, ncol = length(seeds))
# record average correlation across bootstrap samples
averaging = matrix(0, ncol = q, nrow = length(seeds))

for(index in 1:length(seeds)){
  # seed
  seed = seeds[index]
  # import bootstrap result
  setwd("/projects/guo_lab/neurostat/project/dFC/result/decomposition_bootstrap")
  load(file = paste0("result_q", q, "_seed", seed, ".Rdata"))
  
  # extract the S matrix
  S_B = result$S
  
  # calculate their correlation
  dat = data.frame(t(rbind(S_A, S_B)))
  Cor = cor(dat)[1:q, (q+1):(2*q)]
  Cor = abs(Cor)
  
  # sort the correlations
  comb = matrix(0, nrow = q, ncol = 3)
  for(i in 1:q){
    comb[i,] = c(i, which.max(Cor[i,]), max(Cor[i,]))
  }
  
  # matching 
  matching[,index] = comb[,3]
  # index
  indexing[,index] = comb[,2]
  # averaging
  averaging[index,] = apply(Cor, MARGIN = 1, mean)
  

  # save matching, averaging and indexing
  setwd("/projects/guo_lab/neurostat/project/dFC/result/reproducibility")
  save(matching, file = paste0("matching_q", q, ".RData"))
  save(averaging, file = paste0("averaging_q", q, ".RData"))
  save(indexing, file = paste0("indexing_q", q, ".RData"))
}

# calculate the mean
avg.matching  = apply(matching, 1, mean)
avg.averaging = apply(averaging, 2, mean)

# calculate RI
RI = (avg.matching - avg.averaging)/(1-avg.averaging)
setwd("/projects/guo_lab/neurostat/project/dFC/result/reproducibility")
save(RI, file = paste0("RI_q", q, ".RData"))

# calculate reliability
reliability = (avg.matching - avg.averaging)/avg.matching
save(reliability, file = paste0("reliability_q", q, ".RData"))
