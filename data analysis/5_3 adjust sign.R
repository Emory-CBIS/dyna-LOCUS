# ---------------------------
# adjust sign
# ---------------------------

# load results
load("/Users/scarlett/Dropbox/dFC/analysis/result/decomposition/result15TR_q30.RData")

# source 5, 11, 19, 20, 22, 25

# extract S and A
S = result$S # 30 34716
A = result$A # 54484 30
theta = result$theta

# change S
S[c(5, 11, 19, 20, 22, 25), ] = -S[c(5, 11, 19, 20, 22, 25), ]
# change A
A[,c(5, 11, 19, 20, 22, 25)] = -A[, c(5, 11, 19, 20, 22, 25)]
# change theta
for(index in c(5, 11, 19, 20, 22, 25)){
  theta_index = theta[[index]]
  theta_index$lam_l = -theta_index$lam_l
  theta_index$M_l = -theta_index$M_l
  theta[[index]] = theta_index
}

# formulate new result
resultnew = list()
resultnew$Conver = result$Conver
resultnew$S = S
resultnew$A = A
resultnew$theta = theta

for(index in c(5, 11, 19, 20, 22, 25)){
  print(cor(resultnew$S[index,], result$S[index,]))
}

# save results
result = resultnew
save(result, file = "/Users/scarlett/Dropbox/dFC/analysis/result/decomposition/result15TR_q30_adjusted.RData")
