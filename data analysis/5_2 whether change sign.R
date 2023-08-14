#########################################
# whether changing sign
#########################################

# load decomposition result
load("/Users/scarlett/Dropbox/dFC/analysis/result/decomposition/result15TR_q30.RData")
# get S
S = result$S # 30 34716

# load averaged dynamic connectivity matrix
load(file = "/Users/scarlett/Dropbox/dFC/analysis/result/fc\ change\ sign/Y_avg15TR.Rdata")
Y_avg_avg = apply(Y_avg, 2, mean)

# transform a symmetric matrix to a vector
# input x: a symmetric matrix
#       V: number of nodes
#       d: the locus vector include diagonal or not
# output: a symmetric matrix
Ltrans = function(X, d = T){ 
  X[upper.tri(X, d)] 
} 

# ---------------------------------
# consider all edges
# ---------------------------------

# transform to adjacency 'matrix'
S_adj = ifelse(S > 0, 1, -1)
Y_avg_adj = ifelse(Y_avg_avg > 0, 1, -1)

# multiply them together and calculate overall agreement
agreement_overall = numeric(30)
for(index in 1:30){
  multiplication = S_adj[index,] * Y_avg_adj
  agreement_overall[index] = mean(multiplication == 1)
}
round(agreement_overall, 2)

# calculate diagnoal agreement
# create a mask
source("/Users/scarlett/Dropbox/dFC/analysis/code/5_4\ latent\ source\ visualization.R")
M = matrix(0, nrow = 264, ncol = 264)
grid_mode = c(0, grid.mode, 264)
for(index in 1:11){
  start = grid_mode[index] + 1
  end = grid_mode[index + 1]
  # make within module elements to 1
  M[start:end, start:end] = 1
}
M_vector = Ltrans(M, d = F)
# agreement
agreement_diagnoal = numeric(30)
for(index in 1:30){
  multiplication = (S_adj[index,] * Y_avg_adj)[M_vector == 1]
  agreement_diagnoal[index] = mean(multiplication == 1)
}
round(agreement_diagnoal, 2)

# summarize the results
agree = cbind(round(agreement_overall, 2), round(agreement_diagnoal, 2))
colnames(agree) = c("overall", "diagnoal")

# static results (subject-wise agreement, average across subjects)
# load each subjectstatic connectivity matrices
load("/Users/scarlett/Dropbox/dFC/analysis/data/generated\ data/Yraw.Rdata")
Y = Yraw[-135,] # 514 34716
n_subj = dim(Y)[1]

# calculate subject level overall and diagnoal agreement
agreement_overall_mat = matrix(nrow = 0, ncol = 30)
agreement_diagnoal_mat = matrix(nrow = 0, ncol = 30)
for(i in 1:n_subj){
  # transform to adjacency matrix
  Y_subject_adj = ifelse(Y[i,] > 0, 1, -1)
  
  # multiply them together and calculate overall agreement
  multiplication = numeric(30)
  for(index in 1:30){
    multiplication[index] = mean((S_adj[index,] * Y_subject_adj) == 1)
  }
  agreement_overall_mat = rbind(agreement_overall_mat, multiplication)
  
  # agreement diagnoal
  multiplication_diag = numeric(30)
  for(index in 1:30){
    temp = (S_adj[index,] * Y_subject_adj)[M_vector == 1]
    multiplication_diag[index] = mean(temp == 1)
  }
  agreement_diagnoal_mat = rbind(agreement_diagnoal_mat, multiplication_diag)
}

# average across subjects
agree = data.frame(agree)
agree$agree_subj_avg_overall = round(apply(agreement_overall_mat, 2, mean), 2)
agree$agree_subj_avg_diagnoal = round(apply(agreement_diagnoal_mat, 2, mean), 2)

setwd("/Users/scarlett/Dropbox/dFC/analysis/result/fc\ change\ sign")
write.csv(agree, file = "agreement_alledges.csv")

# ---------------------------------
# only consider significant edges
# ---------------------------------

# penalize 0.9 percent edges to 0
S_adj = S
S_adj[abs(S) < quantile(abs(S), probs = 0.98)] = 0
# transform to adjacency 'matrix'
S_adj = ifelse(S_adj > 0, 1, ifelse(S_adj < 0, -1, 0))

# multiply them together and calculate overall agreement
agreement_overall = numeric(30)
for(index in 1:30){
  multiplication = S_adj[index,] * Y_avg_adj
  nonzero_index = which(S_adj[index,] != 0)
  agreement_overall[index] = mean(multiplication[nonzero_index] == 1)
}
round(agreement_overall, 2)
  
# agreement diagnoal
agreement_diagnoal = numeric(30)
for(index in 1:30){
  multiplication = S_adj[index,] * Y_avg_adj
  nonzero_M_index = which(S_adj[index,] != 0 & M_vector == 1)
  agreement_diagnoal[index] = mean(multiplication[nonzero_M_index] == 1)
}
round(agreement_diagnoal, 2)

# summarize the results
agree = cbind(round(agreement_overall, 2), round(agreement_diagnoal, 2))
colnames(agree) = c("overall", "diagnoal")

# calculate subject level overall and diagnoal agreement
agreement_overall_mat = matrix(nrow = 0, ncol = 30)
agreement_diagnoal_mat = matrix(nrow = 0, ncol = 30)
for(i in 1:n_subj){
  # transform to adjacency matrix
  Y_subject_adj = ifelse(Y[i,] > 0, 1, -1)
  
  # multiply them together and calculate overall agreement
  multiplication = numeric(30)
  for(index in 1:30){
    multi = S_adj[index,] * Y_subject_adj
    nonzero_index = which(S_adj[index,] != 0)
    multiplication[index] = mean(multi[nonzero_index] == 1)
  }
  agreement_overall_mat = rbind(agreement_overall_mat, multiplication)
  
  # agreement diagnoal
  multiplication_diag = numeric(30)
  for(index in 1:30){
    multi = S_adj[index,] * Y_subject_adj
    nonzero_M_index = which(S_adj[index,] != 0 & M_vector == 1)
    multiplication_diag[index] = mean(multi[nonzero_M_index] == 1)
  }
  agreement_diagnoal_mat = rbind(agreement_diagnoal_mat, multiplication_diag)
}

# average across subjects
agree = data.frame(agree)
agree$agree_subj_avg_overall = round(apply(agreement_overall_mat, 2, mean), 2)
agree$agree_subj_avg_diagnoal = round(apply(agreement_diagnoal_mat, 2, mean), 2)

setwd("/Users/scarlett/Dropbox/dFC/analysis/result/fc\ change\ sign")
write.csv(agree, file = "agreement_sigedges0.98.csv")
  
# ---------------------------------
# decide whether to change sign
# ---------------------------------

setwd("/Users/scarlett/Dropbox/dFC/analysis/result/fc\ change\ sign")
agreement_alledges = read.csv("agreement_sigedges0.98.csv")
# the bigger the better
ranks = sapply(agreement_alledges, rank, ties.method = "first")
# average rank (the bigger the better)
avg_rank = apply(ranks[,-1], 1, mean)
df = data.frame(cbind(1:30, avg_rank))
colnames(df) = c("index", "avg_rank")
df = df[order(df$avg_rank, decreasing = T),]
write.csv(df, file = "agreement_sigedges0.98_rank.csv")









