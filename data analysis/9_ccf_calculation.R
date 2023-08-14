###################################
# synchrony between latent sources
###################################

# libraries
library(plot.matrix)
library(ggplot2)
library(patchwork)

# load decomposition result
load("/Users/scarlett/Dropbox/dFC/analysis/result/decomposition/result15TR_q30_adjusted.RData")
load("/Users/scarlett/Dropbox/dFC/analysis/result/reproducibility/RI_q30.RData")

# get S
S = result$S # 30 34716
# get loading matrix
A = result$A # 54484  30

# sort latent sources and RI (descending order)
df = as.data.frame(cbind(order(RI, decreasing = T), sort(RI, decreasing = T)))
colnames(df) = c("source", "reproducibility")

# permutate S and A
S = S[df$source, ]
A = A[, df$source]

# split A into 514 parts
matsplitter=function(M, r, c) {
  rg = (row(M)-1)%/%r+1
  cg = (col(M)-1)%/%c+1
  rci = (rg-1)*max(cg) + cg
  N = prod(dim(M))/r/c
  cv = unlist(lapply(1:N, function(x) M[rci==x]))
  dim(cv)=c(r,c,N)
  cv
} 
dat = matsplitter(A, 106, 30) # 106 30 514

# ----------------------------------------------------------------------------------------------
# each subject, calculate the lag which achieves the largest abs ccf between latent source pairs
# ----------------------------------------------------------------------------------------------

ccf_lag_initial = array(NA, c(30, 30, 514))

for(k in 1:514){
  ccf_mat = matrix(NA, 30, 30)
  for(i in 1:29){
    for(j in (i+1):30){
      # calculate ccfs and lags
      fit = ccf(dat[,i,k], dat[,j,k], ylab = "Cross-correlation", lag.max = 40)
      mat = data.frame(acf = fit$acf, lag = fit$lag)
      
      # identify the lag which achieves the largest abs ccf value
      index = which.max(abs(mat$acf))
      # record the lag
      ccf_mat[i,j] = mat[index,]$lag
    }
  }
  ccf_lag_initial[,,k] = ccf_mat
}

# save results
setwd("/projects/guo_lab/neurostat/project/dFC/result/synchrony")
save(ccf_lag_initial, file = "ccf_lag_initial.RData")

# -----------------------------------------------------------------------------------------------------------------------------
# each subject, recalculate the lag which achieves the largest abs ccf lag between latent source pairs using determined sign
# -----------------------------------------------------------------------------------------------------------------------------

# use this info to determine the sign of ccf each latent source 
# -1: ccf smaller than 0; 1: ccf larger than 0
lag_sign = matrix(NA, nrow = 30, ncol = 30)
for(i in 1:29){
  for(j in (i+1):30){
    # calculate ccfs and lags
    prop = c(mean(ccf_lag_initial[i,j,]<0), mean(ccf_lag_initial[i,j,]==0), mean(ccf_lag_initial[i,j,]>0))
    index = which.max(prop)
    
    if(index == 1){
      lag_sign[i, j] = -1
    }else if(index == 2){
      lag_sign[i, j] = 0
    }else{
      lag_sign[i, j] = 1
    }
  }
}

# recalculate the max value lag between latent source pairs
ccf_lag_intermediate = array(NA, c(30, 30, 514))

for(k in 1:514){
  ccf_mat = matrix(0, 30, 30)
  for(i in 1:29){
    for(j in (i+1):30){
      # calculate ccfs and lags
      fit = ccf(dat[,i,k], dat[,j,k], ylab = "Cross-correlation", lag.max = 40)
      mat = data.frame(acf = fit$acf, lag = fit$lag)
      
      # identify the lag which achieves the largest ccf value (with sign information)
      if(lag_sign[i,j]==-1){
        index = which.max(abs(mat[mat$lag<0,]$acf)) 
      }else if(lag_sign[i,j]==1){
        index = which.max(abs(mat[mat$lag>0,]$acf)) + 41
      }else{
        index = 41
      }
      
      # record the lag
      ccf_mat[i,j] = mat[index,]$lag
    }
  }
  ccf_lag_intermediate[,,k] = ccf_mat
}

# save results
setwd("/projects/guo_lab/neurostat/project/dFC/result/synchrony")
save(ccf_lag_intermediate, file = "ccf_lag_intermediate.RData")

# ------------------------------------------------------
# determine the mode of lag for each latent source pair
# ------------------------------------------------------

# get mode function
getmode = function(v) {
  uniqv = unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

lag_mode = matrix(NA, nrow = 30, ncol = 30)
for(i in 1:29){
  for(j in (i+1):30){
    lag_mode[i,j] = getmode(ccf_lag_intermediate[i,j,])
  }
}

# basically all lag mode are either -1 or 1, only 1 pair with lag 0
# lag mode = -1, proportion 56.09%
# lag mode = 1, proportion 43.68% 

# ------------------------------------------------------------------------
# calculate the ultimate ccf value of each latent source pair per subject
# ------------------------------------------------------------------------

ccf_value_final = array(NA, c(30, 30, 514))

for(k in 1:514){
  ccf_mat = matrix(NA, 30, 30)
  for(i in 1:29){
    for(j in (i+1):30){
      # calculate ccfs and lags
      fit = ccf(dat[,i,k], dat[,j,k], ylab = "Cross-correlation", lag.max = 1)
      mat = data.frame(acf = fit$acf, lag = fit$lag)
      # record the value
      ccf_mat[i,j] = mat[mat$lag == lag_mode[i,j],]$acf
    }
  }
  ccf_value_final[,,k] = ccf_mat
}


for(k in 1:514){
  ccf_mat = matrix(NA, 30, 30)
  for(i in 1:29){
    for(j in (i+1):30){
      # calculate ccfs and lags
      fit = ccf(dat[,i,k], dat[,j,k], ylab = "Cross-correlation", lag.max = 1)
      mat = data.frame(acf = fit$acf, lag = fit$lag)
      # record the value
      ccf_mat[i,j] = mat[mat$lag == lag_mode[i,j],]$acf
    }
  }
  ccf_value_final[,,k] = ccf_mat
}

save(ccf_value_final, file = "ccf_value_final.RData")

# ----------------------------------
# visualization of ccf distributions
# ----------------------------------

# transform the matrix to be a vector
Ltrans = function(X, d=T){ 
  X[upper.tri(X,d)]  
} 

ccf_value_mat = matrix(nrow = 514, ncol = choose(30,2))
for(i in 1:514){
  ccf_mat = ccf_value_final[,,i]
  ccf_value_mat[i,] = ccf_mat[upper.tri(ccf_mat,d = F)]
}

# median ccf of latent source pairs across subjects
ccf_median_map = apply(ccf_value_final, MARGIN = c(1,2), median)
ccf_median_map[lower.tri(ccf_median_map)] = 0
diag(ccf_median_map) = 0
# visualization
par(mar=c(5.1, 4.1, 4.1, 4.1))
value.max = max(abs(ccf_median_map))
value.cut = sort(abs(ccf_median_map[upper.tri(ccf_median_map,d = F)]), decreasing = T)[6]
cutoff = c(-value.max, -value.cut, -0.0001, 0.0001, value.cut, value.max)
plot(ccf_median_map, breaks = cutoff, col=c("#F28ED4", "#CFF28E", "white", "#F2EC8E","#F28ED4"), 
     main = "", axis.col=list(side=1, cex.axis=0.7), axis.row=list(cex.axis=0.7), 
     border = NA, xlab = "Latent Source Index", ylab = "Latent Source Index")
# find which part of the matrix is larger than 0.4
which(abs(ccf_median_map) >= 0.40, arr.ind = TRUE)
# row col           lag
# 13  15 [-0.546]   1
# find which part of the matrix is larger than 0.35
which(abs(ccf_median_map) > 0.35, arr.ind = TRUE)
#  row col           lag
#   2   3  [-0.391]  -1
#   6  14  [0.40]    -1
#   13  15 [-0.546]   0
#   14  22 [0.382]   -1
#   20  22 [0.368]   -1
#   22  29 [0.359]   1 
#   20  30 [-0.359]  1

# ------------------------------------------------
# visualization of synchronization (subject level)
# ------------------------------------------------

# visualization of synchronization (subject level)
subj_loading = function(ic1, ic2){
  # determine which subject to display
  vec = ccf_value_final[ic1, ic2, ]
  
  if(median(vec) < 0){
    subj_id = which(vec == min(vec))
  }else{
    subj_id = which(vec == max(vec))
  }

  
  # create a dataframe for plot depends on which ic to display first
  # if(lag_mode[ic1, ic2] == 1){
  #   data = data.frame(x = rep(1:106, 2), y = c(dat[, ic1, subj_id], dat[, ic2, subj_id]), 
  #                     Trait = rep(c(paste0("x: ", ic1), paste0("y: ", ic2)), each = 106))
  # }else{
  #   data = data.frame(x = rep(1:106, 2), y = c(dat[, ic1, subj_id], dat[, ic2, subj_id]), 
  #                     Trait = factor(rep(c(paste0("x: ", ic1), paste0("y: ", ic2)), each = 106),
  #                                    level = c(paste0("y: ", ic2), paste0("x: ", ic1))))
  # }
  
  data = data.frame(x = rep(1:106, 2), y = c(dat[, ic1, subj_id], dat[, ic2, subj_id]), 
                    Trait = rep(c(paste0("x: ", ic1), paste0("y: ", ic2)), each = 106))
  
  p = ggplot(aes(x = x, y = y, group = Trait, col = Trait), data = data) + 
    geom_line(alpha = 0.8)+ theme_few() + 
    ggtitle(bquote(CCF[.(lag_mode[ic1, ic2])](x,y) == .(round(ccf_median_map[ic1, ic2], 2)))) + 
    xlab("Time window") + ylab("Loading") + ylim(c(-0.35, 0.35)) + 
    theme(text = element_text(size=12, family="Times"),
          plot.title = element_text(size = 12),
          legend.position = "none")
  
  return(p)
}

p1 = subj_loading(13, 15) 
p2 = subj_loading(6, 14) 
p3 = subj_loading(2, 3) 
p4 = subj_loading(14, 22) 
p5 = subj_loading(20, 22)
p6 = subj_loading(20, 30) 

(p1|p2)/
  (p3|p4)/
  (p5|p6)

# theme(axis.title.x=element_blank(),
#       axis.text.x=element_blank(),
#       axis.ticks.x=element_blank())

#  row col           lag
#   2   3  [-0.391]  -1
#   6  14  [0.40]    -1
#   13  15 [0.546]   0
#   14  22 [0.382]   -1
#   20  22 [0.368]   -1
#   22  29 [0.359]   1 
#   20  30 [-0.359]  1


# --------------------------------------------------------------
# visualization of synchronization (subject level, powerpoint)
# --------------------------------------------------------------

# visualization of synchronization (subject level)
subj_loading = function(ic1, ic2, subj_id = NULL){
  # determine which subject to display
  vec = ccf_value_final[ic1, ic2, ]
  
  if(is.null(subj_id)){
    if(median(vec) < 0){
      subj_id = which(vec == min(vec))
    }else{
      subj_id = which(vec == max(vec))
    }
  }
  
  data = data.frame(x = rep(1:106, 2), y = c(dat[, ic1, subj_id], dat[, ic2, subj_id]), 
                    Trait = rep(c(paste0("x: ", ic1), paste0("y: ", ic2)), each = 106))
  
  p = ggplot(aes(x = x, y = y, group = Trait, col = Trait), data = data) + 
    geom_line(alpha = 0.8)+ theme_few() + 
    ggtitle(bquote(CCF(x,y)[.(lag_mode[ic1, ic2])] == .(round(ccf_median_map[ic1, ic2], 2)))) + 
    xlab("Time window") + ylab("Loading") + ylim(c(-0.35, 0.35)) + 
    theme(text = element_text(size=12, family="Times"),
          plot.title = element_text(size = 12),
          legend.position = "none")
  
  return(p)
}

p1 = subj_loading(13, 15) 
p2 = subj_loading(6, 14, 296) 
p3 = subj_loading(2, 3) 
p4 = subj_loading(14, 22) 
p5 = subj_loading(20, 22)
p6 = subj_loading(20, 30) 

(p1|p2)/
  (p3|p4)/
  (p5|p6)

# theme(axis.title.x=element_blank(),
#       axis.text.x=element_blank(),
#       axis.ticks.x=element_blank())

#  row col           lag
#   2   3  [-0.391]  -1
#   6  14  [0.40]    -1
#   13  15 [0.546]   0
#   14  22 [0.382]   -1
#   20  22 [0.368]   -1
#   22  29 [0.359]   1 
#   20  30 [-0.359]  1
