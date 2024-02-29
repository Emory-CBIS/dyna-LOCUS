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