########################################################
# sliding window calculation and temporal concatentation
########################################################

# libraries 
library(psych)
library(R.matlab)

# import the data
load(file = "/Users/scarlett/Dropbox/dFC/analysis/data/pnc\ data/PNC515.RData")
# load("/Users/scarlett/Dropbox/LOCUS/data/PNC\ data/PNC515.RData")

# use matlab generate the tapered window
# below is the matlab code
% window length
L = 15
% generate rectangular window
rec = ones(L, 1)/L
% generate gaussian window
sigma = 1
alpha = (L-1)/(2*sigma)
w = gausswin(L,alpha)
% convolve two windows
conv_Signal = conv(w, rec, 'same')
% plot the convolved window
n = (-(L-1)/2):((L-1)/2)
plot(n, conv_Signal, type = "l")

# read the conv signal
conv_Signal = unlist(readMat('/Users/scarlett/Dropbox/dFC/analysis/code/1_conv_signal.mat'))

# calculate the correlation matrix 
bin = window.size 
Yraw.subj = list()
for(subj in 1:length(dat)){
  Corrmat = matrix(0, ncol = 120-bin+1, nrow = V*(V-1)/2)
  for(i in 0:(120-bin)){
    data.sd = dat[[subj]]
    Sigma0 = cor(data.sd[(i+1):(i+bin),]*conv_Signal)
    Corrmat[,i+1] = Sigma0[upper.tri(Sigma0)]
  } 
  # fisher z transformation
  corrmat = fisherz(Corrmat)
  # change to Yraw.subj
  Yraw.subj[[subj]] = t(corrmat)
}

# delete 135th subject
Yraw.subj = Yraw.subj[-135]

# combine all people's 106*p matrices together for decomposition [dimention: (106*514)*p]
Y = do.call(rbind, Yraw.subj)

# store the result
setwd("/projects/guo_lab/neurostat/project/dFC/data/raw_data")
save(Yraw.subj, file = paste0("slidingwindow", bin, "TR.Rdata"))
save(Y,file = paste0("Ycombined", bin, "TR.Rdata"))


