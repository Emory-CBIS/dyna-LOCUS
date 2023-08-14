##############################################
## latent source visualization
##############################################

# libraries
library(ggplot2)
library(ggthemes)
library(gridExtra)
library(lattice)

## visualization settings ------------------------------------------

# node mode(each node's mode, in total 264 nodes)
node.mode = read.table("/Users/scarlett/Dropbox/LOCUS/data/pnc\ data/roi_id.txt")
node.mode = node.mode$V1
modID = data.frame(id = 0:10, modname = c("Other","med vis","op vis","lat vis","DMN","CB","SM","Aud","EC","FPR","FPL"))
# note: UC stands for 'uncertain'

## change order
# 1-3 6 7 4 8 10  9 5  0 orig
# 1-3 4 5 6 7  8  9 10 11 new
modname.new = c("med vis","op vis","lat vis","SM","Aud","DMN","EC","FPL","FPR","CB","Other")
ids = sapply(c(0,6,7,4,8,10,5), match, x = node.mode, nomatch = 0)
for(i in 5:8){
  node.mode[which(ids[,i-3]==1)] = i-1
}
node.mode[which(ids[,1] == 1)] = 11
node.mode[which(ids[,6] == 1)] = 8
node.mode[which(ids[,7] == 1)] = 10
# count the number of nodes <=i (sum up each mode)
k = 0
grid.mode = vector()
count.mode = table(node.mode)
for(i in 1:(length(count.mode)-1)){   
  k = k + count.mode[i]
  grid.mode[i] = k
}
V = 264

loc = (c(grid.mode, V) + c(0, grid.mode))/2

## plot the latent source based on source specific color bar-------------
re.mod = function(Mat){
  Mat = as.matrix(Mat)
  return(Mat[order(node.mode),order(node.mode)])  # determines S's layout
}
Ltrans = function(X, d = T){ X[upper.tri(X,d)]  } 
Ltrinv = function(x, V, d = T){ 
  Y = matrix(0,ncol = V,nrow = V)
  Y[upper.tri(Y,d)]=x
  return(Y + t(Y) -d*diag(diag(Y)))  
}
levelplotX = function(X, maint="", p = 1, colbar=T, li=NULL){
  if(is.null(li)){
    li= quantile(abs(X),probs=p)
  }
  p1 = levelplot(re.mod(Ltrinv(X,V,F)),
                 col.regions = colorRampPalette(c("blue","white", "red"))(1e3),
                 at = do.breaks(c(-li,li), 100), panel = function(...){
                   panel.levelplot(...)
                   panel.abline(v=grid.mode+1, col="black")
                   panel.abline(h=grid.mode+1, col="black")},
                 scales=list(x=list(rot=65,at=loc+1,labels=modname.new,cex=1.05),y=list(at=loc+1,labels=modname.new,cex=1.05)),
                 xlab="",ylab="",main=maint,colorkey=colbar,
                 par.settings=list(axis.text=list(fontfamily="Times"),
                                   par.main.text=list(fontfamily="Times")))
  return(p1)
}


# # --------------------------------------
# # dynamic results (same color bar)
# # --------------------------------------
# 
# load("/Users/scarlett/Dropbox/dFC/analysis/result/decomposition/result15TR_q30_adjusted.RData")
# # produce the picture
# S = result$S
# plots = list()
# for(j in 1:30){
#   plots[[j]] = levelplotX(S[j,], p = 0.99, maint = paste0("Latent Source ",j), li = 2)
# }
# g = arrangeGrob(plots[[1]], plots[[2]], plots[[3]], plots[[4]],plots[[5]],
#                 plots[[6]], plots[[7]], plots[[8]], plots[[9]], plots[[10]],
#                 plots[[11]], plots[[12]], plots[[13]], plots[[14]], plots[[15]],
#                 plots[[16]], plots[[17]], plots[[18]], plots[[19]], plots[[20]],
#                 plots[[21]], plots[[22]], plots[[23]], plots[[24]], plots[[25]],
#                 plots[[26]], plots[[27]], plots[[28]], plots[[29]], plots[[30]],
#                 nrow = 5, ncol = 6)
# # store the plot
# # dynamic
# setwd("/Users/scarlett/Dropbox/dFC/analysis/result/decomposition/")
# ggsave(filename = "q30_15TR_adjusted.png", g, width = 30, height = 20)
