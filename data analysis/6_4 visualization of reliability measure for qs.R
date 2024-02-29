#########################################
# visualization of reliability measure 
# of different number of latent source q
#########################################

# libraries
library(ggplot2)
library(ggthemes)

# load reliability result and construct a dataframe
qs = c()
reliabilitys = c()
for(q in seq(5, 40, 5)){
  # read in reliability result
  load(paste0("/Users/scarlett/Dropbox/dFC/analysis/result/reproducibility/reliability_q", q, ".RData"))
  # reliability
  reliabilitys = c(reliabilitys, reliability)
  qs = c(qs, rep(q, times = q))
}
data = data.frame(reliability = reliabilitys, q = as.factor(qs))

# plot reliability result
# without y labels
ggplot(data, aes(x=q, y=reliability)) + 
  geom_boxplot() +
  xlab("Number of latent source q") +
  ylab("Reproducibility of extracted latent sources") + 
  theme_bw() + theme(text = element_text(size=16, family="serif"), 
                     axis.text.y=element_blank(),
                     legend.position="none")

# with y labels
ggplot(data, aes(x=q, y=reliability)) + 
  geom_boxplot() +
  xlab("Number of latent source q") +
  ylab("Reproducibility of extracted latent sources") + 
  theme_bw() + theme(text = element_text(size=16, family="serif"),
                     legend.position="none")

# bquote("Reliability measure"~gamma[q])
# text family = "Optima"

#########################################
# sources ranked by reproducibility 
#########################################

# q = 30
setwd("/Users/scarlett/Dropbox/dFC/analysis/result/reproducibility")
load("RI_q30.RData")

# sort latent sources and RI (descending order)
df = as.data.frame(cbind(order(RI, decreasing = T), sort(RI, decreasing = T)))
colnames(df) = c("source", "reproducibility")

# visualization based on reproducibility
source("/Users/scarlett/Dropbox/dFC/analysis/code/5_4\ latent\ source\ visualization.R")
load("/Users/scarlett/Dropbox/dFC/analysis/result/decomposition/result15TR_q30_adjusted.RData")

# produce the picture
S = result$S
plots = list()
for(j in 1:30){
  plots[[j]] = levelplotX(S[j,], p = 0.99, li = 2, maint = paste0("Latent source ",j, " (RI = ", round(RI[j], 4), ")"))
}
g = arrangeGrob(plots[[df$source[1]]], plots[[df$source[2]]], plots[[df$source[3]]], plots[[df$source[4]]],plots[[df$source[5]]],
                plots[[df$source[6]]], plots[[df$source[7]]], plots[[df$source[8]]], plots[[df$source[9]]], plots[[df$source[10]]],
                plots[[df$source[11]]], plots[[df$source[12]]], plots[[df$source[13]]], plots[[df$source[14]]], plots[[df$source[15]]],
                plots[[df$source[16]]], plots[[df$source[17]]], plots[[df$source[18]]], plots[[df$source[19]]], plots[[df$source[20]]],
                plots[[df$source[21]]], plots[[df$source[22]]], plots[[df$source[23]]], plots[[df$source[24]]], plots[[df$source[25]]],
                plots[[df$source[26]]], plots[[df$source[27]]], plots[[df$source[28]]], plots[[df$source[29]]], plots[[df$source[30]]],
                nrow = 5, ncol = 6)
# store the plot
setwd("/Users/scarlett/Dropbox/dFC/analysis/result/decomposition") # dynamic
# setwd("/Volumes/LaCie/LOCUS/results/static/latent\ source\ visualization") # static
ggsave(filename = "q30_ri_adjusted.png", g, width = 30, height = 20)
