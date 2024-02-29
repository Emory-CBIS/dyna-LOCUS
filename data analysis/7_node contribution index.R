################################
# identification of driven nodes
################################

# library
library(ggplot2)
library(ggthemes)
library(dplyr)
library(gridExtra)
library(patchwork)
# load decomposition result
load("/Users/scarlett/Dropbox/dFC/analysis/result/decomposition/result15TR_q30_adjusted.RData")
# load node mode information
source("/Users/scarlett/Dropbox/dFC/analysis/code/5_4\ latent\ source\ visualization.R")

# --------------------------
# driven module histogram
# --------------------------

driven_nodes = function(index){
  # latent source
  S_index = result$S[index,]
  # extract its Xl, D
  xl = result$theta[[index]]$X_l
  D = diag(result$theta[[index]]$lam_l)
  Rl = nrow(xl)
  
  # node contribution index
  vec = ((diag(Rl) - ifelse(D>0, 1, 0))*1i + ifelse(D>0, 1, 0)) %*% sqrt(abs(D))%*%xl
  # absolute value of node contribution index
  contribution = abs(Re(apply(vec, 2, function(x){sum(x^2)})))
  
  # arrange the nodes by mode
  position = order(node.mode)
  # reorder the contribution using new position
  contribution_mode = contribution[position]
  
  # plot the contribution of each node
  modname.new = c("med vis","op vis","lat vis","SM","Aud","DMN","EC","FPL","FPR","CB","Other")
  modname.color = c("#ffbe00", "#7ac52e", "#00bfd8", "#ff9100", "#c9dd00", "#3a53bc", "#009988", "#aa14b6", "#0099fa", "#fe0061", "#ff2725")
  
  
  module = factor(rep(modname.new, as.vector(count.mode)), levels = modname.new)
  color = factor(rep(modname.color, as.vector(count.mode)), levels = modname.color)
  color = rep(modname.color, as.vector(count.mode))
  node_contri = data.frame(index = 1:length(contribution_mode), 
                           position = position,
                           contribution = contribution_mode, module = module, color = color)

  p2 = node_contri %>% 
    ggplot(aes(x = module, y = abs(contribution), color = module)) + 
    geom_boxplot(width=0.8,lwd=0.5) +
    scale_color_manual(values=modname.color) +
    geom_jitter(width=0) +
    xlab("") + ylab("") +
    theme_few() + theme(text = element_text(size=18, family="Times New Roman"), legend.position="none")
  
  p2
}

# latent source 24
driven_nodes(index = 24) + scale_x_discrete(guide = guide_axis(angle = 60)) + ylim(c(0, 0.4)) + ylab("") + 
  theme(axis.text.y=element_blank(), axis.ticks.y = element_blank())

# latent source 21
driven_nodes(index = 21) + scale_x_discrete(guide = guide_axis(angle = 60)) + ylim(c(0, 0.4)) + ylab("") + 
  theme(axis.text.y=element_blank(), axis.ticks.y = element_blank())
