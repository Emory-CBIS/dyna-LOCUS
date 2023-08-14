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
  
  # original
  # color = factor(rep(modname.new, as.vector(count.mode)), levels = modname.new)
  # node_contri = data.frame(index = 1:length(contribution_mode), contribution = contribution_mode, module = color)
  
  # p1 = ggplot(node_contri, aes(x=index, y=contribution, color = module)) + geom_point() + theme_tufte() + geom_rangeframe(col = "black") +
  #   xlab("Node")  + ylab("Node Contribution Index") + theme(text = element_text(size=16, family="serif"))
  
  # plot the mode level contribution
  # mode_contri = aggregate(abs(node_contri$contribution), by = list(node_contri$module), median)
  # names(mode_contri) = c("module", "contribution")
  
  # # boxplot by module
  # p2 = ggplot(node_contri, aes(x = factor(module), y = abs(contribution), fill = color)) + 
  #   geom_boxplot() + geom_jitter(width=0.15) + # geom_jitter(size=0.4, alpha=0.9, position=position_jitter(0)) + 
  #   xlab("") + ylab("Node contribution index") + geom_rangeframe() + 
  #   theme_few() + theme(text = element_text(size=16, family="sans")) + ylim(c(0, 1.25)) # serif
  # # show.legend = F, aes(fill=color)
  # p2
  
  # updated
  module = factor(rep(modname.new, as.vector(count.mode)), levels = modname.new)
  color = factor(rep(modname.color, as.vector(count.mode)), levels = modname.color)
  color = rep(modname.color, as.vector(count.mode))
  node_contri = data.frame(index = 1:length(contribution_mode), 
                           position = position,
                           contribution = contribution_mode, module = module, color = color)

  # ggplot(node_contri, aes(x = module, y = abs(contribution), fill = module)) +
  #   geom_boxplot() + 
  #   scale_fill_manual(values=modname.color) + 
  #   theme_few() + theme(text = element_text(size=16, family="serif"), legend.position="none")
  
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


# repro 1 - 12:   1 - 12
# 24(1), 21(2), 11(3), 23(4), 14(5), 2(6)
# 26(7), 19(8), 5(9), 20(10), 28(11), 17(12)

# 7(13), 18(14), 13(15), 10(16), 30(17), 27(18)
# 16(19), 12(20), 4(21), 15(22), 9(23), 22(24)
# 29(25), 3(26), 25(27), 1(28), 6(29), 8(30)

# ---------------------------
# find where trouble node is
# ---------------------------

# load information
node_info = read.csv("/Users/scarlett/Dropbox/locus\ summary/data/pnc\ data/power264_sorted_node_information.csv")

# trait 2 (6)
sub_node_contri = node_contri[node_contri$module == "Other",]
selected = sub_node_contri[sub_node_contri$contribution > median(sub_node_contri$contribution),]
original_roi = selected[order(selected$contribution, decreasing = T),]$position
# original_roi = sub_node_contri[order(sub_node_contri$contribution, decreasing = T),][1:10,]$position
original_roi # 127 184 248  13 247 250  11  81   5 129 124  31 132  85 126  38
# extract node info
table(node_info[node_info$Original_ROI %in% original_roi, ]$V10)
# Default mode Sensory/somatomotor Hand                Uncertain 
#            5                        3                        8 
node_info_sub = node_info[node_info$Original_ROI %in% original_roi, c("ROINAME2", "ROINAME3", "ROINAME4", "ROINAME5", "ROINAME6", "V10")]
table(node_info_sub[node_info_sub$V10 %in% c("Sensory/somatomotor Hand", "Default mode"),]$ROINAME2)

# trait 19 (8)
sub_node_contri = node_contri[node_contri$module == "Other",]
# selected = sub_node_contri[sub_node_contri$contribution > median(sub_node_contri$contribution),]
# original_roi = selected[order(selected$contribution, decreasing = T),]$position
original_roi = sub_node_contri[order(sub_node_contri$contribution, decreasing = T),][1:3,]$position
original_roi # 220 198 207
# extract node info
table(node_info[node_info$Original_ROI %in% original_roi, ]$V10)
# Default mode             Dorsal attention Fronto-parietal Task Control 
# 1                            1                            4 
# Salience                    Uncertain            Ventral attention 
# 2                            1                            1 

# trait 28 (11)
sub_node_contri = node_contri[node_contri$module == "Other",]
# selected = sub_node_contri[sub_node_contri$contribution > median(sub_node_contri$contribution),]
# original_roi = selected[order(selected$contribution, decreasing = T),]$position
original_roi = sub_node_contri[order(sub_node_contri$contribution, decreasing = T),][1:2,]$position
original_roi # 197 100
# extract node info
table(node_info[node_info$Original_ROI %in% original_roi, ]$V10)
# Default mode Fronto-parietal Task Control 
# 1                            1 

# trait 13 (15)
sub_node_contri = node_contri[node_contri$module == "Other",]
# selected = sub_node_contri[sub_node_contri$contribution > median(sub_node_contri$contribution),]
# original_roi = selected[order(selected$contribution, decreasing = T),]$position
original_roi = sub_node_contri[order(sub_node_contri$contribution, decreasing = T),][1:5,]$position
original_roi # 85 124 126   6   4
# extract node info
table(node_info[node_info$Original_ROI %in% original_roi, ]$V10)
# Default mode    Uncertain 
# 2            3 
node_info_sub = node_info[node_info$Original_ROI %in% original_roi, c("ROINAME2", "ROINAME3", "ROINAME4", "ROINAME5", "ROINAME6", "V10")]
table(node_info_sub[node_info_sub$V10 %in% "Default mode",]$ROINAME2)

# trait 30 (17)
sub_node_contri = node_contri[node_contri$module == "Other",]
# selected = sub_node_contri[sub_node_contri$contribution > median(sub_node_contri$contribution),]
# original_roi = selected[order(selected$contribution, decreasing = T),]$position
original_roi = sub_node_contri[order(sub_node_contri$contribution, decreasing = T),][1:10,]$position
original_roi # 124 248 126 184   6  10  81   4  11 249
# extract node info
table(node_info[node_info$Original_ROI %in% original_roi, ]$V10)
# Default mode    Uncertain 
# 3            7 
node_info_sub = node_info[node_info$Original_ROI %in% original_roi, c("ROINAME2", "ROINAME3", "ROINAME4", "ROINAME5", "ROINAME6", "V10")]
table(node_info_sub[node_info_sub$V10 %in% "Default mode",]$ROINAME2)


# trait 29 (25)
sub_node_contri = node_contri[node_contri$module == "Other",]
# selected = sub_node_contri[sub_node_contri$contribution > median(sub_node_contri$contribution),]
# original_roi = selected[order(selected$contribution, decreasing = T),]$position
original_roi = sub_node_contri[order(sub_node_contri$contribution, decreasing = T),][1:10,]$position
original_roi # 126 124   6  83 248 250  10 184   8  81
# extract node info
table(node_info[node_info$Original_ROI %in% original_roi, ]$V10)
# Default mode    Uncertain 
# 4            6 

# trait 25 (27)
sub_node_contri = node_contri[node_contri$module == "Other",]
# selected = sub_node_contri[sub_node_contri$contribution > median(sub_node_contri$contribution),]
# original_roi = selected[order(selected$contribution, decreasing = T),]$position
original_roi = sub_node_contri[order(sub_node_contri$contribution, decreasing = T),][1:10,]$position
original_roi # 37  35   4 128  83  10 250  31 248 129
# extract node info
table(node_info[node_info$Original_ROI %in% original_roi, ]$V10)
# Default mode Sensory/somatomotor Hand                Uncertain 
# 3                        3                        4 

# --------------------------
# skewness 
# --------------------------

library(e1071)
contribution_index = function(index){
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
  # round to 2 decimal places
  contribution = round(contribution, 2)
  # penalize less than third quantile to 0
  # contribution[contribution <= quantile(contribution, 0.9)] = 0
  
  # return 
  return(contribution)
}

# calculate skewness (node level)
skew = numeric(30)
for(index in 1:30){
  # calculate node contribution of a specific latent source
  contribution = contribution_index(index)
  # re-arrange contribution based on node
  contribution = contribution[order(node.mode)]
  # skewness of the distribution
  skew[index] = round(skewness(contribution), 2)
}
names(skew) = 1:30
sort(skew, decreasing = T)

energy_skew_source = data.frame(x = as.factor(1:30), 
                              energy = log(energy_source), 
                              skewness = skew)

ggplot(energy_skew_source, aes(x=energy, y=skewness, color=x)) + 
  geom_point()  + theme_few(base_family = "serif") + 
  theme(legend.position = "none") + 
  xlab("Energy") + ylab("Skewness") +
  # xlim(-1.8, -1) + ylim(-3.65, -3) +
  # xlim(-5, 0) + ylim(-5, 0)  +
  geom_label_repel(aes(x = energy, y = skewness, label = x), size = 3.5, alpha = 1, data=energy_skew_source) 
# geom_abline(intercept = -2.5, slope = 0.6)

cor(energy_source_log, skew)

# calculate skewness (module level)
skew = numeric(30)
for(index in 1:30){
  # calculate node contribution of a specific latent source
  contribution = contribution_index(index)
  # arrange the nodes by mode
  position = order(node.mode)
  # reorder the contribution using new position
  contribution_mode = contribution[position]
  
  # module level contribution
  modname.new = c("med vis","op vis","lat vis","SM","Aud","DMN","EC","FPL","FPR","CB","Other")
  color = factor(rep(modname.new, as.vector(count.mode)), levels = modname.new)
  node_contri = data.frame(contribution = contribution_mode, module = color)
  mode_contri = aggregate(abs(node_contri$contribution), by = list(node_contri$module), median)
  names(mode_contri) = c("module", "contribution")
  
  # skewness of the distribution
  skew[index] = round(skewness(mode_contri$contribution), 2)
  
  print(mode_contri[order(mode_contri$contribution, decreasing = T)[1:2], ]$module)
}
names(skew) = 1:30
sort(skew, decreasing = T)



