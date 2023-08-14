#####################################
# loading analysis: energy, variation
#####################################

# libraries
library(stringr)
library(ggplot2)
library(ggthemes)
library(ggrepel)
library(dplyr)
library(patchwork)
library(DescTools)
library(extrafont)

# load decomposition result
load(file = "/Users/scarlett/Dropbox/dFC/analysis/result/decomposition/result15TR_q30_adjusted.RData")
A = result$A  # (106*514)*30

# Total number of subj
N = 514
# Total number of latent sources
q = 30
# Total length of time
T_len = 106

# split A into 514 parts
matsplitter = function(M, r, c) {
  rg = (row(M)-1)%/%r+1
  cg = (col(M)-1)%/%c+1
  rci = (rg-1)*max(cg) + cg
  N = prod(dim(M))/r/c
  cv = unlist(lapply(1:N, function(x) M[rci==x]))
  dim(cv) = c(r,c,N)
  cv
} 
dat = matsplitter(A, 106, 30) # 106*30*514


# ----------------------------
# energy calculation
# ----------------------------

# input: subj index, source index
# output: the energy of the source of the subj
energy_subj_source = function(subj, source_index){
  begin_index = (subj - 1) * T_len 
  energy = sum(A[(begin_index + 1):(begin_index + T_len), source_index]^2)
  return(energy)
}

# calculate energy matrix
# element i, j represent subj i's j th latent source energy
energy_mat = matrix(nrow = N, ncol = q)
for(i in 1:N){
  for(j in 1:q){
    energy_mat[i, j] = energy_subj_source(i, j)
  }
}

# energy is highly skewed
# use logorithm transformation
energy_mat_log = log(energy_mat)

# source specific energy
# average the energy matrix over subjects
# energy_source = apply(energy_mat, MARGIN = 2, FUN = mean)
energy_source = apply(energy_mat_log, MARGIN = 2, FUN = mean)

# ------------------------------------------------------------------
# sd calculation (using the temporal derivative of the time seires)
# ------------------------------------------------------------------

# input: subj index, source index
# output: the estimated std of the source of the subj
sd_subj_source = function(subj, source_index){
  # the correponding time series
  begin_index = (subj - 1) * T_len 
  vec = A[(begin_index + 1):(begin_index + T_len), source_index]
  # sd of temporal derivate
  standard_d = mean(abs(diff(vec, lag = 1)))
  
  return(standard_d)
}

# calculate standard_d matrix
# element i, j represent subj i's j th latent source standard_d
sd_mat = matrix(nrow = N, ncol = q)
for(i in 1:N){
  for(j in 1:q){
    sd_mat[i, j] = sd_subj_source(i, j)
  }
}

# source specific sd
# average the sd matrix over subjects
sd_source = apply(sd_mat, MARGIN = 2, FUN = mean)

# -----------------------------------------
# subject-wise correlation of energy and sd
# -----------------------------------------

# calculate the correlation of log energy and sd within each subject
cor_energy_sd = c()
for(i in 1:514){
  cor_energy_sd[i] = cor(energy_mat_log[i,], sd_mat[i,])
}

hist(cor_energy_sd, xlab = "correlation", main = "Histogram of energy and sd correlation")
median(cor_energy_sd)

# test if significant larger than 0
# p-value < 2.2e-16
wilcox.test(cor_energy_sd, mu = 0, alternative = "two.sided")

# ----------------------------------------------------
# source-wise relationship between energy and sd
# ----------------------------------------------------

# load reproducibility
load("/Users/scarlett/Dropbox/dFC/analysis/result/reproducibility/RI_q30.RData")
# labels determined by the rank of reproducibility (decreasing)
labels_repro = rank(-RI)
# corresponding colors
colors_repro = c(1, 2, rep(1, 10), 3, rep(1, 10), 4, 5, rep(1, 3), 6, 7)
colors_values = c('#E5E7E9', '#619CFF', '#00BA38', '#F8766D', '#F8766D', '#F8766D', '#F4D03F')
# minor jitter for better plot
energy_source_jitter = energy_source
energy_source_jitter[12] = energy_source_jitter[12] + 0.02
energy_source_jitter[19] = energy_source_jitter[19] + 0.02
energy_source_jitter[27] = energy_source_jitter[27] + 0.02
# scatter plot (each source one point)
# energy_sd_source = data.frame(energy = energy_source_jitter, 
#                               sd = sd_source,
#                               col = colors_repro,
#                               labels = labels_repro)
# 
# ggplot(energy_sd_source, aes(x = energy, y = sd, color = as.factor(col))) + 
#   geom_point(size = 7.5, shape = 18) + 
#   scale_color_manual(values=colors_values) + 
#   geom_text(aes(x = energy, y = sd, label = labels), col = "black", size = 3.2, alpha = 1, data=energy_sd_source) +
#   theme_few() + theme(text = element_text(size=13, family="Times New Roman"), legend.position = "none") + 
#   # xlim(c(-1.8, -1)) + 
#   xlab("Logorithm of Energy") + ylab("Variation")

ggplot(energy_sd_source, aes(x = energy, y = sd, color = as.factor(col))) + 
  geom_point(size = 9, shape = 18, alpha = 0.8) + 
  scale_color_manual(values=colors_values) + 
  geom_text(aes(x = energy, y = sd, label = labels), col = "black", size = 3.5, alpha = 1, fontface = "bold", data=energy_sd_source) +
  theme_few() + theme(text = element_text(size=13, family="Times New Roman"), legend.position = "none") + 
  # xlim(c(-1.8, -1)) + 
  xlab("Logorithm of Energy") + ylab("Variation")

# ----------------------------------------------------
# example of different energy and sd level sources
# ----------------------------------------------------

# plot_loading_ts = function(ic_index, subj_index){
#   # extract the time series
#   ts = dat[, ic_index, subj_index]
#   # create a dataframe
#   ts_dat = data.frame(time = rep(1:106, times = length(ic_index)), 
#                       ts = as.vector(ts), 
#                       ic = rep(ic_index, each = 106))
#   
#   # plot the time series
#   ggplot(data = ts_dat, aes(x = time, y = ts))+
#     geom_line(aes(color = as.factor(ic)), size = 1) +
#     scale_color_manual(values = c("#619CFF", '#00BA38', "#F8766D")) +  # if one wants to add trait 30(17) '#F4D03F'
#     xlab("Time") + ylab("Loadings") +
#     theme_few() + 
#     theme(text = element_text(size=11, family="sans"), legend.position = "none") 
# }
# plot_loading_ts(c(2, 13, 29), 17) + ylim(c(-0.2, 0.3))

# connectivity traits waiting to be plotted
ic_index = c(2, 13, 30, 29)  # 6, 13, 14, 25

# find a representative subject for each index
subj_index = numeric(length(ic_index))
subj_index[1] = order(abs(energy_mat[, ic_index[1]])) # 390 185  86 317 404
subj_index[2] = order(abs(energy_mat[, ic_index[2]])) # 459 224 293  20  74

# 2(6)
# choose a subject with low energy and low variation
match(order(abs(energy_mat[,2]), decreasing = T), order(sd_mat[,2], decreasing = T))
order(sd_mat[,2], decreasing = T)[493] # 390
order(sd_mat[,2], decreasing = T)[452] # 317

# 13(15)
# choose a subject with low energy and median variation
match(order(abs(energy_mat[,13]), decreasing = T), order(sd_mat[,13], decreasing = T))
order(sd_mat[,13], decreasing = T)[84] # 125 *
order(sd_mat[,13], decreasing = T)[159] # 215
order(sd_mat[,13], decreasing = T)[68] # 422
order(sd_mat[,13], decreasing = T)[37] # 86

# 30(17)
# choose a subject with high energy and median variation
match(order(abs(energy_mat[,30]), decreasing = T), order(sd_mat[,30], decreasing = T))
order(sd_mat[,30], decreasing = T)[96] # 179
order(sd_mat[,30], decreasing = T)[175] # 438
order(sd_mat[,30], decreasing = T)[335] # 470

# 29(25)
# choose a subject with high energy and high variation
# the place is the highest energy and its rank in variation
match(order(abs(energy_mat[,29]), decreasing = T), order(sd_mat[,29], decreasing = T))
order(sd_mat[,29], decreasing = T)[2] # 111
order(sd_mat[,29], decreasing = T)[1] # 147
order(sd_mat[,29], decreasing = T)[3] # 38

# extract the time series
ic_index = c(2, 13, 30, 29)
subj_index = c(390, 86, 438, 111)
ts = cbind(dat[,ic_index[1],subj_index[1]], dat[,ic_index[2],subj_index[2]], 
           dat[,ic_index[3],subj_index[3]], dat[,ic_index[4],subj_index[4]])
# create a dataframe
ts_dat = data.frame(time = rep(1:106, times = length(ic_index)), 
                    ts = as.vector(ts), 
                    ic = rep(ic_index, each = 106))
# plot the time series
ggplot(data = ts_dat, aes(x = time, y = ts))+
  geom_line(aes(color = as.factor(ic)), size = 1) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +
  scale_color_manual(values = c("#619CFF", '#00BA38', "#F8766D", '#F4D03F')) +  # if one wants to add trait 30(17) '#F4D03F'
  ylim(c(-0.35, 0.5)) + 
  xlab("Time") + ylab("Loadings") + 
  theme_few() + 
  theme(text = element_text(size=13, family="Times New Roman"), legend.position = "none")

# ----------------------------------------------------
# important gender and sex information
# ----------------------------------------------------

# PNC whole data
PNCinfo = read.table("/Users/scarlett/Dropbox/locus\ summary/data/pnc\ data/PNC_Phenotypes_imaging_subset.txt", 
                     sep="\t",header=T,stringsAsFactors=F) 
PNCinfo = PNCinfo[,c(2,4,757,185)]

# load selected subject list
load(file = "/Users/scarlett/Dropbox/locus\ summary/data/pnc\ data/PNC515.RData")
ID_list = subj.list[-135,]

# PNC selected data
data = PNCinfo[PNCinfo$SUBJID%in%ID_list,]
data = merge(data.frame(SUBJID = ID_list), data, all=T)

# extract gender
gender = as.factor(data$Sex)
# extract birth year
age = data$age_at_cnb

# subjects with missing values
miss_index = which(is.na(gender) | is.na(age))

# ----------------------------------------------------
# linear regression with gender and age group effect
# ----------------------------------------------------

# energy ~ gender*age group
gender_agegroup_energy_inter = list()
for(l in 1:30){
  data = data.frame(energy = energy_mat_log[-miss_index,l], gender = gender[-miss_index], age = age[-miss_index])
  
  # create an age group variable
  data$age_group = factor(NA, levels = c("middle and late childhood ", "adolescence", "early adulthood"))
  data$age_group[data$age >= 8 & data$age <= 11] = "middle and late childhood "
  data$age_group[data$age > 11 & data$age <= 17] = "adolescence"
  data$age_group[data$age >= 18 & data$age <= 21] = "early adulthood"
  
  gender_agegroup_energy_inter[[l]] = anova(lm(energy ~ gender*age_group, data = data))
}

for(i in 1:30){
  print(i)
  print(gender_agegroup_energy_inter[[i]])
}

# significant interaction terms
# 1(0.007), 4(0.004), 13 (0.026), 30 (0.03)


# others consider main effect
# main effect
gender_agegroup_energy_main = list()
for(l in 1:30){
  data = data.frame(energy = energy_mat_log[-miss_index,l], gender = gender[-miss_index], age = age[-miss_index])
  
  # create an age group variable
  data$age_group = factor(NA, levels = c("middle and late childhood ", "adolescence", "early adulthood"))
  data$age_group[data$age >= 8 & data$age <= 11] = "middle and late childhood "
  data$age_group[data$age > 11 & data$age <= 17] = "adolescence"
  data$age_group[data$age >= 18 & data$age <= 21] = "early adulthood"
  
  gender_agegroup_energy_main[[l]] = anova(lm(energy ~ gender+age_group, data = data))
}

for(i in (1:30)[-c(1, 4, 30)]){
  print(i)
  print(gender_agegroup_energy_main[[i]])
}

# gender effect
# 10(0.0004), 13(0.002), 27(0.036)

# age effect
# 2(0.006), 3(0.005), 9(0.02), 14(0.007), 15(0.002), 16(0.023), 27(0.001), 29(0.07)

# sd ~ gender*age group
gender_agegroup_sd_inter = list()
for(l in 1:30){
  data = data.frame(sd = log(sd_mat[-miss_index,l]), gender = gender[-miss_index], age = age[-miss_index])
  
  # create an age group variable
  data$age_group = factor(NA, levels = c("middle and late childhood ", "adolescence", "early adulthood"))
  data$age_group[data$age >= 8 & data$age <= 11] = "middle and late childhood "
  data$age_group[data$age > 11 & data$age <= 17] = "adolescence"
  data$age_group[data$age >= 18 & data$age <= 21] = "early adulthood"
  
  gender_agegroup_sd_inter[[l]] = anova(lm(sd ~ gender*age_group, data = data))
}

for(i in 1:30){
  print(i)
  print(gender_agegroup_sd_inter[[i]])
}

# significant interaction effect
# 4 (0.05), 13 (0.003), 15(0.025), 16(0.008), 17(0.08), 25(0.02), 30(0.0014)

# main effect
gender_agegroup_sd_main = list()
for(l in 1:30){
  data = data.frame(sd = log(sd_mat[-miss_index,l]), gender = gender[-miss_index], age = age[-miss_index])
  
  # create an age group variable
  data$age_group = factor(NA, levels = c("middle and late childhood ", "adolescence", "early adulthood"))
  data$age_group[data$age >= 8 & data$age <= 11] = "middle and late childhood "
  data$age_group[data$age > 11 & data$age <= 17] = "adolescence"
  data$age_group[data$age >= 18 & data$age <= 21] = "early adulthood"
  
  gender_agegroup_sd_main[[l]] = anova(lm(sd ~ gender + age_group, data = data))
}

for(i in (1:30)[-c(4, 13, 15, 16, 25, 30)]){
  print(i)
  print(gender_agegroup_sd_main[[i]])
}


# gender
# 2 (0.098), 10 (1e-6), 22(0.02), 24 (0.023), 26 (0.002), 27 (0.007)
# age
# 1 (0.015), 3 (0.069), 5 (0.066), 7 (0.003), 9 (0.01), 10 (0.03), 18 (0.02), 20 (0.045), 21 (0.095), 29 (0.001)

# ----------------------------------------------------
# gender effect on energy/sd for specific latent source 
# ----------------------------------------------------

gender_effect_plot = function(ic_index){
  # create a dataset
  data = data.frame(
    gender = as.factor(gender[-miss_index]),
    energy = log(energy_mat[-miss_index, ic_index]),
    sd = log(sd_mat[-miss_index, ic_index])
  )
  
  # basic box plot
  p1 = ggplot(data, aes(x=gender, y=energy, fill=gender)) + 
    geom_boxplot() + theme_bw() + 
    xlab("") + ylab("Log(Energy)") + theme_few() +
    theme(legend.position = "none",
          text = element_text(size = 16, family="Times"))
  
  
  p2 = ggplot(data, aes(x=gender, y=sd)) + 
    geom_boxplot() + theme_bw() + 
    xlab("") + ylab("Log(Variation)") + theme_few() +
    theme(legend.position = "none",
          text = element_text(size = 14, family="Times"))
  
  p1
}

# latent source 13
gender_effect_plot(13)

# ----------------------------------------------------
# age effect on energy/sd for specific latent source 
# ----------------------------------------------------

age_effect_plot = function(ic_index){
  # create a dataset
  data = data.frame(
    age = age[-miss_index],
    energy = log(energy_mat[-miss_index, ic_index]),
    sd = log(sd_mat[-miss_index, ic_index])
  )
  
  # create an age group variable
  data$age_group = factor(NA, levels = c("middle childhood", "teens", "adolescence"))
  data$age_group[data$age >= 8 & data$age <= 11] = "middle childhood"
  data$age_group[data$age > 11 & data$age <= 17] = "teens"
  data$age_group[data$age >= 18 & data$age <= 21] = "adolescence"
  
  # side by side boxplot
  # grouped boxplot
  p1 = ggplot(data, aes(x=age_group, y=energy)) + 
    geom_boxplot() + theme_bw() +
    xlab("") + ylab("Log(Energy)") +
    theme(legend.position = "none",
          text = element_text(size = 16, family="serif")) +
    scale_x_discrete(guide = guide_axis(angle = -60))
  
  p2 = ggplot(data, aes(x=age_group, y=sd)) + 
    geom_boxplot() + theme_bw() +
    xlab("") + ylab("Log(Variation)")  + 
    theme(legend.position = "none",
          text = element_text(size = 16, family="serif")) + 
    scale_x_discrete(guide = guide_axis(angle = -60))
  
  p1|p2
} 

# latent source 30
age_effect_plot(30)

# ---------------------------------------------------------
# age*gender effect on energy/sd for specific latent source 
# ---------------------------------------------------------

age_gender_energy_effect_boxplot = function(ic_index){
  # create a dataset
  data = data.frame(
    age = age[-miss_index],
    gender = as.factor(gender[-miss_index]),
    energy = log(energy_mat[-miss_index, ic_index]),
    sd = sd_mat[-miss_index, ic_index]
  )
  
  
  # create an age group variable
  data$age_group = factor(NA, levels = c("middle and late childhood ", "adolescence", "early adulthood"))
  data$age_group[data$age >= 8 & data$age <= 11] = "middle and late childhood "
  data$age_group[data$age > 11 & data$age <= 17] = "adolescence"
  data$age_group[data$age >= 18 & data$age <= 21] = "early adulthood"
  
  # side by side boxplot
  # grouped boxplot
  p1 = ggplot(data, aes(x=age_group, y=energy, fill=gender)) + 
    geom_boxplot() + theme_bw(base_family = "Times") +
    xlab("") + ylab("Log(Energy)") + 
    scale_x_discrete(guide = guide_axis(angle = -60)) + theme_few() +
    theme(legend.position = "none", text = element_text(family="Times")) 
  
  p2 = ggplot(data, aes(x=age_group, y=sd, fill=gender)) + 
    geom_boxplot() + theme_bw(base_family = "Times") +
    xlab("") + ylab("Variation") + theme_few() +
    scale_x_discrete(guide = guide_axis(angle = -60)) +
    theme(text = element_text(family="Times")) 
  
  p1|p2  
}

# trait 15
age_gender_energy_effect_boxplot(13)

# trait 17
age_gender_energy_effect_boxplot(30)

