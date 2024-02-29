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
  vec_diff = diff(vec, lag = 1)
  standard_d = mean(abs(vec_diff/vec[-106]))
  
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

sd_mat_log = log(sd_mat)

# source specific sd
# average the sd matrix over subjects
sd_source = apply(sd_mat, MARGIN = 2, FUN = mean)
# sd_source = apply(sd_mat_log, MARGIN = 2, FUN = mean)

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
# labels determined by the rank of reproducibility (decreasing) 1:30
labels_repro = rank(-RI)
# corresponding colors
# colors_repro = c(1, 2, 1, rep(3, 2), rep(1, 8), 4, rep(1, 10), 5, rep(1, 4), 6)
colors_repro = c(1, 2, rep(1, 2), 3, rep(1, 7), 1, 4, rep(1, 10), 5, rep(1, 4), 6)
colors_values = c('#E5E7E9', '#619CFF', '#BB8FCE', '#00BA38', '#F8766D', '#F4D03F') # white # 40*40

# minor jitter for better plot
energy_source_jitter = energy_source
energy_source_jitter[12] = energy_source_jitter[12] + 0.02
energy_source_jitter[18] = energy_source_jitter[18] - 0.02
energy_source_jitter[27] = energy_source_jitter[27] + 0.02
# scatter plot (each source one point)
energy_sd_source = data.frame(energy = energy_source_jitter,
                              sd = sd_source,
                              col = colors_repro,
                              labels = labels_repro)

ggplot(energy_sd_source, aes(x = energy, y = sd, color = as.factor(col))) + 
  geom_point(size = 9, shape = 18, alpha = 0.8) + 
  scale_color_manual(values=colors_values) + 
  geom_text(aes(x = energy, y = sd, label = labels), col = "black", size = 3.5, alpha = 1, fontface = "bold", data=energy_sd_source) +
  theme_few() + theme(text = element_text(size=13, family="Times New Roman"), legend.position = "none") + 
  # xlim(c(-1.8, -1)) + 
  xlab("Logorithm of Energy") + ylab("Variation")

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

# # create an age group variable
# data$age = data$age_at_cnb
# data$age_group = factor(NA, levels = c("middle and late childhood ", "adolescence", "early adulthood"))
# data$age_group[data$age >= 8 & data$age <= 11] = "middle and late childhood "
# data$age_group[data$age > 11 & data$age <= 17] = "adolescence"
# data$age_group[data$age >= 18 & data$age <= 21] = "early adulthood"
# table(data$age_group)

# middle and late childhood                 adolescence            early adulthood 
#                        69                        291                        143 

# -----------------------------------------------------------
# linear regression with gender and age group effect (energy)
# -----------------------------------------------------------

# energy ~ gender*age group
gender_agegroup_energy_inter = list()
l = 30
data = data.frame(energy = energy_mat_log[-miss_index,l], gender = gender[-miss_index], age = age[-miss_index])

# create an age group variable
data$age_group = factor(NA, levels = c("middle and late childhood ", "adolescence", "early adulthood"))
data$age_group[data$age >= 8 & data$age <= 11] = "middle and late childhood "
data$age_group[data$age > 11 & data$age <= 17] = "adolescence"
data$age_group[data$age >= 18 & data$age <= 21] = "early adulthood"

gender_agegroup_energy_inter[[l]] = anova(lm(energy ~ gender*age_group, data = data))

# test for contrast
l = 30
# create an age group variable
data = data.frame(energy = energy_mat_log[-miss_index,l], gender = gender[-miss_index], age = age[-miss_index])
data$age_group = factor(NA, levels = c("middle and late childhood ", "adolescence", "early adulthood"))
data$age_group[data$age >= 8 & data$age <= 11] = "middle and late childhood "
data$age_group[data$age > 11 & data$age <= 17] = "adolescence"
data$age_group[data$age >= 18 & data$age <= 21] = "early adulthood"
# fit the model
fit = lm(energy ~ gender*age_group, data = data)
summary(fit)
# create the contrast
library(multcomp)
K = matrix(c(0, 1, 0, 1, 0, 1), 1)
t = glht(fit, linfct = K)
summary(t)
# 30 is significant, while 13 is not

# ---------------------------------------------------------
# age*gender effect on energy/sd for specific latent source 
# ---------------------------------------------------------

age_gender_energy_effect_boxplot = function(ic_index){
  # create a dataset
  data = data.frame(
    age = age[-miss_index],
    gender = as.factor(gender[-miss_index]),
    energy = log(energy_mat[-miss_index, ic_index]),
    sd = log(sd_mat[-miss_index, ic_index])
  )
  
  
  # create an age group variable
  data$age_group = factor(NA, levels = c("childhood ", "adolescence", "early adulthood"))
  data$age_group[data$age >= 8 & data$age <= 11] = "childhood "
  data$age_group[data$age > 11 & data$age <= 17] = "adolescence"
  data$age_group[data$age >= 18 & data$age <= 21] = "early adulthood"
  
  # side by side boxplot
  # grouped boxplot
  p1 = ggplot(data, aes(x=age_group, y=energy, fill=gender)) + 
    geom_boxplot() +
    xlab("") + ylab("Log(Energy)") + 
    scale_x_discrete(guide = guide_axis(angle = -60)) + theme_few() +
    theme(text = element_text(size = 18, family="Times")) +
    ylim(c(-5, 2)) # legend.position = "none"
  

  p1
}

# trait 14
age_gender_energy_effect_boxplot(30)

