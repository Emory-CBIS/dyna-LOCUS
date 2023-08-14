###################################
## rho and phi selection
###################################

# libraries
library(ggplot2)

# simulation parameters
parameter_grid = expand.grid(
  # phi
  phi = seq(0.5, 3, 0.5),
  # rho
  rho = seq(0.7, 0.95, 0.05)
)

# merge BIC results
result = matrix(nrow = nrow(parameter_grid), ncol = 9)
for(i in 1:nrow(parameter_grid)){
  bic = get(load(paste0("/Users/scarlett/Dropbox/dFC/analysis/result/bic/raw\ results/bic_phi", 
                        parameter_grid$phi[i], "_rho", parameter_grid$rho[i],".RData")))
  result[i,] = c(parameter_grid$phi[i], parameter_grid$rho[i], bic[-1])
}
colnames(result) = c("phi", "rho", "bic1", "bic2", "bic3", "loglike", "L11", "L12", "L13")

# subset the results
result = as.data.frame(result[,c(1, 2, 3)])
colnames(result) = c("phi", "rho", "bic")

# -------------------------------------------------
# selection of phi: visualize the trend
# -------------------------------------------------
phi_result = result
phi_result$rho = as.factor(phi_result$rho)
phi_result$bic = log(phi_result$bic)
ggplot(phi_result, aes(x=phi, y=bic, group=rho)) + geom_line(aes(color=rho))+
  geom_point(aes(color=rho)) + theme_bw() + 
  xlab(bquote(phi)) + ylab("Log BIC") + scale_color_discrete(bquote(rho)) + 
  theme(text = element_text(size=16, family="serif"))

# -------------------------------------------------
# selection of rho: visualize the trend
# -------------------------------------------------
rho_result = result
rho_result$phi = as.factor(rho_result$phi)
rho_result$bic = log(rho_result$bic)
ggplot(rho_result, aes(x=rho, y=bic, group=phi)) + geom_line(aes(color=phi)) +
  geom_point(aes(color=phi)) + theme_bw() + 
  xlab(bquote(rho)) + ylab("Log BIC") + scale_color_discrete(bquote(phi)) + 
  theme(text = element_text(size=16, family="serif"))

# -------------------------------------------------
# heatmap for visualization
# -------------------------------------------------

ggplot(result, aes(as.factor(rho), as.factor(phi), fill= log(bic))) + 
  geom_tile() + theme_bw() + 
  scale_fill_distiller(palette = "Spectral") +
  xlab(bquote(rho)) + ylab(bquote(phi)) +
  theme(text = element_text(size=16, family="serif")) + 
  guides(fill=guide_legend(title="Log BIC"))
# https://garthtarr.github.io/meatR/ggplot_extensions.html

# range (continuous)

# -------------------------------------------------
# value of Rl when phi = 2 for different rho values
# -------------------------------------------------
setwd("/Users/scarlett/Dropbox/dFC/analysis/result/bic/raw\ results")

# maximum Rl among 30 latent sources
Rl_max = c()
for(rho in seq(0.7, 0.95, 0.05)){
  load(paste0("result_phi2_rho", rho, ".RData"))
  Rl_max = c(Rl_max, max(sapply(result$theta, FUN = function(x){
    length(x$lam_l)
  })))
}
# minimum Rl among 30 latent sources
Rl_min = c()
for(rho in seq(0.7, 0.95, 0.05)){
  load(paste0("result_phi2_rho", rho, ".RData"))
  Rl_min = c(Rl_min, min(sapply(result$theta, FUN = function(x){
    length(x$lam_l)
  })))
}
# display 
mat = rbind(Rl_min, Rl_max)
colnames(mat) = seq(0.7, 0.95, 0.05)
mat

# rho = 0.9
rho = 0.9
load(paste0("result_phi2_rho", rho, ".RData"))
Rl_9 = sapply(result$theta, FUN = function(x){
  length(x$lam_l)
})
# rho = 0.95
rho = 0.95
load(paste0("result_phi2_rho", rho, ".RData"))
Rl_95 = sapply(result$theta, FUN = function(x){
  length(x$lam_l)
})
# display differences
cbind(Rl_9, Rl_95)
