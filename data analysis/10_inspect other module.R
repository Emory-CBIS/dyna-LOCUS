# ------------------------------------------
# inspect the other module node contribution
# ------------------------------------------

# specify a index

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

# add node system V10
power = read_excel("/Users/scarlett/Dropbox/dFC/analysis/data/pnc\ data/node_system.xls")
power = power$V10[position]

# plot the contribution of each node
modname.new = c("med vis","op vis","lat vis","SM","Aud","DMN","EC","FPL","FPR","CB","Other")
modname.color = c("#ffbe00", "#7ac52e", "#00bfd8", "#ff9100", "#c9dd00", "#3a53bc", "#009988", "#aa14b6", "#0099fa", "#fe0061", "#ff2725")

# updated
module = factor(rep(modname.new, as.vector(count.mode)), levels = modname.new)
color = factor(rep(modname.color, as.vector(count.mode)), levels = modname.color)
color = rep(modname.color, as.vector(count.mode))
node_contri = data.frame(index = 1:length(contribution_mode), 
                         position = position,
                         contribution = contribution_mode, module = module, color = color)
node_contri$power = power

# get other model node contribution
other = node_contri[node_contri$module == "Other",]

# order by contribution
other = other[order(other$contribution, decreasing = T),]
other

# get the modules
table(other$power[1:10])
