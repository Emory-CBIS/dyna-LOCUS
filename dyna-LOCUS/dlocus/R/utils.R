#' Helper function to transform the upper triangular of a matrix to a vector
#' 
#' @param X A symmetric matrix
#' @param d If \code{True}, include the diagonal part of the matrix X. Defaults to \code{True}.
Ltrans = function(X, d = T){ 
  X[upper.tri(X, d)]  
} 

#' Helper function to transform a vector to a matrix
#' 
#' @param X A vector representing the upper triangular of a matrix
#' @param V The dimension of the target matrix
#' @param d If \code{True}, the vector X include the diagonal part of the matrix. Defaults to \code{True}.
Ltrinv = function(x, V, d = T){ 
  Y = matrix(0, ncol = V, nrow = V)
  Y[upper.tri(Y, d)] = x
  
  return(Y + t(Y) - d*diag(diag(Y)))  
}

#' Helper function to calculate a perturb inverse of a matrix
#' 
#' @param mat A symmetric matrix 
psdmat_inverse = function(mat){
  # originally defined for psd matrices
  p = dim(mat)[1]
  eigendecomp = eigen(mat)
  
  #if(-min(eigendecomp$values)> 0.0001){print("Error: Matrix has negative eigenvalue!");stop()}
  
  if(min(eigendecomp$values)<=0.0001){
    print("Matrix has nearly zero eigenvalue, perturbation is added.")
    print(round(eigendecomp$values,3))
    perturb = max(max(eigendecomp$values) - p * min(eigendecomp$values), 0)/(p - 1)
  }else{
    perturb = 0
  }
  # mat = mat + diag(p) * perturb
  # eigendecomp = eigen(mat)
  return((eigendecomp$vectors)%*%diag(1/(perturb+eigendecomp$values))%*%t(eigendecomp$vectors))
}

#' Helper function to calculate the scad penalty
#' 
#' @param yvpen A vector upon penalization
#' @param lambda_ch A tunning parameter
#' @param gamm A tunning parameter
SCAD_func = function(yvpen,lambda_ch=0.01,gamma=3){
  if(gamma <= 2){gamma = 2.01; print("Gamma needs > 2!!!")}
  ynew = sign(yvpen)*(abs(yvpen)-lambda_ch)*(abs(yvpen)>=lambda_ch)*(abs(yvpen)<=2*lambda_ch)+ 
    yvpen*(abs(yvpen) > gamma*lambda_ch) + 
    ((gamma-1)*yvpen-sign(yvpen)*gamma*lambda_ch)/(gamma-2)*(abs(yvpen)<=gamma*lambda_ch)*(abs(yvpen)>2*lambda_ch)
  if(sd(ynew) < 0.0000001){print("Parmeters are not correctly specified!");return(ynew)}
  return(ynew)#/sd(ynew)*sd(yvpen))
}

