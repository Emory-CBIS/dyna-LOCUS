#' Node-rotation algorithm (an extremely efficient approximation version)
#' An extremely efficient approximation version with potentially higher performance
#' 
#' @param Y The (preprocessed) concatenated group-level connectivity data.
#' @param A A loading matrix from previous iteration.
#' @param theta A list containing ic-specific parameters and channels.
#' @param q The number of latent sources one would like to extract.
#' @param K The number of edges.
#' @param V The number of nodes.
#' @param n_subject The number of subjects in the study.
#' @param t_length The number of sliding windows.
#' @param penalt_S Penalty added on latent sources. Defaults to \code{TRUE}.
#' @param lambda_ch A tunning parameters with respect to the element-wise penalty on S.
#' @param gamma A tunning parameters with respect to the element-wise penalty on S.
#' @param penalt_A If \code{TRUE}, temporal smooth lasso penalty is added on loadings. Defaults to \code{FALSE}.
#' @param nu A tunning parameter with respect to the temporal smooth lasso penalty on A.
#' @param H_star dewhitening matrix. 
#' @param preprocess If \code{TRUE}, the concatenated group-level connectivity data is preprocessed.
#' @param imput_method A character describing how to input the diagonal values of latent sources before eigen decomposition. 
#' Defaults to \code{"Previous"}.
#' @param silent If \code{TRUE}, print out the penalty added on A and S. Defaults to \code{FALSE}.
#' 
#' @return Updated values of loading matrix \code{A} and latent source \code{S}.
dlocus_update_approx = function(Y, A, theta, 
                                q, K, V, n_subject, t_length, 
                                penalt_S = TRUE, lambda_ch, gamma, 
                                penalt_A = FALSE, nu, H_star,
                                preprocess, imput_method = "Previous", silent = FALSE){
  
  # whether add penalty on S (latent sources)
  if(is.null(penalt_S)){
    if(!silent) print(paste("Locus without penalty on latent sources"))
  }else{
    if(!silent) print(paste("Locus with", penalt_S, "penalty on latent sources"))
  }
  
  # whether add smooth penalty on A (loading matrix)
  if(!penalt_A){
    if(!silent) print(paste("Locus without penalty on loading matrix"))
  }else{
    if(t_length == 1){
      if(!silent) print(paste("Locus without penalty on loading matrix"))
    }else{
      if(!silent) print(paste("Locus with smooth penalty on loading matrix"))
    }
  }
  
  # use theta to store ic-specific parameters and channels
  theta_new = list()

  # use Rls to store the value of Rl in each latent sources
  Rls = numeric(q)
  for(curr_ic in 1:q){
    # current latent source related theta
    theta_ic = theta[[curr_ic]]
    # dimension of Rl in current latent source
    Rls[curr_ic] = dim(theta_ic$X_l)[1]
    # lth row of A_tilde transpose (1*q) = lth column of A_tilde (1*q)
    # times Y (q*p)
    # result in a 1*p vector -> transpose -> p*1 vector
    # the same as Yic in Locus_update
    S_lold = t(theta_ic$M_l%*%Y)
    
    # update V nodes in S simultaneously
    if(is.null(penalt_S)){
      S_new_0 = S_lold
    }else if(penalt_S == "SCAD"){
      if(gamma<=2){ print("Gamma needs to be > 2!"); gamma = 2.01}
  
      S_new_0= SCAD_func(S_lold, lambda_ch = lambda_ch, gamma = gamma)
      S_new_0 = S_new_0 /sd(S_new_0)*sd(S_lold)
      
    }else if(penalt_S == "Hardthreshold"){
      S_new_0 = S_lold*(abs(S_lold)>=lambda_ch)
    }else if(penalt_S == "L1"){
      S_new_0 = sign(S_lold)*(abs(S_lold)-lambda_ch)*(abs(S_lold)>=lambda_ch)
      S_new_0 = S_new_0 /sd(S_new_0)*sd(S_lold)
    }else{
      print("No Penalty available!")
      stop()
    }
    
    if(imput_method == "Previous"){
      Sl = Ltrinv(S_new_0, V, F) + diag(diag(t(theta_ic$X_l)%*%diag(theta_ic$lam_l)%*%theta_ic$X_l ))
    }else if(imput_method == "Average"){ 
      Sl = Ltrinv(S_new_0, V, F) + diag(rep(mean(S_new_0), V)) 
    }else{print("No Imputation available!")
      stop()
    }
    
    # reconstruct Sl based on the number of eigen values
    eigenSl = eigen(Sl)
    orderEigen = order(abs(eigenSl$values),decreasing = T)
    Rl = Rls[curr_ic]
    eigenset = orderEigen[1:Rl]
    
    theta_new[[curr_ic]] = list()
    theta_new[[curr_ic]]$lam_l = eigenSl$values[eigenset]
    if( theta_new[[curr_ic]]$lam_l[1]<0 ){
      theta_new[[curr_ic]]$lam_l = -1*theta_new[[curr_ic]]$lam_l
    }
    for(j in 1:Rl){
      theta_ic$X_l[j,]= eigenSl$vectors[,eigenset[j]]
    }
    theta_new[[curr_ic]]$X_l = theta_ic$X_l
  }
  # reconstruct S based on X_l and lam_l values
  newS = array(dim = c(q, K))
  for(l in 1:q){
    newS[l,] = Ltrans(t(theta_new[[l]]$X_l)%*% diag(theta_new[[l]]$lam_l) %*%theta_new[[l]]$X_l, F)
  }
  
  # preparations to update A
  if(t_length > 1){
    # construct R
    R = matrix(0, nrow = t_length-1, ncol = t_length)
    for(t in 1:(t_length-1)){
      R[t,t] = 1
      R[t, t+1] = -1
    }
    # construct R_star
    R_star = kronecker(diag(1, n_subject), R)
    if(preprocess){
      # construct W = R_star*H_star
      W = R_star %*% H_star
      # calculate eigen values and vectors of W'W
      eigenW = eigen(t(W) %*% W, T)
      # eigen vectors
      Q2 = eigenW$vectors
      # eigen values
      Lmbd2 = eigenW$values
    }else{
      # calculate eigen values and vectors of R_star'R_star
      eigenR_star = eigen(t(R_star) %*% R_star, T)
      # eigen vectors
      Q2 = eigenR_star$vectors
      # eigen values
      Lmbd2 = eigenR_star$values
    }
  }
  
  # update A (useful for preprocessed or original version)
  if((!penalt_A)|(t_length == 1)){
    # if no preprocessing or length = 1
    newA = Y %*% t(newS) %*% solve(newS%*%t(newS)) 
  }else{
    # calculate eigen values and vectors of SS'
    eigenS = eigen(newS %*% t(newS),T)
    # eigen vectors
    Q1 = eigenS$vectors
    # eigen values
    Lmbd1 = eigenS$values
    # construct D
    D = t(Q2) %*% Y %*% t(newS) %*% Q1 
    
    # calculate a scaler for maintaining nu at 0-10 level
    n_digits = nchar(as.integer(floor(max(abs(D)))))
    scaler = 10^n_digits
      
    # calculate A_star
    A_star = matrix(nrow = dim(D)[1], ncol = dim(D)[2])
    for(i in 1:(dim(D)[1])){
      for(j in 1:(dim(D)[2])){
        A_star[i, j] = D[i,j]/(nu*scaler*Lmbd2[i] + Lmbd1[j])
      }
    }
    # calculate A_tilde or newA
    newA = Q2 %*% A_star %*% t(Q1)
  }

  # scale A so that each column has unit variance
  for(i in 1:q){
    ai = sd(newA[,i])
    theta_new[[i]]$lam_l = theta_new[[i]]$lam_l * ai
    newA[,i] = newA[,i] / ai
    newS[i,] = newS[i,] * ai
  }
  
  # save M_l, X_l into theta_new
  # g-inverse of newA
  Mnew = solve(t(newA) %*% newA) %*% t(newA) 
  for(l in 1:q){
    theta_new[[l]]$M_l = Mnew[l,]
  }
  
  return(list(A = newA, theta = theta_new, S = newS))
}
