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
#' @param silent If \code{TRUE}, print out the penalty added on A and S. Defaults to \code{FALSE}.
#' 
#' @return Updated values of loading matrix \code{A} and latent source \code{S}.
dlocus_update = function(Y, A, theta, 
                         q, K, V, n_subject, t_length, 
                         penalt_S = TRUE, lambda_ch, gamma, 
                         penalt_A = FALSE, nu, H_star,
                         preprocess, silent = FALSE){
  
  # whether add penalty on S (latent sources)
  if(is.null(penalt_S)){
    if(!silent) cat(paste("Locus without penalty on latent sources"))
  }else{
    if(!silent) cat(paste("Locus with", penalt_S, "penalty on latent sources"))
  }
  
  # whether add smooth penalty on A (loading matrix)
  if(!penalt_A){
    if(!silent) cat(paste("Locus without penalty on loading matrix"))
  }else{
    if(t_length == 1){
      if(!silent) cat(paste("Locus without penalty on loading matrix"))
    }else{
      if(!silent) cat(paste("Locus with smooth penalty on loading matrix"))
    }
  }

  # use theta to store ic-specific parameters and channels
  theta_new = list()
  
  # node to node edge list 
  rmat = matrix(rep(1:V, V), ncol=V)
  cmat = t(rmat)
  Lcoorsym = matrix(0, ncol=2, nrow=V*(V-1)/2)
  Lcoorsym[,1] = Ltrans(rmat, F)
  Lcoorsym[,2] = Ltrans(cmat, F)
  
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
  
  # for each ic, conditioning on others, estimate latent channels X
  newS = array(dim = c(q, K))
  for(curr_ic in 1:q){
    theta_new[[curr_ic]] = list()
    
    # update Xl: 
    theta_ic = theta[[curr_ic]]
    # D inverse (of length R)
    Dinverse = 1/theta[[curr_ic]]$lam_l 
    # if Xl is NA then set it to be 0
    theta[[curr_ic]]$X_l[is.na(theta[[curr_ic]]$X_l)] = 0
    # lth row of A_tilde transpose (1*q) = lth column of A_tilde (1*q)
    # times Y (q*p)
    # result in a 1*p vector -> transpose -> p*1 vector
    Yic = t(theta[[curr_ic]]$M_l%*%Y)  
    # update starting from node 1
    v = 1 
    while(v <= V){
      # Xl with the vth column removed -> transpose -> (V-1)*Rl
      # *** here, the Xl is t(Xl) in the paper
      # here: Sl = t(Xl)Dl(Xl), Xl(Rl*V)
      Hlv = t(theta[[curr_ic]]$X_l[,-v])   
      # Yic(p*1 vector) with elements relate to node v -> (V-1)*1
      yvpen = Yic[which(Lcoorsym[,1] == v | Lcoorsym[, 2] == v),]  
      # inverse of t(Hlv)*Hlv
      Sigmalv = psdmat_inverse(t(Hlv)%*%Hlv)
      # solution based on different penalization
      if(is.null(penalt_S)){
        beta = yvpen
      }else if(penalt_S == "SCAD"){
        if(gamma<=2){print("Gamma needs to be > 2!"); gamma = 2.01}
        beta = SCAD_func(yvpen, lambda_ch = lambda_ch, gamma = gamma)
      }else if(penalt_S == "L1"){
        beta = sign(yvpen)*(abs(yvpen)-lambda_ch)*(abs(yvpen)>=lambda_ch)
      }else if(penalt_S == "Hardthreshold"){
        beta = yvpen*(abs(yvpen)>=lambda_ch)
      }else{
        print("No Penalty available!")
        stop()
      }
      if(sd(beta) == 0){
        theta[[curr_ic]]$X_l[,v] = 0
      }else{
        theta[[curr_ic]]$X_l[,v] = (Dinverse*(Sigmalv%*%t(Hlv)%*%beta))/sd(beta)*sd(yvpen) 
      }
      v = v+1
    }
    theta_new[[curr_ic]]$X_l = theta[[curr_ic]]$X_l
    
    # update Dl: 
    # calculate Zl with rth column being L(Xl%*%t(Xl)) in the paper
    Xstarstack = t(apply(theta[[curr_ic]]$X_l, 1, function(x){ 
      x = matrix(x, ncol=1)
      return(Ltrans(x%*%t(x),F))
      }))
    
    if(is.null(penalt_S)){
      beta = Yic
    }else if(penalt_S == "SCAD"){
      if(gamma<=2){print("Gamma needs to be > 2!");gamma = 2.01}
      beta = SCAD_func(Yic, lambda_ch = lambda_ch, gamma = gamma)
    }else if(penalt_S == "L1"){
      beta = sign(Yic)*(abs(Yic)-lambda_ch)*(abs(Yic)>=lambda_ch)
    }else if(penalt_S == "Hardthreshold"){
      beta = Yic*(abs(Yic)>=lambda_ch)
    }
    
    theta_new[[curr_ic]]$lam_l = as.numeric(solve(Xstarstack%*%t(Xstarstack))%*%Xstarstack%*%beta /sd(beta)*sd(Yic)) # Rx1
    
    # calculate Sl based on Xl and Dl
    newS[curr_ic,] = Ltrans(t(theta_new[[curr_ic]]$X_l)%*% 
                              diag(theta_new[[curr_ic]]$lam_l) %*%
                              theta_new[[curr_ic]]$X_l,F)  
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
  
  # orthnologize A
  # newA = newA %*% real( solve( t(newA)%*%newA )^(1/2)); # This could lead to numerical problems. 
  
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
