#' A novel decomposition method for brain network dynamic connectivity matrices 
#' using low-rank structure with uniform sparsity and temporal group lasso
#' 
#' @param Y Group-level connectivity data of dimension \code{N} \times \code{K}. 
#' \code{N} is number of subjects, \code{K} is number of edges. 
#' @param q The number of latent sources one would like to extract.
#' @param V The number of nodes.
#' @param n_subject The number of subjects in the study.
#' @param preprocess If \code{TRUE}, the concatenated group-level connectivity data will be preprocessed. 
#' Defaults to \code{TRUE}.
#' @param penalty_S The penalty added on latent sources. Defaults to \code{"SCAD"}.
#' @param phi A tunning parameters with respect to the element-wise penalty on S.
#' @param penalty_A If \code{TRUE}, temporal smooth lasso penalty is added on loadings. Defaults to \code{FALSE}.
#' @param nu A numeric tunning parameter with respect to the temporal smooth lasso penalty on A. Defaults to \code{NULL}.
#' @param maxIteration The number of maximum iterations. Defaults to value 100.
#' @param approximation If \code{TRUE}, use approximated algorithm based on SVD for update.
#' Defaults to \code{TRUE}. Don't change if you are not sure.
#' @param espli1 A number describing toleration for convergence on S. Defaults to value 0.01.
#' @param espli2 A number describing toleration for convergence on A. Defaults to value 0.01.
#' @param rho A tuning parameter for selecting number of ranks in each subnetwork's decomposition. Defaults to value 0.9.
#' @param silent If \code{TRUE}, print out the penalty added on A and S. Defaults to \code{FALSE}.
#' 
#' @export
#' 
#' @return The resulting loading matrix \code{A} and latent sources \code{S} based on the decomposition method
#' @import ica far

dlocus = function(Y, q, V, n_subject, 
                  preprocess = TRUE, penalty_S = "SCAD", phi = 2, penalty_A = FALSE, nu = NULL,
                  maxIteration = 100, approximation = TRUE, 
                  espli1 = 0.01, espli2 = 0.01, rho = 0.9, silent = FALSE){
  # demean Y 
  Y = sweep(Y, 2, apply(Y, 2, mean), "-") 
  
  # total number of observations
  N = dim(Y)[1]
  # number of edges
  K = dim(Y)[2]
  # time length
  t_length = N/n_subject
  # verify number of nodes
  if(V != (sqrt(1+8*K)+1)/2){
    print("V is not correctly specified! Please double check the dimension of your input data.")
    stop()
  }
  
  # preprocess Y based on choice
  if(preprocess){ 
    Yraw = Y
    prep_rslt = dlocus_preprocess(Y,q) 
    Y = prep_rslt$Ynew
    H_star = prep_rslt$H_star
  }else{
    H_star = NULL
  }
  
  # initial estimation based on ica
  theta_ini = dlocus_initial(Y, q, V, rho = rho)
  A = theta_ini$A
  S = theta_ini$S
  theta = theta_ini$theta
  
  # update parameters
  Iter = 1
  while(Iter <= maxIteration){
    if(approximation){
      # update with approximation
      theta_new = dlocus_update_approx(Y = Y, A = A, theta = theta, 
                                       q = q, K = K, V = V, n_subject = n_subject, t_length = t_length,
                                       penalt_S = penalty_S, lambda_ch = phi, gamma = 2.1,
                                       penalt_A = penalty_A, H_star = H_star, nu = nu, 
                                       preprocess = preprocess, imput_method = "Previous", silent = silent)
    }else{
      # update without approximation
      theta_new = dlocus_update(Y = Y, A = A, theta = theta, 
                                q = q, K = K, V = V, n_subject = n_subject, t_length = t_length,
                                penalt_S = penalty_S, lambda_ch = phi, gamma = 2.1,
                                penalt_A = penalty_A, H_star = H_star, nu = nu, 
                                preprocess = preprocess, silent = silent)
    }
    
    # orthogonize A here
    if(preprocess){
      A_new = far::orthonormalization(theta_new$A)
    }else{
      A_new = theta_new$A
    }
    
    # estimation results from one iteration
    A_new = theta_new$A
    S_new = theta_new$S
    theta_new  = theta_new$theta
    
    # calculate error from ica based S
    errS = norm(as.matrix(S_new-S))/norm(as.matrix(S))
    errA = norm(as.matrix(A_new-A))/norm(as.matrix(A))
    
    # if any NA values generated, break the iteration
    if(sum(is.na(c(errS, errA))) > 0){
      return(list(conver = F, A = A, S = S, theta = theta))
    }
    
    if(!silent){
      print(paste("Iter ", Iter, "; Percentage change on S: ", round(errS, 3),
                  "; Percentage change on A: ", round(errA, 3), "." , sep = ""))
    }
    
    A = A_new
    S = S_new
    theta = theta_new
    
    # print(performance_ASq3(S,Struth,Atrue,Atrue)) # this is only for simulation test
    
    if(errA < espli1 & errS < espli2){
      if(!silent){
        print("Converged!")
      }
      if(preprocess){
        # transform back to the original scale
        A = H_star %*% A
      }
      return(list(conver = T, A = A, S = S, theta = theta))
    }
    
    Iter = Iter + 1
  }
  
  # if maximum iteration is reached
  if(!silent){
    print("Failed to converge!")
  }
  if(preprocess){
    # transform back to the original scale
    A = H_star %*% A
  }
  return(list(conver = F, A = A, S = S, theta = theta))
}
