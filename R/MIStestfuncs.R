# multiICSKAT functions part 1

#' Left Survival Function
#'
#' Calculate the left surivival function value.
#' @param l observation value of interest
#' @param d number of Gaussian Quadrature points
#' @param temp_beta vector of parameter estimates
#' @param phen list of data matrices containing both left and right information
#' @param r1 nodes of quadrature points
#' @param k number of outcomes
surv_left <- function(l, d, temp_beta, phen, r1, k){
  
  # get beta values for chosen outcome
  covcol <- ncol(phen$dmats$right_dmat)
  sub_beta <- temp_beta[((l-1)*covcol + 1): (l *covcol)]
  
  # get sigma squared value from parameter list
  sigmasq <- temp_beta[k*covcol + 1]
  
  
  #left and right times + censoring
  lt <- phen$lt
  rt <- phen$rt
  tpos_ind <- as.numeric(lt > 0)
  obs_ind <- as.numeric(rt != Inf)
  
  #left design matrix
  left_dmat <- phen$dmats$left_dmat
  
  #Calculate left survival times
  hl1 <- as.numeric(exp(left_dmat %*% t(matrix(sub_beta, nrow = 1)) + sqrt(2 * sigmasq)* r1[d]))
  sl1 <- ifelse(tpos_ind == 0, 1, exp(-hl1))
  
  #return the survival terms
  return(sl1)
  
}

#' Calculate the right survival function value.
#' @param l observation value of interest
#' @param d number of Gaussian Quadrature points
#' @param temp_beta vector of parameter estimates
#' @param phen list of data matrices containing both left and right information
#' @param r1 nodes of quadrature points
#' @param k number of outcomes
# Right Surivival Function
surv_right <- function(l, d, temp_beta, phen, r1, k){
  
  # get beta values for chosen outcome
  covcol <- ncol(phen$dmats$right_dmat)
  sub_beta <- temp_beta[((l-1)*covcol + 1): (l *covcol)]
  
  # get sigma squared value from parameter list
  sigmasq <- temp_beta[k*covcol + 1]
  
  
  #left and right times + censoring
  lt <- phen$lt
  rt <- phen$rt
  tpos_ind <- as.numeric(lt > 0)
  obs_ind <- as.numeric(rt != Inf)
  
  
  # right design matrix
  right_dmat <- phen$dmats$right_dmat
  
  
  #Calculate right survival times
  hr1 <- as.numeric(exp(right_dmat %*% t(matrix(sub_beta, nrow = 1)) + sqrt(2 * sigmasq)* r1[d]))
  sr1 <- ifelse(obs_ind == 0, 0, exp(-hr1))
  sr1[!is.finite(sr1)] <- 0
  
  return(sr1)
  
}

#' Calculate the left hazard function value.
#' @param l observation value of interest
#' @param d number of Gaussian Quadrature points
#' @param temp_beta vector of parameter estimates
#' @param phen list of data matrices containing both left and right information
#' @param r1 nodes of quadrature points
#' @param k number of outcomes
# Left Hazard Function
haz_left <- function(l, d, temp_beta, phen, r1, k){
  
  # get beta values for chosen outcome
  covcol <- ncol(phen$dmats$right_dmat)
  sub_beta <- temp_beta[((l-1)*covcol + 1): (l *covcol)]
  
  # get sigma squared value from parameter list
  sigmasq <- temp_beta[k*covcol + 1]
  
  
  #left and right times + censoring
  lt <- phen$lt
  rt <- phen$rt
  tpos_ind <- as.numeric(lt > 0)
  obs_ind <- as.numeric(rt != Inf)
  
  # left design matrix
  left_dmat <- phen$dmats$left_dmat
  
  # calculate left hazard terms
  hl1 <- as.numeric(exp(left_dmat %*% t(matrix(sub_beta, nrow = 1)) + sqrt(2 * sigmasq)* r1[d]))
  return(hl1)
}

#' Calculate the right hazard function value.
#' @param l observation value of interest
#' @param d number of Gaussian Quadrature points
#' @param temp_beta vector of parameter estimates
#' @param phen list of data matrices containing both left and right information
#' @param r1 nodes of quadrature points
#' @param k number of outcomes
# Right Hazard Function
haz_right <- function(l, d, temp_beta, phen, r1, k){
  
  # get beta values for chosen outcome
  covcol <- ncol(phen$dmats$right_dmat)
  sub_beta <- temp_beta[((l-1)*covcol + 1): (l *covcol)]
  
  # get sigma squared value from parameter list
  sigmasq <- temp_beta[k*covcol + 1]
  
  
  #left and right times + censoring
  lt <- phen$lt
  rt <- phen$rt
  tpos_ind <- as.numeric(lt > 0)
  obs_ind <- as.numeric(rt != Inf)
  
  #right design matrix
  right_dmat <- phen$dmats$right_dmat
  
  # calculate right hazard times
  hr1 <- as.numeric(exp(right_dmat %*% t(matrix(sub_beta, nrow = 1)) + sqrt(2 * sigmasq)* r1[d]))
  
  return(hr1)
  
}

#' Calculate the difference between left and right survival terms.
#' @param l observation value of interest
#' @param d number of Gaussian Quadrature points
#' @param temp_beta vector of parameter estimates
#' @param phen list of data matrices containing both left and right information
#' @param r1 nodes of quadrature points
#' @param k number of outcomes
# Surv diff
surv_diff <- function(l,d,temp_beta, phen, r1, k){
  
  (surv_left(l, d, temp_beta, phen, r1, k) - surv_right(l, d, temp_beta, phen, r1, k) )
  
}


#' Calculate the product of the difference between survival terms excluding that of the outcome of interest
#' @param l observation value of interest
#' @param k number of outcomes
#' @param store array of difference between left and right survival values
without_one_phen <- function(l, k, store){
  
  #if total two outcomes
  if(k == 2){
    
    out <- store[,-l,]
    
    return(out)
    
    # case for more than two outcomes
  } else {
    
    #make index of outcomes not including the one of interest
    idx <- (1:k)[-l]
    
    # subset the data
    sub <- store[,idx,]
    
    # multiply the survival differences across the other outcomes (idx)
    sub_prod <- apply(sub, c(1,3), prod) 
    
    return(sub_prod)
  }
  
}

#' Calculate the product of the difference between survival terms excluding that of two outcomes of interest
#' @param l observation value of interest
#' @param m second observation of interest
#' @param k number of outcomes
#' @param store array of difference between left and right survival values
#' @param n number of observations
#' @param d number of quadrature nodes
#' @return A list with the elements:
without_two_phen <- function(l,m, k, store, n, d){
  
  # get index of outcomes for all not equal to outcomes l and m
  idx <- (1:k)[-c(l,m)]
  
  # subset the array of differences
  sub <- array(store[,idx,], dim=c (n,length(idx),d))
  
  # multiply the survival differences across the other outcomes (idx)
  sub_prod <- apply(sub, c(1,3), prod) 
  
}

#' Calculate the denominator term of the derivative of the likelihood
#' @param store array of difference between left and right survival values
#' @param weights quadrature weights
#' @param d number of quadrature nodes
#' @param n number of observations
get_A <- function(store, weights, d, n){
  
  #number of outcomes
  k <- ncol(store[,,d])
  
  # multiply the survival differences across number of outcomes
  mult_across_k <- array(apply(store, c(1,3), prod), dim = c(n,1,d))
  
  # add quadrature weights
  add_weights <- t(matrix(mult_across_k, nrow = n)) * weights
  
  # sum terms and divide by sqrt(pi)
  A_i <- colSums(add_weights)/sqrt(pi)
  
  return(A_i)
  
}