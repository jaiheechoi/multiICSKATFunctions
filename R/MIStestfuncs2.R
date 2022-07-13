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


################################################# F3 NR TERMS
fd_term <- function(l, temp_beta, phen,d, apply_diffs, A_i, no_l_all,HL_array, HR_array){

  # quadrature weights and roots
  ghDat <- fastGHQuad::gaussHermiteData(d)
  w1 <- ghDat$w
  r1 <- ghDat$x


  # Get the survival (exp(-exp(eta))) / hazard (exp(eta)) terms
  # for given l, for all 100 weights
  # Result is a list of 100, each n x 1
  hl_d <- HL_array[,,l]
  hr_d <- HR_array[,,l]


  # prod l not equal to k, SL - SR
  mult_k_without_l <- no_l_all[,,l]


  # just the first derivative term times the weight
  # weight_d*(S(L)(-H(L))U - S(R)(-H(R))V)
  first_deriv <- function(l, d, sl_d, sr_d, hl_d, hr_d, phen){


    # get design matrices
    left_dmat <- phen$dmats$left_dmat
    right_dmat <- phen$dmats$right_dmat

    #left and right times + censoring
    lt <- phen$lt
    rt <- phen$rt
    tpos_ind <- as.numeric(lt > 0)
    obs_ind <- as.numeric(rt != Inf)


    # first derivative terms
    U1 <- left_dmat * ifelse(tpos_ind == 0, 0, (exp(-hl_d[,d]) * -hl_d[,d]))
    U2 <- right_dmat * ifelse(obs_ind == 0, 0, (exp(-hr_d[,d]) * -hr_d[,d] ))

    # check to make sure there are no NAs
    U2[is.na(U2)] <- 0

    # the whole term is a difference between the left values and right values
    inside <- U1 - U2

    return(inside)
  }


  #Get apply for all 100 weights
  #Result is a list of 100, each 1000 x 5
  insides <- lapply(1:d, first_deriv, l = l, hl_d = hl_d, hr_d = hr_d, phen = phen)


  #Make function that multiplies the prod of surv-diffs and the inside
  mult_together <- function(d, arrayA, listB, weights){

    # multiply all the terms together
    arrayA[,d]* listB[[d]] * weights[d]
  }


  #make for all 100 weights
  deriv_prod <- lapply(1:d, mult_together, arrayA = mult_k_without_l, listB = insides, weights = w1)

  # Make the list to an array
  dp_array <- simplify2array(deriv_prod)


  # Sum over D
  # output is 1000 x 5
  sum_over_d <- apply(dp_array, c(1,2), sum)


  # Combine with other values and sum over n
  # result is 1 x 5
  fd <- apply(((1/sqrt(pi)) *(sum_over_d/A_i)), 2, sum)


  return(fd)

}



sd_on <- function(l, k, temp_beta, phen, d, apply_diffs, A_i, no_l_all, HL_array, HR_array)
{
  # quadrature weights and roots
  ghDat <- fastGHQuad::gaussHermiteData(d)
  w1 <- ghDat$w
  r1 <- ghDat$x


  # prod l not equal to k, SL - SR
  no_l <- no_l_all[,,l]


  # left and righ design matrices

  left_dmat <- phen$dmats$left_dmat
  right_dmat <- phen$dmats$right_dmat

  #left and right times + censoring terms
  lt <- phen$lt
  rt <- phen$rt
  tpos_ind <- as.numeric(lt > 0)
  obs_ind <- as.numeric(rt != Inf)


  # survival terms
  hl_d <- HL_array[,,l]
  hr_d <- HR_array[,,l]


  # second derivative term
  get_sd <- function(hl_d, hr_d, phen, d, no_l){

    # left and right design matrices
    left_dmat <- phen$dmats$left_dmat
    right_dmat <- phen$dmats$right_dmat


    #left and right times + censoring terms
    lt <- phen$lt
    rt <- phen$rt
    tpos_ind <- as.numeric(lt > 0)
    obs_ind <- as.numeric(rt != Inf)



    # first derivative terms
    ul_1 <- ifelse(tpos_ind == 0, 0, -hl_d[,d] * exp(-hl_d[,d]) + (hl_d[,d]^2 * exp(-hl_d[,d])))
    ur_1 <- ifelse(obs_ind == 0, 0, -hr_d[,d] * exp(-hr_d[,d]) + (hr_d[,d]^2 * exp(-hr_d[,d])))
    ur_1[which(is.na(ur_1))] <- 0

    # second derivative terms
    sd_term1 <- t(left_dmat) %*% ( (no_l[,d] * as.numeric(ul_1/A_i)) * left_dmat)
    sd_term2 <- t(right_dmat) %*% ( (no_l[,d]* as.numeric(ur_1)/A_i) * right_dmat )


    # difference between left and right term multiplied by GQ weights
    sd_5x5 <- (sd_term1 - sd_term2) * w1[d]

    return(sd_5x5)
  }

  # apply the derivatives to each node of the quadrature d
  derivs <- lapply(1:d, get_sd, hl_d = hl_d, hr_d = hr_d, phen = phen, no_l = no_l)

  # Make the list to an array
  derivs_array <- simplify2array(derivs)

  # Sum over n and divide by sqrt(pi)
  term1 <- apply(derivs_array, c(1,2), sum)/sqrt(pi)


  ################################################
  # term 2

  first_deriv <- function(l, d, hl_d, hr_d, phen){

    # left and right design matrices
    left_dmat <- phen$dmats$left_dmat
    right_dmat <- phen$dmats$right_dmat


    #left and right times + censoring terms
    lt <- phen$lt
    rt <- phen$rt
    tpos_ind <- as.numeric(lt > 0)
    obs_ind <- as.numeric(rt != Inf)

    #first derivative terms
    U1 <- left_dmat * ifelse(tpos_ind == 0, 0, (exp(-hl_d[,d]) * -hl_d[,d]/A_i))
    U2 <- right_dmat * ifelse(obs_ind == 0, 0, (exp(-hr_d[,d]) * -hr_d[,d]/A_i ))
    U2[is.na(U2)] <- 0

    # whole term is the difference between left and right terms
    inside <- U1 - U2

    return(inside)

  }

  #Get apply for all 100 weights
  #Result is a list of 100, each 1000 x 5
  insides <- lapply(1:d, first_deriv, l = l, hl_d = hl_d, hr_d = hr_d, phen = phen)


  #Make function that multiplies the prod of surv-diffs and the inside
  mult_together <- function(d, arrayA, listB, weights){
    arrayA[,d]* listB[[d]] * weights[d]
  }


  #make for all 100 weights
  deriv_prod <- lapply(1:d, mult_together, arrayA = no_l, listB = insides, weights = w1)


  # Make the list to an array
  dp_array <- simplify2array(deriv_prod)


  # Sum over D
  # output is 1000 x 5
  term2_a <- apply(dp_array, c(1,2), sum)/sqrt(pi)

  # full term is the term multiplied to itself
  term2 <- (t(term2_a) %*% term2_a)


  #subtract the two terms
  sd_on <- term1 - term2

  return(sd_on)

}



sd_off <- function(l, m, phen_l, phen_m, temp_beta, d, apply_diffs,A_i, HL_array, HR_array, no_l_all, no_two_all, tpos_all, obs_all,k){

  # quadrature weights and roots
  ghDat <- fastGHQuad::gaussHermiteData(d)
  w1 <- ghDat$w
  r1 <- ghDat$x

  # get the left and right hazard for observation l
  hl_d <- HL_array[,,l]
  hr_d <- HR_array[,,l]

  #left and right design matrices for observation l
  ld_l <- phen_l$dmats$left_dmat
  rd_l <- phen_l$dmats$right_dmat

  #left and right design matrices for observation m
  ld_m <- phen_m$dmats$left_dmat
  rd_m <- phen_m$dmats$right_dmat


  #First derivative terms
  U1 <- sweep((exp(-HL_array) * -HL_array), c(1,3), ifelse(tpos_all == 0, 0,1), FUN = "*")
  U2 <- sweep((exp(-HR_array) * -HR_array), c(1,3), ifelse(obs_all == 0, 0,1), FUN = "*")
  U2[is.na(U2)] <- 0

  # product of diff of survival without l
  surv_no_l <- no_l_all

  # term looks different for k = 2 observations vs more than 2 observations
  if(k == 2){
    sg_sd <- function(d, l, m) {

      #second derivative term
      w1[d]* ( t( ld_l/A_i *  U1[,d,l] -  rd_l/A_i* U2[,d,l] ) %*%
                 ( ld_m* U1[,d,m] -  rd_m *U2[,d,m]  ) )

    }

    # apply to all quadrature nodes
    term1 <- apply(simplify2array(lapply(1:d, sg_sd, l = l, m = m)), c(1,2), sum)/sqrt(pi)

  } else {

    # Product of survival terms without two phenotypes
    # get combination of all the indices of outcomes
    combs <- combn(1:k, 2)

    # order the observation indices
    min_k <- min(l,m)
    max_k <- max(l,m)

    # get the product of the differences of survival without observation l and m
    no_l_m <- no_two_all[,,which(combs[1,] == min_k  & combs[2,] == max_k)]

    # Function for the second deriv
    sg_sd <- function(d, l, m) {

      # term for the second derivative
      w1[d]* ( t( ld_l * (no_l_m[,d]/A_i * U1[,d,l]) -  rd_l*(no_l_m[,d]/A_i* U2[,d,l]) ) %*%
                 ( ld_m* U1[,d,m] -  rd_m *U2[,d,m]  ) )

    }

    #apply first derivative term to all
    term1 <- apply(simplify2array(lapply(1:d, sg_sd, l = l, m = m)), c(1,2), sum)/sqrt(pi)
  }

  # the first derivative of observation l
  fd_term_l <- function(d, l){

    w1[d] * (ld_l*U1[,d,l] - rd_l*U2[,d,l]) * (no_l_all[,d,l]/A_i)

  }

  #the first derivative of observation m
  fd_term_m <- function(d, m){

    w1[d] * (ld_m*U1[,d,m] - rd_m*U2[,d,m]) * (no_l_all[,d,m]/A_i)

  }

  # apply the first derivative functions to all the quadrature nodes
  term2_a <- apply(simplify2array(lapply(1:d, fd_term_l, l = l)), c(1,2), sum)/sqrt(pi)
  term2_b <- apply(simplify2array(lapply(1:d, fd_term_m, m = m)), c(1,2), sum)/sqrt(pi)

  # Multiply the terms together, matching the index for weight d
  term2 <- (t(term2_a) %*% term2_b)

  # combine all the terms for the second derivative
  sd_off <- term1 - term2

  return(sd_off)

}


############################################### F4 Sigma function
# SIGMA



ss_fd <- function(l, phen, HL_array, HR_array, tpos_all, obs_all, apply_diffs, temp_beta, A_i, no_l_all, k, d){

  # quadrature weights and roots
  ghDat <- fastGHQuad::gaussHermiteData(d)
  w1 <- ghDat$w
  r1 <- ghDat$x

  # define sigma term
  nocol <- ncol(phen$dmats$right_dmat)
  sigmasq <- temp_beta[k*nocol + 1]

  #calculate the first derivative terms
  U1 <- sweep((exp(-HL_array) * -HL_array), c(1,3), ifelse(tpos_all == 0, 0,1), FUN = "*")
  U2 <- sweep((exp(-HR_array) * -HR_array), c(1,3), ifelse(obs_all == 0, 0,1), FUN = "*")
  U2[is.na(U2)] <- 0

  # get the product of the differences of survival for all k not outcome of interest
  surv_no_l <- no_l_all

  # multiply the terms together and sum across observations
  fd_out <- apply((U1 - U2) * surv_no_l, c(1,2), sum)


  # get the total first derivative term
  deriv <- sum(rowSums(sweep(fd_out, 2, w1 * r1/sqrt(2*sigmasq), FUN = "*"))/A_i)/sqrt(pi)

  return(deriv)

}


ss_sd <- function(HL_array, HR_array, xAll, apply_diffs, temp_beta, A_i, no_l_all, no_two_all, k, d){

  # quadrature weights and roots
  ghDat <- fastGHQuad::gaussHermiteData(d)
  w1 <- ghDat$w
  r1 <- ghDat$x

  # define sigma term
  nocol <- ncol(xAll$xDats[[1]]$dmats$right_dmat)
  sigmasq <- temp_beta[k*nocol + 1]

  # get censoring terms
  tpos_all <- xAll$ts_all
  obs_all <- xAll$ob_all

  # get first derivative terms
  U1 <- sweep((exp(-HL_array) * -HL_array), c(1,3), ifelse(tpos_all == 0, 0,1), FUN = "*")
  U2 <- sweep((exp(-HR_array) * -HR_array), c(1,3), ifelse(obs_all == 0, 0,1), FUN = "*")
  U2[is.na(U2)] <- 0


  # make function for the third term of the derivative
  # Looks different for k = 2 vs more than 2 outcomes
  term_3_rep <- function(l){

    # get index of all outcomes not including outcome of interest
    idx <- (1:k)[-l]

    # get combination of all possible pairs of outcomes
    combs <- combn(1:k, 2)

    # make function to apply to above indices
    to_idx <- function(m){

      # order outcomes
      min_k <- min(l,m)
      max_k <- max(l,m)

      if(k == 2){
        out <- (U1-U2)[,,l] * (U1-U2)[,,m]
      } else {
        out <- (U1-U2)[,,l] * (U1-U2)[,,m] * no_two_all[,,which(combs[1,] == min_k  & combs[2,] == max_k)]
      }
      return(out)
    }


    # apply function to all the index of other outcomes
    apply_idx <- simplify2array(lapply(idx, to_idx))

    #sum across all observations
    return(apply(apply_idx, c(1,2), sum))

  }

  # apply above function to all observations
  term_three <- simplify2array(lapply(1:k, term_3_rep))

  # get product of differences between survival terms for all but observation of interest
  surv_no_l <- no_l_all

  # calculate second derivative terms
  Y1 <- sweep(-HL_array * exp(-HL_array) + (HL_array^2 * exp(-HL_array)), c(1,3), ifelse(tpos_all == 0, 0,1), FUN = "*")
  Y2 <- sweep(-HR_array * exp(-HR_array) + (HR_array^2 * exp(-HR_array)), c(1,3), ifelse(obs_all == 0, 0,1), FUN = "*")
  Y2[is.na(Y2)] <- 0

  # get the first derivative terms
  S1 <- sweep((exp(-HL_array) * -HL_array), c(1,3), ifelse(tpos_all == 0, 0,1), FUN = "*")
  S2 <- sweep((exp(-HR_array) * -HR_array), c(1,3), ifelse(obs_all == 0, 0,1), FUN = "*")
  S2[is.na(S2)] <- 0

  # combine all the terms to get full second derivative
  score <- sweep((U1 - U2) * surv_no_l, 2, -r1/(2*sqrt(2*sigmasq)), FUN = "*") +
    sweep((Y1 - Y2) * surv_no_l, 2, (r1/sqrt(2 * sigmasq))^2, FUN = "*") +
    sweep(term_three,2, (r1/sqrt(2 * sigmasq))^2, FUN = "*")

  # sum across all observations
  out_sumk <- apply(score, c(1,2), sum)

  # sum across node number and multiply by weights
  sum_d <- apply((t(out_sumk) * w1), 2, sum)

  # sum across all subjects and divide by sqrt(pi)
  ss_sd_t1 <- sum(sum_d/A_i)/sqrt(pi)


  # get the second term in the derivative
  ss_sd_t2 <- sum((apply(t(apply((S1 - S2)* surv_no_l, c(1,2), sum))*w1 * r1/sqrt(2*sigmasq),2,sum)/A_i/sqrt(pi))^2)

  # combine the two terms
  return(ss_sd_t1 - ss_sd_t2)
}




st_off<- function(l, HL_array, HR_array, xAll, apply_diffs, temp_beta, A_i, no_l_all, no_two_all, k, d) {

  # quadrature weights and roots
  ghDat <- fastGHQuad::gaussHermiteData(d)
  w1 <- ghDat$w
  r1 <- ghDat$x

  # for phenotype l, get incides of other observations
  idx <- (1:k)[-l]

  # subset the data for just the observation of interest
  phen <- xAll$xDats[[l]]

  # get censoring terms
  tpos_all <- xAll$ts_all
  obs_all <- xAll$ob_all

  # left and right design matrices
  ldm <- phen$dmats$left_dmat
  rdm <- phen$dmats$right_dmat

  # define sigma term
  nocol <- ncol(ldm)
  sigmasq <- temp_beta[k*nocol + 1]

  #calculate the second derivative terms
  Y1 <- sweep(-HL_array * exp(-HL_array) + (HL_array^2 * exp(-HL_array)), c(1,3), ifelse(tpos_all == 0, 0,1), FUN = "*")
  Y2 <- sweep(-HR_array * exp(-HR_array) + (HR_array^2 * exp(-HR_array)), c(1,3), ifelse(obs_all == 0, 0,1), FUN = "*")
  Y2[is.na(Y2)] <- 0

  # calculate the first derivative terms
  U1 <- sweep((exp(-HL_array) * -HL_array), c(1,3), ifelse(tpos_all == 0, 0,1), FUN = "*")
  U2 <- sweep((exp(-HR_array) * -HR_array), c(1,3), ifelse(obs_all == 0, 0,1), FUN = "*")
  U2[is.na(U2)] <- 0

  # product of diff of survival without l
  surv_no_l <- no_l_all

  # second derivative function
  # term is different for k
  st_sd <- function(d) {
    if(k == 2){
      out <- (ldm *(surv_no_l[,d,l] * Y1[,d,l])) - (rdm * (surv_no_l[,d,l]* Y2[,d,l]) ) +
        ( (ldm *U1[,d,l]) - (rdm *U2[,d,l] ) ) * (U1 - U2)[,d,idx]
    } else{
      #Product of survival terms without two phenotypes
      combs <- combn(1:k, 2)
      with_l_idx <- which(apply(combs, 2, function(x) which(x == l)) == 1 | apply(combs, 2, function(x) which(x == l)) == 2)
      surv_diff_no_two <- no_two_all[,,with_l_idx]
      out <- (ldm *(surv_no_l[,d,l] * Y1[,d,l])) - (rdm * (surv_no_l[,d,l]* Y2[,d,l]) ) +
        ( (ldm *U1[,d,l]) - (rdm *U2[,d,l] ) ) * apply((U1 - U2)[,d,idx]* surv_diff_no_two[,d,], 1, sum)

    }
    return(out)
  }


  # Apply to all the node points
  to_d <- apply(sweep(simplify2array(lapply(1:d, st_sd)), 3, w1 * (r1/sqrt(2*sigmasq)), FUN = "*"), c(1,2), sum)


  # first term, sum over d and divide by sqrt(pi)
  st_t1 <- colSums(sweep(to_d, 1, A_i, FUN = "/"))/sqrt(pi)

  # First deriv times the data mat
  term3_func <- function(d){
    (ldm * U1[,d,l] * surv_no_l[,d,l] ) - (rdm * U2[,d,l] * surv_no_l[,d,l] )
  }


  # THe term that is like the first deriv of the theta term
  t2a <- sweep(apply(sweep(simplify2array(lapply(1:d, term3_func)), 3, w1, FUN = "*"), c(1,2), sum),1, A_i, FUN = "/")/sqrt(pi)

  # The term that is like the first deriv of the sigma term
  t2b <- rowSums(sweep(apply((U1 - U2) * surv_no_l/A_i, c(1,2), sum),2, w1*r1/sqrt(2*sigmasq), FUN = "*"))/sqrt(pi)


  # Sum over n after multiplying both terms
  st_t2 <- colSums(t2a *t2b)

  # combine both terms to get the full derivative
  deriv <- st_t1 - st_t2

  return(deriv)
}


################################################### F5 Gamma terms
gamma_fd <- function(l, HL_array, HR_array, tpos_all, obs_all, temp_beta, A_i, no_l_all, gMat, a1, a2, d){

  # quadrature weights and roots
  ghDat <- fastGHQuad::gaussHermiteData(d)
  w1 <- ghDat$w
  r1 <- ghDat$x

  #Get weights from MAF
  MAF <- apply(gMat, 2, function(x) mean(x, na.rm = T)/2)

  #Weight MAF using beta distribution
  beta <- dbeta(MAF, a1, a2, ncp = 0, log = FALSE)

  #Apply weights to matrix of genotypes
  Z_w = t(t(gMat) * (beta))

  # calculate first derivative term
  U1 <- sweep((exp(-HL_array) * -HL_array), c(1,3), ifelse(tpos_all == 0, 0,1), FUN = "*")
  U2 <- sweep((exp(-HR_array) * -HR_array), c(1,3), ifelse(obs_all == 0, 0,1), FUN = "*")
  U2[is.na(U2)] <- 0

  # difference between survival terms multiplied across observations not including k
  no_l <- no_l_all[,,l]

  # mutliply elements together to get the first term of the derivative
  t1 <- t(Z_w) %*% ((U1 - U2)[,,l] * no_l/A_i)

  #multiply quadrature weights and sum across subjects
  deriv <- apply(sweep(t1, 2, w1, FUN = "*"), 1, sum)/sqrt(pi)

  return(deriv)
}


gamma_on <- function(l, HL_array, HR_array, tpos_all, obs_all, temp_beta, A_i, no_l_all, gMat, a1, a2, d){

  # quadrature weights and roots
  ghDat <- fastGHQuad::gaussHermiteData(d)
  w1 <- ghDat$w
  r1 <- ghDat$x

  #Get weights from MAF
  MAF <- apply(gMat, 2, function(x) mean(x, na.rm = T)/2)

  #Weight MAF using beta distribution
  beta <- dbeta(MAF, a1, a2, ncp = 0, log = FALSE)

  #Apply weights to matrix of genotypes
  Z_w = t(t(gMat) * (beta))

  # calculate second derivative terms
  Y1 <- sweep(-HL_array * exp(-HL_array) + (HL_array^2 * exp(-HL_array)), c(1,3), ifelse(tpos_all == 0, 0,1), FUN = "*")
  Y2 <- sweep(-HR_array * exp(-HR_array) + (HR_array^2 * exp(-HR_array)), c(1,3), ifelse(obs_all == 0, 0,1), FUN = "*")
  Y2[is.na(Y2)] <- 0

  # function for the first term
  term1 <- function(d,l){
    (t(Z_w * (Y1 - Y2)[,d,l] *no_l_all[,d,l]/A_i)) %*% Z_w
  }

  # apply to the quadrature notes
  gt_t1 <- apply(sweep(simplify2array(lapply(1:d, term1, l = l)),3, w1, FUN = "*"), c(1,2), sum)/sqrt(pi)

  # calculate the first derivative term
  U1 <- sweep((exp(-HL_array) * -HL_array), c(1,3), ifelse(tpos_all == 0, 0,1), FUN = "*")
  U2 <- sweep((exp(-HR_array) * -HR_array), c(1,3), ifelse(obs_all == 0, 0,1), FUN = "*")
  U2[is.na(U2)] <- 0

  # term multiplying the difference of the survival terms for all observations other than k
  no_l <- no_l_all[,,l]

  # function combining all the terms
  fd_t <- function(l, d){
    out <- sweep(Z_w, 1, (U1 - U2)[,d,l] * no_l[,d]/A_i, FUN = "*")
    return(out)
  }

  # apply to all the quadrature nodes
  to_d <- apply(sweep(simplify2array(lapply(1:d, fd_t, l = l)), 3, w1, FUN = "*"), c(1,2), sum)/sqrt(pi)

  #multiply term to itself
  gt_t2 <- t(t(to_d) %*% to_d)

  return(gt_t1 - gt_t2)

}



# dl/dgamma_k d_gamma_j
# 50 x 50
gamma_off <- function(l, m, HL_array, HR_array, tpos_all, obs_all, temp_beta, A_i, no_l_all, no_two_all, gMat, a1, a2, k, d){
  #Get weights from MAF
  MAF <- apply(gMat, 2, function(x) mean(x, na.rm = T)/2)

  # quadrature weights and roots
  ghDat <- fastGHQuad::gaussHermiteData(d)
  w1 <- ghDat$w
  r1 <- ghDat$x

  #Weight MAF using beta distribution
  beta <- dbeta(MAF, a1, a2, ncp = 0, log = FALSE)

  #Apply weights to matrix of genotypes
  Z_w = t(t(gMat) * (beta))

  # Product of survival terms without two phenotypes

  # first derivative term
  U1 <- sweep((exp(-HL_array) * -HL_array), c(1,3), ifelse(tpos_all == 0, 0,1), FUN = "*")
  U2 <- sweep((exp(-HR_array) * -HR_array), c(1,3), ifelse(obs_all == 0, 0,1), FUN = "*")
  U2[is.na(U2)] <- 0

  # function for the first term combining all the elements
  # different for k = 2 vs more than 2 observations
  term1 <- function(d,l,m){
    if(k == 2){
      out <- (t(Z_w * (U1 - U2)[,d,l] *(U1 - U2)[,d,m]/A_i)) %*% Z_w
    } else{
      # combination of all the outcome indices
      combs <- combn(1:k, 2)

      # order the observation indices
      min_k <- min(l,m)
      max_k <- max(l,m)

      # product of survival differences without observation l and m
      no_l_m <- no_two_all[,,which(combs[1,] == min_k  & combs[2,] == max_k)]
      out <- (t(Z_w * (U1 - U2)[,d,l] *(U1 - U2)[,d,m] * no_l_m[,d]/A_i)) %*% Z_w
    }

    return(out)
  }

  # apply to all quadrature nodes
  gt_t1 <- apply(sweep(simplify2array(lapply(1:d, term1, l = l, m = m)),3, w1, FUN = "*"), c(1,2), sum)/sqrt(pi)


  # function same as first derivative term
  fd_t <- function(l, d){
    out <- sweep(Z_w, 1, (U1 - U2)[,d,l] * no_l_all[,d,l]/A_i, FUN = "*")
    return(out)
  }

  # apply to all quadrature nodes, multiply by weight, sum over all d
  # apply for both observation l and observation m
  to_d_l <- apply(sweep(simplify2array(lapply(1:d, fd_t, l = l)), 3, w1, FUN = "*"), c(1,2), sum)/sqrt(pi)
  to_d_m <- apply(sweep(simplify2array(lapply(1:d, fd_t, l = m)), 3, w1, FUN = "*"), c(1,2), sum)/sqrt(pi)

  # multiply together
  gt_t2 <- t(t(to_d_l) %*% to_d_m)

  return(gt_t1 - gt_t2)
}


gammatheta <- function(l, HL_array, HR_array, tpos_all, obs_all, apply_diffs, temp_beta, A_i, xDats, no_l_all, gMat, a1, a2, d){

  #Get weights from MAF
  MAF <- apply(gMat, 2, function(x) mean(x, na.rm = T)/2)

  #Weight MAF using beta distribution
  beta <- dbeta(MAF, a1, a2, ncp = 0, log = FALSE)

  #Apply weights to matrix of genotypes
  Z_w = t(t(gMat) * (beta))

  # product of differences of survival terms excluding observation l
  no_l <- no_l_all[,,l]

  # design matrices
  phen = xDats[[l]]
  left_dmat <- phen$dmats$left_dmat
  right_dmat <- phen$dmats$right_dmat

  # left and right survival terms
  lt <- phen$lt
  rt <- phen$rt

  # quadrature weights and roots
  ghDat <- fastGHQuad::gaussHermiteData(d)
  w1 <- ghDat$w
  r1 <- ghDat$x

  # calculate the second derivative terms
  Y1 <- sweep(-HL_array * exp(-HL_array) + (HL_array^2 * exp(-HL_array)), c(1,3), ifelse(tpos_all == 0, 0,1), FUN = "*")
  Y2 <- sweep(-HR_array * exp(-HR_array) + (HR_array^2 * exp(-HR_array)), c(1,3), ifelse(obs_all == 0, 0,1), FUN = "*")
  Y2[is.na(Y2)] <- 0

  # function for total second derivative term
  get_sg <- function(l, HL_array, HR_array, phen, d, no_l_all, Z_w){
    left_dmat <- phen$dmats$left_dmat
    right_dmat <- phen$dmats$right_dmat
    lt <- phen$lt
    rt <- phen$rt

    #left term
    sg_t1 <-  left_dmat * (no_l_all[,d,l]/A_i * as.numeric(Y1[,d,l]))

    #right term
    sg_t2 <-  right_dmat * (no_l_all[,d,l]/A_i * as.numeric(Y2[,d,l]))

    # multiply difference of left and right term times weighted genetic matrix
    out <- t(Z_w) %*% (sg_t1 - sg_t2)

    return(out)
  }

  # apply to all quadrature nodes
  derivs <- sweep(simplify2array(lapply(1:d, get_sg, l = l, HL_array = HL_array, HR_array = HR_array, phen = phen, no_l_all = no_l_all, Z_w = Z_w)), 3, w1, FUN = "*")

  # Sum over D
  term1 <- apply(derivs, c(1,2), sum)/sqrt(pi)

  #### now calculate term 2

  # calculate first derivative terms
  U1 <- sweep((exp(-HL_array) * -HL_array), c(1,3), ifelse(tpos_all == 0, 0,1), FUN = "*")
  U2 <- sweep((exp(-HR_array) * -HR_array), c(1,3), ifelse(obs_all == 0, 0,1), FUN = "*")
  U2[is.na(U2)] <- 0

  # calculate the first part of the second term
  term2_a <- function(d,l){
    sweep(Z_w, 1, (U1 - U2)[,d,l] * no_l_all[,d,l]/A_i, FUN = "*")
  }

  # apply to all quadrature nodes, multiply weights, and sum over d
  t2a_tod <- apply(sweep(simplify2array(lapply(1:d, term2_a, l = l)),3,w1,FUN = "*"), c(1,2), sum)/sqrt(pi)

  # function for the third term of the derivative
  term3_func <- function(d,l){
    phen = xDats[[l]]
    left_dmat <- phen$dmats$left_dmat
    right_dmat <- phen$dmats$right_dmat

    # combine left and right terms
    (left_dmat * U1[,d,l] * no_l_all[,d,l]) - (right_dmat * U2[,d,l] * no_l_all[,d,l] )
  }

  # The term that is like the first deriv of the theta term
  t2b_tod <- sweep(apply(sweep(simplify2array(lapply(1:d, term3_func, l = l)), 3, w1, FUN = "*"), c(1,2), sum), 1, A_i, FUN = "/")/sqrt(pi)

  # mulply the two parts together
  term2 <- (t(t2a_tod) %*% t2b_tod)

  return(t(term1 - term2))

}


gammatheta_off <- function(l,m, HL_array, HR_array, xAll, apply_diffs, temp_beta, A_i, no_l_all, no_two_all, gMat, a1, a2, k, d){

  # quadrature weights and roots
  ghDat <- fastGHQuad::gaussHermiteData(d)
  w1 <- ghDat$w
  r1 <- ghDat$x

  # design matrices
  phen = xAll$xDats[[l]]
  left_dmat <- phen$dmats$left_dmat
  right_dmat <- phen$dmats$right_dmat

  # censor terms
  tpos_all <- xAll$ts_all
  obs_all <- xAll$ob_all

  # first derivative term
  U1 <- sweep((exp(-HL_array) * -HL_array), c(1,3), ifelse(tpos_all == 0, 0,1), FUN = "*")
  U2 <- sweep((exp(-HR_array) * -HR_array), c(1,3), ifelse(obs_all == 0, 0,1), FUN = "*")
  U2[is.na(U2)] <- 0

  #Get weights from MAF
  MAF <- apply(gMat, 2, function(x) mean(x, na.rm = T)/2)

  #Weight MAF using beta distribution
  beta <- dbeta(MAF, a1, a2, ncp = 0, log = FALSE)

  #Apply weights to matrix of genotypes
  Z_w = t(t(gMat) * (beta))

  # product of differences of left and right survival excluding observation l
  no_l <- no_l_all[,,l]

  # second derivative term
  # different for k = 2 vs more than k outcomes
  get_sg_off <- function(l, m, HL_array, HR_array, phen, d, Z_w){

    # design matrices and left and right times
    left_dmat <- phen$dmats$left_dmat
    right_dmat <- phen$dmats$right_dmat
    lt <- phen$lt
    rt <- phen$rt

    if(k == 2) {
      # left_term of the first term
      sg_t1a <-  left_dmat * (as.numeric(U1[,d,l])/A_i)

      # right term of the first term
      sg_t1b <-  right_dmat * (as.numeric(U2[,d,l])/A_i)

    } else{
      combs <- combn(1:k, 2)
      min_k <- min(l,m)
      max_k <- max(l,m)
      no_l_m <- no_two_all[,,which(combs[1,] == min_k  & combs[2,] == max_k)]

      #left_term of first term
      sg_t1a <-  left_dmat * (no_l_m[,d]/A_i * as.numeric(U1[,d,l]))

      #right term of first term
      sg_t1b <-  right_dmat * (no_l_m[,d]/A_i * as.numeric(U2[,d,l]))
    }

    # second term
    sg_t2 <- Z_w *(U1[,d,m] - U2[,d,m])

    # combine all the terms
    out <-t(sg_t2)%*% (sg_t1a - sg_t1b)

    return(out)

  }

  # apply first derivative function to all quadrature nodes
  derivs <- simplify2array(lapply(1:d, get_sg_off, l = l, m = m, HL_array = HL_array, HR_array = HR_array, phen = phen, Z_w = Z_w))

  # Sum over D
  term1 <- apply(sweep(derivs, 3, w1, FUN = '*'), c(1,2), sum)/sqrt(pi)

  #### now term 2

  # first derivative terms
  U1 <- sweep((exp(-HL_array) * -HL_array), c(1,3), ifelse(tpos_all == 0, 0,1), FUN = "*")
  U2 <- sweep((exp(-HR_array) * -HR_array), c(1,3), ifelse(obs_all == 0, 0,1), FUN = "*")
  U2[is.na(U2)] <- 0


  # first part of second term
  term2_a <- function(d,l){

    # design matrices
    left_dmat <- phen$dmats$left_dmat
    right_dmat <- phen$dmats$right_dmat

    # combine left and right terms
    sweep((left_dmat * U1[,d,l] - right_dmat*U2[,d,l]), 1, no_l_all[,d,l]/A_i, FUN = "*")
  }

  # apply to all quadrature nodes, sum over d, multiply by weights
  t2a_tod <- apply(sweep(simplify2array(lapply(1:d, term2_a, l = l)),3,w1,FUN = "*"), c(1,2), sum)/sqrt(pi)

  # function for second part of second term
  term2_b <- function(d,m){
    sweep(Z_w, 1, (U1 - U2)[,d,m] * no_l_all[,d,m]/A_i, FUN = "*")
  }

  # apply to all quadrature nodes, sum over d, multiply by weights
  t2b_tod <- apply(sweep(simplify2array(lapply(1:d, term2_b, m = m)),3,w1,FUN = "*"), c(1,2), sum)/sqrt(pi)

  # combine both parts to get second term
  term2 <- t(t(t2a_tod) %*% t2b_tod)

  return(t(term1 - term2))
}


# dl/dgamma_k dsigma
# 50 x 1
gammasigma <- function(l, HL_array, HR_array, tpos_all, obs_all, apply_diffs, temp_beta, A_i, xDats, no_l_all, no_two_all, gMat, a1, a2, k, d) {

  # quadrature weights and roots
  ghDat <- fastGHQuad::gaussHermiteData(d)
  w1 <- ghDat$w
  r1 <- ghDat$x

  #Get weights from MAF
  MAF <- apply(gMat, 2, function(x) mean(x, na.rm = T)/2)

  #Weight MAF using beta distribution
  beta <- dbeta(MAF, a1, a2, ncp = 0, log = FALSE)

  #Apply weights to matrix of genotypes
  Z_w = t(t(gMat) * (beta))

  # for phenotype l
  idx <- (1:k)[-l]

  # design matrices
  phen <- xDats[[l]]
  left_dmat <- phen$dmats$left_dmat
  right_dmat <- phen$dmats$right_dmat

  # get sigma term
  nocol <- ncol(left_dmat)
  sigmasq <- temp_beta[k*nocol + 1]

  # calculate second derivative terms
  Y1 <- sweep(-HL_array * exp(-HL_array) + (HL_array^2 * exp(-HL_array)), c(1,3), ifelse(tpos_all == 0, 0,1), FUN = "*")
  Y2 <- sweep(-HR_array * exp(-HR_array) + (HR_array^2 * exp(-HR_array)), c(1,3), ifelse(obs_all == 0, 0,1), FUN = "*")
  Y2[is.na(Y2)] <- 0

  # calculate first derivative terms
  U1 <- sweep((exp(-HL_array) * -HL_array), c(1,3), ifelse(tpos_all == 0, 0,1), FUN = "*")
  U2 <- sweep((exp(-HR_array) * -HR_array), c(1,3), ifelse(obs_all == 0, 0,1), FUN = "*")
  U2[is.na(U2)] <- 0

  # product of diff of survival without l
  surv_no_l <- no_l_all

  # Function for the second deriv
  sg_sd <- function(d, l) {

    if(k == 2){
      out <- Z_w*( (surv_no_l[,d,l] * Y1[,d,l]) - (surv_no_l[,d,l]* Y2[,d,l]) ) +
        ( Z_w* (U1[,d,l] - U2[,d,l] ) ) * (U1 - U2)[,d,idx]
    } else {
      # Product of survival terms without two phenotypes
      combs <- combn(1:k, 2)
      no_l_m <- no_two_all[,,-which(combs[1,] != l & combs[2,] != l)]
      out <- Z_w*( (surv_no_l[,d,l] * Y1[,d,l]) - (surv_no_l[,d,l]* Y2[,d,l]) ) +
        ( Z_w* (U1[,d,l] - U2[,d,l] ) ) * apply((U1 - U2)[,d,idx]* no_l_m[,d,], 1, sum)
    }
    return(out)
  }

  # Apply to all the node points
  to_d <- apply(sweep(simplify2array(lapply(1:d, sg_sd, l = l)), 3, w1 * (r1/sqrt(2*sigmasq)), FUN = "*"), c(1,2), sum)

  # sum over d to get first term
  st_t1 <- colSums(sweep(to_d, 1, A_i, FUN = "/"))/sqrt(pi)

  # function for first part of second term
  term2_a <- function(d,l){
    sweep(Z_w, 1, (U1 - U2)[,d,l] * no_l_all[,d,l]/A_i, FUN = "*")
  }

  # apply to all quadrature nodes and sum
  t2a <- apply(sweep(simplify2array(lapply(1:d, term2_a, l = l)),3,w1,FUN = "*"), c(1,2), sum)/sqrt(pi)

  # The term that is like the first deriv of the sigma term
  t2b <- rowSums(sweep(apply((U1 - U2) * surv_no_l/A_i, c(1,2), sum),2, w1*r1/sqrt(2*sigmasq), FUN = "*"))/sqrt(pi)

  # Sum over n after multiplying both terms
  st_t2 <- colSums(t2a *t2b)

  # combine the two terms
  deriv <- st_t1 - st_t2

  return(matrix(deriv, nrow = 1))

}

