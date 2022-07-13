# newton raphson, get p-value, generate data

#' Newton Raphson Function
#'
#' @param init_beta starting values for NR
#' @param epsilon stopping criterion for NR
#' @param xDats list of design matrices
#' @param lt_all n x k matrix of left times
#' @param rt_all n x k matrix of right times
#' @param k number of outcomes
#' @param d number of quadrature points
#' @export
#' fit_null_general()
fit_null_general <- function(init_beta, epsilon, xDats, lt_all, rt_all, k, d) {

  # number of observations
  n = nrow(xDats[[1]]$dmats$right_dmat)

  # number of covariates
  nocol <- ncol(xDats[[1]]$dmats$right_dmat)

  # create matrix of censoring
  tpos_all <- matrix(NA, nrow = n, ncol = k)
  obs_all <- matrix(NA, nrow = n, ncol = k)

  # get censoring terms
  for(j in 1:k){
    tpos_all[,j] <- as.numeric(lt_all[,j] > 0)
    obs_all[,j] <- as.numeric(rt_all[,j] != Inf)
  }

  # make data in correct form
  threedmat <- list()
  for(i in 1:k){
    threedmat[[i]] <- list(dmats = xDats[[i]]$dmats, lt = lt_all[,i], rt = rt_all[,i])
  }

  # compile all elements
  xAll <- list(xDats = threedmat, ts_all = tpos_all, ob_all = obs_all)

  # get values for format of function input
  xDats <- xAll$xDats
  t_all <- xAll$ts_all
  o_all <- xAll$ob_all

  # quadrature weights and roots
  ghDat <- fastGHQuad::gaussHermiteData(d)
  w1 <- ghDat$w
  r1 <- ghDat$x

  # conditions for convergence
  iter = 1
  diff = 1

  # inital beta values
  temp_beta <- matrix(init_beta, nrow = 1)

  # loop for newton raphson
  while(diff > epsilon & iter < 200){

    # make left and right survival and hazard values
    HL_array <- array (c (NA, NA), dim=c (n,d,k))
    HR_array <- array (c (NA, NA), dim=c (n,d,k))

    # calculate the values of the array
    for(i in 1:k){
      hl_d <- lapply(1:d, haz_left, l = i, temp_beta = temp_beta, phen = xDats[[i]], r1 = r1, k =  k)
      hr_d <- lapply(1:d, haz_right, l = i, temp_beta = temp_beta, phen = xDats[[i]], r1 = r1, k = k)
      HL_array[,,i] <- simplify2array(hl_d)
      HR_array[,,i] <- simplify2array(hr_d)
    }

    # calculate the differences between survival values
    apply_diffs <- array (c (NA, NA), dim=c (n,k,d))

    for(i in 1:k){

      # subset data for observation number
      phen = xDats[[i]]

      # cacluate the values and fill in the array
      apply_diffs[,i,] <- simplify2array(lapply(1:d,surv_diff, l = i, temp_beta = temp_beta, phen= phen, r1 = r1, k = k))
    }

    # get denominator of values
    A_i <- get_A(apply_diffs, w1, 100, n)

    # make sure there are no zeroes
    A_i[which(A_i == 0)]<- min(A_i[which(A_i!= 0)])

    # calculate the product of the difference between the survival terms excluding outcome l
    no_l_all <- simplify2array(lapply(1:k, without_one_phen, k = k, store=apply_diffs))

    # get the combination of all outcome indices
    combs <- combn(1:k, 2)

    # calculate the product of the difference between the survival terms excluding two outcomes
    if(k == 2){
      no_two_all <- 1
    } else {
      no_two_all <- array(data = NA, dim = c(n,d, choose(k,2)))
      for(i in 1:choose(k,2)){
        no_two_all[,,i] <- without_two_phen(combs[1,i], combs[2,i], k, apply_diffs, n, d)

      }
    }

    # generalizable way to gt the gradient
    grad <- c()

    # loop through all outcomes to get the full gradient
    for(i in 1:k){
      temp_grad <- fd_term(i, temp_beta, xDats[[i]], d, apply_diffs = apply_diffs, A_i =A_i, no_l_all = no_l_all, HL_array = HL_array, HR_array = HR_array)

      # combine all the gradients
      grad <- c(grad, temp_grad)

    }

    # calculate the gradient of the sigma squared term
    grad_ss <- ss_fd(1, xDats[[1]], HL_array, HR_array, t_all, o_all, apply_diffs = apply_diffs, temp_beta = temp_beta, A_i =A_i, no_l_all = no_l_all, k = k, d = d)

    #combine all the gradient terms
    grad <- matrix(c(grad, grad_ss), ncol = 1)

    # total number of rows and columns for the information matrix
    totd <- (nocol*k) + 1

    # start building the information matrix
    jmat <- matrix(NA, nrow = totd, ncol = totd)

    # start filling in the matrix
    for(i in 1:k){

      # get the correct indices
      d1 <- (nocol*(i-1))+ 1
      d2 <- nocol* i

      # calculate the second derivatives
      jmat[d1:d2, d1:d2] <- sd_on(i, k, temp_beta, xDats[[i]], d, apply_diffs = apply_diffs, A_i =A_i, no_l_all = no_l_all, HL_array = HL_array, HR_array = HR_array)
      jmat[d1:d2, totd] <- st_off(i, HL_array, HR_array, xAll, apply_diffs = apply_diffs, temp_beta = temp_beta, A_i =A_i, no_l_all = no_l_all, no_two_all = no_two_all, k = k, d = d)
      jmat[totd, d1:d2] <- t(st_off(i, HL_array, HR_array, xAll, apply_diffs = apply_diffs, temp_beta = temp_beta, A_i =A_i, no_l_all = no_l_all, no_two_all = no_two_all, k = k, d = d))

      # get the off diagonal terms
      idx <- (1:k)[-i]

      # fill in the matrix
      for(j in 1:length(idx))
      {
        od1 <- (nocol*(idx[j]-1))+ 1
        od2 <- nocol* idx[j]

        # second derivative terms
        jmat[d1:d2, od1:od2] <- sd_off(i,idx[j], phen_l = xDats[[i]], phen_m = xDats[[idx[j]]], temp_beta, d = d, apply_diffs = apply_diffs, A_i =A_i, HL_array = HL_array, HR_array = HR_array, no_l_all = no_l_all, no_two_all = no_two_all, tpos_all = t_all, obs_all = o_all, k = k)
        jmat[od1:od2, d1:d2] <- sd_off(idx[j],i, phen_l = xDats[[idx[j]]], phen_m = xDats[[i]], temp_beta, d = d, apply_diffs = apply_diffs, A_i =A_i, HL_array = HL_array, HR_array = HR_array, no_l_all = no_l_all, no_two_all = no_two_all, tpos_all = t_all, obs_all = o_all, k = k)

      }
    }

    # second derivative of the sigma term
    jmat[totd, totd] <- ss_sd(HL_array, HR_array, xAll, apply_diffs = apply_diffs, temp_beta = temp_beta, A_i =A_i, no_l_all = no_l_all, no_two_all = no_two_all, k = k, d = d)

    # newton raphson
    beta_new <- temp_beta - t(grad) %*% solve(jmat)

    # calculate the difference
    diff = (-t(grad) %*% solve(jmat)) %*% t(-t(grad) %*% solve(jmat))

    # update the iteration
    iter = iter + 1

    # update the new beta
    temp_beta <- matrix(beta_new, nrow = 1)

    print(iter)
    print(diff)
  }

  return(list(temp_beta = temp_beta, iter = iter, diff = diff, jmat = jmat, grad = grad))

}


############### get p value
#' test statistic and pvalue calculation
#'
#' @param nullFit beta values fitted from NR
#' @param xDats list of design matrices
#' @param lt_all n x k matrix of left times
#' @param rt_all n x k matrix of right times
#' @param Itt information of theta terms from NR
#' @param a1 first beta distribution shape parameter for weights
#' @param a2 second beta distribution shape parameter for weights
#' @param G n x q matrix of genetic information
#' @param k number of outcomes
#' @param d number of quadrature points
#' @export
#' multiICSKAT_p_general()
multiICSKAT_p_general <- function(nullFit, xDats, lt_all, rt_all, Itt, a1, a2, G, k, d){

  # quadrature weights and roots
  ghDat <- fastGHQuad::gaussHermiteData(d)
  w1 <- ghDat$w
  r1 <- ghDat$x

  # number of observations
  n = nrow(xDats[[1]]$dmats$right_dmat)

  # number of covariates
  nocol <- ncol(xDats[[1]]$dmats$right_dmat)

  # create matrix of censoring
  tpos_all <- matrix(NA, nrow = n, ncol = k)
  obs_all <- matrix(NA, nrow = n, ncol = k)

  # get censoring terms
  for(j in 1:k){
    tpos_all[,j] <- as.numeric(lt_all[,j] > 0)
    obs_all[,j] <- as.numeric(rt_all[,j] != Inf)
  }

  # make data in correct form
  threedmat <- list()
  for(i in 1:k){
    threedmat[[i]] <- list(dmats = xDats[[i]]$dmats, lt = lt_all[,i], rt = rt_all[,i])
  }

  # compile all elements
  xAll <- list(xDats = threedmat, ts_all = tpos_all, ob_all = obs_all)

  # get values for format of function input
  xDats <- xAll$xDats
  tpos_all <- xAll$ts_all
  obs_all <- xAll$ob_all

  # get null fit from NR
  temp_beta <- nullFit

  # make left and right survival and hazard values
  HL_array <- array (c (NA, NA), dim=c (n,d,k))
  HR_array <- array (c (NA, NA), dim=c (n,d,k))

  # calculate the values of the array
  for(i in 1:k){
    hl_d <- lapply(1:d, haz_left, l = i, temp_beta = temp_beta, phen = xDats[[i]], r1 = r1, k =  k)
    hr_d <- lapply(1:d, haz_right, l = i, temp_beta = temp_beta, phen = xDats[[i]], r1 = r1, k = k)
    HL_array[,,i] <- simplify2array(hl_d)
    HR_array[,,i] <- simplify2array(hr_d)
  }

  # calculate the differences between survival values
  apply_diffs <- array (c (NA, NA), dim=c (n,k,d))

  for(i in 1:k){

    # subset data for observation number
    phen = xDats[[i]]

    # cacluate the values and fill in the array
    apply_diffs[,i,] <- simplify2array(lapply(1:d,surv_diff, l = i, temp_beta = temp_beta, phen= phen, r1 = r1, k = k))
  }

  # get denominator of values
  A_i <- get_A(apply_diffs, w1, 100, n)

  # make sure there are no zeroes
  A_i[which(A_i == 0)]<- min(A_i[which(A_i!= 0)])

  # calculate the product of the difference between the survival terms excluding outcome l
  no_l_all <- simplify2array(lapply(1:k, without_one_phen, k = k, store=apply_diffs))

  # get the combination of all outcome indices
  combs <- combn(1:k, 2)

  # calculate the product of the difference between the survival terms excluding two outcomes
  if(k == 2){
    no_two_all <- 1
  } else {
    no_two_all <- array(data = NA, dim = c(n,d, choose(k,2)))
    for(i in 1:choose(k,2)){
      no_two_all[,,i] <- without_two_phen(combs[1,i], combs[2,i], k, apply_diffs, n, d)

    }
  }

  # build the test statistic

  # get correct dimension values
  dmatdim <- ncol(xDats[[1]]$dmats$right_dmat)
  gdim <- ncol(G)

  # information matrix of gamma and theta terms
  Igt <- matrix(NA, nrow = (dmatdim * k + 1),ncol = k*gdim)

  # loop through observations and run derivative functions
  for(i in 1:k){

    # gamma theta
    Igt[(((i - 1)*dmatdim) + 1): (i*dmatdim), (((i - 1)*gdim) + 1): (i*gdim)] <- gammatheta(i, HL_array, HR_array, tpos_all, obs_all, apply_diffs, temp_beta, A_i, xDats, no_l_all, G, a1, a2, d)
    Igt[(dmatdim * k + 1), (((i - 1)*gdim) + 1): (i*gdim)] <- gammasigma(i, HL_array, HR_array, tpos_all, obs_all, apply_diffs, temp_beta, A_i, xDats, no_l_all, no_two_all, G, a1, a2, k, d)

    # get indices for off diagonal terms
    idx <- (1:k)[-i]

    # loop through index
    for(j in 1:length(idx)){

      Igt[(((idx[j] - 1)*dmatdim) + 1): (idx[j]*dmatdim), (((i - 1)*gdim) + 1): (i*gdim)] <- gammatheta_off(idx[j],i, HL_array, HR_array, xAll, apply_diffs, temp_beta, A_i, no_l_all, no_two_all, G, a1, a2, k, d)
      Igt[(dmatdim * k + 1), (((idx[j] - 1)*gdim) + 1): (idx[j]*gdim)] <- gammasigma(idx[j], HL_array, HR_array, tpos_all, obs_all, apply_diffs, temp_beta, A_i, xDats, no_l_all, no_two_all, G, a1, a2, k, d)

    }
  }

  # information matrix of gamma gamma
  Igg <- matrix(NA, nrow = k * gdim, ncol = k*gdim)

  for(i in 1:k){

    # get correct dimension values
    d1 <- (gdim*(i-1))+ 1
    d2 <- gdim* i

    # second derivative of gamma term on the diagonal
    Igg[d1:d2, d1:d2] <- gamma_on(i, HL_array, HR_array, tpos_all, obs_all, temp_beta, A_i, no_l_all, G, a1, a2, d)

    # get indices for off diagonals
    idx <- (1:k)[-i]

    # loop through indices
    for(j in 1:length(idx))
    {
      od1 <- (ncol(G)*(idx[j]-1))+ 1
      od2 <- ncol(G)* idx[j]

      # second derivative for gamma terms on the off diagonals
      Igg[d1:d2, od1:od2] <- gamma_off(i, idx[j], HL_array, HR_array, tpos_all, obs_all, temp_beta, A_i, no_l_all, no_two_all, G, a1, a2, k, d)
      Igg[od1:od2, d1:d2] <- gamma_off(idx[j], i, HL_array, HR_array, tpos_all, obs_all, temp_beta, A_i, no_l_all, no_two_all, G, a1, a2, k, d)
    }
  }

  # get the full information matrix
  sigmat <- (-Igg) - t(-Igt) %*% solve(-Itt) %*% (-Igt)

  # get the gradient term
  U_g <- c()

  # loop through all the observation values
  for(i in 1:k){

    temp_grad <- gamma_fd(i, HL_array, HR_array, tpos_all, obs_all, temp_beta, A_i, no_l_all, G, a1, a2, d)
    U_g <- c(U_g, temp_grad)
  }

  # compute the score statistic
  gamma_score <- t(U_g) %*% U_g

  # get the eigen values
  lams <- eigen(sigmat)$values

  # compute the p-value using davies method
  pval <- CompQuadForm::davies(q=gamma_score, lambda=lams)


  return(list(Q = gamma_score, pval = pval$Qq))
}


################################################# F6
#' Sample generation
#'
#' @param bhFunInv A function, the inverse of the baseline hazard function.
#' @param obsTimes Vector of the intended observation times.
#' @param windowHalf The amount of time before or after the intended obsTimes that a visit might take place.
#' @param n number of observations
#' @param k number of outcomes
#' @param tausq variance of subject specific random effect
#' @param gMatCausal matrix of subsetted genetic information for only a select causal SNPs
#' @param effectSizes vector of genetic effects
#' @export
#' gen_mICSKAT_dat()
gen_mICSKAT_dat <- function(bhFunInv, obsTimes = 1:3, windowHalf = 0.5, n, k, tauSq, gMatCausal, effectSizes) {

  # true model has nothing
  fixedMat <- matrix(data=0, nrow=n, ncol=K)

  # get genetic effects
  geneticVec <- c()

  # calculate the effect size
  for (outcomes in 1:k)
  {
    # create vector for each SNP
    Bk <- rep(NA, ncol(gMatCausal))

    # loop through all SNPs
    for(j in 1:length(Bk)){

      # calculate minor allele frequency
      MAF <- apply(gMatCausal, 2, function(x) mean(x)/2)

      # multply effect size and genetic matrix
      Bk[j] = effectSizes[outcomes]* abs(log10(MAF[j]))
    }
    geneticVec <- c(geneticVec, (gMatCausal %*% Bk))
  }

  # random unif vector for PIT
  unifVec <- runif(n=n * k)

  # vectorize fixed effects, all first phenotype, then all second phenotype, etc.
  fixedVec <- c(fixedMat)

  # random intercept
  randomInt <- rnorm(n=n, sd = sqrt(tauSq))

  # get full term
  randomIntRep <- rep(randomInt, k)

  # add random effect
  etaVec <- fixedVec + randomIntRep + geneticVec

  # probability integral transform - assumes PH model
  toBeInv <- -log(1 - unifVec) / exp(etaVec)

  # all n*K exact failure times
  exactTimesVec <- bhFunInv(toBeInv)

  # all exact failure times for all K phenotypes
  exactTimesMat <- matrix(data=exactTimesVec, nrow=n, ncol=k, byrow = FALSE)

  # hold left and right intervals data for all K phenotypes
  leftTimesMat <- matrix(data=NA, nrow=n, ncol=k)
  rightTimesMat <- matrix(data=NA, nrow=n, ncol=k)
  obsInd <- matrix(data=NA, nrow=n, ncol=k)
  tposInd <- matrix(data=NA, nrow=n, ncol=k)

  # do visits separately for each phenotype
  nVisits <- length(obsTimes)

  for (pheno_it in 1:k) {

    # your visit is uniformly distributed around the intended obsTime, windowHalf on each side
    allVisits <- sweep(matrix(data=runif(n=n*nVisits, min=-windowHalf, max=windowHalf), nrow=n, ncol=nVisits),
                       MARGIN=2, STATS=obsTimes, FUN="+")

    # make the interval for each subject
    allInts <- t(mapply(FUN=createInt, obsTimes = data.frame(t(allVisits)), eventTime=exactTimesMat[, pheno_it]))

    leftTimesMat[, pheno_it] <- allInts[, 1]
    rightTimesMat[, pheno_it] <- allInts[, 2]

    # event time indicators
    obsInd[, pheno_it] <- ifelse(rightTimesMat[, pheno_it] == Inf, 0, 1)
    tposInd[, pheno_it] <- ifelse(leftTimesMat[, pheno_it] == 0, 0, 1)
  }

  # return
  return(list(exactTimesMat = exactTimesMat, leftTimesMat = leftTimesMat,
              rightTimesMat = rightTimesMat, obsInd = obsInd, tposInd = tposInd))

}

# function needed to generate data
createInt <- function(obsTimes, eventTime) {

  # order the times in case the random portion causes them to go out of order
  orderedTimes <- sort(obsTimes)

  # left end of interval
  minIdx <- which(orderedTimes < eventTime)

  if (length(minIdx) == 0) {
    minTime <- 0
  } else {
    minTime <- orderedTimes[max(minIdx)]
  }

  # right end of interval
  maxIdx <- which(orderedTimes >= eventTime)

  if (length(maxIdx) == 0) {
    maxTime <- Inf
  } else {
    maxTime <- orderedTimes[min(maxIdx)]

  }
  return(c(minTime, maxTime))

}
