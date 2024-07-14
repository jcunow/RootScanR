#### sourced from: https://rdrr.io/github/diffCircadian/diffCircadian/
## author: Zhiguang Huo, Haocheng Ding
# Feb. 24, 2023, 9:07 a.m.; Version 0.0.0
# install.packages("remotes")
# remotes::install_github("diffCircadian/diffCircadian")




#' Fit sin function
##'
##' Fit a sine curve where tt is time, and yy is expression value.
##' @title Fit Data Based on Sine Curve
##' @param tt Time vector.
##' @param yy Expression vector.
##' @param period Period of the sine curve. Default is 24.
##' @param parStart Initial value for optimization purpose.
##' @return A list of amp, phase, offset, peak, A, B, SST, SSE, R2.
##' Formula 1: \eqn{yy = amp * sin(2\pi/period * (phase + tt)) + offset}
##' Formula 2: \eqn{yy = A * sin(2\pi/period \times tt) + B * cos(2*\pi/period * tt) + offset}
##' \item{amp}{Amplitude based on formula 1.}
##' \item{phase}{Phase based on formula 1, phase is restricted within (0, period).}
##' \item{offset}{Basal level (vertical shift) based on formula 1 or on formula 2.}
##' \item{A}{A based on formula 2.}
##' \item{B}{B based on formula 2.}
##' \item{tss}{Total sum of square.}
##' \item{rss}{Residual sum of square, SSE/n is the MLE of the variance sigma square.}
##' \item{R2}{Pseudo R2 defined as (tss - rss)/tss.}
##' @author Caleb
##' @import minpack.lm
##' @export
##' @examples
##' set.seed(32608)
##' n <- 50
##' tt <- runif(n,0,24)
##' Amp <- 2
##' Phase <- 6
##' Offset <- 3
##' yy <- Amp * sin(2*pi/24 * (tt + Phase)) + Offset + rnorm(n,0,1)
##' fitSinCurve(tt, yy)

fitSinCurve <- function(tt, yy, period = 24, parStart = list(amp=3,phase=0, offset=0)){

	getPred <- function(parS, tt) {
		parS$amp * sin(2*pi/period * (tt + parS$phase)) + parS$offset
	}

	residFun <- function(p, yy, tt) yy - getPred(p,tt)

	nls.out <- minpack.lm::nls.lm(par=parStart, fn = residFun, yy = yy,	tt = tt)

	apar <- nls.out$par

	amp0 <- apar$amp
	asign <- sign(amp0)
	## restrict amp > 0
	amp <- amp0 * asign

	phase0 <- apar$phase
	#phase <- (round(apar$phase) + ifelse(asign==1,0,12)) %% period
	phase <- (phase0 + ifelse(asign==1,0,period/2)) %% period
	offset <- apar$offset

	peak <- (period/2 * sign(amp0) - period/4 - phase) %%period
	if(peak > period/4*3) peak = peak - period

	A <- amp0 * cos(2*pi/period * phase0)
	B <- amp0 * sin(2*pi/period * phase0)

	rss <- sum(nls.out$fvec^2)
	tss <- sum((yy - mean(yy))^2)
	R2 <- 1 - rss/tss

	if(F){
		amp <- apar$amp
		phase <- apar$phase
		offset <- apar$offset
	}

	res <- list(amp=amp, phase=phase, offset=offset, peak=peak, A=A, B=B, tss=tss, rss=rss, R2=R2)
	res
}



##' Fisher information matrix when two conditions exist
##'
##' Obtain the Fisher information matrix when two conditions exist
##' @title Fisher information matrix when two conditions exist
##' @param beta parameter vector of 8 with the following order: amp_1, phase_1, offset_1, theta_1, amp_2, phase_2, offset_2, theta_2
##' @param tt1 time vector of condition 1
##' @param yy1 expression vector of condition 1
##' @param tt2 time vector of condition 2
##' @param yy2 expression vector of condition 2
##' @param period Period of the since curve. Default is 24.
##' @return The Fisher information matrix, this is a 8*8 matrix, with the same order as the input beta parameter.
##' @author Caleb
##' @noRd
##' @examples
##' set.seed(32608)
##' n <- 50
##' tt1 <- runif(n,0,24)
##' Amp1 <- 2
##' Phase1 <- 6
##' Offset1 <- 3
##' yy1 <- Amp1 * sin(2*pi/24 * (tt1 + Phase1)) + Offset1 + rnorm(n,0,1)
##' tt2 <- runif(n,0,24)
##' Amp2 <- 3
##' Phase2 <- 15
##' Offset2 <- 2
##' yy2 <- Amp2 * sin(2*pi/24 * (tt2 + Phase2)) + Offset2 + rnorm(n,0,1)
##' beta <- c(Amp1,Phase1,Offset1,1,Amp2,Phase2,Offset2,2)
##' fisherInformation2(beta, tt1, yy1, tt2, yy2)

fisherInformation2 <- function(beta, tt1, yy1, tt2, yy2, period = 24){

  n1 <- length(tt1)
  stopifnot(n1 == length(yy1))
  n2 <- length(tt2)
  stopifnot(length(tt2) == length(yy2))

  stopifnot(length(beta) == 8)

  w <- 2*pi/period

  amp_1 <- beta[1]
  phase_1 <- beta[2]
  offset_1 <- beta[3]
  theta_1 <- beta[4]
  amp_2 <- beta[5]
  phase_2 <- beta[6]
  offset_2 <- beta[7]
  theta_2 <- beta[8]

  asin_1 <- sin(w * (tt1 + phase_1))
  asin_2 <- sin(w * (tt2 + phase_2))

  acos_1 <- cos(w * (tt1 + phase_1))
  acos_2 <- cos(w * (tt2 + phase_2))

  yhat_1 <- amp_1 * asin_1 + offset_1
  yhat_2 <- amp_2 * asin_2 + offset_2

  diffy_1 <- yy1-yhat_1
  diffy_2 <- yy2-yhat_2

  h_amp1_amp1 <- theta_1 * sum(asin_1^2)
  h_amp1_phase1 <- theta_1 * w * sum((amp_1 * asin_1 - diffy_1) * acos_1)
  h_amp1_offset1 <- theta_1 * sum(asin_1)
  h_amp1_theta1 <- - sum(diffy_1 * asin_1)

  h_phase1_phase1 <- theta_1 * amp_1 * w^2 * sum(amp_1 * acos_1^2 + diffy_1 * asin_1)
  h_phase1_offset1 <- theta_1 * amp_1 * w * sum(acos_1)
  h_phase1_theta1 <- - amp_1 * w * sum(diffy_1 * acos_1)

  h_offset1_offset1 <- theta_1 * n1
  h_offset1_theta1 <- - sum(diffy_1)

  h_theta1_theta1 <- n1 / 2 / theta_1^2

  h_amp2_amp2 <- theta_2 * sum(asin_2^2)
  h_amp2_phase2 <- theta_2 * w * sum((amp_2 * asin_2 - diffy_2) * acos_2)
  h_amp2_offset2 <- theta_2 * sum(asin_2)
  h_amp2_theta2 <- - sum(diffy_2 * asin_2)

  h_phase2_phase2 <- theta_2 * amp_2 * w^2 * sum(amp_2 * acos_2^2 + diffy_2 * asin_2)
  h_phase2_offset2 <- theta_2 * amp_2 * w * sum(acos_2)
  h_phase2_theta2 <- - amp_2 * w * sum(diffy_2 * acos_2)

  h_offset2_offset2 <- theta_2 * n2
  h_offset2_theta2 <- - sum(diffy_2)

  h_theta2_theta2 <- n2 / 2 / theta_2^2

  hmatrix <- matrix(0,nrow=8,ncol = 8)

  hmatrix[1,1] <- h_amp1_amp1
  hmatrix[1,2] <- hmatrix[2,1] <- h_amp1_phase1
  hmatrix[1,3] <- hmatrix[3,1] <- h_amp1_offset1
  hmatrix[1,4] <- hmatrix[4,1] <- h_amp1_theta1
  hmatrix[2,2] <- h_phase1_phase1
  hmatrix[2,3] <- hmatrix[3,2] <- h_phase1_offset1
  hmatrix[2,4] <- hmatrix[4,2] <- h_phase1_theta1
  hmatrix[3,3] <- h_offset1_offset1
  hmatrix[3,4] <- hmatrix[4,3] <- h_offset1_theta1
  hmatrix[4,4] <- h_theta1_theta1

  hmatrix[5,5] <- h_amp2_amp2
  hmatrix[5,6] <- hmatrix[6,5] <- h_amp2_phase2
  hmatrix[5,7] <- hmatrix[7,5] <- h_amp2_offset2
  hmatrix[5,8] <- hmatrix[8,5] <- h_amp2_theta2
  hmatrix[6,6] <- h_phase2_phase2
  hmatrix[6,7] <- hmatrix[7,6] <- h_phase2_offset2
  hmatrix[6,8] <- hmatrix[8,6] <- h_phase2_theta2
  hmatrix[7,7] <- h_offset2_offset2
  hmatrix[7,8] <- hmatrix[8,7] <- h_offset2_theta2
  hmatrix[8,8] <- h_theta2_theta2

  hmatrix
}






#' Finite sample/Large sample Wald test for differential amplitude.
##'
##' Test differential amplitude of circadian curve fitting using Wald test
##' @title Wald test for detecting differential amplitude
##' @param tt1 Time vector of condition 1
##' @param yy1 Expression vector of condition 1
##' @param tt2 Time vector of condition 2
##' @param yy2 Expression vector of condition 2
##' @param period Period of the since curve. Default is 24.
##' @param type Type of Wald test to use, "FN" or "LS". Default is finite sample.
##' @return A list, see details below.
##' Formula 1: \eqn{yy = amp \times sin(2\pi/period \times (phase + tt)) + offset}
##' Formula 2: \eqn{yy = A \times sin(2\pi/period \times tt) + B * cos(2*pi/period * tt) + offset}
##' \item{amp_1}{Amplitude estimate of the 1st data}
##' \item{amp_2}{Amplitude estimate of the 2nd data}
##' \item{amp_c}{Amplitude estimate pooling all data together}
##' \item{df}{Degree of freedom for the Wald test}
##' \item{stat}{Wald statistics}
##' \item{pvalue}{P-value from the Wald test}
##' @author Caleb
##' @noRd
##' @examples
##' set.seed(32608)
##' n <- 50
##' tt1 <- runif(n,0,24)
##' Amp1 <- 2
##' Phase1 <- 6
##' Offset1 <- 3
##' yy1 <- Amp1 * sin(2*pi/24 * (tt1 + Phase1)) + Offset1 + rnorm(n,0,1)
##' tt2 <- runif(n,0,24)
##' Amp2 <- 3
##' Phase2 <- 5
##' Offset2 <- 2
##' yy2 <- Amp2 * sin(2*pi/24 * (tt2 + Phase2)) + Offset2 + rnorm(n,0,1)
##' WaldTest_diff_amp(tt1, yy1, tt2, yy2)


WaldTest_diff_amp <- function(tt1, yy1, tt2, yy2, period = 24,FN=TRUE){
  n1 <- length(tt1)
  stopifnot(n1 == length(yy1))
  n2 <- length(tt2)
  stopifnot(length(tt2) == length(yy2))

  #period <- 24
  w <- 2*pi/period

  fit1 <- fitSinCurve(tt1, yy1, period = period)
  fit2 <- fitSinCurve(tt2, yy2, period = period)

  A1 <- fit1$amp
  A2 <- fit2$amp

  phase1 <- fit1$phase
  phase2 <- fit2$phase

  basal1 <- fit1$offset
  basal2 <- fit2$offset

  sigma2_1 <- 1/n1 * fit1$rss
  sigma2_2 <- 1/n2 * fit2$rss

  theta1 <- 1/sigma2_1
  theta2 <- 1/sigma2_2

  p1 <- c(A1, phase1, basal1, theta1)
  p2 <- c(A2, phase2, basal2, theta2)

  x_Ha <- c(p1, p2)

  eval_f_list <- function(x) {
    p1 <- x[1:4]
    p2 <- x[5:8]

    A1 <- p1[1]
    phase1 <- p1[2]
    basel1 <- p1[3]
    theta1 <- p1[4]
    asin1 <- sin(w * (tt1 + phase1) )
    acos1 <- cos(w * (tt1 + phase1) )
    yhat1 <- A1 * asin1 + basel1

    A2 <- p2[1]
    phase2 <- p2[2]
    basel2 <- p2[3]
    theta2 <- p2[4]
    asin2 <- sin(w * (tt2 + phase2) )
    acos2 <- cos(w * (tt2 + phase2) )
    yhat2 <- A2 * asin2 + basel2

    ll1_a <- log(theta1)/2
    ll1_b <- (yy1 - yhat1)^2 * theta1 / 2
    ll1 <- ll1_a - ll1_b

    ll2_a <- log(theta2)/2
    ll2_b <- (yy2 - yhat2)^2 * theta2 / 2
    ll2 <- ll2_a - ll2_b

    partial_A1 <- - theta1 * sum((yy1 - yhat1) * asin1)
    partial_phase1 <- - theta1 * A1 * w * sum((yy1 - yhat1) * acos1)
    partial_C1 <- - theta1 * sum(yy1 - yhat1)
    partial_theta1 <-  sum((yy1 - yhat1)^2)/2 - n1/2/theta1

    partial_A2 <- - theta2 * sum((yy2 - yhat2) * asin2)
    partial_phase2 <- - theta2 * A2 * w * sum((yy2 - yhat2) * acos2)
    partial_C2 <- - theta2 * sum(yy2 - yhat2)
    partial_theta2 <-  sum((yy2 - yhat2)^2)/2 - n2/2/theta2


    return( list( "objective" = -sum(ll1) - sum(ll2),
                  "gradient"  = c(partial_A1, partial_phase1, partial_C1, partial_theta1,
                                  partial_A2, partial_phase2, partial_C2, partial_theta2)
    )
    )
  }

  # Equality constraints
  eval_g_eq <- function(x)
  {
    p1 <- x[1:4]
    p2 <- x[5:8]

    A1 <- p1[1]
    A2 <- p2[1]

    A1 - A2
  }

  # Equality constraints
  eval_g_eq_jac <- function(x)
  {
    c(1, 0, 0, 0,
      -1, 0, 0, 0)
  }


  # Lower and upper bounds
  lb <- c(0,-Inf,-Inf,0, 0, -Inf,-Inf, 0)
  ub <- c(Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf)
  #initial values

  ## Error in is.nloptr(ret) :
  #  If you want to use equality constraints, then you should use one of these algorithms NLOPT_LD_AUGLAG, NLOPT_LN_AUGLAG, NLOPT_LD_AUGLAG_EQ, NLOPT_LN_AUGLAG_EQ, NLOPT_GN_ISRES, NLOPT_LD_SLSQP

  local_opts <- list( "algorithm" = "NLOPT_LD_MMA", "xtol_rel" = 1.0e-15 )
  "local_opts" = local_opts
  opts <- list( "algorithm"= "NLOPT_LD_SLSQP",
                "xtol_rel"= 1.0e-15,
                "maxeval"= 160000,
                "local_opts" = local_opts,
                "print_level" = 0
                #"check_derivatives"=TRUE
  )

  res <- nloptr::nloptr( x0 = x_Ha,
                  eval_f = eval_f_list,
                  #eval_grad_f=eval_g,
                  lb = lb,
                  ub = ub,
                  #eval_g_ineq = eval_g_ineq,
                  eval_g_eq = eval_g_eq,
                  eval_jac_g_eq = eval_g_eq_jac,
                  opts = opts)

  #
  #x_Ha
  x_H0 <- res$solution

  beta_ha <- x_Ha
  beta_h0 <- x_H0
  #
  I8 <- fisherInformation2(beta_ha, tt1, yy1, tt2, yy2, period=period)
  beta_diff <- matrix(beta_ha - beta_h0,ncol=1)
  #stat <- t(beta_diff) %*% I8 %*% beta_diff

  beta_diff2 <- beta_diff[c(1,5)]
  I2 <- solve(solve(I8)[c(1,5),c(1,5)])
  stat <- as.numeric(t(beta_diff2) %*% I2 %*% beta_diff2)

  dfdiff <- 1

  if(FN==FALSE){
    pvalue <- pchisq(stat,dfdiff,lower.tail = F)
  }
  else if(FN==TRUE){
    r <- 1
    k <- 6
    n <- n1+n2
    Fstat <- stat*(n-k)/n/r
    pvalue <- pf(Fstat,df1=r,df2=n-k,lower.tail = F)
  }

  amp_c <- x_H0[1]
  amp_c2 <- x_H0[5]

  res <- list(amp_1=A1, amp_2=A2, amp_c=amp_c,
              #df = dfdiff,
              stat = stat,
              pvalue = pvalue)
  return(res)
}


##' Finite sample/Large sample Wald test for differential basal level (vertical shift).
##'
##' Test differential basal level of circadian curve fitting using Wald test
##' @title Wald test for detecting differential offset
##' @param tt1 Time vector of condition 1
##' @param yy1 Expression vector of condition 1
##' @param tt2 Time vector of condition 2
##' @param yy2 Expression vector of condition 2
##' @param period Period of the since curve. Default is 24.
##' @param type Type of Wald test to use, "FN" or "LS". Default is finite sample.
##' @return A list, see details below.
##' Formula 1: \eqn{yy = amp \times sin(2\pi/period \times (phase + tt)) + offset}
##' Formula 2: \eqn{yy = A \times sin(2\pi/period \times tt) + B * cos(2*pi/period * tt) + offset}
##' \item{offset_1}{Basal level estimate of the 1st data}
##' \item{offset_2}{Basal level estimate of the 2nd data}
##' \item{offset_c}{Basal level estimate pooling all data together}
##' \item{df}{Degree of freedom for the Wald test}
##' \item{stat}{Wald statistics}
##' \item{pvalue}{P-value from the Wald test}
##' @author Caleb
##' @noRd
##' @examples
##' set.seed(32608)
##' n <- 50
##' tt1 <- runif(n,0,24)
##' Amp1 <- 2
##' Phase1 <- 6
##' Offset1 <- 3
##' yy1 <- Amp1 * sin(2*pi/24 * (tt1 + Phase1)) + Offset1 + rnorm(n,0,1)
##' tt2 <- runif(n,0,24)
##' Amp2 <- 3
##' Phase2 <- 5
##' Offset2 <- 2
##' yy2 <- Amp2 * sin(2*pi/24 * (tt2 + Phase2)) + Offset2 + rnorm(n,0,1)
##' WaldTest_diff_offset(tt1, yy1, tt2, yy2)



WaldTest_diff_offset <- function(tt1, yy1, tt2, yy2, period = 24,FN=TRUE){
  n1 <- length(tt1)
  stopifnot(n1 == length(yy1))
  n2 <- length(tt2)
  stopifnot(length(tt2) == length(yy2))

  #period <- 24
  w <- 2*pi/period

  fit1 <- fitSinCurve(tt1, yy1, period = period)
  fit2 <- fitSinCurve(tt2, yy2, period = period)

  A1 <- fit1$amp
  A2 <- fit2$amp

  phase1 <- fit1$phase
  phase2 <- fit2$phase

  basal1 <- fit1$offset
  basal2 <- fit2$offset

  sigma2_1 <- 1/n1 * fit1$rss
  sigma2_2 <- 1/n2 * fit2$rss

  theta1 <- 1/sigma2_1
  theta2 <- 1/sigma2_2

  p1 <- c(A1, phase1, basal1, theta1)
  p2 <- c(A2, phase2, basal2, theta2)

  x_Ha <- c(p1, p2)

  eval_f_list <- function(x) {
    p1 <- x[1:4]
    p2 <- x[5:8]

    A1 <- p1[1]
    phase1 <- p1[2]
    basel1 <- p1[3]
    theta1 <- p1[4]
    asin1 <- sin(w * (tt1 + phase1) )
    acos1 <- cos(w * (tt1 + phase1) )
    yhat1 <- A1 * asin1 + basel1

    A2 <- p2[1]
    phase2 <- p2[2]
    basel2 <- p2[3]
    theta2 <- p2[4]
    asin2 <- sin(w * (tt2 + phase2) )
    acos2 <- cos(w * (tt2 + phase2) )
    yhat2 <- A2 * asin2 + basel2

    ll1_a <- log(theta1)/2
    ll1_b <- (yy1 - yhat1)^2 * theta1 / 2
    ll1 <- ll1_a - ll1_b

    ll2_a <- log(theta2)/2
    ll2_b <- (yy2 - yhat2)^2 * theta2 / 2
    ll2 <- ll2_a - ll2_b

    partial_A1 <- - theta1 * sum((yy1 - yhat1) * asin1)
    partial_phase1 <- - theta1 * A1 * w * sum((yy1 - yhat1) * acos1)
    partial_C1 <- - theta1 * sum(yy1 - yhat1)
    partial_theta1 <-  sum((yy1 - yhat1)^2)/2 - n1/2/theta1

    partial_A2 <- - theta2 * sum((yy2 - yhat2) * asin2)
    partial_phase2 <- - theta2 * A2 * w * sum((yy2 - yhat2) * acos2)
    partial_C2 <- - theta2 * sum(yy2 - yhat2)
    partial_theta2 <-  sum((yy2 - yhat2)^2)/2 - n2/2/theta2


    return( list( "objective" = -sum(ll1) - sum(ll2),
                  "gradient"  = c(partial_A1, partial_phase1, partial_C1, partial_theta1,
                                  partial_A2, partial_phase2, partial_C2, partial_theta2)
    )
    )
  }

  # Equality constraints
  eval_g_eq <- function(x)
  {
    p1 <- x[1:4]
    p2 <- x[5:8]

    basal1 <- p1[3]
    basal2 <- p2[3]

    basal1 - basal2
  }

  # Equality constraints
  eval_g_eq_jac <- function(x)
  {
    c(0, 0, 1, 0,
      0, 0, -1, 0)
  }


  # Lower and upper bounds
  lb <- c(0,-Inf,-Inf,0, 0, -Inf,-Inf, 0)
  ub <- c(Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf)
  #initial values

  ## Error in is.nloptr(ret) :
  #  If you want to use equality constraints, then you should use one of these algorithms NLOPT_LD_AUGLAG, NLOPT_LN_AUGLAG, NLOPT_LD_AUGLAG_EQ, NLOPT_LN_AUGLAG_EQ, NLOPT_GN_ISRES, NLOPT_LD_SLSQP

  local_opts <- list( "algorithm" = "NLOPT_LD_MMA", "xtol_rel" = 1.0e-15 )
  "local_opts" = local_opts
  opts <- list( "algorithm"= "NLOPT_LD_SLSQP",
                "xtol_rel"= 1.0e-15,
                "maxeval"= 160000,
                "local_opts" = local_opts,
                "print_level" = 0
                #"check_derivatives"=TRUE
  )

  res <- nloptr::nloptr( x0 = x_Ha,
                  eval_f = eval_f_list,
                  #eval_grad_f=eval_g,
                  lb = lb,
                  ub = ub,
                  #eval_g_ineq = eval_g_ineq,
                  eval_g_eq = eval_g_eq,
                  eval_jac_g_eq = eval_g_eq_jac,
                  opts = opts)

  #
  #x_Ha
  x_H0 <- res$solution

  beta_ha <- x_Ha
  beta_h0 <- x_H0

  #
  I8 <- fisherInformation2(beta_ha, tt1, yy1, tt2, yy2, period=period)
  beta_diff <- matrix(beta_ha - beta_h0,ncol=1)
  #stat <- t(beta_diff) %*% I8 %*% beta_diff

  beta_diff2 <- beta_diff[c(3,7)]
  I2 <- solve(solve(I8)[c(3,7),c(3,7)])
  stat <- as.numeric(t(beta_diff2) %*% I2 %*% beta_diff2)

  dfdiff <- 1

  if(FN==FALSE){
    pvalue <- pchisq(stat,dfdiff,lower.tail = F)
  }
  else if(FN==TRUE){
    r <- 1
    k <- 6
    n <- n1+n2
    Fstat <- stat*(n-k)/n/r
    pvalue <- pf(Fstat,df1=r,df2=n-k,lower.tail = F)
  }


  offset_c <- x_H0[3]
  offset_c2 <- x_H0[7]

  res <- list(offset_1=basal1, offset_2=basal2, offset_c=offset_c,
              #df = dfdiff,
              stat = stat,
              pvalue = pvalue)
  return(res)
}


##' Finite sample/Large sample Wald test for differential phase.
##'
##' Test differential phase of circadian curve fitting using Wald test
##' @title Wald test for detecting differential phase
##' @param tt1 Time vector of condition 1
##' @param yy1 Expression vector of condition 1
##' @param tt2 Time vector of condition 2
##' @param yy2 Expression vector of condition 2
##' @param period Period of the since curve. Default is 24.
##' @param type Type of Wald test to use, "FN" or "LS". Default is finite sample.
##' @return A list, see details below.
##' Formula 1: \eqn{yy = amp \times sin(2\pi/period \times (phase + tt)) + offset}
##' Formula 2: \eqn{yy = A \times sin(2\pi/period \times tt) + B * cos(2*pi/period * tt) + offset}
##' \item{phase_1}{Phase estimate of the 1st data, phase is restricted in (0, period)}
##' \item{phase_2}{Phase estimate of the 2nd data, phase is restricted in (0, period)}
##' \item{phase_c}{Phase estimate pooling all data together, phase is restricted in (0, period)}
##' \item{df}{Degree of freedom for the Wald test}
##' \item{stat}{Wald statistics}
##' \item{pvalue}{P-value from the Wald test}
##' @author Caleb
##' @noRd
##' @examples
##' set.seed(32608)
##' n <- 50
##' tt1 <- runif(n,0,24)
##' Amp1 <- 2
##' Phase1 <- 6
##' Offset1 <- 3
##' yy1 <- Amp1 * sin(2*pi/24 * (tt1 + Phase1)) + Offset1 + rnorm(n,0,1)
##' tt2 <- runif(n,0,24)
##' Amp2 <- 3
##' Phase2 <- 5
##' Offset2 <- 2
##' yy2 <- Amp2 * sin(2*pi/24 * (tt2 + Phase2)) + Offset2 + rnorm(n,0,1)
##' WaldTest_diff_phase(tt1, yy1, tt2, yy2)


WaldTest_diff_phase <- function(tt1, yy1, tt2, yy2, period = 24,FN=TRUE){
  n1 <- length(tt1)
  stopifnot(n1 == length(yy1))
  n2 <- length(tt2)
  stopifnot(length(tt2) == length(yy2))

  #period <- 24
  w <- 2*pi/period

  fit1 <- fitSinCurve(tt1, yy1, period = period)
  fit2 <- fitSinCurve(tt2, yy2, period = period)

  A1 <- fit1$amp
  A2 <- fit2$amp

  phase1 <- fit1$phase
  phase2 <- fit2$phase

  if(phase2 - phase1 > period/2){
    phase2 <- phase2 - period
  } else if(phase1 - phase2 > period/2){
    phase1 <- phase1 - period
  }

  basal1 <- fit1$offset
  basal2 <- fit2$offset

  sigma2_1 <- 1/n1 * fit1$rss
  sigma2_2 <- 1/n2 * fit2$rss

  theta1 <- 1/sigma2_1
  theta2 <- 1/sigma2_2

  p1 <- c(A1, phase1, basal1, theta1)
  p2 <- c(A2, phase2, basal2, theta2)

  x_Ha <- c(p1, p2)

  eval_f_list <- function(x) {
    p1 <- x[1:4]
    p2 <- x[5:8]

    A1 <- p1[1]
    phase1 <- p1[2]
    basel1 <- p1[3]
    theta1 <- p1[4]
    asin1 <- sin(w * (tt1 + phase1) )
    acos1 <- cos(w * (tt1 + phase1) )
    yhat1 <- A1 * asin1 + basel1

    A2 <- p2[1]
    phase2 <- p2[2]
    basel2 <- p2[3]
    theta2 <- p2[4]
    asin2 <- sin(w * (tt2 + phase2) )
    acos2 <- cos(w * (tt2 + phase2) )
    yhat2 <- A2 * asin2 + basel2

    ll1_a <- log(theta1)/2
    ll1_b <- (yy1 - yhat1)^2 * theta1 / 2
    ll1 <- ll1_a - ll1_b

    ll2_a <- log(theta2)/2
    ll2_b <- (yy2 - yhat2)^2 * theta2 / 2
    ll2 <- ll2_a - ll2_b

    partial_A1 <- - theta1 * sum((yy1 - yhat1) * asin1)
    partial_phase1 <- - theta1 * A1 * w * sum((yy1 - yhat1) * acos1)
    partial_C1 <- - theta1 * sum(yy1 - yhat1)
    partial_theta1 <-  sum((yy1 - yhat1)^2)/2 - n1/2/theta1

    partial_A2 <- - theta2 * sum((yy2 - yhat2) * asin2)
    partial_phase2 <- - theta2 * A2 * w * sum((yy2 - yhat2) * acos2)
    partial_C2 <- - theta2 * sum(yy2 - yhat2)
    partial_theta2 <-  sum((yy2 - yhat2)^2)/2 - n2/2/theta2


    return( list( "objective" = -sum(ll1) - sum(ll2),
                  "gradient"  = c(partial_A1, partial_phase1, partial_C1, partial_theta1,
                                  partial_A2, partial_phase2, partial_C2, partial_theta2)
    )
    )
  }

  # Equality constraints
  eval_g_eq <- function(x)
  {
    p1 <- x[1:4]
    p2 <- x[5:8]

    phase1 <- p1[2]
    phase2 <- p2[2]

    phase1 - phase2
  }

  # Equality constraints
  eval_g_eq_jac <- function(x)
  {
    c(0, 1, 0, 0,
      0, -1, 0, 0)
  }


  # Lower and upper bounds
  lb <- c(0,-Inf,-Inf,0, 0, -Inf,-Inf, 0)
  ub <- c(Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf)
  #initial values

  ## Error in is.nloptr(ret) :
  #  If you want to use equality constraints, then you should use one of these algorithms NLOPT_LD_AUGLAG, NLOPT_LN_AUGLAG, NLOPT_LD_AUGLAG_EQ, NLOPT_LN_AUGLAG_EQ, NLOPT_GN_ISRES, NLOPT_LD_SLSQP

  local_opts <- list( "algorithm" = "NLOPT_LD_MMA", "xtol_rel" = 1.0e-15 )
  "local_opts" = local_opts
  opts <- list( "algorithm"= "NLOPT_LD_SLSQP",
                "xtol_rel"= 1.0e-15,
                "maxeval"= 160000,
                "local_opts" = local_opts,
                "print_level" = 0
                #"check_derivatives"=TRUE
  )

  res <- nloptr::nloptr( x0 = x_Ha,
                  eval_f = eval_f_list,
                  #eval_grad_f=eval_g,
                  lb = lb,
                  ub = ub,
                  #eval_g_ineq = eval_g_ineq,
                  eval_g_eq = eval_g_eq,
                  eval_jac_g_eq = eval_g_eq_jac,
                  opts = opts)

  #
  #x_Ha
  x_H0 <- res$solution

  beta_ha <- x_Ha
  beta_h0 <- x_H0

  #
  I8 <- fisherInformation2(beta_ha, tt1, yy1, tt2, yy2, period=period)
  beta_diff <- matrix(beta_ha - beta_h0,ncol=1)
  #stat <- t(beta_diff) %*% I8 %*% beta_diff

  beta_diff2 <- beta_diff[c(2,6)]
  I2 <- solve(solve(I8)[c(2,6),c(2,6)])
  stat <- as.numeric(t(beta_diff2) %*% I2 %*% beta_diff2)

  dfdiff <- 1



  if(FN==FALSE){
    pvalue <- pchisq(stat,dfdiff,lower.tail = F)
  }
  else if(FN==TRUE){
    r <- 1
    k <- 6
    n <- n1+n2
    Fstat <- stat*(n-k)/n/r
    pvalue <- pf(Fstat,df1=r,df2=n-k,lower.tail = F)
  }


  phase_c <- x_H0[2]
  phase_c2 <- x_H0[6]

  res <- list(phase_1=phase1, phase_2=phase2, phase_c=phase_c,
              #df = dfdiff,
              stat = stat,
              pvalue = pvalue)
  return(res)
}



##' Finite sample/Large sample Wald test for differential sigma square.
##'
##' Test differential sigma square of circadian curve fitting using Wald test
##' @title Wald test for detecting differential sigma square
##' @param tt1 Time vector of condition 1
##' @param yy1 Expression vector of condition 1
##' @param tt2 Time vector of condition 2
##' @param yy2 Expression vector of condition 2
##' @param period Period of the since curve. Default is 24.
##' @param type Type of Wald test to use, "FN" or "LS". Default is finite sample.
##' @return A list, see details below.
##' Formula 1: \eqn{yy = amp \times sin(2\pi/period \times (phase + tt)) + offset}
##' Formula 2: \eqn{yy = A \times sin(2\pi/period \times tt) + B * cos(2*pi/period * tt) + offset}
##' \item{sigma2_1}{Variance estimate of the 1st data}
##' \item{sigma2_2}{Variance estimate of the 2nd data}
##' \item{sigma2_C}{Variance estimate pooling all data together}
##' \item{df}{Degree of freedom for the Wald test}
##' \item{stat}{Wald statistics}
##' \item{pvalue}{P-value from the Wald test}
##' @author Caleb
##' @noRd
##' @examples
##' set.seed(32608)
##' n <- 50
##' tt1 <- runif(n,0,24)
##' Amp1 <- 2
##' Phase1 <- 6
##' Offset1 <- 3
##' yy1 <- Amp1 * sin(2*pi/24 * (tt1 + Phase1)) + Offset1 + rnorm(n,0,1)
##' tt2 <- runif(n,0,24)
##' Amp2 <- 3
##' Phase2 <- 5
##' Offset2 <- 2
##' yy2 <- Amp2 * sin(2*pi/24 * (tt2 + Phase2)) + Offset2 + rnorm(n,0,1)
##' WaldTest_diff_sigma2(tt1, yy1, tt2, yy2)


WaldTest_diff_sigma2 <- function(tt1, yy1, tt2, yy2, period = 24,FN=TRUE){
  n1 <- length(tt1)
  stopifnot(n1 == length(yy1))
  n2 <- length(tt2)
  stopifnot(length(tt2) == length(yy2))

  #period <- 24
  w <- 2*pi/period

  fit1 <- fitSinCurve(tt1, yy1, period = period)
  fit2 <- fitSinCurve(tt2, yy2, period = period)

  A1 <- fit1$amp
  A2 <- fit2$amp

  phase1 <- fit1$phase
  phase2 <- fit2$phase

  basal1 <- fit1$offset
  basal2 <- fit2$offset

  sigma2_1 <- 1/n1 * fit1$rss
  sigma2_2 <- 1/n2 * fit2$rss
  sigma2_c <- 1/(n1 + n2) * (sigma2_1 * n1 + sigma2_2 * n2)

  theta1 <- 1/sigma2_1
  theta2 <- 1/sigma2_2
  thetaC <- 1/sigma2_c

  p1 <- c(A1, phase1, basal1, theta1)
  p2 <- c(A2, phase2, basal2, theta2)

  x_Ha <- c(p1, p2)
  x_H0 <- c(A1, phase1, basal1, thetaC, A2, phase2, basal2, thetaC)

  beta_ha <- x_Ha
  beta_h0 <- x_H0


  I8 <- fisherInformation2(beta_ha, tt1, yy1, tt2, yy2, period=period)
  beta_diff <- matrix(beta_ha - beta_h0,ncol=1)
  #stat <- t(beta_diff) %*% I8 %*% beta_diff

  beta_diff2 <- beta_diff[c(4,8)]
  I2 <- solve(solve(I8)[c(4,8),c(4,8)])
  stat <- as.numeric(t(beta_diff2) %*% I2 %*% beta_diff2)

  dfdiff <- 1

  if(FN==FALSE){
    pvalue <- pchisq(stat,dfdiff,lower.tail = F)
  }
  else if(FN==TRUE){
    r <- 1
    k <- 6
    n <- n1+n2
    Fstat <- stat*(n-k)/n/r
    pvalue <- pf(Fstat,df1=r,df2=n-k,lower.tail = F)
  }


  res <- list(sigma2_1=sigma2_1, sigma2_2=sigma2_2, sigma2_c=sigma2_c,
              stat=stat,
              pvalue=pvalue)
  return(res)
}



##' Finite sample/ Large sample Wald test for circadian pattern detection
##'
##' Test the significance of circadian curve fitting using finite sample Wald test
##' @title Wald test for detecting circadian pattern.
##' @param tt Time vector
##' @param yy Expression vector
##' @param period Period of the since curve. Default is 24.
##' @param type Type of Test, finite sample "FN" or large sample "LS", default is "FN".
##' @return A list of A, B, offset, df, stat, and pvalue
##' Formula 1: \eqn{yy = amp \times sin(2\pi/period \times (phase + tt)) + offset}
##' Formula 2: \eqn{yy = A \times sin(2\pi/period \times tt) + B * cos(2*pi/period * tt) + offset}
##' \item{A}{A based on formula 2}
##' \item{B}{B based on formula 2}
##' \item{amp}{Amplitude based on formula 1}
##' \item{phase}{Phase based on formula 1, phase is restricted within (0, period)}
##' \item{peakTime}{Phase based on formula 1, peakTime is restricted within (0, period). phase + peakTime = period/4}
##' \item{offset}{Basal level based on formula 1 or on formula 2}
##' \item{df}{Degree of freedom for the Wald test}
##' \item{stat}{Wald statistics}
##' \item{pvalue}{P-value from the Wald test}
##' \item{R2}{Pseudo R2 defined as (tss - rss)/tss}
##' @author Caleb
##' @noRd
##' @examples
##' set.seed(32608)
##' n <- 50
##' tt <- runif(n,0,24)
##' Amp <- 2
##' Phase <- 6
##' Offset <- 3
##' yy <- Amp * sin(2*pi/24 * (tt + Phase)) + Offset + rnorm(n,0,1)
##' WaldTest(tt, yy)

WaldTest <- function(tt, yy, period = 24, type=TRUE){
  afit <- fitSinCurve(tt, yy, period=period)
  n <- length(tt)
  rss <- afit$rss
  tss <- afit$tss
  amp <- afit$amp
  phase <- afit$phase

  ## HA:
  A <- afit$A
  B <- afit$B
  offset <- afit$offset

  asin <- sin(2*pi/period * tt)
  acos <- cos(2*pi/period * tt)

  yhat <- A * asin + B * acos + offset
  residual <- yy - yhat

  sigmaA2 <- 1/n * sum(residual^2)
  invSigmaA2 <- 1/sigmaA2


  det_A_A <- - invSigmaA2 * sum(asin^2)
  det_A_B <- - invSigmaA2 * sum(asin * acos)
  det_A_offset <- - invSigmaA2 * sum(asin)
  det_A_sigma2 <- - invSigmaA2^2 * sum(residual * asin)

  det_B_B <- - invSigmaA2 * sum(acos^2)
  det_B_offset <- - invSigmaA2 * sum(acos)
  det_B_sigma2 <- - invSigmaA2^2 * sum(residual * acos)

  det_offset_offset <- - n * invSigmaA2
  det_offset_sigma2 <- - invSigmaA2^2 * sum(residual)

  det_sigma2_sigma2 <- n/2*invSigmaA2^2 - invSigmaA2^3 * sum(residual^2)

  r1 <- c(det_A_A,det_A_B,det_A_offset,det_A_sigma2)
  r2 <- c(det_A_B,det_B_B,det_B_offset,det_B_sigma2)
  r3 <- c(det_A_offset,det_B_offset,det_offset_offset,det_offset_sigma2)
  r4 <- c(det_A_sigma2,det_B_sigma2,det_offset_sigma2,det_sigma2_sigma2)
  I <- - rbind(r1,r2,r3,r4)

  I_test <- solve(solve(I)[1:2,1:2])
  Waldstat <- matrix(c(A, B), nrow = 1, ncol = 2) %*%
    I_test %*% matrix(c(A, B), nrow = 2, ncol = 1)

  df <- 2
  stat <- as.numeric(Waldstat)
  if(type==FALSE){
    pvalue <- pchisq(stat,2,lower.tail = F)
  }

  else if(type==TRUE){
    r <- 2
    k <- 3
    stat <- stat * (n-k)/n/r
    pvalue <- pf(stat,df1 = r, df2 = n-k, lower.tail = F)
  }

  R2 <- 1-rss/tss
  res <- list(
    amp=amp,phase=phase,
    peakTime = (6 - phase) %% period,
    offset=offset,
    stat=stat,
    pvalue=pvalue,R2=R2
  )

  return(res)

}




##' Finite sample/Large sample Likelihood ratio test for differential basal level.
##'
##' Test differential offset of circadian curve fitting using likelihood ratio test
##' @title Likelihood ratio test for detecting differential basal level.
##' @param tt1 Time vector of condition 1
##' @param yy1 Expression vector of condition 1
##' @param tt2 Time vector of condition 2
##' @param yy2 Expression vector of condition 2
##' @param period Period of the since curve. Default is 24.
##' @param type Type of likelihood ratio test to use, "FN" or "LS". Default is finite sample.
##' @return A list, see details below.
##' Formula 1: \eqn{yy = amp \times sin(2\pi/period \times (phase + tt)) + offset}
##' Formula 2: \eqn{yy = A \times sin(2\pi/period \times tt) + B * cos(2*pi/period * tt) + offset}
##' \item{offset_1}{Basal level estimate of the 1st data}
##' \item{offset_2}{Basal level estimate of the 2nd data}
##' \item{offset_c}{Basal level estimate pooling all data together}
##' \item{l0}{Log likelihood under the null (same variance between the two groups)}
##' \item{l1}{Log likelihood under the alternative (different variance between the two groups)}
##' \item{df}{Degree of freedom for the LR test}
##' \item{stat}{LR statistics}
##' \item{pvalue}{P-value from the LR test}
##' @author Caleb
##' @noRd
##' @examples
##' set.seed(32608)
##' n <- 50
##' tt1 <- runif(n,0,24)
##' Amp1 <- 2
##' Phase1 <- 6
##' Offset1 <- 3
##' yy1 <- Amp1 * sin(2*pi/24 * (tt1 + Phase1)) + Offset1 + rnorm(n,0,1)
##' tt2 <- runif(n,0,24)
##' Amp2 <- 3
##' Phase2 <- 5
##' Offset2 <- 2
##' yy2 <- Amp2 * sin(2*pi/24 * (tt2 + Phase2)) + Offset2 + rnorm(n,0,1)
##' LRTest_diff_offset(tt1, yy1, tt2, yy2)


LRTest_diff_offset <- function(tt1, yy1, tt2, yy2, period = 24,FN=TRUE){
  n1 <- length(tt1)
  stopifnot(n1 == length(yy1))
  n2 <- length(tt2)
  stopifnot(length(tt2) == length(yy2))

  #period <- 24
  w <- 2*pi/period

  fit1 <- fitSinCurve(tt1, yy1, period = period)
  fit2 <- fitSinCurve(tt2, yy2, period = period)

  A1 <- fit1$amp
  A2 <- fit2$amp

  phase1 <- fit1$phase
  phase2 <- fit2$phase

  E1 <- A1 * cos(w * phase1)
  F1 <- A1 * sin(w * phase1)

  E2 <- A2 * cos(w * phase2)
  F2 <- A2 * sin(w * phase2)

  basal1 <- fit1$offset
  basal2 <- fit2$offset

  sigma2_1 <- 1/n1 * fit1$rss
  sigma2_2 <- 1/n2 * fit2$rss

  theta1 <- 1/sigma2_1
  theta2 <- 1/sigma2_2

  p1 <- c(E1, F1, basal1, theta1)
  p2 <- c(E2, F2, basal2, theta2)

  x_Ha <- c(p1, p2)

  asin1 <- sin(w * tt1)
  acos1 <- cos(w * tt1)
  asin2 <- sin(w * tt2)
  acos2 <- cos(w * tt2)

  eval_f_list <- function(x,asin1,acos1,asin2,acos2) {
    p1 <- x[1:4]
    p2 <- x[5:8]

    E1 <- p1[1]
    F1 <- p1[2]
    basel1 <- p1[3]
    theta1 <- p1[4]
    yhat1 <- E1 * asin1 + F1 * acos1 + basel1

    E2 <- p2[1]
    F2 <- p2[2]
    basel2 <- p2[3]
    theta2 <- p2[4]
    yhat2 <- E2 * asin2 + F2 * acos2 + basel2

    ll1_a <- log(theta1)/2
    ll1_b <- (yy1 - yhat1)^2 * theta1 / 2
    ll1 <- ll1_a - ll1_b

    ll2_a <- log(theta2)/2
    ll2_b <- (yy2 - yhat2)^2 * theta2 / 2
    ll2 <- ll2_a - ll2_b

    partial_E1 <- - theta1 * sum((yy1 - yhat1) * asin1)
    partial_F1 <- - theta1 * sum((yy1 - yhat1) * acos1)
    partial_C1 <- - theta1 * sum(yy1 - yhat1)
    partial_theta1 <-  sum((yy1 - yhat1)^2)/2 - n1/2/theta1

    partial_E2 <- - theta2 * sum((yy2 - yhat2) * asin2)
    partial_F2 <- - theta2 * sum((yy2 - yhat2) * acos2)
    partial_C2 <- - theta2 * sum(yy2 - yhat2)
    partial_theta2 <-  sum((yy2 - yhat2)^2)/2 - n2/2/theta2


    return( list( "objective" = -sum(ll1) - sum(ll2),
                  "gradient"  = c(partial_E1, partial_F1, partial_C1, partial_theta1,
                                  partial_E2, partial_F2, partial_C2, partial_theta2)
    )
    )
  }

  # Equality constraints
  eval_g_eq <- function(x,asin1,acos1,asin2,acos2)
  {
    p1 <- x[1:4]
    p2 <- x[5:8]

    E1 <- p1[1]
    F1 <- p1[2]
    basel1 <- p1[3]
    theta1 <- p1[4]
    #yhat1 <- E1 * asin1 + F1 * acos1 + basel1

    E2 <- p2[1]
    F2 <- p2[2]
    basel2 <- p2[3]
    theta2 <- p2[4]
    #yhat2 <- E2 * asin2 + F2 * acos2 + basel2

    basel1 - basel2
  }

  # Equality constraints
  eval_g_eq_jac <- function(x,asin1,acos1,asin2,acos2)
  {
    p1 <- x[1:4]
    p2 <- x[5:8]

    E1 <- p1[1]
    F1 <- p1[2]
    #basel1 <- p1[3]
    theta1 <- p1[4]
    #yhat1 <- E1 * asin1 + F1 * acos1 + basel1

    E2 <- p2[1]
    F2 <- p2[2]
    #basel2 <- p2[3]
    theta2 <- p2[4]
    #yhat2 <- E2 * asin2 + F2 * acos2 + basel2

    A2_1 <- (E1^2 + F1^2)
    A2_2 <- (E2^2 + F2^2)

    c(0, 0, 1, 0,
      0, 0, -1, 0)
  }


  # Lower and upper bounds
  lb <- c(-Inf,-Inf,-Inf,0, -Inf, -Inf,-Inf, 0)
  ub <- c(Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf)
  #initial values

  ## Error in is.nloptr(ret) :
  #  If you want to use equality constraints, then you should use one of these algorithms NLOPT_LD_AUGLAG, NLOPT_LN_AUGLAG, NLOPT_LD_AUGLAG_EQ, NLOPT_LN_AUGLAG_EQ, NLOPT_GN_ISRES, NLOPT_LD_SLSQP

  local_opts <- list( "algorithm" = "NLOPT_LD_MMA", "xtol_rel" = 1.0e-15 )
  "local_opts" = local_opts
  opts <- list( "algorithm"= "NLOPT_LD_SLSQP",
                "xtol_rel"= 1.0e-15,
                "maxeval"= 160000,
                "local_opts" = local_opts,
                "print_level" = 0
                #"check_derivatives"=TRUE
  )

  res <- nloptr::nloptr( x0 = x_Ha,
                  eval_f = eval_f_list,
                  #eval_grad_f=eval_g,
                  lb = lb,
                  ub = ub,
                  #eval_g_ineq = eval_g_ineq,
                  eval_g_eq = eval_g_eq,
                  eval_jac_g_eq = eval_g_eq_jac,
                  opts = opts,
                  asin1=asin1,
                  acos1=acos1,
                  asin2=asin2,
                  acos2=acos2)

  #
  #x_Ha
  x_H0 <- res$solution

  l0 <- - eval_f_list(x_H0,asin1,acos1,asin2,acos2)$objective
  la <- - eval_f_list(x_Ha,asin1,acos1,asin2,acos2)$objective

  LR_stat <- -2*(l0-la)

  dfdiff <- 1
  if(!FN){
    pvalue <- pchisq(LR_stat,dfdiff,lower.tail = F)
  } else if(FN){
    r <- 1
    k <- 6
    n <- n1+n2
    Fstat <- (exp(LR_stat/n) - 1) * (n-k) / r
    pvalue <- pf(Fstat,df1 = r, df2 = n-k, lower.tail = F)
  } else{
    stop("FN has to be TRUE or FALSE")
  }

  offset_c <- x_H0[3]
  offset_c2 <- x_H0[7]


  res <- list(offset_1=basal1, offset_2=basal2, offset_c=offset_c,
              l0=l0,
              la=la,
              #df = dfdiff,
              stat=-2*(l0-la),
              pvalue=pvalue)
  return(res)
}


##' Finite sample/Large sample Likelihood ratio test for differential amplitude.
##'
##' Test differential amplitude of circadian curve fitting using likelihood ratio test
##' @title Likelihood ratio test for detecting differential amplitudes.
##' @param tt1 Time vector of condition 1
##' @param yy1 Expression vector of condition 1
##' @param tt2 Time vector of condition 2
##' @param yy2 Expression vector of condition 2
##' @param period Period of the since curve. Default is 24.
##' @param type Type of likelihood ratio test to use, "FN" or "LS". Default is finite sample.
##' @return A list, see details below.
##' Formula 1: \eqn{yy = amp \times sin(2\pi/period \times (phase + tt)) + offset}
##' Formula 2: \eqn{yy = A \times sin(2\pi/period \times tt) + B * cos(2*pi/period * tt) + offset}
##' \item{amp_1}{Amplitude estimate of the 1st data}
##' \item{amp_2}{Amplitude estimate of the 2nd data}
##' \item{amp_c}{Amplitude estimate pooling all data together}
##' \item{l0}{Log likelihood under the null (same variance between the two groups)}
##' \item{l1}{Log likelihood under the alternative (different variance between the two groups)}
##' \item{df}{Degree of freedom for the LR test}
##' \item{stat}{LR statistics}
##' \item{pvalue}{P-value from the LR test}
##' @author Caleb
##' @noRd
##' @examples
##' set.seed(32608)
##' n <- 50
##' tt1 <- runif(n,0,24)
##' Amp1 <- 2
##' Phase1 <- 6
##' Offset1 <- 3
##' yy1 <- Amp1 * sin(2*pi/24 * (tt1 + Phase1)) + Offset1 + rnorm(n,0,1)
##' tt2 <- runif(n,0,24)
##' Amp2 <- 3
##' Phase2 <- 5
##' Offset2 <- 2
##' yy2 <- Amp2 * sin(2*pi/24 * (tt2 + Phase2)) + Offset2 + rnorm(n,0,1)
##' LRTest_diff_amp(tt1, yy1, tt2, yy2)


LRTest_diff_amp<- function(tt1, yy1, tt2, yy2, period = 24,FN=TRUE){

  n1 <- length(tt1)
  stopifnot(n1 == length(yy1))
  n2 <- length(tt2)
  stopifnot(length(tt2) == length(yy2))

  #period <- 24
  w <- 2*pi/period

  fit1 <- fitSinCurve(tt1, yy1, period = period)
  fit2 <- fitSinCurve(tt2, yy2, period = period)

  A1 <- fit1$amp
  A2 <- fit2$amp

  phase1 <- fit1$phase
  phase2 <- fit2$phase

  E1 <- A1 * cos(w * phase1)
  F1 <- A1 * sin(w * phase1)

  E2 <- A2 * cos(w * phase2)
  F2 <- A2 * sin(w * phase2)

  basal1 <- fit1$offset
  basal2 <- fit2$offset

  sigma2_1 <- 1/n1 * fit1$rss
  sigma2_2 <- 1/n2 * fit2$rss

  theta1 <- 1/sigma2_1
  theta2 <- 1/sigma2_2

  p1 <- c(E1, F1, basal1, theta1)
  p2 <- c(E2, F2, basal2, theta2)

  x_Ha <- c(p1, p2)

  asin1 <- sin(w * tt1)
  acos1 <- cos(w * tt1)
  asin2 <- sin(w * tt2)
  acos2 <- cos(w * tt2)

  eval_f_list <- function(x,asin1,acos1,asin2,acos2) {
    p1 <- x[1:4]
    p2 <- x[5:8]

    E1 <- p1[1]
    F1 <- p1[2]
    basel1 <- p1[3]
    theta1 <- p1[4]
    yhat1 <- E1 * asin1 + F1 * acos1 + basel1

    E2 <- p2[1]
    F2 <- p2[2]
    basel2 <- p2[3]
    theta2 <- p2[4]
    yhat2 <- E2 * asin2 + F2 * acos2 + basel2

    ll1_a <- log(theta1)/2
    ll1_b <- (yy1 - yhat1)^2 * theta1 / 2
    ll1 <- ll1_a - ll1_b

    ll2_a <- log(theta2)/2
    ll2_b <- (yy2 - yhat2)^2 * theta2 / 2
    ll2 <- ll2_a - ll2_b

    partial_E1 <- - theta1 * sum((yy1 - yhat1) * asin1)
    partial_F1 <- - theta1 * sum((yy1 - yhat1) * acos1)
    partial_C1 <- - theta1 * sum(yy1 - yhat1)
    partial_theta1 <-  sum((yy1 - yhat1)^2)/2 - n1/2/theta1

    partial_E2 <- - theta2 * sum((yy2 - yhat2) * asin2)
    partial_F2 <- - theta2 * sum((yy2 - yhat2) * acos2)
    partial_C2 <- - theta2 * sum(yy2 - yhat2)
    partial_theta2 <-  sum((yy2 - yhat2)^2)/2 - n2/2/theta2


    return( list( "objective" = -sum(ll1) - sum(ll2),
                  "gradient"  = c(partial_E1, partial_F1, partial_C1, partial_theta1,
                                  partial_E2, partial_F2, partial_C2, partial_theta2)
    )
    )
  }

  # Equality constraints
  eval_g_eq <- function(x,asin1,acos1,asin2,acos2)
  {
    p1 <- x[1:4]
    p2 <- x[5:8]

    E1 <- p1[1]
    F1 <- p1[2]
    #basel1 <- p1[3]
    theta1 <- p1[4]
    #yhat1 <- E1 * asin1 + F1 * acos1 + basel1

    E2 <- p2[1]
    F2 <- p2[2]
    #basel2 <- p2[3]
    theta2 <- p2[4]
    #yhat2 <- E2 * asin2 + F2 * acos2 + basel2

    A2_1 <- (E1^2 + F1^2)
    A2_2 <- (E2^2 + F2^2)
    A2_1 - A2_2
  }

  # Equality constraints
  eval_g_eq_jac <- function(x,asin1,acos1,asin2,acos2)
  {
    p1 <- x[1:4]
    p2 <- x[5:8]

    E1 <- p1[1]
    F1 <- p1[2]
    #basel1 <- p1[3]
    theta1 <- p1[4]
    #yhat1 <- E1 * asin1 + F1 * acos1 + basel1

    E2 <- p2[1]
    F2 <- p2[2]
    #basel2 <- p2[3]
    theta2 <- p2[4]
    #yhat2 <- E2 * asin2 + F2 * acos2 + basel2

    A2_1 <- (E1^2 + F1^2)
    A2_2 <- (E2^2 + F2^2)
    A2_1 * theta1 - A2_2 * theta2

    c(2 * E1, 2 * F1, 0, 0,
      - 2 * E2, - 2 * F2, 0, 0)
  }


  # Lower and upper bounds
  lb <- c(-Inf,-Inf,-Inf,0, -Inf, -Inf,-Inf, 0)
  ub <- c(Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf)
  #initial values

  ## Error in is.nloptr(ret) :
  #  If you want to use equality constraints, then you should use one of these algorithms NLOPT_LD_AUGLAG, NLOPT_LN_AUGLAG, NLOPT_LD_AUGLAG_EQ, NLOPT_LN_AUGLAG_EQ, NLOPT_GN_ISRES, NLOPT_LD_SLSQP

  local_opts <- list( "algorithm" = "NLOPT_LD_MMA", "xtol_rel" = 1.0e-15 )
  "local_opts" = local_opts
  opts <- list( "algorithm"= "NLOPT_LD_SLSQP",
                "xtol_rel"= 1.0e-15,
                "maxeval"= 160000,
                "local_opts" = local_opts,
                "print_level" = 0
                #"check_derivatives"=TRUE
  )

  res <- nloptr:nloptr( x0 = x_Ha,
                  eval_f = eval_f_list,
                  #eval_grad_f=eval_g,
                  lb = lb,
                  ub = ub,
                  #eval_g_ineq = eval_g_ineq,
                  eval_g_eq = eval_g_eq,
                  eval_jac_g_eq = eval_g_eq_jac,
                  opts = opts,
                  asin1=asin1,
                  acos1=acos1,
                  asin2=asin2,
                  acos2=acos2)

  #
  #x_Ha
  x_H0 <- res$solution

  l0 <- - eval_f_list(x_H0,asin1,acos1,asin2,acos2)$objective
  la <- - eval_f_list(x_Ha,asin1,acos1,asin2,acos2)$objective

  LR_stat <- -2*(l0-la)

  dfdiff <- 1
  if(!FN){
    pvalue <- pchisq(LR_stat,dfdiff,lower.tail = F)
  } else if(FN){
    r <- 1
    k <- 6
    n <- n1+n2
    Fstat <- (exp(LR_stat/n) - 1) * (n-k) / r
    pvalue <- pf(Fstat,df1 = r, df2 = n-k, lower.tail = F)
  } else{
    stop("FN has to be TRUE or FALSE")
  }

  amp_c <- sqrt(x_H0[1]^2 + x_H0[2]^2)
  amp_c2 <- sqrt(x_H0[5]^2 + x_H0[6]^2)

  res <- list(amp_1=A1, amp_2=A2, amp_c=amp_c,
              l0=l0,
              la=la,
              #df = dfdiff,
              stat=-2*(l0-la),
              pvalue=pvalue)
  return(res)
}


##' Finite sample/Large sample Likelihood ratio test for differential phase.
##'
##' Test differential phase of circadian curve fitting using likelihood ratio test
##' @title Likelihood ratio test for detecting differential phase.
##' @param tt1 Time vector of condition 1
##' @param yy1 Expression vector of condition 1
##' @param tt2 Time vector of condition 2
##' @param yy2 Expression vector of condition 2
##' @param period Period of the since curve. Default is 24.
##' @param type Type of likelihood ratio test to use, "FN" or "LS". Default is finite sample.
##' @return A list, see details below.
##' Formula 1: \eqn{yy = amp \times sin(2\pi/period \times (phase + tt)) + offset}
##' Formula 2: \eqn{yy = A \times sin(2\pi/period \times tt) + B * cos(2*pi/period * tt) + offset}
##' \item{phase_1}{Phase estimate of the 1st data, phase is restricted in (0, period)}
##' \item{phase_2}{Phase estimate of the 2nd data, phase is restricted in (0, period)}
##' \item{phase_c}{Phase estimate pooling all data together, phase is restricted in (0, period)}
##' \item{l0}{Log likelihood under the null (same variance between the two groups)}
##' \item{l1}{Log likelihood under the alternative (different variance between the two groups)}
##' \item{df}{Degree of freedom for the LR test}
##' \item{stat}{LR statistics}
##' \item{pvalue}{P-value from the LR test}
##' @author Caleb
##' @noRd
##' @examples
##' set.seed(32608)
##' n <- 50
##' tt1 <- runif(n,0,24)
##' Amp1 <- 2
##' Phase1 <- 6
##' Offset1 <- 3
##' yy1 <- Amp1 * sin(2*pi/24 * (tt1 + Phase1)) + Offset1 + rnorm(n,0,1)
##' tt2 <- runif(n,0,24)
##' Amp2 <- 3
##' Phase2 <- 5
##' Offset2 <- 2
##' yy2 <- Amp2 * sin(2*pi/24 * (tt2 + Phase2)) + Offset2 + rnorm(n,0,1)
##' LRTest_diff_phase(tt1, yy1, tt2, yy2)


LRTest_diff_phase <- function(tt1, yy1, tt2, yy2, period = 24,FN=TRUE){

  n1 <- length(tt1)
  stopifnot(n1 == length(yy1))
  n2 <- length(tt2)
  stopifnot(length(tt2) == length(yy2))

  #period <- 24
  w <- 2*pi/period

  fit1 <- fitSinCurve(tt1, yy1, period = period)
  fit2 <- fitSinCurve(tt2, yy2, period = period)

  A1 <- fit1$amp
  A2 <- fit2$amp

  phase1 <- fit1$phase
  phase2 <- fit2$phase

  if(phase2 - phase1 > period/2){
    phase2 <- phase2 - period
  } else if(phase1 - phase2 > period/2){
    phase1 <- phase1 - period
  }

  basal1 <- fit1$offset
  basal2 <- fit2$offset

  sigma2_1 <- 1/n1 * fit1$rss
  sigma2_2 <- 1/n2 * fit2$rss

  theta1 <- 1/sigma2_1
  theta2 <- 1/sigma2_2

  p1 <- c(A1, phase1, basal1, theta1)
  p2 <- c(A2, phase2, basal2, theta2)

  x_Ha <- c(p1, p2)

  eval_f_list <- function(x) {
    p1 <- x[1:4]
    p2 <- x[5:8]

    A1 <- p1[1]
    phase1 <- p1[2]
    basel1 <- p1[3]
    theta1 <- p1[4]
    asin1 <- sin(w * (tt1 + phase1) )
    acos1 <- cos(w * (tt1 + phase1) )
    yhat1 <- A1 * asin1 + basel1

    A2 <- p2[1]
    phase2 <- p2[2]
    basel2 <- p2[3]
    theta2 <- p2[4]
    asin2 <- sin(w * (tt2 + phase2) )
    acos2 <- cos(w * (tt2 + phase2) )
    yhat2 <- A2 * asin2 + basel2

    ll1_a <- log(theta1)/2
    ll1_b <- (yy1 - yhat1)^2 * theta1 / 2
    ll1 <- ll1_a - ll1_b

    ll2_a <- log(theta2)/2
    ll2_b <- (yy2 - yhat2)^2 * theta2 / 2
    ll2 <- ll2_a - ll2_b

    partial_A1 <- - theta1 * sum((yy1 - yhat1) * asin1)
    partial_phase1 <- - theta1 * A1 * w * sum((yy1 - yhat1) * acos1)
    partial_C1 <- - theta1 * sum(yy1 - yhat1)
    partial_theta1 <-  sum((yy1 - yhat1)^2)/2 - n1/2/theta1

    partial_A2 <- - theta2 * sum((yy2 - yhat2) * asin2)
    partial_phase2 <- - theta2 * A2 * w * sum((yy2 - yhat2) * acos2)
    partial_C2 <- - theta2 * sum(yy2 - yhat2)
    partial_theta2 <-  sum((yy2 - yhat2)^2)/2 - n2/2/theta2


    return( list( "objective" = -sum(ll1) - sum(ll2),
                  "gradient"  = c(partial_A1, partial_phase1, partial_C1, partial_theta1,
                                  partial_A2, partial_phase2, partial_C2, partial_theta2)
    )
    )
  }

  # Equality constraints
  eval_g_eq <- function(x)
  {
    p1 <- x[1:4]
    p2 <- x[5:8]

    phase1 <- p1[2]
    phase2 <- p2[2]

    phase1 - phase2
  }

  # Equality constraints
  eval_g_eq_jac <- function(x)
  {
    c(0, 1, 0, 0,
      0, -1, 0, 0)
  }


  # Lower and upper bounds
  lb <- c(0,-Inf,-Inf,0, 0, -Inf,-Inf, 0)
  ub <- c(Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf)
  #initial values

  ## Error in is.nloptr(ret) :
  #  If you want to use equality constraints, then you should use one of these algorithms NLOPT_LD_AUGLAG, NLOPT_LN_AUGLAG, NLOPT_LD_AUGLAG_EQ, NLOPT_LN_AUGLAG_EQ, NLOPT_GN_ISRES, NLOPT_LD_SLSQP

  local_opts <- list( "algorithm" = "NLOPT_LD_MMA", "xtol_rel" = 1.0e-15 )
  "local_opts" = local_opts
  opts <- list( "algorithm"= "NLOPT_LD_SLSQP",
                "xtol_rel"= 1.0e-15,
                "maxeval"= 160000,
                "local_opts" = local_opts,
                "print_level" = 0
                #"check_derivatives"=TRUE
  )

  res <- nloptr::nloptr( x0 = x_Ha,
                  eval_f = eval_f_list,
                  #eval_grad_f=eval_g,
                  lb = lb,
                  ub = ub,
                  #eval_g_ineq = eval_g_ineq,
                  eval_g_eq = eval_g_eq,
                  eval_jac_g_eq = eval_g_eq_jac,
                  opts = opts)

  #
  #x_Ha
  x_H0 <- res$solution

  l0 <- - eval_f_list(x_H0)$objective
  la <- - eval_f_list(x_Ha)$objective

  LR_stat <- -2*(l0-la)

  dfdiff <- 1
  if(!FN){
    pvalue <- pchisq(LR_stat,dfdiff,lower.tail = F)
  } else if(FN){
    r <- 1
    k <- 6
    n <- n1+n2
    Fstat <- (exp(LR_stat/n) - 1) * (n-k) / r
    pvalue <- pf(Fstat,df1 = r, df2 = n-k, lower.tail = F)
  } else{
    stop("FN has to be TRUE or FALSE")
  }

  phase_c <- x_H0[2]
  phase_c2 <- x_H0[6]

  res <- list(phase_1=phase1, phase_2=phase2, phase_c=phase_c,
              l0=l0,
              la=la,
              #df = dfdiff,
              stat=-2*(l0-la),
              pvalue=pvalue)
  return(res)
}


##' Finite sample/Large sample Likelihood ratio test for differential sigma square.
##'
##' Test differential sigma2 of circadian curve fitting using likelihood ratio test
##' @title Likelihood ratio test for detecting differential sigma square.
##' @param tt1 Time vector of condition 1
##' @param yy1 Expression vector of condition 1
##' @param tt2 Time vector of condition 2
##' @param yy2 Expression vector of condition 2
##' @param period Period of the since curve. Default is 24.
##' @param type Type of likelihood ratio test to use, "FN" or "LS". Default is finite sample.
##' @return A list, see details below.
##' Formula 1: \eqn{yy = amp \times sin(2\pi/period \times (phase + tt)) + offset}
##' Formula 2: \eqn{yy = A \times sin(2\pi/period \times tt) + B * cos(2*pi/period * tt) + offset}
##' \item{sigma2_1}{Variance estimate of the 1st data}
##' \item{sigma2_2}{Variance estimate of the 2nd data}
##' \item{sigma2_C}{Variance estimate pooling all data together}
##' \item{l0}{Log likelihood under the null (same variance between the two groups)}
##' \item{l1}{Log likelihood under the alternative (different variance between the two groups)}
##' \item{df}{Degree of freedom for the LR test}
##' \item{stat}{LR statistics}
##' \item{pvalue}{P-value from the LR test}
##' @author Caleb
##' @noRd
##' @examples
##' set.seed(32608)
##' n <- 50
##' tt1 <- runif(n,0,24)
##' Amp1 <- 2
##' Phase1 <- 6
##' Offset1 <- 3
##' yy1 <- Amp1 * sin(2*pi/24 * (tt1 + Phase1)) + Offset1 + rnorm(n,0,1)
##' tt2 <- runif(n,0,24)
##' Amp2 <- 3
##' Phase2 <- 5
##' Offset2 <- 2
##' yy2 <- Amp2 * sin(2*pi/24 * (tt2 + Phase2)) + Offset2 + rnorm(n,0,1)
##' LRTest_diff_sigma2(tt1, yy1, tt2, yy2)


LRTest_diff_sigma2 <- function(tt1, yy1, tt2, yy2, period = 24,FN=TRUE){
  n1 <- length(tt1)
  stopifnot(n1 == length(yy1))
  n2 <- length(tt2)
  stopifnot(length(tt2) == length(yy2))

  fit1 <- fitSinCurve(tt1, yy1, period = period)
  fit2 <- fitSinCurve(tt2, yy2, period = period)
  sigma2_1 <- 1/n1 * fit1$rss
  sigma2_2 <- 1/n2 * fit2$rss
  sigma2_C <- 1/(n1 + n2) * (sigma2_1 * n1 + sigma2_2 * n2)

  ## H0 (sigma2_C)
  l0_1 <- -n1/2*log(2*pi*sigma2_C) # - 1/(2*sigma2_C)*sum(residual_1^2)
  l0_2 <- -n2/2*log(2*pi*sigma2_C) # - 1/(2*sigma2_C)*sum(residual_2^2)
  l0 <- l0_1 + l0_2

  ## H1 (sigma2_1, sigma2_2)
  la_1 <- -n1/2*log(2*pi*sigma2_1) # - 1/(2*sigma2_1)*sum(residual_1^2)
  la_2 <- -n2/2*log(2*pi*sigma2_2) # - 1/(2*sigma2_2)*sum(residual_2^2)
  la <- la_1 + la_2

  dfdiff <- 1
  if(FN==FALSE){
    pvalue <- pchisq(-2*(l0-la),dfdiff,lower.tail = F)
  } else if(FN==TRUE){
    LR_stat <- -2*(l0-la)
    r <- 1
    k <- 6
    n <- n1+n2
    Fstat <- (exp(LR_stat/n) - 1) * (n-k) / r
    pvalue <- pf(Fstat,df1 = r, df2 = n-k, lower.tail = F)
  }


  res <- list(sigma2_1=sigma2_1, sigma2_2=sigma2_2, sigma2_c=sigma2_C,
              l0=l0,
              la=la,
              #df = dfdiff,
              stat=-2*(l0-la),
              pvalue=pvalue)
  return(res)
}



##' Finite sample/ large sample Likelihood ratio test for circadian pattern detection
##'
##' Test the significance of circadian curve fitting using finite sample likelihood ratio test
##' @title LR Test for detecting circadian pattern.
##' @param tt Time vector
##' @param yy Expression vector
##' @param period Period of the since curve. Default is 24.
##' @param type Type of Test, finite sample "FN" or large sample "LS", default is "FN".
##' @return A list of amp, phase, offset, sigma02, sigmaA2, l0, l1, df, stat, and pvalue.
##' Formula 1: \eqn{yy = amp \times sin(2\pi/period \times (phase + tt)) + offset}
##' Formula 2: \eqn{yy = A \times sin(2\pi/period \times tt) + B * cos(2*pi/period * tt) + offset}
##' \item{amp}{Amplitude based on formula 1}
##' \item{phase}{Phase based on formula 1, phase is restricted within (0, period)}
##' \item{peakTime}{Phase based on formula 1, peakTime is restricted within (0, period). phase + peakTime = period/4}
##' \item{offset}{Basal level(vertical shift) based on formula 1 or on formula 2}
##' \item{sigma02}{Variance estimate under the null (intercept only)}
##' \item{sigmaA2}{Variance estimate under the alternative (since curve fitting)}
##' \item{l0}{Log likelihood under the null (intercept only)}
##' \item{l1}{Log likelihood under the alternative (since curve fitting)}
##' \item{df}{Degree of freedom for the LR test}
##' \item{stat}{LR statistics}
##' \item{pvalue}{P-value from the LR test}
##' \item{R2}{Pseudo R2 defined as (tss - rss)/tss}
##' @author Caleb
##' @noRd
##' @examples
##' set.seed(32608)
##' n <- 50
##' tt <- runif(n,0,24)
##' Amp <- 2
##' Phase <- 6
##' Offset <- 3
##' yy <- Amp * sin(2*pi/24 * (tt + Phase)) + Offset + rnorm(n,0,1)
##' LRTest(tt, yy)

#model: y=A*sin(2*pi*x+B)+C
#y: a 1*n vector of data y
#A: estimated A^hat from fitCurve
#B: estimated B^hat from fitCurve
#C: estimated C^hat from fitCurve
#sigma0: sigma0^hat under H0
#sigmaA: sigmaA^hat under H1
#n: length of data y
#df0: df under H0
#df1: df under H1

LRTest <- function(tt,yy, period = 24,type=TRUE){
  fitCurveOut <- fitSinCurve(tt,yy,period=period)
  n <- length(yy)
  rss <- fitCurveOut$rss
  tss <- fitCurveOut$tss

  amp <- fitCurveOut$amp
  phase <- fitCurveOut$phase
  offset <- fitCurveOut$offset

  sigma02 <- 1/(n)*sum((yy-mean(yy))^2)
  sigmaA2 <- 1/(n)*sum((yy-amp*sin(2*pi/period*(tt+phase))-offset)^2)

  l0 <- -n/2*log(2*pi*sigma02)-1/(2*sigma02)*sum((yy-mean(yy))^2)
  l1 <- -n/2*log(2*pi*sigmaA2)-1/(2*sigmaA2)*sum((yy-amp*sin(2*pi/period*(tt+phase))-offset)^2)

  dfdiff <- (n-1)-(n-3)
  LR_stat <- -2*(l0-l1)

  if(type==FALSE){
    pvalue <- pchisq(LR_stat,dfdiff,lower.tail = F)
  }
  else if(type==TRUE){
    r <- 2
    k <- 3
    LR_stat <- (exp(LR_stat/n) - 1) * (n-k) / r
    pvalue <- pf(LR_stat,df1 = r, df2 = n-k, lower.tail = F)
  }
  R2 <- 1-rss/tss
  res <- list(
    amp = amp,
    phase = phase,
    peakTime = (6 - phase) %% period,
    offset = offset,
    sigma02=sigma02, sigmaA2=sigmaA2,
    l0=l0,
    l1=l1,
    stat=LR_stat,
    pvalue=pvalue,R2=R2)
  return(res)
}


##' Likelihood-based tests for circadian pattern detection.
##'
##' Test the significance of circadian curve fitting using likelihood-based tests.
##' @title Likelihood-based Tests for Detecting Circadian Pattern.
##' @param tt Time vector
##' @param yy Expression vector
##' @param period Period of the since curve. Default is 24.
##' @param method Testing methods can be "Wald" or "LR". Default is "LR".
##' @param FN Type of Test, finite sample if TRUE or large sample if FALSE. Default is TRUE.
##' @return A list of amp, phase, offset, sigma02, sigmaA2, l0, l1, df, stat, and pvalue.
##' Formula 1: \eqn{yy = amp * sin(2\pi/period * (phase + tt)) + offset}
##' Formula 2: \eqn{yy = A * sin(2\pi/period * tt) + B * cos(2*\pi/period * tt) + offset}
##' \item{amp}{Amplitude based on formula 1.}
##' \item{phase}{Phase based on formula 1, phase is restricted within (0, period).}
##' \item{peakTime}{Phase based on formula 1, peakTime is restricted within (0, period). phase + peakTime = period/4}
##' \item{offset}{Basal level (vertical shift) based on formula 1 or on formula 2.}
##' \item{sigma02}{Variance estimate under the null (intercept only).}
##' \item{sigmaA2}{Variance estimate under the alternative (since curve fitting).}
##' \item{l0}{Log likelihood under the null (intercept only).}
##' \item{l1}{Log likelihood under the alternative (since curve fitting).}
##' \item{stat}{Test statistic.}
##' \item{pvalue}{P-value from the test.}
##' \item{R2}{Pseudo R2 defined as (tss - rss)/tss.}
##' @author Zhiguang Huo, Haocheng Ding
##' @export
##' @examples
##' set.seed(32608)
##' n <- 50
##' tt <- runif(n,0,24)
##' Amp <- 2
##' Phase <- 6
##' Offset <- 3
##' yy <- Amp * sin(2*pi/24 * (tt + Phase)) + Offset + rnorm(n,0,1)
##' LR_rhythmicity(tt, yy, period=24, method="LR", FN=TRUE)
LR_rhythmicity <- function(tt,yy,period=24,method="LR",FN=TRUE){
  if(method=="Wald"){
    WaldTest(tt=tt, yy=yy, period = period, type=FN)
  }
  else if(method=="LR"){
    LRTest(tt=tt, yy=yy, period = period, type=FN)
  }
  else(("Please check your input! Method only supports 'Wald','LR','F' or 'Permutation'."))
}


##' F test for sin function.
##'
##' Test the significance of circadian curve fitting using F test.
##' @title F test for detecting circadian pattern
##' @param tt Time vector
##' @param yy Expression vector
##' @param period Period of the since curve. Default is 24.
##' @return A list of amp, phase, offset, peak, SST, SSE, R2.
##' Formula 1: \eqn{yy = amp \times sin(2\pi/period \times (phase + tt)) + offset}
##' Formula 2: \eqn{yy = A \times sin(2\pi/period \times tt) + B * cos(2*pi/period * tt) + offset}
##' \item{amp}{Amplitude based on formula 1}
##' \item{phase}{Phase based on formula 1, phase is restricted within (0, period)}
##' \item{offset}{Basal level(vertical shift) based on formula 1 or on formula 2}
##' \item{tss}{Total sum of square}
##' \item{rss}{Residual sum of square, SSE/n is the MLE of the variance sigma2}
##' \item{R2}{Pseudo R2 defined as (tss - rss)/tss}
##' @author Caleb
##' @noRd
##' @examples
##' set.seed(32608)
##' n <- 50
##' tt <- runif(n,0,24)
##' Amp <- 2
##' Phase <- 6
##' Offset <- 3
##' yy <- Amp * sin(2*pi/24 * (tt + Phase)) + Offset + rnorm(n,0,1)
##' FTest(tt, yy)

FTest <- function(tt,yy, period = 24){
  fitCurveOut <- fitSinCurve(tt,yy,period=period)
  n <- length(yy)
  rss <- fitCurveOut$rss
  tss <- fitCurveOut$tss
  amp <- fitCurveOut$amp
  phase <- fitCurveOut$phase
  offset <- fitCurveOut$offset

  fss <-tss-rss

  df1 <- 2
  df2 <- n-3
  fvalue <- fss/df1/(rss/df2)
  R2 <- 1-rss/tss
  pvalue <- pf(fvalue,df1,df2,lower.tail = F)

  return(list(amp=amp, phase=phase, offset=offset, rss=rss, tss=tss,R2=R2 ,df1 = df1, df2= df2, stat=fvalue,pvalue=pvalue))
}
