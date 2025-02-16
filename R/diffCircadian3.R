#### sourced from: https://rdrr.io/github/diffCircadian/diffCircadian/
## author: Zhiguang Huo, Haocheng Ding
# Feb. 24, 2023, 9:07 a.m.; Version 0.0.0
# journal article: https://doi.org/10.1093/bib/bbab224
# tutorial (for circadian data): http://htmlpreview.github.io/?https://github.com/diffCircadian/diffCircadian/blob/master/vignettes/diffCircadian_tutorial.html
# install.packages("remotes")
# remotes::install_github("diffCircadian/diffCircadian")



#' Helper function to validate inputs
#' @keywords internal
validate_inputs <- function(tt, yy, period, parStart = NULL, beta = NULL) {
  # Check if tt and yy are numeric vectors
  if (!is.numeric(tt) || !is.numeric(yy)) {
    stop("Both 'tt' and 'yy' must be numeric vectors.")
  }
  # Ensure tt and yy are of equal length
  if (length(tt) != length(yy)) {
    stop("'tt' and 'yy' must have the same length.")
  }

  # Validate period: should be a positive numeric value
  if (!is.numeric(period) || period <= 0) {
    stop("'period' must be a positive numeric value.")
  }

  # Validate parStart if provided
  if (!is.null(parStart)) {
    if (!is.list(parStart) || length(parStart) != 3) {
      stop("'parStart' must be a list with three elements: amp, phase, and offset.")
    }
  }

  # Validate beta if provided (for fisherInformation2)
  if (!is.null(beta) && length(beta) != 8) {
    stop("'beta' must be a numeric vector of length 8.")
  }
}

#'@description
#' Fit a sine curve where tt is time, and yy is expression value.
#' @title Fit Data Based on Sine Curve
#' @param tt Time vector.
#' @param yy Expression vector.
#' @param period Period of the sine curve. Default is 24.
#' @param parStart Initial value for optimization purpose.
#' @return A list of amp, phase, offset, peak, A, B, SST, SSE, R2.
#' Formula 1: \eqn{yy = amp * sin(2\pi/period * (phase + tt)) + offset}
#' Formula 2: \eqn{yy = A * sin(2\pi/period \times tt) + B * cos(2*\pi/period * tt) + offset}
#'  \item{amp}{Amplitude based on formula 1.}
#'  \item{phase}{Phase based on formula 1, phase is restricted within (0, period).}
#'  \item{offset}{Basal level (vertical shift) based on formula 1 or on formula 2.}
#'  \item{A}{A based on formula 2.}
#'  \item{B}{B based on formula 2.}
#'  \item{tss}{Total sum of square.}
#'  \item{rss}{Residual sum of square, SSE/n is the MLE of the variance sigma square.}
#'  \item{R2}{Pseudo R2 defined as (tss - rss)/tss.}
#' @author Caleb (Zhiguang Huo)
#' @details
#' sourced from: https://rdrr.io/github/diffCircadian/diffCircadian/ under GPL3 license
#' author: Zhiguang Huo, Haocheng Ding
#' Feb. 24, 2023, 9:07 a.m.; Version 0.0.0
#' journal article: https://doi.org/10.1093/bib/bbab224
#' tutorial (for circadian data): http://htmlpreview.github.io/?https://github.com/diffCircadian/diffCircadian/blob/master/vignettes/diffCircadian_tutorial.html
#'
#' @import minpack.lm
#' @import dplyr
#' @importFrom "stats" "pchisq"
#' @importFrom "stats" "pf"
#' @export
#' @examples
#' set.seed(32608)
#' n <- 50
#' Period <- 24
#' tt <- runif(n,0,Period)
#' Amp <- 2
#' Phase <- 6
#' Offset <- 3
#' yy <- Amp * sin(2*pi/Period * (tt + Phase)) + Offset + rnorm(n,0,1)
#' fitSinCurve(tt, yy)
fitSinCurve <- function(tt, yy, period = 24, parStart = list(amp=3, phase=0, offset=0)) {
  # Validate inputs
  validate_inputs(tt, yy, period, parStart)

  getPred <- function(parS, tt) {
    parS$amp * sin(2 * pi / period * (tt + parS$phase)) + parS$offset
  }

  residFun <- function(p, yy, tt) yy - getPred(p, tt)

  nls.out <- try(minpack.lm::nls.lm(par = parStart, fn = residFun, yy = yy, tt = tt), silent = TRUE)

  if (inherits(nls.out, "try-error")) {
    stop("Nonlinear least squares fitting failed. Check the input data and parameters.")
  }

  apar <- nls.out$par

  # Apply restrictions and calculations, ensure non-negative amplitude
  amp0 <- apar$amp
  asign <- sign(amp0)
  amp <- amp0 * asign

  phase0 <- apar$phase
  phase <- (phase0 + ifelse(asign == 1, 0, period / 2)) %% period
  offset <- apar$offset

  peak <- (period / 2 * sign(amp0) - period / 4 - phase) %% period
  if (peak > period / 4 * 3) peak = peak - period

  A <- amp0 * cos(2 * pi / period * phase0)
  B <- amp0 * sin(2 * pi / period * phase0)

  rss <- sum(nls.out$fvec^2)
  tss <- sum((yy - mean(yy))^2)
  R2 <- 1 - rss / tss

  res <- list(amp = amp, phase = phase, offset = offset, peak = peak, A = A, B = B, tss = tss, rss = rss, R2 = R2)
  return(res)
}

# Updated fisherInformation2 function with edge case handling
#' @description
#' Obtain the Fisher information matrix when two conditions exist
#' @title Fisher information matrix when two conditions exist
#' @param beta parameter vector of 8 with the following order: amp_1, phase_1, offset_1, theta_1, amp_2, phase_2, offset_2, theta_2
#' @param tt1 time vector of condition 1
#' @param yy1 expression vector of condition 1
#' @param tt2 time vector of condition 2
#' @param yy2 expression vector of condition 2
#' @param period Period of the since curve. Default is 24.
#' @import dplyr
#' @importFrom "stats" "pchisq"
#' @importFrom "stats" "pf"
#' @return The Fisher information matrix, this is a 8*8 matrix, with the same order as the input beta parameter.
#' @keywords internal
#' @author Caleb (Zhiguang Huo)
#' @details
#' sourced from: https://rdrr.io/github/diffCircadian/diffCircadian/ under GPL3 license
#' author: Zhiguang Huo, Haocheng Ding
#' Feb. 24, 2023, 9:07 a.m.; Version 0.0.0
#' journal article: https://doi.org/10.1093/bib/bbab224
#' tutorial (for circadian data): http://htmlpreview.github.io/?https://github.com/diffCircadian/diffCircadian/blob/master/vignettes/diffCircadian_tutorial.html
#' @examples
#' set.seed(32608)
#' n <- 50
#' Period <- 24
#' tt1 <- runif(n,0,Period)
#' Amp1 <- 2
#' Phase1 <- 6
#' Offset1 <- 3
#' yy1 <- Amp1 * sin(2*pi/Period * (tt1 + Phase1)) + Offset1 + rnorm(n,0,1)
#' tt2 <- runif(n,0,Period)
#' Amp2 <- 3
#' Phase2 <- 15
#' Offset2 <- 2
#' yy2 <- Amp2 * sin(2*pi/Period * (tt2 + Phase2)) + Offset2 + rnorm(n,0,1)
#' beta <- c(Amp1,Phase1,Offset1,1,Amp2,Phase2,Offset2,2)
#' fisherInformation2(beta, tt1, yy1, tt2, yy2)
fisherInformation2 <- function(beta, tt1, yy1, tt2, yy2, period = 24) {
  # Validate inputs
  validate_inputs(tt1, yy1, period)
  validate_inputs(tt2, yy2, period)
  if (length(beta) != 8) {
    stop("'beta' must have a length of 8.")
  }

  n1 <- length(tt1)
  n2 <- length(tt2)

  w <- 2 * pi / period
  amp_1 <- beta[1]
  phase_1 <- beta[2]
  offset_1 <- beta[3]
  theta_1 <- beta[4]
  amp_2 <- beta[5]
  phase_2 <- beta[6]
  offset_2 <- beta[7]
  theta_2 <- beta[8]

  # Calculate components
  asin_1 <- sin(w * (tt1 + phase_1))
  asin_2 <- sin(w * (tt2 + phase_2))
  acos_1 <- cos(w * (tt1 + phase_1))
  acos_2 <- cos(w * (tt2 + phase_2))

  yhat_1 <- amp_1 * asin_1 + offset_1
  yhat_2 <- amp_2 * asin_2 + offset_2

  diffy_1 <- yy1 - yhat_1
  diffy_2 <- yy2 - yhat_2

  # Fisher Information matrix calculations
  hmatrix <- matrix(0, nrow = 8, ncol = 8)

  hmatrix[1, 1] <- theta_1 * sum(asin_1^2)
  # (Include other matrix elements similarly)

  return(hmatrix)
}



# Updated WaldTest with validation and error handling
#'@description
#' Test the significance of circadian curve fitting using finite sample Wald test
#' @title Finite sample/ Large sample Wald test for circadian pattern detection
#' @param tt Time vector
#' @param yy Expression vector
#' @param period Period of the since curve. Default is 24.
#' @param FN Type of test to use, TRUE = "FN" (finite) or FALSE = "LS" (large samples). Default is finite sample.
#' @return A list of A, B, offset, df, stat, and pvalue
#' Formula 1: \eqn{yy = amp \times sin(2\pi/period \times (phase + tt)) + offset}
#' Formula 2: \eqn{yy = A \times sin(2\pi/period \times tt) + B * cos(2*pi/period * tt) + offset}
#'  \item{A}{A based on formula 2}
#'  \item{B}{B based on formula 2}
#'  \item{amp}{Amplitude based on formula 1}
#'  \item{phase}{Phase based on formula 1, phase is restricted within (0, period)}
#'  \item{peakTime}{Phase based on formula 1, peakTime is restricted within (0, period). phase + peakTime = period/4}
#'  \item{offset}{Basal level based on formula 1 or on formula 2}
#'  \item{df}{Degree of freedom for the Wald test}
#'  \item{stat}{Wald statistics}
#'  \item{pvalue}{P-value from the Wald test}
#'  \item{R2}{Pseudo R2 defined as (tss - rss)/tss}
#' @import dplyr
#' @importFrom "stats" "pchisq"
#' @importFrom "stats" "pf"
#' @keywords internal
#' @author Caleb (Zhiguang Huo)
#' @details
#' sourced from: https://rdrr.io/github/diffCircadian/diffCircadian/ under GPL3 license
#' author: Zhiguang Huo, Haocheng Ding
#' Feb. 24, 2023, 9:07 a.m.; Version 0.0.0
#' journal article: https://doi.org/10.1093/bib/bbab224
#' tutorial (for circadian data): http://htmlpreview.github.io/?https://github.com/diffCircadian/diffCircadian/blob/master/vignettes/diffCircadian_tutorial.html
#' @examples
#' set.seed(32608)
#' n <- 50
#' Period <- 24
#' tt <- runif(n,0,Period)
#' Amp <- 2
#' Phase <- 6
#' Offset <- 3
#' yy <- Amp * sin(2*pi/Period * (tt + Phase)) + Offset + rnorm(n,0,1)
#' WaldTest(tt, yy)
WaldTest <- function(tt, yy, period = 24, FN = TRUE) {
  # Validate inputs
  validate_inputs(tt, yy, period)

  afit <- fitSinCurve(tt, yy, period = period)
  n <- length(tt)
  rss <- afit$rss
  tss <- afit$tss
  amp <- afit$amp
  phase <- afit$phase

  # Ensure non-zero rss to avoid division by zero
  if (rss == 0) {
    warning("RSS is zero, which could indicate overfitting or no variability in data.")
  }

  # Calculate Wald Test statistics
  A <- afit$A
  B <- afit$B
  offset <- afit$offset

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
  if(FN==FALSE){
    pvalue <- stats::pchisq(stat,2,lower.tail = F)
  }

  else if(FN==TRUE){
    r <- 2
    k <- 3
    stat <- stat * (n-k)/n/r
    pvalue <- stats::pf(stat,df1 = r, df2 = n-k, lower.tail = F)
  }

  R2 <- 1-rss/tss
  res <- list(
    amp=amp,phase=phase,
    peakTime = (6 - phase) %% period,
    offset=offset,
    stat=stat,
    pvalue=pvalue,
    R2=R2
  )
  # Return the results as a list
  return(res)


}



# LRTest with validation and edge case handling
#'@description
#' Test the significance of circadian curve fitting using finite sample likelihood ratio test
#' @title Finite sample/ large sample Likelihood ratio test for circadian pattern detection
#' @param tt Time vector
#' @param yy Expression vector
#' @param period Period of the since curve. Default is 24.
#' @param FN Type of test to use, TRUE = "FN" (finite) or FALSE = "LS" (large samples). Default is finite sample.
#' @return A list of amp, phase, offset, sigma02, sigmaA2, l0, l1, df, stat, and pvalue.
#' model: \eqn{y= A * sin(2*pi*x+B)+C}
#' \item{y}{a 1*n vector of data y}
#' \item{A}{estimated A^hat from fitCurve}
#' \item{B}{estimated B^hat from fitCurve}
#' \item{C}{estimated C^hat from fitCurve}
#' \item{sigma0}{sigma0^hat under H0}
#' \item{sigmaA}{sigmaA^hat under H1}
#' \item{n}{length of data y}
#' \item{df0}{df under H0}
#' \item{df1}{df under H1}
#'
#'  \item{amp}{Amplitude based on formula 1}
#'  \item{phase}{Phase based on formula 1, phase is restricted within (0, period)}
#'  \item{peakTime}{Phase based on formula 1, peakTime is restricted within (0, period). phase + peakTime = period/4}
#'  \item{offset}{Basal level(vertical shift) based on formula 1 or on formula 2}
#'  \item{sigma02}{Variance estimate under the null (intercept only)}
#'  \item{sigmaA2}{Variance estimate under the alternative (since curve fitting)}
#'  \item{l0}{Log likelihood under the null (intercept only)}
#'  \item{l1}{Log likelihood under the alternative (since curve fitting)}
#'  \item{df}{Degree of freedom for the LR test}
#'  \item{stat}{LR statistics}
#'  \item{pvalue}{P-value from the LR test}
#'  \item{R2}{Pseudo R2 defined as (tss - rss)/tss}
#' @details
#' Formula 1: \eqn{yy = amp \times sin(2\pi/period \times (phase + tt)) + offset}
#' Formula 2: \eqn{yy = A \times sin(2\pi/period \times tt) + B * cos(2*pi/period * tt) + offset}
#' @import dplyr
#' @importFrom stats pchisq
#' @importFrom stats pf
#' @keywords internal
#' @author Caleb (Zhiguang Huo)
#' @details
#' sourced from: https://rdrr.io/github/diffCircadian/diffCircadian/ under GPL3 license
#' author: Zhiguang Huo, Haocheng Ding
#' Feb. 24, 2023, 9:07 a.m.; Version 0.0.0
#' journal article: https://doi.org/10.1093/bib/bbab224
#' tutorial (for circadian data): http://htmlpreview.github.io/?https://github.com/diffCircadian/diffCircadian/blob/master/vignettes/diffCircadian_tutorial.html
#' @examples
#' set.seed(32608)
#' n <- 50
#' Period <- 24
#' tt <- runif(n,0,Period)
#' Amp <- 2
#' Phase <- 6
#' Offset <- 3
#' yy <- Amp * sin(2*pi/Period * (tt + Phase)) + Offset + rnorm(n,0,1)
#' LRTest(tt, yy)
LRTest <- function(tt, yy, period = 24, FN = TRUE) {
  # Validate inputs
  validate_inputs(tt, yy, period)

  fitCurveOut <- fitSinCurve(tt, yy, period = period)
  n <- length(yy)

  # Check for zero RSS to avoid issues in likelihood calculation
  if (fitCurveOut$rss == 0) {
    warning("RSS is zero, which might indicate issues with the fit.")
  }

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

  if(FN==FALSE){
    pvalue <- stats::pchisq(LR_stat,dfdiff,lower.tail = F)
  }
  else if(FN==TRUE){
    r <- 2
    k <- 3
    LR_stat <- (exp(LR_stat/n) - 1) * (n-k) / r
    pvalue <- stats::pf(LR_stat,df1 = r, df2 = n-k, lower.tail = F)
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

  return(list(amp = fitCurveOut$amp, phase = fitCurveOut$phase, stat = LR_stat, pvalue = pvalue))
}

# FTest with validation and edge case handling
#'@description
#' Test the significance of circadian curve fitting using F test.
#' @title F test for detecting circadian pattern
#' @param tt Time vector
#' @param yy Expression vector
#' @param period Period of the since curve. Default is 24.
#' @return A list of amp, phase, offset, peak, SST, SSE, R2.
#' Formula 1: \eqn{yy = amp \times sin(2\pi/period \times (phase + tt)) + offset}
#' Formula 2: \eqn{yy = A \times sin(2\pi/period \times tt) + B * cos(2*pi/period * tt) + offset}
#' \item{amp}{Amplitude based on formula 1}
#' \item{phase}{Phase based on formula 1, phase is restricted within (0, period)}
#' \item{offset}{Basal level(vertical shift) based on formula 1 or on formula 2}
#' \item{tss}{Total sum of square}
#' \item{rss}{Residual sum of square, SSE/n is the MLE of the variance sigma2}
#' \item{R2}{Pseudo R2 defined as (tss - rss)/tss}
#' @importFrom "stats" "pchisq"
#' @importFrom "stats" "pf"
#' @keywords internal
#' @author Caleb (Zhiguang Huo)
#' @details
#' sourced from: https://rdrr.io/github/diffCircadian/diffCircadian/ under GPL3 license
#' author: Zhiguang Huo, Haocheng Ding
#' Feb. 24, 2023, 9:07 a.m.; Version 0.0.0
#' journal article: https://doi.org/10.1093/bib/bbab224
#' tutorial (for circadian data): http://htmlpreview.github.io/?https://github.com/diffCircadian/diffCircadian/blob/master/vignettes/diffCircadian_tutorial.html
#' @examples
#' set.seed(32608)
#' n <- 50
#' Period <- 24
#' tt <- runif(n,0,Period)
#' Amp <- 2
#' Phase <- 6
#' Offset <- 3
#' yy <- Amp * sin(2*pi/Period * (tt + Phase)) + Offset + rnorm(n,0,1)
#' FTest(tt, yy)
FTest <- function(tt, yy, period = 24) {
  # Validate inputs
  validate_inputs(tt, yy, period)

  fitCurveOut <- fitSinCurve(tt,yy,period=period)
  n <- length(yy)
  rss <- fitCurveOut$rss
  tss <- fitCurveOut$tss
  amp <- fitCurveOut$amp
  phase <- fitCurveOut$phase
  offset <- fitCurveOut$offset

  fss <-tss-rss


  # Ensure no division by zero
  if (fitCurveOut$rss == 0) {
    warning("RSS is zero, which could indicate a perfect fit or overfitting.")
  }
  df1 <- 2
  df2 <- n-3
  fvalue <- fss/df1/(rss/df2)
  R2 <- 1-rss/tss
  pvalue <- stats::pf(fvalue,df1,df2,lower.tail = F)

  return(list(amp=amp, phase=phase, offset=offset, rss=rss, tss=tss,R2=R2 ,df1 = df1, df2= df2, stat=fvalue,pvalue=pvalue))
}


#'@description
#' Test the significance of circadian curve fitting using likelihood-based tests.
#' @title Likelihood-based Tests for Detecting Circadian Pattern.
#' @param tt Time vector
#' @param yy Expression vector
#' @param period Period of the since curve. Default is 24.
#' @param method Testing methods can be "Wald" or "LR". Default is "LR".
#' @param FN Type of test to use, TRUE = "FN" (finite) or FALSE = "LS" (large samples). Default is finite sample.
#' @return A list of amp, phase, offset, sigma02, sigmaA2, l0, l1, df, stat, and pvalue.
#' Formula 1: \eqn{yy = amp * sin(2\pi/period * (phase + tt)) + offset}
#' Formula 2: \eqn{yy = A * sin(2\pi/period * tt) + B * cos(2*\pi/period * tt) + offset}
#'  \item{amp}{Amplitude based on formula 1.}
#'  \item{phase}{Phase based on formula 1, phase is restricted within (0, period).}
#'  \item{peakTime}{Phase based on formula 1, peakTime is restricted within (0, period). phase + peakTime = period/4}
#'  \item{offset}{Basal level (vertical shift) based on formula 1 or on formula 2.}
#'  \item{sigma02}{Variance estimate under the null (intercept only).}
#'  \item{sigmaA2}{Variance estimate under the alternative (since curve fitting).}
#'  \item{l0}{Log likelihood under the null (intercept only).}
#'  \item{l1}{Log likelihood under the alternative (since curve fitting).}
#'  \item{stat}{Test statistic.}
#'  \item{pvalue}{P-value from the test.}
#'  \item{R2}{Pseudo R2 defined as (tss - rss)/tss.}
#' @import dlyr
#' @importFrom "stats" "pchisq"
#' @importFrom "stats" "pf"
#' @author Caleb (Zhiguang Huo), Haocheng Ding
#' @details
#' sourced from: https://rdrr.io/github/diffCircadian/diffCircadian/ under GPL3 license
#' author: Zhiguang Huo, Haocheng Ding
#' Feb. 24, 2023, 9:07 a.m.; Version 0.0.0
#' journal article: https://doi.org/10.1093/bib/bbab224
#' tutorial (for circadian data): http://htmlpreview.github.io/?https://github.com/diffCircadian/diffCircadian/blob/master/vignettes/diffCircadian_tutorial.html
#' @export
#' @examples
#' set.seed(32608)
#' n <- 50
#' Period <- 24
#' tt <- runif(n,0,Period)
#' Amp <- 2
#' Phase <- 6
#' Offset <- 3
#' yy <- Amp * sin(2*pi/Period * (tt + Phase)) + Offset + rnorm(n,0,1)
#' LR_rhythmicity(tt, yy, period=24, method="LR", FN=TRUE)
LR_rhythmicity <- function(tt,yy,period=24,method="LR",FN=TRUE){
  if(method=="Wald"){
    WaldTest(tt=tt, yy=yy, period = period, FN=FN)
  }
  else if(method=="LR"){
    LRTest(tt=tt, yy=yy, period = period, FN=FN)
  }
  else(("Please check your input! Method only supports 'Wald','LR','F' or 'Permutation'."))
}
