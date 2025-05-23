% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/diffCircadian3.R
\name{WaldTest}
\alias{WaldTest}
\title{Finite sample/ Large sample Wald test for circadian pattern detection}
\usage{
WaldTest(tt, yy, period = 24, FN = TRUE)
}
\arguments{
\item{tt}{Time vector}

\item{yy}{Expression vector}

\item{period}{Period of the since curve. Default is 24.}

\item{FN}{Type of test to use, TRUE = "FN" (finite) or FALSE = "LS" (large samples). Default is finite sample.}
}
\value{
A list of A, B, offset, df, stat, and pvalue
Formula 1: \eqn{yy = amp \times sin(2\pi/period \times (phase + tt)) + offset}
Formula 2: \eqn{yy = A \times sin(2\pi/period \times tt) + B * cos(2*pi/period * tt) + offset}
 \item{A}{A based on formula 2}
 \item{B}{B based on formula 2}
 \item{amp}{Amplitude based on formula 1}
 \item{phase}{Phase based on formula 1, phase is restricted within (0, period)}
 \item{peakTime}{Phase based on formula 1, peakTime is restricted within (0, period). phase + peakTime = period/4}
 \item{offset}{Basal level based on formula 1 or on formula 2}
 \item{df}{Degree of freedom for the Wald test}
 \item{stat}{Wald statistics}
 \item{pvalue}{P-value from the Wald test}
 \item{R2}{Pseudo R2 defined as (tss - rss)/tss}
}
\description{
Test the significance of circadian curve fitting using finite sample Wald test
}
\details{
sourced from: https://rdrr.io/github/diffCircadian/diffCircadian/ under GPL3 license
author: Zhiguang Huo, Haocheng Ding
Feb. 24, 2023, 9:07 a.m.; Version 0.0.0
journal article: https://doi.org/10.1093/bib/bbab224
tutorial (for circadian data): http://htmlpreview.github.io/?https://github.com/diffCircadian/diffCircadian/blob/master/vignettes/diffCircadian_tutorial.html
}
\examples{
set.seed(32608)
n <- 50
Period <- 24
tt <- runif(n,0,Period)
Amp <- 2
Phase <- 6
Offset <- 3
yy <- Amp * sin(2*pi/Period * (tt + Phase)) + Offset + rnorm(n,0,1)
WaldTest(tt, yy)
}
\author{
Caleb (Zhiguang Huo)
}
\keyword{internal}
