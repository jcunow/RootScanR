% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/diffCircadian3.R
\name{FTest}
\alias{FTest}
\title{F test for detecting circadian pattern}
\usage{
FTest(tt, yy, period = 24)
}
\arguments{
\item{tt}{Time vector}

\item{yy}{Expression vector}

\item{period}{Period of the since curve. Default is 24.}
}
\value{
A list of amp, phase, offset, peak, SST, SSE, R2.
Formula 1: \eqn{yy = amp \times sin(2\pi/period \times (phase + tt)) + offset}
Formula 2: \eqn{yy = A \times sin(2\pi/period \times tt) + B * cos(2*pi/period * tt) + offset}
\item{amp}{Amplitude based on formula 1}
\item{phase}{Phase based on formula 1, phase is restricted within (0, period)}
\item{offset}{Basal level(vertical shift) based on formula 1 or on formula 2}
\item{tss}{Total sum of square}
\item{rss}{Residual sum of square, SSE/n is the MLE of the variance sigma2}
\item{R2}{Pseudo R2 defined as (tss - rss)/tss}
}
\description{
Test the significance of circadian curve fitting using F test.
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
FTest(tt, yy)
}
\author{
Caleb (Zhiguang Huo)
}
\keyword{internal}
