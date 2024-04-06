#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector TubeColorationCpp(List img, double r = 0.21, double g = 0.72, double b = 0.07) {
  NumericVector vr = img[0];
  NumericVector vg = img[1];
  NumericVector vb = img[2];

  NumericVector bright = vr + vg + vb;
  NumericVector lumGray = vr * r + vg * g + vb * b;

  double meanBright = mean(bright);
  double meanLum = mean(lumGray);
  double rcc = mean(vr / bright);
  double gcc = mean(vg / bright);
  double bcc = mean(vb / bright);

  return NumericVector::create(rcc, gcc, bcc, meanBright, meanLum);
}

/*** R
Tube_coloration <- function(img, r = 0.21, g = 0.72, b = 0.07) {
  vr <- values(img[[1]])
  vg <- values(img[[2]])
  vb <- values(img[[3]])

  bright <- vr + vg + vb
  lum.gray <- vr * r + vg * g + vb * b

  mean.bright <- mean(bright, na.rm = TRUE)
  mean.lum <- mean(lum.gray, na.rm = TRUE)

  rcc <- mean(vr / bright, na.rm = TRUE)
  gcc <- mean(vg / bright, na.rm = TRUE)
  bcc <- mean(vb / bright, na.rm = TRUE)

  return(c(rcc, gcc, bcc, mean.bright, mean.lum))
}

# Example usage:
img <- list(vr = c(0.1, 0.2, 0.3), vg = c(0.2, 0.3, 0.4), vb = c(0.3, 0.4, 0.5))
TubeColorationCpp(img)
*/
