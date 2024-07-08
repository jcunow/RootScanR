## performs an affine transform on a set of images
affine <- function (x, m, filter = c("bilinear", "none"), output.dim, bg.col = "black", antialias = TRUE) {
  ## check arguments
  validImage(x)
  if ( !is.matrix(m) || !identical(dim(m), c(3L, 2L)) ) stop("'m' must be a 3x2 matrix")
  if ( any(is.na(m)) ) stop("'m' shouldn't contain any NAs")

  filter <- switch(match.arg(filter), none=0L, bilinear=1L)

  # dimensions of output image
  d = dim(x)
  if (!missing(output.dim)) {
    if( length(output.dim)!=2L || !is.numeric(output.dim) ) stop("'output.dim' must be a numeric vector of length 2")
    d[1:2] = as.integer(round(output.dim))
  }

  ## inverse of the transformation matrix
  ## this is needed because in C we iterate over (x',y') for which we calculate (x,y)
  m <- solve(cbind(m, c(0, 0, 1)))

  ## image background
  nf = numberOfFrames(x, "render")

  bg <- Image(rep_len(bg.col, nf), c(1L, 1L, d[-c(1,2)]), colorMode(x))

  .Call(C_affine, castImage(x), d, castImage(bg), m, filter, isTRUE(antialias))
}


resize <- function(x, w, h, output.dim = c(w, h), output.origin = c(0, 0), antialias = FALSE, ...) {
  ## check arguments
  if ( missing(h) && missing(w) ) stop("either 'w' or 'h' must be specified")
  if ( length(output.origin)!=2L || !is.numeric(output.origin) ) stop("'output.origin' must be a numeric vector of length 2")

  d = dim(x)[1:2]
  if ( missing(w) ) w <- round(h*d[1]/d[2])
  if ( missing(h) ) h <- round(w*d[2]/d[1])
  if ( missing(output.dim) ) output.dim <- c(w, h)
  ratio <- c(w, h)/d

  m <- matrix(c(ratio[1], 0, (1-ratio[1]) * output.origin[1],
                0, ratio[2], (1-ratio[2]) * output.origin[2]), 3L, 2L)

  affine(x = x, m = m, output.dim = output.dim, antialias = antialias, ...)
}
