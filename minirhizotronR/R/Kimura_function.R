#' Kimura Root Length
#'
#' Estimates the root length. Adapted from Kimura, K., Kikuchi, S., & Yamasaki, S. I. (1999). Accurate root length measurement by image analysis. Plant and Soil, 216(1), 117-127.
#' @param image binarized image with root pixels = 1 and background = 0
#'
#' @return root length in px
#' @export
#'
#' @examples kimura_length(x) = 20
kimura_length = function(image){
  compute_orthogonal_connections = function(image){
    # kernel definition
    k0 = matrix(c(1,1,NA), nrow = 1, ncol = 3)
    k1 = matrix(c(1,1,NA), nrow = 3, ncol = 1)
    r0 = terra::focal(image,w = k0, pad =F,fun = function(x)sum(x,na.rm=T))
    r1 = terra::focal(image,w = k1, pad =F,fun = function(x)sum(x,na.rm=T))
    return(sum((r0 == 2) | (r1 == 2)))
  }
  compute_diagonal_connections = function(image){
    # kernel definition
    k0 = matrix(c(0,1,0,1,0,0,0,0,0), nrow = 3, ncol = 3)
    k1 = matrix(c(1,0,0,0,1,0,0,0,0), nrow = 3, ncol = 3)
    r0 = terra::focal(image ,w = k0, fun = function(x)sum(x,na.rm=T))
    r1 = terra::focal(image ,w = k1,fun = function(x)sum(x,na.rm=T))
    return(sum((r0 == 2) | (r1 == 2)))
  }

  N_o =  compute_orthogonal_connections(image) %>% cellStats("sum",na.rm = T)
  N_d = compute_diagonal_connections(image) %>% cellStats("sum",na.rm = T)
  return(( N_d**2 + (N_d + N_o/2)**2 )**0.5   + N_o/2)
}
#
## kimura sum
#' Kimura sum
#'
#' This function takes the orthogonal pixels and diagonal pixels to derive total root length
#' @param orth.img orthogonal pixels
#' @param diag.img diagonal pixels
#'
#' @return Kimura Root Length in px
#' @export
#'
#' @examples kimura.sum(orth,diag) = 20
kimura.sum = function(orth.img = orth.img,diag.img = diag.img){
  ( diag.img**2 + (diag.img + orth.img/2)**2 )**0.5   + orth.img/2
}
