#### Image Texture

# uses the glcm package
#' Texture corresponds to glmc::glmc()
#'
#' @param img.color image with three color chanels. Will be converted in grayscale.
#' @param grays number of gray shades. Documentation is lacking, see: ?glmc::glmc()
#' @param window convolution window size e.g., c(3,3)
#' @param metrics texture metrics based on illumination differences. see:: ?glmc::glmc() for available methods
#'
#' @return raster with each layer correspond to on metric
#' @export
#'
#' @examples texture(img.color, 7, c(9,9), metrics = "second_moment")
texture = function(img.color,grays = 7, window = c(9,9), metrics = c("variance","second_moment") ){
  img.gray = (img.color[[1]]*0.21 + img.color[[2]]*0.72 + img.color[[3]]*0.07)/255
  tx.im  =  glcm::glcm(img.gray,
                       n_grey = grays,
                       statistics = metrics)
  }
