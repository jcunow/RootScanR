#' Rotate simple features for 3D layers
#' Rotates a simple features layer using a shear matrix transformation on the
#' \code{geometry} column. This can get nice for visualisation and works with
#' points, lines and polygons. Sourced from: https://www.mzes.uni-mannheim.de/socialsciencedatalab/article/geospatial-data/
#'
#' @param data an object of class \code{sf}
#' @param x_add integer; x value to move geometry in space
#' @param y_add integer; x value to move geometry in space
#' @param rotation.angle the angle of rotation
#'
#' @importFrom magrittr %>%
#'
#' @return raster output with rotation
#' @export
#'
#' @examples rotated.raster = RotateRaster(raster,rotation.angle = 45)
RotateRaster <- function(data, x_add = 0, y_add = 0, rotation.angle = 45) {

  if (!any(class(data) %in% c("sf", "sfg"))) {
    data <- stars::st_as_stars(data)
    data <- sf::st_as_sf(data)
  }

  shear_matrix <- function (x) {
    matrix(c(1, 1, 1, 1), 2, 2)
  }

  rotate_matrix <- function(x) {
    matrix(c(cos(x), sin(x), -sin(x), cos(x)), 2, 2)
  }

  data2 = data %>%
    dplyr::mutate(
      geometry =
        .$geometry *
        #shear_matrix() *
        rotate_matrix((pi/rotation.angle)) + c(x_add, y_add)
    )
   tilt.r = stars::st_rasterize(data2)
   tilt.r = as(tilt.r,"Raster")

  return(tilt.r)
}
