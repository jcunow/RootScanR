test_convert_function <- function() {
  library(testthat)
  library(imager)
  library(magick)
  library(terra)
  library(raster)

  # Define supported input and output formats
  input_formats <- list(
    "matrix" = matrix(runif(100), 10, 10),
    "array" = array(runif(1000), dim = c(10, 10, 10)),
    "cimg" = imager::as.cimg(array(runif(1000), dim = c(10, 10, 1, 10))),
    "magick-image" = magick::image_blank(10, 10),
    "RasterLayer" = raster::raster(matrix(runif(100), 10, 10)),
    "RasterBrick" = raster::brick(raster(matrix(runif(100), 10, 10)), raster(matrix(runif(100), 10, 10))),
    "SpatRaster" = terra::rast(matrix(runif(100), 10, 10))
  )

  output_formats <- c("SpatRaster", "RasterLayer", "RasterBrick", "cimg", "magick-image", "matrix", "array")

  for (input_name in names(input_formats)) {
    for (output_format in output_formats) {
      test_that(paste("Convert", input_name, "to", output_format), {

        select.layer = NULL
        expect_error(
          result <- load_flexible_image(input_formats[[input_name]], output_format = output_format, select.layer = select.layer, normalize = TRUE, binarize = FALSE),
          NA  # Expect no error
        )
      })
    }
  }
}

test_convert_function()














