#' Estimate Root Diameters
#'
#' This function estimates root diameters from an input image using skeletonization and distance transform methods.
#' The input can be a file path, raster, image object, or array, which is converted to a binary image before processing.
#'
#' @param image_input A character string (file path), `SpatRaster`, `RasterBrick`, `RasterLayer`, `cimg`, `magick-image`, or array.
#'   The input image to process.
#' @param diagnostics Logical. If `TRUE`, enables diagnostic plots and logging. Default is `FALSE`.
#' @param skeleton_method Character. The method to use for skeletonization. Default is `"Guo-Hall"`.
#' @param layer Integer. Specifies which layer to use if the input is a multi-band image. Default is `2`.
#'
#' @details
#' The function works as follows:
#' - Converts the input image to a binary format (`cimg`).
#' - Applies a distance transform to compute the Euclidean distance for the foreground (root) pixels.
#' - Skeletonizes the binary image to identify root centerlines.
#' - Filters distance values to retain only those corresponding to the skeletonized regions.
#' - Computes diameter statistics, including quantiles, mean, and median diameters.
#'
#' The function supports various input formats and normalizes image values to the range [0, 1] if needed. It uses the `terra` package for raster operations and the `imager` package for image processing.
#'
#' @return A list containing:
#' \describe{
#'   \item{quantiles}{Numeric vector of diameter quantiles (10th to 100th percentile).}
#'   \item{mean_diameter}{Numeric. The mean root diameter.}
#'   \item{median_diameter}{Numeric. The median root diameter.}
#'   \item{diameters}{Numeric vector of all diameter values in the skeletonized regions.}
#'   \item{skeleton_rast}{`SpatRaster`. Binary raster mask of skeletonized regions.}
#'   \item{diameter_rast}{`SpatRaster`. Raster showing diameters in the skeletonized regions.}
#'   \item{distance_map_rast}{`SpatRaster`. Raster showing the distance transform values.}
#' }
#'
#' @examples
#' # Example usage:
#' image_path <- "path/to/image.tif"
#' result <- estimate_root_diameters(image_input = image_path, diagnostics = TRUE)
#'
#' # Access results:
#' print(result$mean_diameter)
#' terra::plot(result$skeleton_rast)
#'
#' @import terra
#' @import magick
#' @import imager
#' @import stats
#' @export
estimate_root_diameters <- function(image_input, diagnostics = FALSE,skeleton_method = "Guo-Hall",layer = 2) {
  # Helper function to convert input to cimg
  convert_to_cimg <- function(input) {
    if (is.character(input)) {
      if (!file.exists(input)) {
        stop("The provided image path does not exist.")
      }

      img0 <- terra::rast(input)
      if(dim(img0)[3] > 1){
        img0 <- img0[[layer]]
      }
      if(terra::global(img0,"max") > 1){
        img0 <- img0 / 255
      }

      img1 <- magick::image_read(terra::as.array(img0))
      return(imager::magick2cimg(img1))
    } else if (inherits(input, "SpatRaster")) {
      img0 <- input
      if(dim(img0)[3] > 1){
        img0 <- img0[[layer]]
      }
      if(terra::global(img0,"max") > 1){
        img0 <- img0 / 255
      }
      img1 <- magick::image_read(terra::as.array(img0))
      return(imager::magick2cimg(img1))
    } else if (inherits(input, "RasterBrick") || inherits(input, "RasterLayer")) {
      img0 <- terra::rast(input)
      if(dim(img0)[3] > 1){
        img0 <- img0[[layer]]
      }
      if(terra::global(img0,"max") > 1){
        img0 <- img0 / 255
      }
      img1 <- magick::image_read(terra::as.array(img0))
      return(imager::magick2cimg(img1))
    } else if (inherits(input, "cimg")) {
      return(input)
    } else if (inherits(input, "magick-image")) {
      return(imager::magick2cimg(input))
    } else if (is.array(input)) {
      return(imager::as.cimg(input))
    } else {
      stop("Unsupported input type. Provide a path, SpatRaster, RasterBrick, RasterLayer, cimg, magick image, or array.")
    }
  }

  # Convert input to cimg
  img_binary <- convert_to_cimg(image_input)


  # Distance transform function
  distance_transform <- function(binary_img) {
    # Compute the Euclidean distance transform
    dt <- imager::as.cimg(imager::distance_transform(imager::as.cimg(binary_img), value = 0))
    return(dt)
  }

  # Compute distance map and diameters
  distance_map <- distance_transform(binary_img = img_binary)
  diameters <- distance_map * 2

  # Convert to SpatRast
  Ds <- terra::rast(as.array(diameters))
  IM <- terra::rast(as.array(img_binary))


  # Skeletonize the binary image
  IMS <- skeletonize_image(IM,methods = skeleton_method)

  # Filter only root regions
  DsSKL <- Ds
  DsSKL[IMS == 0] <- NA

  # Threshold to create a binary skeleton mask
  skl <- DsSKL
  skl[skl > 0] <- 1

  # Compute statistics
  quantile_diameter <- stats::quantile(terra::values(DsSKL), probs = seq(0, 1, 0.1), na.rm = TRUE)
  mean_diameter <- terra::global(DsSKL, fun = function(x) stats::mean(x, na.rm = TRUE))$global
  median_diameter <- terra::global(DsSKL, fun = function(x) stats::median(x, na.rm = TRUE))$global


  # Visualization
  # par(mfrow = c(1, 3))
  # terra::plot(Ds, main = "Distance Map",col = map.pal("inferno",100)[6:100])
  # terra::plot(skl, main = "Skeletonized Diameters", add = FALSE, legend = FALSE, col = "black")

  # Histogram of diameters
  diameters <- terra::values(DsSKL, na.rm = TRUE)
  #hist(diameters, breaks = 50, main = "Histogram of Root Diameters", xlab = "Diameter", col = "skyblue", border = "white")

  return(list(
    quantiles = quantile_diameter,
    mean_diameter = mean_diameter,
    median_diameter = median_diameter,
    diameters <- diameters,
    skeleton_rast = skl,
    diameter_rast = DsSKL,
    distance_map_rast = Ds
  ))
}

