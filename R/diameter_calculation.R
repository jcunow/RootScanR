#' Estimate Root Diameters
#'
#' This function estimates root diameters from an input image using skeletonization and distance transform methods.
#' The input can be a file path, raster, image object, or array, which is converted to a binary image before processing.
#'
#' @param img A character string (file path), `SpatRaster`, `RasterBrick`, `RasterLayer`, `cimg`, `magick-image`, or array.
#'   The input image to process.
#' @param diagnostics Logical. If `TRUE`, enables diagnostic plots and logging. Default is `FALSE`.
#' @param skeleton_method Character. The method to use for skeletonization. Default is `"Guo-Hall"`.
#' @param select.layer Integer. Specifies which layer to use if the input is a multi-band image. Default is `2`.
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
#' data(seg_Oulanka2023_Session01_T067)
#' result <- root.diameters(
#' img = seg_Oulanka2023_Session01_T067,
#' skeleton_method = "MAT", select.layer = 2,
#' diagnostics = TRUE)
#'
#' # Access results:
#' print(result$mean_diameter)
#' terra::plot(result$skeleton_rast)
#'
#' @export
root.diameters <- function(img, diagnostics = FALSE,skeleton_method = "GuoHall",select.layer = 2) {

  # flexible input
  img <- load_flexible_image(img, select.layer = select.layer, output_format = "cimg", normalize = TRUE  )


  # Distance transform function
  distance_transform <- function(img) {
    # Compute the Euclidean distance transform
    dt <- imager::as.cimg(imager::distance_transform(imager::as.cimg(img), value = 0))
    return(dt)
  }

  # Compute distance map and diameters
  distance_map <- distance_transform(img = img)
  diameters <- distance_map * 2

  # Convert to SpatRast
  Ds <- terra::rast(as.array(diameters))
  IM <- terra::rast(as.array(img))


  # Skeletonize the binary image
  IMS <- convert_to_spatrast( skeletonize_image(IM,methods = skeleton_method))

  # Filter only root regions
  DsSKL <- Ds
  DsSKL[IMS == 0] <- NA

  # Threshold to create a binary skeleton mask
  skl <- DsSKL
  skl[skl > 0] <- 1

  # Compute statistics
  quantile_diameter <- stats::quantile(terra::values(DsSKL), probs = seq(0, 1, 0.1), na.rm = TRUE)
  mean_diameter <- terra::global(DsSKL, fun = function(x) base::mean(x, na.rm = TRUE))$global
  median_diameter <- terra::global(DsSKL, fun = function(x) stats::median(x, na.rm = TRUE))$global


  # Histogram of diameters
  diameters <- terra::values(DsSKL, na.rm = TRUE)


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

