#' Filling the gaps in root pics
#'
#' @param img input Image
#' @param max_hole_size by default Inf, but you probably want to adjust it
#' @param select.layer Integer. Specifies which layer to use if the input is a multi-band image. Default is `2`.
#' @param plotting should the result be plotted?
#'
#' @return SpatRast input image but with filled holes
#' @export
#'
#' @examples
#' data(seg_Oulanka2023_Session01_T067)
#' fill_holes(terra::rast(seg_Oulanka2023_Session01_T067),
#' select.layer = 2, max_hole_size = 200, plotting = FALSE)
fill_holes <- function(img, max_hole_size = Inf,select.layer = 2, plotting = TRUE) {


  # flexible input
  img <- load_flexible_image(img, select.layer = select.layer, output_format = "spatrast", normalize = TRUE  )



  # Identify holes using clump (connected components)
  # First invert the image so holes become objects
  inverted <- !img

  # Label connected components
  clumps <- terra::patches(inverted, directions=8, allowGaps = T, zeroAsNA = T)

  # Find border components (background)
  border_vals <- unique(c(
    clumps[1,],              # top border
    clumps[nrow(clumps),],   # bottom border
    clumps[,1],              # left border
    clumps[,ncol(clumps)]    # right border
  ))
  border_vals <- border_vals[!is.na(border_vals)]

  # Create a mask of non-border components (holes)
  holes <- !is.element(clumps, border_vals) & clumps > 0

  if (max_hole_size < Inf) {
    # Calculate size of each component
    freq <- terra::freq(clumps)
    small_holes <- freq$value[freq$count <= max_hole_size]
    holes <- holes & is.element(clumps, small_holes)
  }



  # Fill the holes
  result <- img | holes

  if(plotting){
    graphics::par(mfrow = c(1,1))
    terra::plot(img, main = "Original Image")
    terra::plot(result-img, main = "Only the filled holes")
  }

  return(result)
}


#' Analyse hole characteristica
#'
#' @param img Segmented Root Image Input
#' @param select.layer Integer. Specifies which layer to use if the input is a multi-band image. Default is `2`.
#'
#' @return a summary statistics
#' @export
#'
#' @examples
#' data(seg_Oulanka2023_Session01_T067)
#' im = terra::rast(seg_Oulanka2023_Session01_T067)
#' analyze_holes(im)
analyze_holes <- function(img, select.layer  =2) {

  # flexible input
  img <- load_flexible_image(img, select.layer = select.layer, output_format = "spatrast", normalize = TRUE  )


  # Identify holes
  inverted <- !img
  clumps <- terra::patches(inverted, directions=8, allowGaps = T, zeroAsNA = T)

  # Find border components
  border_vals <- unique(c(
    clumps[1,], clumps[nrow(clumps),],
    clumps[,1], clumps[,ncol(clumps)]
  ))
  border_vals <- border_vals[!is.na(border_vals)]

  # Get frequencies of non-border components
  freq <- terra::freq(clumps)
  hole_sizes <- freq$count[!freq$value %in% border_vals & freq$value > 0]

  if (length(hole_sizes) == 0) {
    return(list(
      n_holes = 0,
      sizes = numeric(0),
      summary = data.frame(
        min_size = NA,
        max_size = NA,
        mean_size = NA,
        median_size = NA
      )
    ))
  }



  summary_stats <- data.frame(
    min_size = min(hole_sizes),
    max_size = max(hole_sizes),
    mean_size = mean(hole_sizes),
    median_size = stats::median(hole_sizes),
    percentile_size = stats::quantile(hole_sizes)
  )

  return(list(
    n_holes = length(hole_sizes),
    sizes = hole_sizes,
    summary = summary_stats
  ))
}


# # Example usage
# example_usage <- function() {
#   # Create sample data
#   m <- matrix(1, nrow=400, ncol=400)
#   m[4:350, 4:350] <- 1
#   m[6:7,5] <- 0         # 1 px
#   m[8:9, 8:9] <- 0     # small hole (4 pixels)
#   m[20:25,20:25] <- 0  # 25 px
#   m[50:100, 50:100] <- 0 # medium hole (16 pixels)
#   m[200:330, 200:330] <- 0 # large hole (36 pixels)
#
#   r <- rast(m)
#
#   # Analyze holes
#   analysis <- analyze_holes(r)
#   print("Hole Analysis:")
#   print(analysis$summary)
#
#   # Fill holes with different thresholds
#   filled_ss <- fill_holes(r, max_hole_size = 2)
#   filled_small <- fill_holes(r, max_hole_size = 50)
#   filled_ms <- fill_holes(r, max_hole_size = 2601)
#   filled_medium <- fill_holes(r, max_hole_size = 5000)
#   filled_all <- fill_holes(r,max_hole_size = 17161)
#
#   # Plot results
#   par(mfrow=c(2,3))
#   plot(r, main="Original")
#   plot(filled_ss, main = "1px hole filled")
#   plot(filled_small, main="Small holes filled (<4 pixels)")
#   plot(filled_ms, main = "medium hole fileld 24 px")
#   plot(filled_medium, main="Medium holes filled (<25 pixels)")
#   plot(filled_all, main="Big holes filled (<500)")
#
#   return(list(
#     original = r,
#     filled_small = filled_small,
#     filled_medium = filled_medium,
#     filled_all = filled_all,
#     analysis = analysis
#   ))
# }

