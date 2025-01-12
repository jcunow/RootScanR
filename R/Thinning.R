### thinning algorithm

#' Prepare Image for Processing (Internal)
#'
#' This internal function ensures that the input image is in binary format and a matrix structure.
#' It accepts matrices, data frames, or `SpatRaster` objects and standardizes the format for further processing.
#'
#' @param img A matrix, data frame, or `SpatRaster` object representing the image to be processed.
#' @param layer Integer indicating the layer to use if `img` is a multi-layer `SpatRaster`. Default is 2.
#'
#' @details
#' - If `img` is a `SpatRaster`, it will be converted to a matrix. For multi-layer rasters, the specified `layer` is used.
#' - If `img` is a data frame, it will be converted to a matrix.
#' - Ensures the matrix contains binary values (`1` for foreground and `0` for background).
#' - Validates that the matrix contains at least one foreground pixel (`1`).
#'
#' @return A binary matrix with dimensions corresponding to the input image.
#'
#' @keywords internal
#'
#' @examples
#' \dontrun{
#' # Example usage
#' raster <- terra::rast(matrix(c(0, 1, 0, 1), nrow = 2))
#' binary_image <- prepare_image(raster)
#' }
prepare_image <- function(img,layer = 2) {
  if (inherits(img, "SpatRaster")) {
    # Convert SpatRaster to matrix
    if(dim(img)[3]>1){
      img = img[[layer]]
    }
    img <- terra::as.array(img)
  } else if (is.data.frame(img)) {
    # Convert data.frame to matrix
    img <- as.matrix(img)
  } else if (!is.matrix(img)) {
    stop("Unsupported input type. Please provide a matrix, data.frame, or SpatRaster.")
  }

  # Ensure numeric type
  img <- matrix(as.numeric(img), nrow = nrow(img))

  # Force binary values (1 for foreground, 0 for background)
  img[abs(img - 1) < 1e-10] <- 1
  img[abs(img) < 1e-10] <- 0

  # Check if the image contains any foreground pixels
  if (sum(img == 1) == 0) {
    stop("No foreground pixels found after type conversion. Ensure the input is binary.")
  }

  return(img)
}


#' Thin Binary Image using Zhang-Suen Algorithm (Internal)
#'
#' This internal function performs image thinning using the Zhang-Suen thinning algorithm.
#' It reduces binary images to their skeleton while preserving the structure and connectivity of the foreground pixels.
#'
#' @param img A matrix, data frame, or `SpatRaster` object representing the binary image to be thinned.
#' @param verbose Logical. If `TRUE`, prints diagnostic information such as iteration progress and pixel removal counts. Default is `TRUE`.
#' @param layer Integer indicating the layer to use if `img` is a multi-layer `SpatRaster`. Default is 2.
#'
#' @details
#' - The function first prepares the image using the \code{\link{prepare_image}} function, ensuring binary matrix format.
#' - Thinning is performed iteratively in two subiterations per cycle:
#'   1. Identifying pixels to be removed based on Zhang-Suen conditions (first subiteration).
#'   2. Refining removal decisions in the second subiteration.
#' - The algorithm continues until no pixels are removed in an iteration or a maximum number of iterations is reached (default: 1000).
#'
#' @return A binary matrix representing the thinned image (skeleton).
#'
#' @keywords internal
#'
#' @examples
#' \dontrun{
#' # Example usage
#' raster <- terra::rast(matrix(c(0, 1, 1, 0, 0, 1, 1, 0), nrow = 4))
#' thinned_image <- thin_image_zhangsuen(raster, verbose = TRUE)
#' }
thin_image_zhangsuen <- function(img, verbose = TRUE,layer = 2) {
    # flexible input
  img = prepare_image(img,layer = layer)

  if(verbose){
    # Initial diagnostics
    cat("Image dimensions:", nrow(img), "x", ncol(img), "\n")
    cat("Initial foreground pixels:", sum(img == 1), "\n")
  }


  # Force binary values
  img[abs(img - 1) < 1e-10] <- 1
  img[abs(img) < 1e-10] <- 0

  if(sum(img == 1) == 0) {
    stop("No foreground pixels found after type conversion")
  }

  count_transitions <- function(p) {
    neighbors <- c(p[2], p[3], p[4], p[5], p[6], p[7], p[8], p[9], p[2])
    transitions <- 0
    for(i in 1:8) {
      if(neighbors[i] == 0 && neighbors[i + 1] == 1) transitions <- transitions + 1
    }
    return(transitions)
  }

  get_neighbors <- function(img, i, j) {
    p <- rep(0, 9)
    p[1] <- img[i, j]  # P1: Center point
    if (i-1 >= 1) p[2] <- img[i-1, j]    # P2: North
    if (i-1 >= 1 && j+1 <= ncol(img)) p[3] <- img[i-1, j+1]  # P3: Northeast
    if (j+1 <= ncol(img)) p[4] <- img[i, j+1]  # P4: East
    if (i+1 <= nrow(img) && j+1 <= ncol(img)) p[5] <- img[i+1, j+1]  # P5: Southeast
    if (i+1 <= nrow(img)) p[6] <- img[i+1, j]  # P6: South
    if (i+1 <= nrow(img) && j-1 >= 1) p[7] <- img[i+1, j-1]  # P7: Southwest
    if (j-1 >= 1) p[8] <- img[i, j-1]  # P8: West
    if (i-1 >= 1 && j-1 >= 1) p[9] <- img[i-1, j-1]  # P9: Northwest
    return(p)
  }

  max_iterations <- 1000
  iteration_count <- 0
  any_changes <- TRUE

  while(any_changes && iteration_count < max_iterations) {
    iteration_count <- iteration_count + 1
    any_changes <- FALSE
    pixels_removed <- 0

    # First subiteration
    to_delete <- matrix(FALSE, nrow=nrow(img), ncol=ncol(img))

    for(i in 2:(nrow(img)-1)) {
      for(j in 2:(ncol(img)-1)) {
        if(img[i,j] == 1) {
          p <- get_neighbors(img, i, j)
          B <- sum(p[-1])
          A <- count_transitions(p)

          if(B >= 2 && B <= 6 && A == 1) {
            if(p[2] * p[4] * p[6] == 0 && p[4] * p[6] * p[8] == 0) {
              to_delete[i,j] <- TRUE
              any_changes <- TRUE
            }
          }
        }
      }
    }

    pixels_removed <- pixels_removed + sum(to_delete)
    if(any(to_delete)) img[to_delete] <- 0

    # Second subiteration
    to_delete <- matrix(FALSE, nrow=nrow(img), ncol=ncol(img))

    for(i in 2:(nrow(img)-1)) {
      for(j in 2:(ncol(img)-1)) {
        if(img[i,j] == 1) {
          p <- get_neighbors(img, i, j)
          B <- sum(p[-1])
          A <- count_transitions(p)

          if(B >= 2 && B <= 6 && A == 1) {
            if(p[2] * p[4] * p[8] == 0 && p[2] * p[6] * p[8] == 0) {
              to_delete[i,j] <- TRUE
              any_changes <- TRUE
            }
          }
        }
      }
    }

    pixels_removed <- pixels_removed + sum(to_delete)
    if(any(to_delete)) img[to_delete] <- 0

    if(verbose & pixels_removed > 0) {
      cat("Iteration", iteration_count, ": Removed", pixels_removed, "pixels\n")
    }
  }

  if(verbose){
    cat("Final foreground pixels:", sum(img == 1), "\n")
    cat("Total iterations:", iteration_count, "\n")
  }


  return(img)
}

#' Thin Binary Image using Guo-Hall Algorithm (Internal)
#'
#' This internal function applies the Guo-Hall thinning algorithm to reduce binary images to their skeletons while preserving connectivity and structure.
#'
#' @param img A matrix, data frame, or `SpatRaster` object representing the binary image to be thinned.
#' @param verbose Logical. If `TRUE`, outputs diagnostic information such as image dimensions, pixel removal counts, and iteration progress. Default is `FALSE`.
#' @param layer Integer indicating the layer to use if `img` is a multi-layer `SpatRaster`. Default is 2.
#'
#' @details
#' - The input image is first processed using \code{\link{prepare_image}} to ensure it is a binary matrix.
#' - Thinning is performed in an iterative process consisting of two subiterations per cycle:
#'   1. In the first subiteration, pixels are marked for removal based on specific Guo-Hall conditions.
#'   2. In the second subiteration, a different set of conditions is applied to mark additional pixels for removal.
#' - The process continues until no pixels are removed in an iteration or the maximum number of iterations (default: 1000) is reached.
#' - The Guo-Hall algorithm ensures that the skeleton of the image is preserved.
#'
#' @return A binary matrix representing the thinned image (skeleton).
#'
#' @keywords internal
#'
#' @examples
#' \dontrun{
#' # Example usage
#' raster <- terra::rast(matrix(c(0, 1, 1, 0, 0, 1, 1, 0), nrow = 4))
#' thinned_image <- thin_image_guohall(raster, verbose = TRUE)
#' }
thin_image_guohall <- function(img, verbose = FALSE,layer = 2) {
  # flexible input
  img = prepare_image(img,layer = layer)

  if(verbose){
    cat("Image dimensions:", nrow(img), "x", ncol(img), "\n")
    cat("Initial foreground pixels:", sum(img == 1), "\n")
  }

  # Force binary values
  img[abs(img - 1) < 1e-10] <- 1
  img[abs(img) < 1e-10] <- 0

  get_neighbors <- function(img, i, j) {
    p <- rep(0, 9)
    p[1] <- img[i, j]      # P1: Center point
    if (i-1 >= 1) p[2] <- img[i-1, j]    # P2: North
    if (i-1 >= 1 && j+1 <= ncol(img)) p[3] <- img[i-1, j+1]  # P3: Northeast
    if (j+1 <= ncol(img)) p[4] <- img[i, j+1]  # P4: East
    if (i+1 <= nrow(img) && j+1 <= ncol(img)) p[5] <- img[i+1, j+1]  # P5: Southeast
    if (i+1 <= nrow(img)) p[6] <- img[i+1, j]  # P6: South
    if (i+1 <= nrow(img) && j-1 >= 1) p[7] <- img[i+1, j-1]  # P7: Southwest
    if (j-1 >= 1) p[8] <- img[i, j-1]  # P8: West
    if (i-1 >= 1 && j-1 >= 1) p[9] <- img[i-1, j-1]  # P9: Northwest
    return(p)
  }

  # Guo-Hall specific function to check if pixel can be removed
  check_pixel <- function(p) {
    # C1(p) condition
    C1 <- function(p) {
      neighbors <- c(p[2], p[3], p[4], p[5], p[6], p[7], p[8], p[9])
      min_val <- 2
      max_val <- 6
      sum(neighbors) >= min_val && sum(neighbors) <= max_val
    }

    # C2(p) condition - specific to Guo-Hall
    C2 <- function(p) {
      Xi <- 0
      if (p[2] == 0 && p[3] == 1) Xi <- Xi + 1
      if (p[3] == 0 && p[4] == 1) Xi <- Xi + 1
      if (p[4] == 0 && p[5] == 1) Xi <- Xi + 1
      if (p[5] == 0 && p[6] == 1) Xi <- Xi + 1
      if (p[6] == 0 && p[7] == 1) Xi <- Xi + 1
      if (p[7] == 0 && p[8] == 1) Xi <- Xi + 1
      if (p[8] == 0 && p[9] == 1) Xi <- Xi + 1
      if (p[9] == 0 && p[2] == 1) Xi <- Xi + 1
      Xi == 1
    }

    # C3(p) and C4(p) conditions - different for each subiteration
    C3_C4_first <- function(p) {
      ((p[2] * p[4] * p[6] == 0) || (p[4] * p[6] * p[8] == 0)) &&
        ((p[2] * p[4] * p[8] == 0) || (p[2] * p[6] * p[8] == 0))
    }

    C3_C4_second <- function(p) {
      ((p[2] * p[4] * p[6] == 0) || (p[2] * p[4] * p[8] == 0)) &&
        ((p[2] * p[6] * p[8] == 0) || (p[4] * p[6] * p[8] == 0))
    }

    list(
      first = C1(p) && C2(p) && C3_C4_first(p),
      second = C1(p) && C2(p) && C3_C4_second(p)
    )
  }

  max_iterations <- 1000
  iteration_count <- 0
  any_changes <- TRUE

  while(any_changes && iteration_count < max_iterations) {
    iteration_count <- iteration_count + 1
    any_changes <- FALSE
    pixels_removed <- 0

    # First subiteration
    to_delete <- matrix(FALSE, nrow=nrow(img), ncol=ncol(img))

    for(i in 2:(nrow(img)-1)) {
      for(j in 2:(ncol(img)-1)) {
        if(img[i,j] == 1) {
          p <- get_neighbors(img, i, j)
          if(check_pixel(p)$first) {
            to_delete[i,j] <- TRUE
            any_changes <- TRUE
          }
        }
      }
    }

    pixels_removed <- pixels_removed + sum(to_delete)
    if(any(to_delete)) img[to_delete] <- 0

    # Second subiteration
    to_delete <- matrix(FALSE, nrow=nrow(img), ncol=ncol(img))

    for(i in 2:(nrow(img)-1)) {
      for(j in 2:(ncol(img)-1)) {
        if(img[i,j] == 1) {
          p <- get_neighbors(img, i, j)
          if(check_pixel(p)$second) {
            to_delete[i,j] <- TRUE
            any_changes <- TRUE
          }
        }
      }
    }

    pixels_removed <- pixels_removed + sum(to_delete)
    if(any(to_delete)) img[to_delete] <- 0

    if(verbose && pixels_removed > 0) {
      cat("Iteration", iteration_count, ": Removed", pixels_removed, "pixels\n")
    }
  }

  if(verbose){
    cat("Final foreground pixels:", sum(img == 1), "\n")
    cat("Total iterations:", iteration_count, "\n")
  }


  return(img)
}

#' Medial Axis Transform (Internal)
#'
#' This internal function computes the medial axis transform of a binary image, identifying the set of skeleton points equidistant to the object's boundaries.
#'
#' @param img A matrix, data frame, or `SpatRaster` object representing the binary image for transformation.
#' @param verbose Logical. If `TRUE`, outputs diagnostic information such as image dimensions, progress of computation, and final skeleton size. Default is `FALSE`.
#' @param layer Integer indicating the layer to use if `img` is a multi-layer `SpatRaster`. Default is 2.
#'
#' @details
#' - The input image is first processed using \code{\link{prepare_image}} to ensure it is a binary matrix.
#' - The algorithm proceeds through the following steps:
#'   1. **Distance Transform**: Computes the distance of each foreground pixel to the nearest background pixel using a two-pass algorithm.
#'   2. **Local Maxima Detection**: Identifies local maxima in the distance transform to mark potential skeleton points.
#'   3. **Skeleton Refinement**: Ensures connectivity by connecting skeleton points within an 8-neighborhood.
#' - The result is a binary image representing the medial axis of the input object.
#'
#' @return A binary matrix where `1` represents skeleton pixels and `0` represents the background.
#'
#' @keywords internal
#'
#' @examples
#' \dontrun{
#' # Example usage
#' raster <- terra::rast(matrix(c(0, 1, 1, 0, 0, 1, 1, 0), nrow = 4))
#' skeleton <- medial_axis_transform(raster, verbose = TRUE)
#' }
medial_axis_transform <- function(img, verbose = FALSE,layer = 2) {

  # flexible input
  img = prepare_image(img,layer = layer)


  # Force binary values
  img[abs(img - 1) < 1e-10] <- 1
  img[abs(img) < 1e-10] <- 0

  if(verbose) {
    cat("Image dimensions:", nrow(img), "x", ncol(img), "\n")
    cat("Initial foreground pixels:", sum(img == 1), "\n")
  }

  # Helper function to compute distance transform
  compute_distance_transform <- function(binary_img) {
    dist_map <- matrix(0, nrow = nrow(binary_img), ncol = ncol(binary_img))

    # First pass - forward scan
    for(i in 2:nrow(binary_img)) {
      for(j in 2:ncol(binary_img)) {
        if(binary_img[i,j] == 1) {
          dist_map[i,j] <- min(
            dist_map[i-1,j],
            dist_map[i,j-1],
            dist_map[i-1,j-1]
          ) + 1
        }
      }
    }

    # Second pass - backward scan
    for(i in (nrow(binary_img)-1):1) {
      for(j in (ncol(binary_img)-1):1) {
        if(binary_img[i,j] == 1) {
          dist_map[i,j] <- min(
            dist_map[i,j],
            dist_map[i+1,j] + 1,
            dist_map[i,j+1] + 1,
            dist_map[i+1,j+1] + 1
          )
        }
      }
    }

    return(dist_map)
  }

  # Helper function to find local maxima
  find_local_maxima <- function(dist_map) {
    maxima <- matrix(FALSE, nrow = nrow(dist_map), ncol = ncol(dist_map))

    for(i in 2:(nrow(dist_map)-1)) {
      for(j in 2:(ncol(dist_map)-1)) {
        if(dist_map[i,j] > 0) {
          # Check 8-neighborhood
          neighborhood <- dist_map[max(1,i-1):min(nrow(dist_map),i+1),
                                   max(1,j-1):min(ncol(dist_map),j+1)]
          if(dist_map[i,j] >= max(neighborhood)) {
            maxima[i,j] <- TRUE
          }
        }
      }
    }

    return(maxima)
  }

  # Compute distance transform
  dist_transform <- compute_distance_transform(img)

  if(verbose) {
    cat("\nDistance transform computed\n")
  }

  # Find local maxima in distance transform
  skeleton <- find_local_maxima(dist_transform)

  # Connect local maxima using path following
  result <- matrix(0, nrow = nrow(img), ncol = ncol(img))
  result[skeleton] <- 1

  # Post-process to ensure connectivity
  for(i in 2:(nrow(img)-1)) {
    for(j in 2:(ncol(img)-1)) {
      if(skeleton[i,j]) {
        # Check 8-neighborhood for other skeleton points
        neighborhood <- skeleton[max(1,i-1):min(nrow(img),i+1),
                                 max(1,j-1):min(ncol(img),j+1)]
        if(sum(neighborhood) > 1) {
          # Connect to nearest skeleton point
          max_dist <- dist_transform[i,j]
          for(di in -1:1) {
            for(dj in -1:1) {
              if(di != 0 || dj != 0) {
                ni <- i + di
                nj <- j + dj
                if(ni >= 1 && ni <= nrow(img) &&
                   nj >= 1 && nj <= ncol(img) &&
                   skeleton[ni,nj]) {
                  result[i,j] <- 1
                  break
                }
              }
            }
          }
        }
      }
    }
  }

  if(verbose) {
    cat("\nFinal skeleton pixels:", sum(result), "\n")
  }

  return(result)
}

#' Skeletonization Wrapper Function
#'
#' This function serves as a wrapper for applying different skeletonization methods to a binary image, including the Zhang-Suen, Guo-Hall, and Medial Axis Transform (MAT) algorithms.
#'
#' @param img A matrix, data frame, or `SpatRaster` object representing the binary image to be skeletonized.
#' @param methods A character vector specifying the skeletonization methods to apply. Valid options are \code{"ZhangSuen"}, \code{"GuoHall"}, and \code{"MAT"}. Defaults to all three methods.
#' @param verbose Logical. If \code{TRUE}, displays progress and diagnostic messages during processing. Defaults to \code{TRUE}.
#' @param layer Integer specifying the layer to use if \code{img} is a multi-layer `SpatRaster`. Defaults to 2.
#'
#' @details
#' This function allows for flexible and streamlined skeletonization of binary images using one or more supported algorithms:
#' \itemize{
#'   \item \code{"ZhangSuen"}: Implements the Zhang-Suen thinning algorithm.
#'   \item \code{"GuoHall"}: Implements the Guo-Hall thinning algorithm.
#'   \item \code{"MAT"}: Computes the Medial Axis Transform to extract the skeleton.
#' }
#'
#' The function processes the input image with the specified methods and returns the results. If multiple methods are chosen, the results are returned as a named list, with each element corresponding to a method.
#'
#' @return If a single method is selected, the function returns a `SpatRaster` object representing the skeletonized image. If multiple methods are selected, a named list of `SpatRaster` objects is returned.
#'
#' @examples
#' \dontrun{
#' # Load a binary image as a SpatRaster
#' binary_image <- terra::rast(matrix(c(0, 1, 1, 0, 0, 1, 1, 0), nrow = 4))
#'
#' # Apply all skeletonization methods
#' skeletons <- skeletonize_image(binary_image, verbose = TRUE)
#'
#' # Apply a single method
#' zhang_skeleton <- skeletonize_image(binary_image, methods = "ZhangSuen")
#'
#' # Apply multiple methods
#' selected_skeletons <- skeletonize_image(binary_image, methods = c("GuoHall", "MAT"))
#' }
#'
#' @seealso \code{\link{thin_image_zhangsuen}}, \code{\link{thin_image_guohall}}, \code{\link{medial_axis_transform}}
#' @export
skeletonize_image <- function(img, methods = c("ZhangSuen", "GuoHall", "MAT"), verbose = TRUE, layer = 2) {
  # Ensure methods are valid
  valid_methods <- c("ZhangSuen", "GuoHall", "MAT")
  methods <- intersect(methods, valid_methods)
  if (length(methods) == 0) {
    stop("No valid methods specified. Choose from: 'ZhangSuen', 'GuoHall', 'MAT'.")
  }

  # # input flexibility
  # img = prepare_image(img = img)

  # Process each method
  results <- list()
  for (method in methods) {
    if (verbose) cat("\nApplying method:", method, "\n")
    result <- switch(
      method,
      "ZhangSuen" =
        terra::rast(thin_image_zhangsuen(img, verbose = verbose)),
      "GuoHall" = terra::rast(thin_image_guohall(img, verbose = verbose)),
      "MAT" = terra::rast(medial_axis_transform(img, verbose = verbose)),
      stop(paste("Unsupported method:", method))
    )
    results[[method]] <- result
  }

  # Return results as a list
  if (length(results) == 1) {
    return(results[[1]])  # Single method result
  } else {
    return(results)  # Multiple method results
  }
}


#' Detect Skeleton Points: Branching Points and Endpoints
#'
#' Identifies the branching points and endpoints of a skeletonized binary image.
#'
#' @param img A matrix, data frame, or `SpatRaster` object representing the skeletonized binary image.
#'
#' @details
#' This function detects key points in a skeletonized binary image:
#' \itemize{
#'   \item \strong{Endpoints}: Pixels with exactly one neighbor in the skeleton.
#'   \item \strong{Branching Points}: Pixels with more than two neighbors in the skeleton.
#' }
#'
#' The function uses a 3x3 neighborhood kernel to count the number of neighbors for each foreground pixel (\code{1}) in the image. Based on the neighbor count, points are classified as endpoints or branching points.
#'
#' The input image should be skeletonized (thin and connected) before using this function. If not already binary, the input image will be binarized internally.
#'
#' @return A named list containing two `SpatRaster` objects:
#' \itemize{
#'   \item \code{endpoints}: A binary raster where endpoints are marked as \code{1}.
#'   \item \code{branching_points}: A binary raster where branching points are marked as \code{1}.
#' }
#'
#' @examples
#' \dontrun{
#' library(terra)
#'
#' # Example skeletonized image
#' skeleton <- rast(matrix(c(0, 1, 1, 0, 0, 1, 1, 0), nrow = 4))
#'
#' # Detect endpoints and branching points
#' points <- detect_skeleton_points(skeleton)
#'
#' # Access results
#' endpoints <- points$endpoints
#' branching_points <- points$branching_points
#' }
#'
#' @seealso \code{\link{skeletonize_image}}, \code{\link{thin_image_zhangsuen}}, \code{\link{thin_image_guohall}}
#' @export
detect_skeleton_points <- function(img) {
  img = prepare_image(img)
  # Count neighbors using a kernel method
  count_neighbors <- function(img) {
    # Ensure the input is binary
    img[img > 0] <- 1

    # Define the kernel (3x3 neighborhood)
    kernel <- matrix(1, nrow = 3, ncol = 3)
    kernel[2, 2] <- 0  # Exclude the center pixel

    # Pad the image with zeros to handle edges
    padded_img <- matrix(0, nrow = nrow(img) + 2, ncol = ncol(img) + 2)
    padded_img[2:(nrow(img) + 1), 2:(ncol(img) + 1)] <- img

    # Initialize a result matrix
    neighbor_count <- matrix(0, nrow = nrow(img), ncol = ncol(img))

    # Perform the convolution manually
    for (i in 2:(nrow(padded_img) - 1)) {
      for (j in 2:(ncol(padded_img) - 1)) {
        region <- padded_img[(i - 1):(i + 1), (j - 1):(j + 1)]
        neighbor_count[i - 1, j - 1] <- sum(region * kernel)
      }
    }

    return(neighbor_count)
  }
  # Count neighbors
  neighbor_count <- count_neighbors(img)

  # Classify points
  endpoints <- (img == 1) & (neighbor_count == 1)
  branching_points <- (img == 1) & (neighbor_count > 2)

  # SpatRast Output
  endpoints <- terra::rast(endpoints)
  branching_points <- terra::rast(branching_points)

  return(list(
    endpoints = endpoints,
    branching_points = branching_points
  ))
}



