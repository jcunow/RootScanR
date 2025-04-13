
#######################

# Example usage
# t0 = Sys.time()
# path = "C:/Users/jocu0013/Desktop/Oulanka/Scan_Raw/Oulanka2023_June/"
# files = list.files(path)
# image_paths <- c(paste0(path,files[269]), paste0(path,files[270]),paste0(path,files[271]),paste0(path,files[272]),paste0(path,files[273]),paste0(path,files[274]))
# stitched <- stitch_sequential_images(image_paths, overlap_px = 200, side1 = "bottom",side2 = "top")
# difftime(Sys.time(),t0)
# dev.off()
# #save.image(stitched, "stitched_output.jpg")
# OpenImageR::imageShow(array(stitched,dim(stitched)[c(1:2,4)]))




##############################################
#' Stitch Sequential Images
#'
#' This function stitches a sequence of images together by aligning overlapping regions based on cross-correlation or phase-correlation.
#'
#' @param imgs A character vector of file paths to the images to be stitched.
#' @param overlap_px An integer specifying the number of pixels to use to search for overlap between images. Default is 200. The CI-600 scanner produces 100-200 px overlap at 300dpi.
#' @param side1 A string specifying the side of the previous image used for overlap (e.g., "bottom", "top", "left", "right"). Default is "bottom".
#' @param method A string specifying the method for alignment, either "crosscorr" or "phasecorr". Default is "crosscorr".
#'
#' @return A stitched image as a cimg object.
#'
#' @examples
#' \dontrun{
#' imgs <- c("image1.jpg", "image2.jpg", "image3.jpg")
#' stitched <- stitch_sequential_images(imgs = imgs,
#'                     overlap_px = 250, side1 = "bottom")
#'}
#' @export
stitch_sequential_images <- function(imgs, overlap_px = 200, side1="bottom", method="crosscorr") {
  # Input validation module
  tryCatch({
    # Validate input parameters
    if (!is.character(imgs) || length(imgs) < 2)
      stop("imgs must be a character vector with at least 2 image paths")
    if (!all(file.exists(imgs)))
      stop("One or more image files do not exist")
    if (!is.numeric(overlap_px) || overlap_px <= 0)
      stop("overlap_px must be a positive number")
    if (!side1 %in% c("bottom", "top", "left", "right"))
      stop("side1 must be one of: 'bottom', 'top', 'left', 'right'")
    if (!method %in% c("crosscorr", "phasecorr"))
      stop("method must be either 'crosscorr' or 'phasecorr'")

    # Automatically determine the opposite side
    opposite_sides <- list("bottom" = "top", "top" = "bottom", "left" = "right", "right" = "left")
    side2 <- opposite_sides[[side1]]

    # Load images with error handling
    images <- lapply(imgs, function(img_path) {
      img <- tryCatch({
        load_flexible_image(img_path, select.layer = NULL, output_format = "cimg", normalize = FALSE)
      }, error = function(e) stop(sprintf("Failed to load image: %s", img_path)))

      if (is.null(img) || any(dim(img) <= 0)) stop("Invalid image dimensions")
      return(img)
    })

    # Validate image dimensions consistency
    dims <- lapply(images, dim)
    if (!all(sapply(dims[-1], function(x) all(x[3:4] == dims[[1]][3:4]))))
      stop("All images must have the same number of channels and frames")

    # Initialize result
    result <- images[[1]]
    total_width <- dim(result)[2]

    # Process images sequentially with progress tracking
    for (i in 2:length(images)) {
      message(sprintf("Processing image pair %d/%d", i-1, length(images)-1))

      prev_img <- result
      curr_img <- images[[i]]

      # Validate overlap size
      if (overlap_px >= dim(prev_img)[2] || overlap_px >= dim(curr_img)[2])
        stop(sprintf("Overlap size (%d) is larger than image width", overlap_px))

      # Extract and validate overlap regions
      prev_overlap <- extract_overlap_region(prev_img, overlap_px, side=side1)
      curr_overlap <- extract_overlap_region(curr_img, overlap_px, side=side2)

      if (is.null(prev_overlap) || is.null(curr_overlap))
        stop("Failed to extract overlap regions")

      # Estimate translation with validation
      translation <- estimate_translation(prev_overlap, curr_overlap, method=method)
      if (any(is.na(translation)) || length(translation) != 2)
        stop("Invalid translation estimation")

      # Validate translation magnitude
      max_allowed_translation <- dim(curr_img)[1] * 0.5
      if (any(abs(translation) > max_allowed_translation))
        warning(sprintf("Large translation detected (%d, %d). Results may be unreliable.", translation[1], translation[2]))

      # Calculate new dimensions
      new_width <- total_width + dim(curr_img)[2] - overlap_px
      R <- new_width
      C <- max(dim(curr_img)[1] + translation[2], dim(result)[1])

      if (new_width <= 0 || C <= 0)
        stop("Invalid resulting image dimensions")

      # Safe matrix operations
      M_curr <- matrix(c(1, 0, total_width - overlap_px, 0, 1, translation[2]), nrow=2, byrow=TRUE)
      M_prev <- diag(3)[1:2,]

      # Safe warping operations
      warped_curr <- tryCatch({
        OpenImageR::warpAffine(as.array(curr_img)[,,,1:3], M=M_curr, R=R, C=C)
      }, error = function(e) stop("Warping failed for current image"))

      warped_prev <- tryCatch({
        OpenImageR::warpAffine(as.array(prev_img)[,,,1:3], M=M_prev, R=R, C=C)
      }, error = function(e) stop("Warping failed for previous image"))

      # Convert warped images back to cimg format
      warped_curr <- tryCatch(imager::as.cimg(warped_curr), error = function(e) stop("Failed to convert warped image"))
      warped_prev <- tryCatch(imager::as.cimg(warped_prev), error = function(e) stop("Failed to convert warped image"))

      # Combine images
      result <- combine_images(warped_prev, warped_curr, overlap_start=total_width - overlap_px, overlap_end=total_width)
      if (is.null(result))
        stop("Failed to combine images")

      total_width <- new_width
    }

    return(result)

  }, error = function(e) {
    stop(paste("Error in stitch_sequential_images:", e$message))
  })
}




#' Extract Overlap Region from an Image
#'
#' Extracts a region of an image corresponding to a specified overlap along a given side.
#'
#' @param img A cimg object representing the image.
#' @param overlap_px An integer specifying the size of the overlap region in pixels. Default is 200.
#' @param side A string specifying the side of the image to extract the overlap from. One of "top", "bottom", "left", or "right".
#'
#' @return A cimg object representing the extracted overlap region.
#'
#' @keywords internal
extract_overlap_region <- function(img, overlap_px = 250, side) {
  tryCatch({
    if (!inherits(img, "cimg"))
      stop("Input must be a cimg object")
    if (!is.numeric(overlap_px) || overlap_px <= 0)
      stop("overlap_px must be a positive number")
    if (!side %in% c("top", "bottom", "left", "right"))
      stop("Invalid side specified")

    w <- imager::width(img)
    h <- imager::height(img)

    if (overlap_px >= w || overlap_px >= h)
      stop("Overlap size exceeds image dimensions")

    # Extract region based on side
    region <- switch(side,
                     "right" = if (w > overlap_px) img[(w - overlap_px+1):w, , 1, ] else NULL,
                     "left" = if (w > overlap_px) img[1:overlap_px, , 1, ] else NULL,
                     "bottom" = if (h > overlap_px) img[,(h - overlap_px+1):h, 1, ] else NULL,
                     "top" = if (h > overlap_px) img[,1:overlap_px, 1, ] else NULL,
                     NULL)

    if (is.null(region))
      stop("Failed to extract region")

    return(region)
  }, error = function(e) {
    stop(paste("Error in extract_overlap_region:", e$message))
  })
}



#' Estimate Translation Between Overlap Regions
#'
#' Computes the translation required to align two overlapping regions using either cross-correlation or phase-correlation.
#'
#' @param prev_region A cimg object representing the overlap region from the previous image.
#' @param curr_region A cimg object representing the overlap region from the current image.
#' @param method A string specifying the method for estimating translation. One of "crosscorr" or "phasecorr". Default is "crosscorr".
#'
#' @return A numeric vector of length two containing the translation offsets (dx, dy).
#'
#' @keywords internal
estimate_translation <- function(prev_region, curr_region, method = "crosscorr") {
  tryCatch({
    if (!all(dim(prev_region) == dim(curr_region)))
      stop("Overlap regions must have the same dimensions")

    # Validate input types
    if (!inherits(prev_region, "cimg") || !inherits(curr_region, "cimg"))
      stop("Input regions must be cimg objects")

    # Safe grayscale conversion
    prev_region_gray <- prev_region[,,1]*0.21 + prev_region[,,2]*0.72 + prev_region[,,3]*0.07
    curr_region_gray <- curr_region[,,1]*0.21 + curr_region[,,2]*0.72 + curr_region[,,3]*0.07

    # Compute correlation based on method
    result <- switch(method,
                     "crosscorr" = {
                       cc <- imagefx::xcorr3d(prev_region_gray, curr_region_gray)
                       if (cc$max.corr < 0.1)
                         warning("Very low correlation detected, results may be unreliable")
                       cat("max. cross-correlation:", round(cc$max.corr, 3), "\n")
                       cc$max.shifts
                     },
                     "phasecorr" = {
                       pp <- imagefx::pcorr3d(prev_region_gray, curr_region_gray)
                       if (pp$max.corr < 0.1)
                         warning("Very low correlation detected, results may be unreliable")
                       cat("max. phase-correlation:", round(pp$max.corr, 3), "\n")
                       pp$max.shifts
                     },
                     stop("Invalid method specified"))

    xy <- as.vector(result)
    if (length(xy) != 2)
      stop("Invalid correlation result")

    return(c(xy[2], xy[1]))

  }, error = function(e) {
    stop(paste("Error in estimate_translation:", e$message))
  })
}


#' Combine Two Warped Images
#'
#' Combines two warped images by blending their overlapping regions and concatenating the non-overlapping parts.
#'
#' @param warped_prev A cimg object representing the first warped image.
#' @param warped_curr A cimg object representing the second warped image.
#' @param overlap_start An integer specifying the starting column index of the overlap region.
#' @param overlap_end An integer specifying the ending column index of the overlap region.
#'
#' @return A cimg object representing the blended and combined image.
#'
#' @keywords internal
combine_images <- function(warped_prev, warped_curr, overlap_start, overlap_end) {
  tryCatch({
    # Validate inputs
    if (!inherits(warped_prev, "cimg") || !inherits(warped_curr, "cimg"))
      stop("Input images must be cimg objects")
    if (!is.numeric(overlap_start) || !is.numeric(overlap_end))
      stop("overlap_start and overlap_end must be numeric")
    if (overlap_start >= overlap_end)
      stop("overlap_start must be less than overlap_end")
    if (overlap_end > dim(warped_curr)[2])
      stop("overlap_end exceeds image dimensions")

    # Initialize result
    blended <- warped_prev

    # Safe image combination
    tryCatch({
      # Add non-overlapping part
      blended[, overlap_end:dim(warped_curr)[2], 1, ] <-
        warped_curr[, overlap_end:dim(warped_curr)[2], 1, ]

      # Blend overlap zone
      for(c in 1:3) {
        blended[, overlap_start:overlap_end, 1, c] <-
          0.5 * warped_prev[, overlap_start:overlap_end, 1, c] +
          0.5 * warped_curr[, overlap_start:overlap_end, 1, c]
      }
    }, error = function(e) {
      stop("Failed to blend images: ", e$message)
    })

    return(blended)

  }, error = function(e) {
    stop(paste("Error in combine_images:", e$message))
  })
}
