#' Stitch Sequential Images
#'
#' This function stitches a sequence of images together by aligning overlapping regions based on cross-correlation or phase-correlation.
#'
#' @param imgs A character vector of file paths to the images to be stitched.
#' @param overlap_px An integer specifying the number of pixels to use for overlap between images. Default is 200.
#' @param side1 A string specifying the side of the previous image used for overlap (e.g., "bottom", "top", "left", "right"). Default is "bottom".
#' @param side2 A string specifying the side of the current image used for overlap. Default is "top".
#' @param method A string specifying the method for alignment, either "crosscorr" or "phasecorr". Default is "crosscorr".
#'
#' @return A stitched image as a cimg object.
#'
#' @examples
#' \dontrun{
#' imgs <- c("image1.jpg", "image2.jpg", "image3.jpg")
#' stitched <- stitch_sequential_images(imgs = imgs,
#'                     overlap_px = 250, side1 = "bottom", side2 = "top")
#'}
#' @export
stitch_sequential_images <- function(imgs, overlap_px = 200, side1="bottom", side2="top", method="crosscorr") {
  images <- lapply(imgs,
                   load_flexible_image( select.layer = NULL, output_format = "cimg", normalize = FALSE  )
  )
  result <- images[[1]]
  total_width <- dim(result)[2]  # Track total canvas width

  for (i in 2:length(images)) {
    # define stitching image pairs and carry-over the already stichted ones
    prev_img <- result
    curr_img <- images[[i]]

    # Find alignment using overlap regions
    prev_overlap <- extract_overlap_region(prev_img, overlap_px, side=side1)
    curr_overlap <- extract_overlap_region(curr_img, overlap_px, side=side2)

    translation <- estimate_translation(prev_overlap, curr_overlap, method=method)


    # Calculate new dimensions
    new_width <- total_width + dim(curr_img)[2] - overlap_px
    R <- new_width
    C <- max(dim(curr_img)[1] + translation[2], dim(result)[1])

    # Position translation matrices
    M_curr <- matrix(c(
      1, 0, total_width - overlap_px,
      0, 1, translation[2]
    ), nrow=2, byrow=TRUE)

    M_prev <- matrix(c(
      1, 0, 0,
      0, 1, 0
    ), nrow=2, byrow=TRUE)

    # convert to OpenImageR:: format (array)
    curr_warp_img <- array(curr_img, dim(curr_img))[,,,1:3]
    prev_warp_img <- array(prev_img, dim(prev_img))[,,,1:3]
    # affine translation warp
    warped_curr <- OpenImageR::warpAffine(curr_warp_img, M=M_curr, R=R, C=C)
    warped_prev <- OpenImageR::warpAffine(prev_warp_img, M=M_prev, R=R, C=C)
    # convert to imager:: format (cimg)
    warped_curr <- imager::as.cimg(warped_curr)
    warped_prev <- imager::as.cimg(warped_prev)

    # blend images and update
    result <- combine_images(warped_prev, warped_curr,
                             overlap_start=total_width - overlap_px,
                             overlap_end=total_width)

    total_width <- new_width
  }

  return(result)
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
  w <- imager::width(img)
  h <- imager::height(img)
  overlap_width <- overlap_px
  overlap_height <- overlap_px

  if (side == "right") {
    region <- img[ (w - overlap_width+1):w, , 1, ]
  }
  if (side == "left") {
    region <- img[ 1:overlap_width, , 1, ]
  }
  if (side == "bottom") {
    region <- img[,(h - overlap_height+1):h  , 1, ]
  }
  if (side == "top") {
    region <- img[,1:overlap_height  , 1, ]
  }
  return(region)
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
  # Grayscale conversion
  prev_region_gray <- prev_region[,,1]*0.21 + prev_region[,,2]*0.72 + prev_region[,,3]*0.07
  curr_region_gray <- curr_region[,,1]*0.21 + curr_region[,,2]*0.72 + curr_region[,,3]*0.07


  if (method == "crosscorr" ) {
    # Compute normalized cross-correlation
    cc <- imagefx::xcorr3d(prev_region_gray, curr_region_gray)
    xy <- as.vector(cc$max.shifts)
    cat("max. cross-correlation:", round(cc$max.corr, 3), "\n")
  } else if (method == "phasecorr") {
    # Compute phase-correlation
    pp <- imagefx::pcorr3d(prev_region_gray, curr_region_gray)
    xy <- as.vector(pp$max.shifts)
    cat("max. phase-correlation:", round(pp$max.corr, 3), "\n")
  } else {
    cat("No appropriate method or overlap detected. Ensure input regions are valid.\n")
  }

  # Estimate translation
  dx <- xy[2]
  dy <- xy[1]


  return(c(dx, dy))
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
  # Initialize result as warped_prev
  blended <- warped_prev

  # Add non-overlapping part of warped_curr
  blended[, overlap_end:dim(warped_curr)[2], 1, ] <- warped_curr[, overlap_end:dim(warped_curr)[2], 1, ]

  # Blend the overlap zone with 0.5 weights
  for(c in 1:3) {
    blended[, overlap_start:overlap_end, 1, c] <-
      0.5 * warped_prev[, overlap_start:overlap_end, 1, c] +
      0.5 * warped_curr[, overlap_start:overlap_end, 1, c]
  }

  return(blended)
}

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

