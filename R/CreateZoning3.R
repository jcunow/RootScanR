
#' Create a buffer halo) around non-zero pixels
#'
#' @param img SpatRaster/matrix/array - segmented image
#' @param width numeric - buffer width in pixels (default: 2)
#' @param halo.only logical - if TRUE, returns only the buffer zone (default: TRUE)
#' @param kernel character - shape of the thickening kernel: "circle" or "diamond"
#'
#' @return SpatRast - buffer zone around non-zero pixels
#' @export
#'
#'
#' @examples
#' data(seg_Oulanka2023_Session03_T067)
#' img <- terra::rast(seg_Oulanka2023_Session03_T067)
#' Halo(img, width = 2)
Halo = function(img, width=2, halo.only=TRUE, kernel="circle") {

  # Validation module
  tryCatch({
    # Check required input
    if (is.null(img)) {
      stop("Input image is required")
    }

    # Parameter validation
    if (!is.numeric(width) || width < 1) {
      stop("width must be a positive integer")
    }
    if (!is.logical(halo.only)) {
      stop("halo.only must be logical")
    }
    if (!kernel %in% c("circle", "diamond")) {
      stop("kernel must be either 'circle' or 'diamond'")
    }

    # Load and validate image
    im <- tryCatch({
      load_flexible_image(img, output_format="SpatRaster", normalize=TRUE, binarize = TRUE)
    }, error = function(e) {
      stop("Failed to load image: ", e$message)
    })

    if (terra::ncell(im) == 0) {
      stop("Input image has no valid cells")
    }

    # Create kernel based on type
    k0 = if (kernel == "circle") {
      matrix(c(1,1,1,1,0,1,1,1,1), nrow=3, ncol=3)
    } else {
      matrix(c(0,1,0,1,0,1,0,1,0), nrow=3, ncol=3)
    }

    # Apply focal operation with error handling
    im2 = im
    for (itr in 1:width) {
      im2 <- tryCatch({
        terra::focal(im2, w=k0, fun="sum", na.policy="omit")
      }, error = function(e) {
        stop(sprintf("Focal operation failed at iteration %d: %s", itr, e$message))
      })
    }

    # Create output image
    out.im = (im2 >= 1) * 1

    if (halo.only) {
      out.im = out.im - im
      out.im = terra::subst(out.im, from=-1, to=0)
    }

    # Validate output
    if (terra::ncell(out.im) == 0) {
      warning("Output image has no valid cells")
    }

    return(out.im)

  }, error = function(e) {
    stop(paste("Error in Halo:", e$message))
  }, warning = function(w) {
    warning(paste("Warning in Halo:", w$message))
  })
}

#' Bin continuous depth values into discrete intervals
#'
#' @param depthmap SpatRaster/matrix/array - continuous depth values
#' @param nn numeric - bin width
#' @param round.option character - binning method: "rounding", "ceiling", or "floor"
#' @return SpatRaster - binned depth values
#' @export
#'
#' @examples
#' data(seg_Oulanka2023_Session01_T067)
#' img = terra::rast(seg_Oulanka2023_Session01_T067)
#' mask = img[[1]] - img[[2]]
#' mask[mask == 255] <- NA
#' img = img
#' depthmap = create.depthmap(img,mask,start.soil = 2.9, select.layer = 2 )
#' binned.map = binning(depthmap,nn = 5)
binning = function(depthmap, nn, round.option="rounding") {

  # Validation module
  tryCatch({
    # Check required inputs
    if (is.null(depthmap)) {
      stop("depthmap is required")
    }
    if (missing(nn)) {
      stop("bin width (nn) is required")
    }

    # Parameter validation
    if (!is.numeric(nn) || nn <= 0) {
      stop("bin width (nn) must be positive numeric")
    }
    if (!round.option %in% c("rounding", "ceiling", "floor")) {
      stop("round.option must be one of: 'rounding', 'ceiling', 'floor'")
    }

    # Load and validate image
    img <- tryCatch({
      load_flexible_image(depthmap, output_format="spatrast", normalize=FALSE)
    }, error = function(e) {
      stop("Failed to load depthmap: ", e$message)
    })

    if (terra::ncell(img) == 0) {
      stop("Depthmap has no valid cells")
    }

    # Check for infinite or NA values
    if (any(is.infinite(terra::values(img)), na.rm=TRUE)) {
      warning("Infinite values detected in depthmap")
    }

    # Perform binning based on selected method
    im = tryCatch({
      if (round.option == "rounding") {
        nn * round(depthmap / nn, 0)
      } else if (round.option == "ceiling") {
        nn * ceiling(depthmap / nn)
      } else {
        nn * floor(depthmap / nn)
      }
    }, error = function(e) {
      stop("Binning operation failed: ", e$message)
    })

    # Validate output
    if (terra::ncell(im) == 0) {
      warning("Output has no valid cells")
    }

    return(im)

  }, error = function(e) {
    stop(paste("Error in binning:", e$message))
  }, warning = function(w) {
    warning(paste("Warning in binning:", w$message))
  })
}


#' Extract a zone of interest based on depth values
#'
#' @param rootpic SpatRaster/matrix/array - source image
#' @param binned.map SpatRaster/matrix/array - binned depth map
#' @param indexD numeric - depth index to extract. Usually used to index in a loop to iterate over depths.
#' @param nn numeric - bin width that has been used in binning
#' @param select.layer Integer. Specifies which layer to use if the input is a multi-band image. Default is `2`.
#' @param silent logical - suppress messages
#' @return SpatRaster - extracted zone
#' @export
#'
#' @examples
#' data(seg_Oulanka2023_Session01_T067)
#' img = terra::rast(seg_Oulanka2023_Session01_T067)
#' mask = img[[1]] - img[[2]]
#' mask[mask == 255] <- NA
#' img = img
#' depthmap = create.depthmap(img,mask,start.soil = 2.9, select.layer = 2 )
#' binned.map = binning(depthmap,nn = 5)
#'
#' data(seg_Oulanka2023_Session01_T067)
#' img = terra::rast(seg_Oulanka2023_Session01_T067)
#' mask = img[[1]] - img[[2]]
#' mask[mask == 255] <- NA
#' depthmap = create.depthmap(img,mask,start.soil = 0, tube.thicc = 7,dpi = 150, select.layer = 2 )
#' binned.map = binning(depthmap,nn = 5)
#' image.zone = zone.fun(img, binned.map, indexD = 10, nn=5, silent = FALSE,select.layer = NULL)
zone.fun = function(rootpic, binned.map, indexD=0, nn=5, silent=FALSE, select.layer=NULL) {

  # Validation module
  tryCatch({
    # Check required inputs
    if (is.null(rootpic) || is.null(binned.map)) {
      stop("Both rootpic and binned.map are required")
    }

    # Parameter validation
    if (!is.numeric(indexD)) {
      stop("indexD must be numeric")
    }
    if (!is.numeric(nn) || nn <= 0) {
      stop("nn must be positive numeric")
    }
    if (!is.logical(silent)) {
      stop("silent must be logical")
    }
    if (!is.null(select.layer) && (!is.numeric(select.layer) || select.layer < 1)) {
      stop("select.layer must be NULL or positive integer")
    }

    # Load and validate images
    rootpic <- tryCatch({
      load_flexible_image(rootpic, select.layer=select.layer,
                          output_format="spatrast", normalize=FALSE)
    }, error = function(e) {
      stop("Failed to load rootpic: ", e$message)
    })

    binned.map <- tryCatch({
      load_flexible_image(binned.map, select.layer=select.layer,
                          output_format="spatrast", normalize=FALSE)
    }, error = function(e) {
      stop("Failed to load binned.map: ", e$message)
    })

    # Match dimensions with validation
    ex.bm = c(0, dim(binned.map)[1], 0, dim(binned.map)[2])
    ex.rp = c(0, dim(rootpic)[1], 0, dim(rootpic)[2])

    if (any(ex.bm != ex.rp)) {
      rootpic = tryCatch({
        terra::t(rootpic)
      }, error = function(e) {
        stop("Failed to transpose rootpic: ", e$message)
      })
    }

    # Update extents after potential transformation
    ex.bm = c(0, dim(binned.map)[1], 0, dim(binned.map)[2])
    ex.rp = c(0, dim(rootpic)[1], 0, dim(rootpic)[2])

    tryCatch({
      terra::ext(binned.map) = ex.bm
      terra::ext(rootpic) = ex.rp
    }, error = function(e) {
      stop("Failed to set extents: ", e$message)
    })

    # Process image
    r = rootpic
    if (terra::ext(rootpic)[2] != terra::ext(binned.map)[2] ||
        terra::ext(rootpic)[4] != terra::ext(binned.map)[4]) {
      r = tryCatch({
        terra::crop(r, ex.bm)
      }, error = function(e) {
        stop("Failed to crop image: ", e$message)
      })
    }

    # Set values with validation
    tryCatch({
      terra::values(r)[terra::values(binned.map != indexD)] <- NA
      terra::values(r)[is.na(terra::values(binned.map))] <- NA
    }, error = function(e) {
      stop("Failed to set NA values: ", e$message)
    })

    # Calculate coverage
    coverNA = round(sum(is.na(terra::values(r[[1]]))) / terra::ncell(r), 2)

    if (coverNA < 0.95) {
      r = tryCatch({
        terra::trim(r)
      }, error = function(e) {
        stop("Failed to trim image: ", e$message)
      })
    } else {
      terra::values(r) = NA
      if (!silent) {
        cat(sprintf("Depth: %d cm. Not enough informative pixels. %.1f%% are NAs in this Depth Slice.\n",
                    indexD, coverNA*100))
      }
    }

    return(r)

  }, error = function(e) {
    stop(paste("Error in zone.fun:", e$message))
  }, warning = function(w) {
    warning(paste("Warning in zone.fun:", w$message))
  })
}



#' Extract zones based on rotation axis
#'
#' @param img array/matrix/SpatRaster - source image
#' @param k numeric vector - which cuts to keep (length 2)
#' @param kk numeric - total number of cuts
#' @param mm numeric vector - region limits (length 2)
#' @param select.layer Integer. Specifies which layer to use if the input is a multi-band image. Default is `NULL`.
#' @return SpatRaster - extracted rotation zone
#' @export
#'
#' @examples
#' data(seg_Oulanka2023_Session01_T067)
#' img = terra::rast(seg_Oulanka2023_Session01_T067)
#' rotationZone = zone.rotation.fun(img, k = c(1,2), kk = 7, mm = c(1500,3000))
zone.rotation.fun = function(img, k=c(3,4), kk=5, mm=c(2000,4000), select.layer=NULL) {
  # Input validation module
  validate_inputs <- function() {
    # Check rootpic
    if (is.null(img)) {
      stop("Input image 'rootpic' cannot be NULL")
    }

    # Validate k parameter
    if (!is.numeric(k) || length(k) != 2) {
      stop("'k' must be a numeric vector of length 2")
    }
    if (k[1] >= k[2]) {
      stop("First value of 'k' must be less than second value")
    }

    # Validate kk parameter
    if (!is.numeric(kk) || length(kk) != 1 || kk <= 0) {
      stop("'kk' must be a positive numeric value")
    }
    if (any(k > kk)) {
      stop("Values in 'k' cannot be greater than 'kk'")
    }

    # Validate mm parameter
    if (!is.null(mm)) {
      if (!is.numeric(mm) || length(mm) != 2) {
        stop("'mm' must be a numeric vector of length 2")
      }
      if (mm[1] >= mm[2]) {
        stop("First value of 'mm' must be less than second value")
      }
    }

    # Validate select.layer
    if (!is.null(select.layer)) {
      if (!is.numeric(select.layer) || length(select.layer) != 1 || select.layer < 1) {
        stop("'select.layer' must be a positive integer")
      }
    }
  }

  # Edge case handling module
  handle_edge_cases <- function(img0) {
    # Check for empty or all-NA image
    if (all(is.na(img0))) {
      warning("Input image contains only NA values")
      return(NULL)
    }

    # Check dimensions
    if (length(dim(img0)) < 2) {
      stop("Input image must have at least 2 dimensions")
    }

    # Check if mm range is within image dimensions
    if (!is.null(mm)) {
      if (mm[2] > dim(img0)[2]) {
        warning("'mm' range exceeds image dimensions, truncating to image size")
        mm[2] <- dim(img0)[2]
      }
    }

    return(img0)
  }

  # Main function execution with try-catch
  tryCatch({
    # Validate inputs
    validate_inputs()


    img0 <- load_flexible_image(img, select.layer = select.layer, output_format = "array", normalize = FALSE  )

    # Handle edge cases
    img0 <- handle_edge_cases(img0)
    if (is.null(img0)) return(NULL)


    if(exists("mm")){
      m = mm[1]:mm[2]
      ## rotation slice test
      img0 = img0[,m,]
    }

    # lower row for the kth bin
    q = floor((k[1]*(nrow(img0)/kk)))
    q = q +1
    p = floor((k[2]*(nrow(img0)/kk)))

    img00 = img0


    if(is.na(dim(img0)[3] > 1 )){
      img00[1:dim(img00)[1],1:dim(img00)[2]] <- NA
    }else{
      img00[1:dim(img00)[1],1:dim(img00)[2],] <- NA
    }

    if(is.na(dim(img0)[3] > 1 )){
      img00[c(p:q),] = img0[c(p:q),]
    }else{
      img00[c(p:q),,] = img0[c(p:q),,]
    }


    img00 = terra::rast(img00)
    img00 = terra::trim(img00)



    return(img00)

  }, error = function(e) {
    stop(paste("Error in zone.rotation.fun:", e$message))
  }, warning = function(w) {
    warning(paste("Warning in zone.rotation.fun:", w$message))
  })
}
