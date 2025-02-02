

#' Calculate a circular mean to determine average Directionality
#'
#' @param angles Numeric vector of input angles
#' @param input_units Character string specifying input units ("radians" or "degrees")
#' @param output_units Character string specifying output units ("radians" or "degrees")
#'
#' @return Numeric value representing the average angle
#' @export
#'
#' @examples circular_mean(angles = c(360,90,0), input_units = "degrees", output_units = "degrees")
circular_mean <- function(angles, input_units = "degrees", output_units = "degrees") {
  # Convert angles to radians if they are in degrees
  if (input_units == "degrees") {
    angles <- angles * pi / 180
  } else if (input_units != "radians") {
    stop("Invalid input unit specified. Use 'radians' or 'degrees'.")
  }

  # Calculate the sine and cosine of the angles
  sin_sum <- sum(sin(angles))
  cos_sum <- sum(cos(angles))

  # Calculate the circular mean
  mean_angle <- atan2(sin_sum, cos_sum)

  # Ensure the mean angle is in the range [0, 2*pi)
  if (mean_angle < 0) {
    mean_angle <- mean_angle + 2 * pi
  }

  # Convert the mean angle to the desired output unit
  if (output_units == "degrees") {
    mean_angle <- mean_angle * 180 / pi
  } else if (output_units != "radians") {
    stop("Invalid output unit specified. Use 'radians' or 'degrees'.")
  }

  return(mean_angle)
}






#' Threshold an image to rebinarize a blurred image
#'
#' @param img SpatRaster object
#' @param threshold Numeric value between 0 and 1 for thresholding
#' @param focus.layer Integer specifying which layer should be used to capture blur
#' @param mask.layer Integer specifying which layer preserves mask sections
#'
#' @return SpatRaster object
#' @export
#'
#' @examples
#' blurred.img = terra::rast(seg_Oulanka2023_Session03_T067)
#' img = blur.correction(blurred.img, 0.3)
blur.correction = function(img,threshold = 0.4,focus.layer = 2,mask.layer = 1){
  img2 = img

  if(length(dim(img)) == 3){
    mx = terra::global(img2[[focus.layer]], "max")[[1]]
    # isolate mask part
    img2[[mask.layer]] = img2[[mask.layer]]-img2[[focus.layer]]
    img2[[mask.layer]] = (img2[[mask.layer]] >= (mx * threshold)) * mx
    # unblurr focus layer
    img2[[focus.layer]] = (img2[[focus.layer]] >= (mx * threshold)) * mx
    # unblurr remaining layer
    other.layer = which(!1:3 %in% focus.layer & !1:3 %in% mask.layer)
    img2[[other.layer]] = (img2[[focus.layer]] >= (mx * threshold)) * mx

  }else{
    mx = terra::global(img2, "max")[[1]]
    # unblurr focus layer
    img2 = (img2 >= (mx * threshold)) * mx
  }

  return(img2)
}






#' Calculate root accumulation
#'
#' @param x Data frame containing group, depth, and variable columns
#' @param group Character vector specifying grouping variable(s)
#' @param depth Character string specifying depth column name
#' @param variable Character string specifying accumulating values column
#' @param stdrz Character string specifying standardization method
#'
#' @return Numeric vector of accumulated values
#'
#' @examples
#'df = data.frame(depth = c(seq(0,80,20),seq(0,80,20)),
#'                Plot = c(rep("a",5),rep("b",5)), rootpx = c(5,50,20,15,5,10,40,30,10,5) )
#' accum_root = root.accumulation(df,group = "Plot", depth = "depth", variable = "rootpx")
root.accumulation = function(x,group,depth,variable,stdrz = "counts"){
  # Split data by group
  split_df <- split(x, x[,group])

  # Initialize an empty list to store results
  result_list <- list()

  # Loop over each group
  for (group in names(split_df)) {
    # Sort the data within the group by depth
    sorted_group <- split_df[[group]][order(split_df[[group]][[depth]]), ]

    # Compute cumulative sum of variable
    #sorted_group$cs <- cumsum(sorted_group[[variable]])
    cs <- cumsum(dplyr::coalesce(sorted_group[[variable]], 0)) + sorted_group[[variable]]*0
    # root accumulation distribution
    if(stdrz == "counts"){
      cs <- cs
    }
    if(stdrz == "additive"){
      mx.roots = max(cs,na.rm=TRUE)
      cs <- cs / mx.roots
    }
    if(stdrz == "relative"){
      sm.roots = sum(cs,na.rm=TRUE)
      cs <- cs / sm.roots
    }
    sorted_group$cs <- cs

    # Append the result to the list
    result_list[[group]] <- sorted_group
  }

  # Combine the list back into a single data frame
  result_df <- do.call(rbind, result_list)

  # Reorder the result to match the original row order
  result_df <- result_df[order(rownames(result_df)), ]

  # only the accumulation values
  out = result_df$cs

  return(out)

}



#' Convert RGB image to grayscale with optimized memory management and parallel processing
#'
#' @param img SpatRaster RGB image
#' @param r Weight for red channel
#' @param g Weight for green channel
#' @param b Weight for blue channel
#'
#' @examples
#' data(seg_Oulanka2023_Session01_T067)
#' img = seg_Oulanka2023_Session01_T067
#' gray.raster = rgb2gray(img)
rgb2gray = function(img, r=0.21,g=0.72,b=0.07){
  gray.im = img[[1]] * r + img[[2]] * g + img[[3]] * b
  return(gray.im)

}








#' Combine multi-dimensional arrays
#'
#' @description
#'Combine multi-dimensional arrays.  This is a
#'generalization of cbind and rbind.  Takes a sequence of
#'vectors, matrices, or arrays and produces a single array of
#'the same or higher dimension.
#'
#' @param ...  Any number of vectors, matrices, arrays, or data frames.
#'The dimensions of all the arrays must match, except on one dimension
#'(specified by \code{along=}).  If these arguments are named, the name
#'will be used for the name of the dimension along which the arrays are
#'joined.  Vectors are treated as having a dim attribute of length one.
#'
#'Alternatively, there can be one (and only one) list argument supplied,
#'whose components are the objects to be bound together.  Names of the
#'list components are treated in the same way as argument names.
#' @param along (optional) The dimension along which to bind the arrays.
#'The default is the last dimension, i.e., the maximum length of the dim
#'attribute of the supplied arrays.  \code{along=} can take any
#'non-negative value up to the minimum length of the dim attribute of
#'supplied arrays plus one.  When \code{along=} has a fractional value, a
#'value less than 1, or a value greater than N (N is the maximum of the
#'lengths of the dim attribute of the objects to be bound together), a new
#'dimension is created in the result.  In these cases, the dimensions of
#'all arguments must be identical.
#' @param rev.along (optional)
#'Alternate way to specify the dimension along which to bind the arrays:
#'  \code{along = N + 1 - rev.along}.  This is provided mainly to allow easy
#'specification of \code{along = N + 1} (by supplying
#'                                       \code{rev.along=0}).  If both \code{along} and \code{rev.along} are
#'supplied, the supplied value of \code{along} is ignored.
#' @param new.names (optional)
#'If new.names is a list, it is the first choice for the
#'dimnames attribute of the result.  It should have the same
#'structure as a dimnames attribute.  If the names for a
#'particular dimension are \code{NULL}, names for this dimension are
#'constructed in other ways.
#'
#'If \code{new.names} is a character vector, it is used for dimension
#'names in the same way as argument names are used.  Zero
#'length ("") names are ignored.
#' @param force.array (optional) If \code{FALSE}, rbind or cbind are
#'called when possible, i.e., when the arguments are all vectors, and
#'along is not 1, or when the arguments are vectors or matrices or data
#'frames and along is 1 or 2.  If rbind or cbind are used, they will
#'preserve the data.frame classes (or any other class that r/cbind
#'preserve).  Otherwise, abind will convert objects to class array.  Thus,
#'to guarantee that an array object is returned, supply the argument
#'\code{force.array=TRUE}.  Note that the use of rbind or cbind introduces
#'some subtle changes in the way default dimension names are constructed:
#'  see the examples below.
#' @param make.names (optional)
#'If \code{TRUE}, the last resort for dimnames for the along
#'dimension will be the deparsed versions of anonymous
#'arguments.  This can result in cumbersome names when
#'arguments are expressions.
#'
#'<p>The default is \code{FALSE}.
#' @param use.anon.names (optional)
#'\code{use.anon.names}
#'is a deprecated synonym for \code{make.names}.
#' @param use.first.dimnames (optional)
#'When dimension names are present on more than one
#'argument, should dimension names for the result be take from
#'the first available (the default is to take them from the
#'                     last available, which is the same behavior as
#'                     \code{rbind} and \code{cbind}.)
#' @param hier.names (optional)
#'If \code{TRUE}, dimension names on the concatenated dimension will be
#'composed of the argument name and the dimension names of the objects
#'being bound.  If a single list argument is supplied, then the names of
#'the components serve as the argument names.  \code{hier.names} can
#'also have values \code{"before"} or \code{"after"}; these determine
#'the order in which the argument name and the dimension name are put
#'together (\code{TRUE} has the same effect as \code{"before"}).
#' @param use.dnns (default \code{FALSE}) Use names on dimensions, e.g.,
#'so that \code{names(dimnames(x))} is non-empty.  When there are
#'multiple possible sources for names of dimnames, the value of
#'\code{use.first.dimnames} determines the result.
#'
#' @details
#' The dimensions of the supplied vectors or arrays do not need
#'to be identical, e.g., arguments can be a mixture of vectors
#'and matrices.  \code{abind} coerces arguments by the addition
#'of one dimension in order to make them consistent with other
#'arguments and \code{along=}.  The extra dimension is
#'added in the place specified by \code{along=}.
#'
#'The default action of abind is to concatenate on the last
#'dimension, rather than increase the number of dimensions.
#'For example, the result of calling abind with vectors is a
#'longer vector (see first example below).  This differs from
#'the action of \code{rbind} and cbind which is to return a matrix when
#'called with vectors.  abind can be made to behave like cbind
#'on vectors by specifying \code{along=2}, and like rbind by
#'specifying \code{along=0}.
#'
#'The dimnames of the returned object are pieced together
#'from the dimnames of the arguments, and the names of the
#'arguments.  Names for each dimension are searched for in the
#'following order: new.names, argument name, dimnames (or
#'names) attribute of last argument, dimnames (or names)
#'attribute of second last argument, etc.  (Supplying the
#'                                          argument \code{use.first.dimnames=TRUE} changes this to
#'                                          cause \code{abind} to use dimnames or names from the
#'                                          first argument first.  The default behavior is the same as
#'                                          for \code{rbind} and \code{cbind}: use dimnames
#'                                          from later arguments.)  If some names are supplied for the
#'along dimension (either as argument names or dimnames in
#'                 arguments), names are constructed for anonymous arguments
#'unless \code{use.anon.names=FALSE}.
#'
#' @author Tony Plate \email{tplate@acm.org} and Richard Heiberger
#' @import utils
#' @return merged multidimensional arrays
#' @keywords internal


abind2 = function (..., along = N, rev.along = NULL, new.names = NULL,
                   force.array = TRUE, make.names = use.anon.names, use.anon.names = FALSE,
                   use.first.dimnames = FALSE, hier.names = FALSE, use.dnns = FALSE)
{
  if (is.character(hier.names))
    hier.names <- match.arg(hier.names, c("before", "after",
                                          "none"))
  else hier.names <- if (hier.names)
    "before"
  else "no"
  arg.list <- list(...)
  if (is.list(arg.list[[1]]) && !is.data.frame(arg.list[[1]])) {
    if (length(arg.list) != 1)
      stop("can only supply one list-valued argument for ...")
    if (make.names)
      stop("cannot have make.names=TRUE with a list argument")
    arg.list <- arg.list[[1]]
    have.list.arg <- TRUE
  }
  else {
    N <- max(1, sapply(list(...), function(x) length(dim(x))))
    have.list.arg <- FALSE
  }
  if (any(discard <- sapply(arg.list, is.null)))
    arg.list <- arg.list[!discard]
  if (length(arg.list) == 0)
    return(NULL)
  N <- max(1, sapply(arg.list, function(x) length(dim(x))))
  if (!is.null(rev.along))
    along <- N + 1 - rev.along
  if (along < 1 || along > N || (along > floor(along) && along <
                                 ceiling(along))) {
    N <- N + 1
    along <- max(1, min(N + 1, ceiling(along)))
  }
  if (length(along) > 1 || along < 1 || along > N + 1)
    stop(paste("\"along\" must specify one dimension of the array,",
               "or interpolate between two dimensions of the array",
               sep = "\n"))
  if (!force.array && N == 2) {
    if (!have.list.arg) {
      if (along == 2)
        return(cbind(...))
      if (along == 1)
        return(rbind(...))
    }
    else {
      if (along == 2)
        return(do.call("cbind", arg.list))
      if (along == 1)
        return(do.call("rbind", arg.list))
    }
  }
  if (along > N || along < 0)
    stop("along must be between 0 and ", N)
  pre <- seq(from = 1, len = along - 1)
  post <- seq(to = N - 1, len = N - along)
  perm <- c(seq(len = N)[-along], along)
  arg.names <- names(arg.list)
  if (is.null(arg.names))
    arg.names <- rep("", length(arg.list))
  if (is.character(new.names)) {
    arg.names[seq(along = new.names)[nchar(new.names) > 0]] <- new.names[nchar(new.names) >
                                                                           0]
    new.names <- NULL
  }
  if (any(arg.names == "")) {
    if (make.names) {
      dot.args <- match.call(expand.dots = FALSE)$...
      if (is.call(dot.args) && identical(dot.args[[1]],
                                         as.name("list")))
        dot.args <- dot.args[-1]
      arg.alt.names <- arg.names
      for (i in seq(along = arg.names)) {
        if (arg.alt.names[i] == "") {
          if (utils::object.size(dot.args[[i]]) < 1000) {
            arg.alt.names[i] <- paste(deparse(dot.args[[i]],
                                              40), collapse = ";")
          }
          else {
            arg.alt.names[i] <- paste("X", i, sep = "")
          }
          arg.names[i] <- arg.alt.names[i]
        }
      }
    }
    else {
      arg.alt.names <- arg.names
      arg.alt.names[arg.names == ""] <- paste("X", seq(along = arg.names),
                                              sep = "")[arg.names == ""]
    }
  }
  else {
    arg.alt.names <- arg.names
  }
  use.along.names <- any(arg.names != "")
  names(arg.list) <- arg.names
  arg.dimnames <- matrix(vector("list", N * length(arg.names)),
                         nrow = N, ncol = length(arg.names))
  dimnames(arg.dimnames) <- list(NULL, arg.names)
  arg.dnns <- matrix(vector("list", N * length(arg.names)),
                     nrow = N, ncol = length(arg.names))
  dimnames(arg.dnns) <- list(NULL, arg.names)
  dimnames.new <- vector("list", N)
  arg.dim <- matrix(integer(1), nrow = N, ncol = length(arg.names))
  for (i in seq(len = length(arg.list))) {
    m <- arg.list[[i]]
    m.changed <- FALSE
    if (is.data.frame(m)) {
      m <- as.matrix(m)
      m.changed <- TRUE
    }
    else if (!is.array(m) && !is.null(m)) {
      if (!is.atomic(m))
        stop("arg '", arg.alt.names[i], "' is non-atomic")
      dn <- names(m)
      m <- as.array(m)
      if (length(dim(m)) == 1 && !is.null(dn))
        dimnames(m) <- list(dn)
      m.changed <- TRUE
    }
    new.dim <- dim(m)
    if (length(new.dim) == N) {
      if (!is.null(dimnames(m))) {
        arg.dimnames[, i] <- dimnames(m)
        if (use.dnns && !is.null(names(dimnames(m))))
          arg.dnns[, i] <- as.list(names(dimnames(m)))
      }
      arg.dim[, i] <- new.dim
    }
    else if (length(new.dim) == N - 1) {
      if (!is.null(dimnames(m))) {
        arg.dimnames[-along, i] <- dimnames(m)
        if (use.dnns && !is.null(names(dimnames(m))))
          arg.dnns[-along, i] <- as.list(names(dimnames(m)))
        dimnames(m) <- NULL
      }
      arg.dim[, i] <- c(new.dim[pre], 1, new.dim[post])
      if (any(perm != seq(along = perm))) {
        dim(m) <- c(new.dim[pre], 1, new.dim[post])
        m.changed <- TRUE
      }
    }
    else {
      stop("'", arg.alt.names[i], "' does not fit: should have `length(dim())'=",
           N, " or ", N - 1)
    }
    if (any(perm != seq(along = perm)))
      arg.list[[i]] <- aperm(m, perm)
    else if (m.changed)
      arg.list[[i]] <- m
  }
  conform.dim <- arg.dim[, 1]
  for (i in seq(len = ncol(arg.dim))) {
    if (any((conform.dim != arg.dim[, i])[-along])) {
      stop("arg '", arg.alt.names[i], "' has dims=", paste(arg.dim[,
                                                                   i], collapse = ", "), "; but need dims=", paste(replace(conform.dim,
                                                                                                                           along, "X"), collapse = ", "))
    }
  }
  if (N > 1)
    for (dd in seq(len = N)[-along]) {
      for (i in (if (use.first.dimnames)
        seq(along = arg.names)
        else rev(seq(along = arg.names)))) {
        if (length(arg.dimnames[[dd, i]]) > 0) {
          dimnames.new[[dd]] <- arg.dimnames[[dd, i]]
          if (use.dnns && !is.null(arg.dnns[[dd, i]]))
            names(dimnames.new)[dd] <- arg.dnns[[dd,
                                                 i]]
          break
        }
      }
    }
  for (i in seq(len = length(arg.names))) {
    if (arg.dim[along, i] > 0) {
      dnm.along <- arg.dimnames[[along, i]]
      if (length(dnm.along) == arg.dim[along, i]) {
        use.along.names <- TRUE
        if (hier.names == "before" && arg.names[i] !=
            "")
          dnm.along <- paste(arg.names[i], dnm.along,
                             sep = ".")
        else if (hier.names == "after" && arg.names[i] !=
                 "")
          dnm.along <- paste(dnm.along, arg.names[i],
                             sep = ".")
      }
      else {
        if (arg.dim[along, i] == 1)
          dnm.along <- arg.names[i]
        else if (arg.names[i] == "")
          dnm.along <- rep("", arg.dim[along, i])
        else dnm.along <- paste(arg.names[i], seq(length = arg.dim[along,
                                                                   i]), sep = "")
      }
      dimnames.new[[along]] <- c(dimnames.new[[along]],
                                 dnm.along)
    }
    if (use.dnns) {
      dnn <- unlist(arg.dnns[along, ])
      if (length(dnn)) {
        if (!use.first.dimnames)
          dnn <- rev(dnn)
        names(dimnames.new)[along] <- dnn[1]
      }
    }
  }
  if (!use.along.names)
    dimnames.new[along] <- list(NULL)
  out <- array(unlist(arg.list, use.names = FALSE), dim = c(arg.dim[-along,
                                                                    1], sum(arg.dim[along, ])), dimnames = dimnames.new[perm])
  if (any(order(perm) != seq(along = perm)))
    out <- aperm(out, order(perm))
  if (!is.null(new.names) && is.list(new.names)) {
    for (dd in seq(len = N)) {
      if (!is.null(new.names[[dd]])) {
        if (length(new.names[[dd]]) == dim(out)[dd])
          dimnames(out)[[dd]] <- new.names[[dd]]
        else if (length(new.names[[dd]]))
          warning(paste("Component ", dd, " of new.names ignored: has length ",
                        length(new.names[[dd]]), ", should be ",
                        dim(out)[dd], sep = ""))
      }
      if (use.dnns && !is.null(names(new.names)) && names(new.names)[dd] !=
          "")
        names(dimnames(out))[dd] <- names(new.names)[dd]
    }
  }
  if (use.dnns && !is.null(names(dimnames(out))) && any(i <- is.na(names(dimnames(out)))))
    names(dimnames(out))[i] <- ""
  out
}







    #' Get list of supported formats as a string
    #' @keywords internal
    supported_formats_string <- function() {
      paste(
        "Supported formats:",
        "- File path (.jpg, .jpeg, .png, .tif, .tiff, .bmp)",
        "- SpatRaster",
        "- RasterLayer",
        "- RasterBrick",
        "- RasterStack",
        "- cimg",
        "- magick-image",
        "- matrix",
        "- array",
        sep = "\n"
      )
    }

    #' Helper function to convert various input formats to terra SpatRaster
    #' Helper function to convert various input formats to terra SpatRaster
    #' @param input Any supported image format
    #' @param normalize Logical, whether to normalize values to 0-1 range if they're in 0-255
    #' @param select.layer Numeric, which layer to select if input has multiple layers
    #' @keywords internal
    convert_to_spatrast <- function(input, normalize = TRUE, select.layer = NULL) {
      if (inherits(input, "SpatRaster")) {
        out <- input
      }
      else if (is.character(input) && file.exists(input)) {
        ext <- tolower(tools::file_ext(input))
        if (ext %in% c("tif", "tiff")) {
          out <- terra::rast(input)
        } else {
          arr <- as.array(imager::load.image(input))
          # Handle 4D arrays from imager
          if (length(dim(arr)) == 4) {
            arr <- arr[,,,1]  # Take first frame if multiple frames
          }
          out <- terra::rast(arr)
        }
      }
      else if (is.matrix(input)) {
        out <- terra::rast(array(input, dim = c(dim(input), 1)))
      }
      else if (is.array(input)) {
        dims <- dim(input)
        if (length(dims) == 2) {
          out <- terra::rast(array(input, dim = c(dims, 1)))
        } else if (length(dims) == 3) {
          out <- terra::rast(input)
        } else if (length(dims) == 4) {
          out <- terra::rast(input[,,,1])  # Take first frame if multiple frames
        }
      }
      else if (inherits(input, c("RasterLayer", "RasterBrick", "RasterStack"))) {
        out <- terra::rast(input)
      }
      else if (inherits(input, "magick-image")) {
        arr <- as.array(input[[1]])
        if (length(dim(arr)) == 4) {
          arr <- arr[,,,1]  # Take first frame if multiple frames
        }
        out <- terra::rast(arr)
      }
      else if (inherits(input, "cimg")) {
        arr <- as.array(input)
        if (length(dim(arr)) == 4) {
          arr <- arr[,,,1]  # Take first frame if multiple frames
        }
        out <- terra::rast(arr)
      }
      else {
        stop("Unsupported input format.\n", supported_formats_string())
      }

      if (!is.null(select.layer) && terra::nlyr(out) > 1) {
        out <- out[[select.layer]]
      }

      if (normalize && terra::global(out, "max", na.rm = TRUE)$max[1] > 1) {
        out <- out / 255
      }

      return(out)
    }

    #' Helper function to convert various input formats to RasterBrick
    #' @param input Any supported image format
    #' @param normalize Logical, whether to normalize values to 0-1 range if they're in 0-255
    #' @param select.layer Numeric, which layer to select if input has multiple layers
    #' @keywords internal
    convert_to_brick <- function(input, normalize = TRUE, select.layer = NULL) {
      if (inherits(input, "RasterBrick")) {
        out <- input
      }
      else if (is.character(input) && file.exists(input)) {
        ext <- tolower(tools::file_ext(input))
        if (ext %in% c("tif", "tiff")) {
          out <- raster::brick(input)
        } else {
          arr <- as.array(imager::load.image(input))
          if (length(dim(arr)) == 4) {
            arr <- arr[,,,1]  # Take first frame if multiple frames
          }
          out <- raster::brick(arr)
        }
      }
      else if (is.matrix(input)) {
        out <- raster::brick(array(input, dim = c(dim(input), 1)))
      }
      else if (is.array(input)) {
        dims <- dim(input)
        if (length(dims) == 2) {
          out <- raster::brick(array(input, dim = c(dims, 1)))
        } else if (length(dims) == 3) {
          out <- raster::brick(input)
        } else if (length(dims) == 4) {
          out <- raster::brick(input[,,,1])  # Take first frame if multiple frames
        }
      }
      else if (inherits(input, c("RasterLayer", "RasterStack"))) {
        out <- raster::brick(input)
      }
      else if (inherits(input, "SpatRaster")) {
        out <- raster::brick(raster::raster(input))
      }
      else if (inherits(input, "magick-image")) {
        arr <- as.array(input[[1]])
        if (length(dim(arr)) == 4) {
          arr <- arr[,,,1]  # Take first frame if multiple frames
        }
        out <- raster::brick(arr)
      }
      else if (inherits(input, "cimg")) {
        arr <- as.array(input)
        if (length(dim(arr)) == 4) {
          arr <- arr[,,,1]  # Take first frame if multiple frames
        }
        out <- raster::brick(arr)
      }
      else {
        stop("Unsupported input format.\n", supported_formats_string())
      }

      if (!is.null(select.layer) && raster::nlayers(out) > 1) {
        out <- raster::brick(out[[select.layer]])
      }

      if (normalize && max(raster::maxValue(out)) > 1) {
        out <- out / 255
      }

      return(out)
    }

    #' Helper function to convert various input formats to RasterLayer
    #' @param input Any supported image format
    #' @param normalize Logical, whether to normalize values to 0-1 range if they're in 0-255
    #' @param select.layer Numeric, which layer to select if input has multiple layers
    #' @keywords internal
    convert_to_raster <- function(input, normalize = TRUE, select.layer = NULL) {
      if (inherits(input, "RasterLayer")) {
        out <- input
      }
      else if (is.character(input) && file.exists(input)) {
        ext <- tolower(tools::file_ext(input))
        if (ext %in% c("tif", "tiff")) {
          out <- raster::raster(input)
        } else {
          arr <- as.array(imager::load.image(input))
          layer_idx <- if(is.null(select.layer)) 1 else select.layer
          if (length(dim(arr)) == 4) {
            arr <- arr[,,,1]  # Take first frame if multiple frames
          }
          out <- raster::raster(arr[,,layer_idx])
        }
      }
      else if (is.matrix(input)) {
        out <- raster::raster(input)
      }
      else if (is.array(input)) {
        dims <- dim(input)
        layer_idx <- if(is.null(select.layer)) 1 else select.layer
        if (length(dims) == 2) {
          out <- raster::raster(input)
        } else if (length(dims) == 3) {
          out <- raster::raster(input[,,layer_idx])
        } else if (length(dims) == 4) {
          out <- raster::raster(input[,,layer_idx,1])  # Take first frame
        }
      }
      else if (inherits(input, c("RasterBrick", "RasterStack"))) {
        layer_idx <- if(is.null(select.layer)) 1 else select.layer
        out <- input[[layer_idx]]
      }
      else if (inherits(input, "SpatRaster")) {
        layer_idx <- if(is.null(select.layer)) 1 else select.layer
        out <- raster::raster(input[[layer_idx]])
      }
      else if (inherits(input, "magick-image")) {
        arr <- as.array(input[[1]])
        layer_idx <- if(is.null(select.layer)) 1 else select.layer
        if (length(dim(arr)) == 4) {
          arr <- arr[,,,1]  # Take first frame if multiple frames
        }
        out <- raster::raster(arr[,,layer_idx])
      }
      else if (inherits(input, "cimg")) {
        arr <- as.array(input)
        layer_idx <- if(is.null(select.layer)) 1 else select.layer
        if (length(dim(arr)) == 4) {
          arr <- arr[,,,1]  # Take first frame if multiple frames
        }
        out <- raster::raster(arr[,,layer_idx])
      }
      else {
        stop("Unsupported input format.\n", supported_formats_string())
      }

      if (normalize && max(raster::maxValue(out)) > 1) {
        out <- out / 255
      }

      return(out)
    }

    #' Helper function to convert various input formats to array
    #' @param input Any supported image format
    #' @param normalize Logical, whether to normalize values to 0-1 range if they're in 0-255
    #' @param select.layer Numeric, which layer to select if input has multiple layers
    #' @keywords internal
    convert_to_array <- function(input, normalize = TRUE, select.layer = NULL) {
      if (is.array(input)) {
        dims <- dim(input)
        if (length(dims) == 4) {
          out <- input[,,,1]  # Take first frame from 4D array
        } else {
          out <- input  # Keep 2D or 3D as is
        }
      }
      else if (is.matrix(input)) {
        out <- array(input, dim = c(dim(input), 1))  # Convert to 3D
      }
      else if (is.character(input) && file.exists(input)) {
        ext <- tolower(tools::file_ext(input))
        if (ext %in% c("tif", "tiff")) {
          out <- as.array(tiff::readTIFF(input))  # Already 3D
        } else {
          arr <- as.array(imager::load.image(input))  # 4D array
          out <- arr[,,,1]  # Convert to 3D by taking first frame
        }
      }
      else if (inherits(input, "SpatRaster")) {
        out <- terra::as.array(input)  # Already 3D
      }
      else if (inherits(input, c("RasterLayer", "RasterBrick", "RasterStack"))) {
        out <- as.array(terra::rast(input))  # Already 3D
      }
      else if (inherits(input, "magick-image")) {
        arr <- as.array(input[[1]])  # 4D array
        out <- arr[,,,1]  # Convert to 3D
      }
      else if (inherits(input, "cimg")) {
        arr <- as.array(input)  # 4D array
        out <- arr[,,,1]  # Convert to 3D
      }
      else {
        stop("Unsupported input format.\n", supported_formats_string())
      }

      # Ensure 2D becomes 3D
      if (length(dim(out)) == 2) {
        out <- array(out, dim = c(dim(out), 1))
      }

      # Handle layer selection for 3D arrays
      if (!is.null(select.layer) && dim(out)[3] > 1) {
        out <- out[,,select.layer,drop=FALSE]
      }

      if (normalize && max(out, na.rm = TRUE) > 1) {
        out <- out / 255
      }

      return(out)
    }

    #' Helper function to convert various input formats to cimg
    #' @param input Any supported image format
    #' @param normalize Logical, whether to normalize values to 0-1 range if they're in 0-255
    #' @param select.layer Numeric, which layer to select if input has multiple layers
    #' @keywords internal
    convert_to_cimg <- function(input, normalize = TRUE, select.layer = NULL) {
      if (inherits(input, "cimg")) {
        out <- input
      }
      else if (is.character(input) && file.exists(input)) {
        ext <- tolower(tools::file_ext(input))
        if (ext %in% c("tif", "tiff")) {
          arr <- tiff::readTIFF(input)  # 3D array
          out <- imager::as.cimg(arr)  # Converts to 4D (adds dim=1)
        } else {
          out <- imager::load.image(input)  # Already 4D
        }
      }
      else if (is.matrix(input)) {
        # Convert 2D to 4D (adding color channel and frame dimensions)
        out <- imager::as.cimg(array(input, dim = c(dim(input), 1, 1)))
      }
      else if (is.array(input)) {
        dims <- dim(input)
        if (length(dims) == 2) {
          # 2D to 4D (adding color channel and frame dimensions)
          out <- imager::as.cimg(array(input, dim = c(dims, 1, 1)))
        } else if (length(dims) == 3) {
          # 3D to 4D (adding frame dimension)
          out <- imager::as.cimg(array(input, dim = c(dims, 1)))
        } else if (length(dims) == 4) {
          # Take first frame if multiple frames
          out <- imager::as.cimg(input[,,,1,drop=FALSE])
        }
      }
      else if (inherits(input, "SpatRaster")) {
        arr <- as.array(input)  # Get 3D array
        out <- imager::as.cimg(array(arr, dim = c(dim(arr), 1)))  # Convert to 4D
      }
      else if (inherits(input, c("RasterLayer", "RasterBrick", "RasterStack"))) {
        arr <- as.array(terra::rast(input))  # Get 3D array
        out <- imager::as.cimg(array(arr, dim = c(dim(arr), 1)))  # Convert to 4D
      }
      else if (inherits(input, "magick-image")) {
        arr <- as.array(input[[1]])  # Get 4D array
        out <- imager::as.cimg(arr[,,,1,drop=FALSE])  # Keep as 4D with single frame
      }
      else {
        stop("Unsupported input format.\n", supported_formats_string())
      }

      # Handle layer selection (for color channels)
      if (!is.null(select.layer) && dim(out)[3] > 1) {
        out <- out[,,select.layer,,drop=FALSE]  # Keep 4D structure
      }

      if (normalize && max(out, na.rm = TRUE) > 1) {
        out <- out / 255
      }

      return(out)
    }

#' Load an image flexibly from file or convert from memory
#' @param input File path or image object
#' @param normalize Logical, whether to normalize values to 0-1 range if they're in 0-255
#' @param output_format Character, one of "cimg", "spatrast", "array", "brick", "raster"
#' @param select.layer Numeric, which layer to select if input has multiple layers
#' @keywords internal
load_flexible_image <- function(input, normalize = TRUE, output_format = "cimg", select.layer = NULL) {
  # And in load_flexible_image():
  converter <- switch(output_format,
                      "cimg" = convert_to_cimg,
                      "spatrast" = convert_to_spatrast,
                      "array" = convert_to_array,
                      "brick" = convert_to_brick,
                      "raster" = convert_to_raster,
                      stop("Unsupported output format.\n", supported_formats_string()))

  return(converter(input, normalize = normalize, select.layer = select.layer))
}
