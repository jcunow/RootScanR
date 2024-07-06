## helper functions

#' Skeletonize a segmented image
#'
#' @param img binary image, can be path or magick object
#' @param itr iterations of thinning. see magick::image_morphology
#'
#' @return image in as magick object
#' @export
#'
#' @examples img.skeleton = skeletonize("/path/image.png",itr = 2)
skeletonize = function(img,itr = 2){
  if(!is.character(img)){
    img_array <- as.array(img)
    if(!(max(img) > 1)){
      img_array <- img_array * 255
    }

    # Convert to raw array (required by magick)
    magick_img <- imager::cimg2magick(imager::as.cimg( img_array))
  }else{
    magick_img = magick::image_read(img)
  }

  grayscale_img = magick::image_channel(magick_img,"lightness")

  # Perform morphological thinning
  thinned_img <- image_morphology(grayscale_img, 'Thinning', 'Skeleton:3',iterations = itr)

  print(thinned_img)

  return(thinned_img)
}

#' Root accumulation Curve
#'
#' @param data dataframe must include group,depth, and variable columns
#' @param group specify the grouping variable e.g., Plot
#' @param depth specifiy column name which includes depth values
#' @param variable accumulating values
#'
#' @return dataframe with one added column "cs" containing the cummulated values
#' @export
#'
#' @examples data1 = root.accumulation(data,group = Plot, depth = depth, variable = rootpx)
root.accumulation = function(data,group,depth,variable){
  pdf = data %>% group_by(group) %>% arrange(depth) %>% mutate(cs = cumsum(rootpx))
  return(pdf)
}



#' Converts RGB to Grayscale
#'
#' @param img rgb raster
#' @param r weight for red color
#' @param g weight for green color
#' @param b weight for blue color
#'
#' @return a single layer gray scale raster
#' @export
#'
#' @examples gray.raster = rgb2gray(img)
rgb2gray = function(img, r=0.21,g=0.72,b=0.07){
  gray.im = img[[1]] * r + img[[2]] * g + img[[3]] * b
  return(gray.im)

}



#' Resize the Image
#'
#' useful if some scanning campaigns have used different dpi
#' @param import.path Input Path
#' @param output.path Output Path
#' @param height target height
#' @param width target width
#'
#' @return a image with new dimensions and resolution
#' @export
#'
#' @examples resize.image(import.path, output.path) = resized.image
resize.image = function(import.path,output.path,height = 2550,width= 2273){
  img.list = list.files(import.path, pattern = ".tiff")

  # sequentially read images, lower resolution and restore them
  t1 = Sys.time()
  for (i in 1:length(img.list)) {
    # read
    temp.im = OpenImageR::readImage(paste0(import.path,img.list[i]))
    # resample bilinear
    temp.im = OpenImageR::resizeImage(temp.im,
                                      method = "bilinear",
                                      height = 2550, width = 2273)
    # transformation to 0-1 range
    temp.im = imagefx::range01(temp.im)
    # write as tiff
    OpenImageR::writeImage(temp.im,file_name = paste0(output.path,img.list[i]))

  }
  t2 = Sys.time();floor(t2-t1)

}




#' read rename write wrapper
#'
#' @param img image name
#' @param pattern pattern
#' @param replace replace
#' @param dir where to find the image
#' @param dir.out where to write the image
#'
#' @export
#'
#' @return image output
#'
#' @examples rrwr(Oulanka2023_T001_L001.tiff,pattern = "_L001", replace = "", dir = getwd(),dir.out = paste0(getwd(),"/renamed/"))
rrwr = function(img,pattern,replace,dir,dir.out){
  im = tiff::readTIFF(paste0(dir,img))
  name = stringr::str_replace(img, pattern = "Ecfg",replacement = "Oulanka_2020")
  tiff::writeTIFF(im,paste0(dir.out,"/",name))
}




#' Split Images
#'
#' Use if too many roots are present for the Root Detector AI image segmentation
#' @param path Input path
#' @param dir.out Output path
#' @param pattern only include images with this pattern
#' @param ratio split point 0-1
#'
#' @return two sets of images
#' @export
#'
#' @examples split_im(path, dir.out, ".tiff") = c(top.im, bot.im)
split_im = function(path,dir.out,pattern = ".tiff",ratio = 0.5){
  file.ls = list.files(path = path, pattern = pattern)
  ## split the image
  for (i in file.ls) {
    im = tiff::readTIFF(paste0(path,i))
    im.top = im[,1:(dim(im)[2]*ratio),]
    im.dwn = im[,(dim(im)[2]*ratio+1):dim(im)[2],]
    tiff::writeTIFF(im.top,where = paste0(dir.out,"Split_top_",i))
    tiff::writeTIFF(im.dwn,where = paste0(dir.out,"Split_dwn_",i))
  }
}



#' Fuse two Images
#'
#' this function is intended to to be used after using the root detector
#' @param path Input path
#' @param dir.out Output path
#' @param pattern only include images with this pattern
#'
#' @return one complete image
#' @export
#'
#' @examples join_im(path,out.path) = img
join_im = function(path,dir.out,pattern = "skeleton"){
  ### segmented
  dir.ls = paste0(list.dirs(path = path)[-1],"/")
  file.ls = list.files(path = dir.ls, pattern = pattern)
  top.files =  file.ls[stringr::str_detect(file.ls,pattern = "Split_top")]
  dwn.files =  file.ls[stringr::str_detect(file.ls,pattern = "Split_dwn")]

  for (i in 1:length(top.files)) {
    im.top = png::readPNG(paste0(dir.ls[i+36],top.files[i]))
    im.dwn = png::readPNG(paste0(dir.ls[i],dwn.files[i]))
    im.all = abind::abind(im.top,im.dwn, along = 2 )
    file.name= top.files[i] %>% stringr::str_remove(pattern="Split_top_") %>% stringr::str_remove(pattern=".tiff.segmentation.png")
    tiff::writeTIFF(im.all,where = paste0(dir.out,"FullSegmented_",file.name,".tiff"))

  }
}


# #' Skeletonize a segmented image
# #'
# #' @param img binary iamge
# #' @param method which approach should be used
# #'
# #' @return binary image with area turned into linear features
# #' @export
# #'
# #' @examples skeletonize(img,"gonzales") = skeleton_image
# skeletonize_function = function(img,method = "gonzales"){
#   if(method == "gonzales"){
#     skeleton = dipr::thinning(img)
#   }
#   if(method == "erode"){
#     skeleton = dipr::skeletonize(img)
#   }
#   if(method == "lantuejoul"){
#     skeleton = mmand::skeletonise(method = method )
#   }
#   if(method == "beucher"){
#     skeleton = mmand::skeletonise(method = method )
#   }
#   if(method == "hitormiss"){
#     skeleton = mmand::skeletonise(method = method )
#   }
#   if(method %in% any(c("gonzales","erode","lantuejoul","beucher","hitormiss"))  ){
#     print("method not appropriate")
#   }
# }

#


#' abind
#'
#' @param ... d
#' @param along dim
#' @param rev.along d
#' @param new.names d
#' @param force.array d
#' @param make.names d
#' @param use.anon.names d
#' @param use.first.dimnames d
#' @param hier.names d
#' @param use.dnns d
#' @name name d
#'
#' @return merged multidimensional arrays
#'
#' @examples marray = abind2()
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
          if (object.size(dot.args[[i]]) < 1000) {
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


## coulumn align

#' To cut top or bottom parts of a raster
#'
#' @param im image input
#' @param fixed.depth how much should be retained
#' @param top.cut discard the top part and retain the botom?
#'
#' @return a shorter raster
#' @export
#'
#' @examples im2 = DepthLimiter(im = image, fixed.depth = 7200, top.cut = T)
DepthLimiter = function(im,fixed.depth = 7200,top.cut = T){
  depth = dim(im)[2]

    im2 = raster::as.array(im)


  if(top.cut == T){
    im3 = im2[,(depth - fixed.depth-1) : depth,]
  }else{
    im3 = im2[,1 : (fixed.depth+1), ]
  }
  if(dim(im)[3]>1){
    im3 = raster::brick(im3)
    im3 = raster::raster(im3)
  }else{
    im3 = raster::raster(im3)
  }


  return(im3)
}




