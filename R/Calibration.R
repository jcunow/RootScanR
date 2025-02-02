#' Estimates rotation from tape coverage
#'
#' This function analyzes image data to determine rotation based on tape coverage,
#' assuming more tape is present on the upper side of the tube.
#'
#' @param img Input image as raster, file name, or array
#' @param tape.brightness Brightness threshold for tape detection (0-1)
#' @param tape.quantile Quantile used to align brightness with tape (0-1)
#' @param extra.rows Additional rows to add for analysis
#' @param search.area Proportion of image to analyze (0-1)
#' @param nclasses Number of classes for pixel clustering
#' @param select.layer Integer. Specifies which layer to use if the input is a multi-band image. Default is `NULL`.
#' @return numeric Position of the center of extruding tape
#' @export
#'
#' @examples
#' img = seg_Oulanka2023_Session01_T067
#' r0 = RotationE(img)
RotationE = function(img,tape.brightness = 0.66,extra.rows = 100,search.area = 0.45,tape.quantile = 0.98, nclasses = 3,select.layer = NULL){


  im <- load_flexible_image(img,select.layer = select.layer, output_format = "array", normalize = TRUE  )


  ## add one row of red tape pixel
  red.line = array(dim = c(dim(im)[1],extra.rows,dim(im)[3]))
  red.line[,,1:dim(im)[3]] <- stats::quantile(im[,,1],tape.quantile,na.rm=TRUE)
  img1 = abind2(red.line,im,along = 2) # adapted function from package:abind

  r.img1 = terra::rast(img1)
  ### make the search are smaller
  r.img1 = terra::crop(r.img1,raster::extent(0,search.area*terra::ext(r.img1)[2],0,terra::ext(r.img1)[4]))

  # identify distinct pixel groups
  r1 = RStoolbox::unsuperClass(r.img1,nClasses = nclasses)
  # determine average group
  clust.center = apply(r1$model$centers,1,mean)
  # silver tape should have highest luminance across clusters -> select max lum.cluster (but not close to == 1 [pure white?])
  clust= which(clust.center ==max(clust.center[clust.center > tape.brightness*terra::global(r.img1,"max")[[1]]]))
  # identify the end of tape by rowsum threshold[]
  rr1 = r1$map == clust
  rr1 = rr1 * 1


  ## determine 0cm offset from the middle

  rsums = rowSums(terra::as.array(rr1),na.rm = TRUE)
  # bin into two parts & take the middle one (my images have a complete tape part and a partial tape part.
  # Partial tape indicates the tube surface. Here, we assume that the partial tape is well placed !! (?)
  # needs manual calibration -- work in progress)
  bin = dplyr::ntile(rsums,2)
  zero.rotation.center = stats::median(which(bin == 2 ))
  return(zero.rotation.center)

}


### Rotation Censoring

# Session = sampling campaign with scans taken in short succession
# the CI-600 doesn't have a full 360 degree rotation and thus sequential images have an region with no overlap
# these regions needs to be censored (next implementation)
# keep in mind that inner and outer tube diameter differ and a structural underestimation of rootlength needs a resize coefficient!!!

## input rotation matters for th conclusion!

#' Detect rotation shift between two images
#'
#' Calculates the rotation shift between two sequential images using either
#' cross-correlation or phase correlation methods.
#'
#' @param img1 Reference image (3-channel RGB)
#' @param img2 Subsequent image to compare (3-channel RGB)
#' @param cor.type Correlation type: "ccf" (cross) or "phase" (frequency domain)
#' @param fixed.depth.pixel Depth range to analyze c(start, end)
#' @param fixed.width Width of analysis region in pixels
#' @param select.layer Integer. Specifies which layer to use if the input is a multi-band image. Default is `NULL`.
#' @return Vector of shifts (x,y) in pixels
#' @export
#'
#' @examples
#' img1 = seg_Oulanka2023_Session01_T067
#' img2 = seg_Oulanka2023_Session03_T067
#' y.lag = RotShiftDet(img1,img2,"phase")
RotShiftDet = function(img1,img2,cor.type = "phase",fixed.depth.pixel  = c(1000,4000),fixed.width  =NULL,select.layer = NULL){

  im1 <- load_flexible_image(img1,select.layer = select.layer, output_format = "array", normalize = FALSE  )
  im2 <- load_flexible_image(img2, select.layer = select.layer, output_format = "array", normalize = FALSE  )

  if(is.null(fixed.width)){
    fixed.width = min(dim(img1)[1],dim(img2)[2])
  }

  # grayscale conversion
  img11 <- im1[,,1]*0.21 + im1[,,2]*0.72 + im1[,,3] * 0.07
  img22 <- im2[,,1]*0.21 + im2[,,2]*0.72 + im2[,,3] * 0.07



  ## resize gray image dimensions between sessions; just to be sure
  dif.dim1 = nrow(img11) - nrow(img22)
  dif.dim2 = ncol(img11) - ncol(img22)
  if(dif.dim1 != 0 | dif.dim2 != 0 ){
    print(paste0("Subsequent images differ in size. If the difference is large, reconsider if this operation is valid! Difference is ",dif.dim1," px & ",dif.dim2," px."))
    #img22=resize(img22,output.dim = dim(img11)[1:2],w=dim(img11)[1],h= dim(img11)[2]) # Used to be EBImage function

    ## fix rotation dimension
    rot.dim = round((nrow(img11)/2)-(fixed.width/2)) : round((nrow(img11)/2)+(fixed.width/2)-1)
    img11 = img11[rot.dim,  ]    #RotCensor(img11, fixed.width  = fixed.width)
    img22 = img22[rot.dim,  ]
    ## fix length dimension
    img11 = img11[,fixed.depth.pixel]
    img22 = img22[,fixed.depth.pixel]

    #img22=OpenImageR::resizeImage(img22,width=dim(img11)[1],height= dim(img11)[2], method = "bilinear")
  }

  img11 = terra::as.array(terra::rast(img11))
  img22 = terra::as.array(terra::rast(img22))

  ########### shift detection part on grayscale images  ###
  if(cor.type == "phase"){
    ## phase correlation which should be faster than cross correlation
    a=imagefx::pcorr3d(img11[,,1],img22[,,1])
  }else{
    if(cor.type == "ccf"){
      ## phase correlation which should be faster than cross correlation
      a=imagefx::xcorr3d(img11[,,1],img22[,,1])
    }
  }
  y.lag = c(a$max.shifts[1])# ["row"]
  x.lag = c(a$max.shifts[2])# ["col"]


  if(x.lag > 10){
    print(paste0("Careful, also depth shift present!\n About ",x.lag," pixel." ))
  }

  return(c(x.lag,y.lag))
}



### Cuts images along rotational center
#' Censor image edges based on rotation
#'
#' Crops image edges to handle non-overlapping regions between sequential scans.
#'
#' @param img Input image to censor
#' @param center.offset Rotation shift in rows (from RotShiftDet)
#' @param cut.buffer Proportion of image to cut when fixed_rotation=FALSE
#' @param fixed.rotation Use fixed output dimensions
#' @param fixed.width Output width when fixed_rotation=TRUE
#' @param select.layer Integer. Specifies which layer to use if the input is a multi-band image. Default is `NULL`.
#' @return Cropped raster image
#' @export
#'
#' @examples
#' data(seg_Oulanka2023_Session01_T067)
#' img = terra::rast(seg_Oulanka2023_Session01_T067)
#' censored.raster = RotCensor(img,center.offset = 120, cut.buffer = 0.02)
RotCensor = function(img, center.offset=0,  cut.buffer = 0.02, fixed.rotation = TRUE,
                     fixed.width  =1000, select.layer = NULL ){

  ###### uses the detected shift to whiten the region that doesnt appear in another image (y-dimension) or slide them up and down the tube (x-dimension)
  img.c <- load_flexible_image(img, select.layer = select.layer, output_format = "spatrast", normalize = FALSE  )

  offset.ratio = (abs(center.offset))/  dim(img.c)[1]
  cut.ratio = offset.ratio + cut.buffer

if(fixed.rotation == FALSE){
  if(center.offset >0){
    # positive lag (upper rows cut)
    ex = terra::ext(img.c)
    ex[4] = (1-cut.ratio) * ex[4]
    img.cc = terra::crop(img.c, ex)
  }

  if(center.offset <0){
    # negative lag (lower rows cut)
    ex = terra::ext(img.c)
    ex[3] = cut.ratio * ex[3]
    img.cc = terra::crop(img.c, ex)
  }

  if(center.offset == 0){
    ex = terra::ext(img.c)
    ex[3] = (cut.ratio/2 ) * ex[3]
    ex[4] = (1-cut.ratio/2) * ex[4]
    img.cc = terra::crop(img.c, ex)
  }
}

if(fixed.rotation == TRUE){
  center.row = center.offset

  ex = terra::ext(img.c)
  ex[2] = dim(img.c)[2]
  ex[4] = dim(img.c)[1]
  terra::ext(img.c) <- ex
  ex[3] = ((dim(img.c)[1]/2)+center.row) - fixed.width/2
  ex[4] = ((dim(img.c)[1]/2)+center.row) + fixed.width/2
  img.cc = terra::crop(img.c, ex)
}

  return(img.cc)
}




#' Estimate soil surface position using tape markers
#'
#' Detects the soil surface by analyzing tape coverage patterns in the image.
#'
#' @param img Input image (raster, filename, or array)
#' @param search.area Proportion of image to analyze
#' @param tape.tresh Minimum tape coverage ratio
#' @param dpi Image resolution
#' @param nclasses Number of clustering classes
#' @param inverse Invert detection for dark markers
#' @param tape.overlap Safety margin for tape (cm)
#' @param tape.brightness Brightness threshold for tape
#' @param extra.rows Additional analysis rows
#' @param select.layer Integer. Specifies which layer to use if the input is a multi-band image. Default is `NULL`.
#' @param tape.quantile Brightness alignment quantile
#' @return data.frame with soil surface and tape end positions
#' @export
#'
#' @examples
#' img = seg_Oulanka2023_Session01_T067
#' Soil0Estimates = SoilSurfE(img)
SoilSurfE = function(img,search.area = 0.45, tape.tresh = 0.33,dpi = 300, nclasses = 3, inverse = FALSE,
                     tape.overlap = 0.5,tape.brightness = 0.6,extra.rows = 100,tape.quantile = 0.98, select.layer = NULL ){

  im <- load_flexible_image(img, select.layer = select.layer, output_format = "array", normalize = FALSE  )

  if(inverse == TRUE){
    tape.quantile = 1- tape.quantile
    tape.brightness = 1 / tape.brightness
  }


  ## add one row of red tape pixel
  red.line = array(dim = c(dim(im)[1],extra.rows,dim(im)[3]))
  red.line[,,1:dim(im)[3]] <- stats::quantile(im[,,1],tape.quantile,na.rm=TRUE)
  img1 = abind2(red.line,im,along = 2)

  r.img1 = terra::rast(img1)
  ### make the search are smaller
  r.img1 = terra::crop(r.img1,c(0,search.area*terra::ext(r.img1)[2],0,terra::ext(r.img1)[4]))

  # identify distinct pixel groups
  r1 = RStoolbox::unsuperClass(r.img1,nClasses = nclasses)
  # determine average group
  clust.center = apply(r1$model$centers,1,mean)

  if(inverse == TRUE){
    clust= which(clust.center ==min(clust.center[clust.center < tape.brightness*terra::global(r.img1,"min")[[1]][1:3]]))
  }else{
    # silver tape should have highest luminance across clusters -> select max lum.cluster (but not close to == 1 [pure white?])
    clust= which(clust.center ==max(clust.center[clust.center > tape.brightness*terra::global(r.img1,"max")[[1]][1:3]]))
  }

  # identify the end of tape by rowsum threshold[]
  rr1 = r1$map == clust
  rr1 = rr1*1



  # iterate over depth and check proportion covered by tape
  # checks also 0.1 cm (+12 rows) and 0.2 cm (+24rows) and 0.3 cm (+36) down if the tape reappears  (in case a row impurities or other reasons for missclassification)
  i=1
  while (
    sum(rr1[,i,])/ncol(rr1) > tape.tresh |
    sum(rr1[,i+24,])/ncol(rr1) > tape.tresh |
    sum(rr1[,i+48,])/ncol(rr1) > tape.tresh) {
    i=i+1
    # when loop exceeds image limits, we assume there is no tape = no offset
    if(i > dim(rr1)[2]){
      i = 1
    }
  }

  # store row index to determine soil 0cm offset
  rw.ind = i
  # remove a fixed amount of offset which corresponds to the tape that was installed deeper than 0 cm to prevent light intrusion
  rw.ind = rw.ind - extra.rows - round(tape.overlap*dpi/2.54) # cm tape in the ground converted as pixel offset from soil 0
  rw.ind = round(rw.ind)
  #start.soil = rw.ind # in case
  tape.end = rw.ind + round(tape.overlap*dpi/2.54)
  tape.end = round(tape.end)

  out = data.frame(soil0 = rw.ind, tape.end = tape.end  )
  return(out)
}






