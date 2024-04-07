## Rotation Estimate


### read image
#img1 <- readTIFF(paste(path,Session,file.name,sep=''))

#' Estimates Rotation from Tape Coverage
#'
#' @param im image.tiff
#' @param tape.brightness used for clustering. Tape appears bright e.g., 0.66
#' @param extra.rows In case no tape is present. Best leave unchanged - some extra.rows are recommended and will be substracted from the output anyway.
#'
#' @return numeric value corresponding to the center of extruding tape
#' @export
#'
#' @examples value = RotationE(img)
RotationE = function(im,tape.brightness = 0.66,extra.rows = 100){

  if(class(im) != "array"){
    im = as.array(im)
  }

  ## add one row of red tape pixel
  red.line = array(dim = c(dim(im)[1],extra.rows,dim(im)[3]))
  red.line[,,1:dim(im)[3]] <- max(im[,,1])
  img1 = abind::abind(red.line,im,along = 2)


  r.img1 = brick(img1)
  ### make the search are smaller
  r.img1 = crop(r.img1,extent(0,search.area,0,1))

  # identify distinct pixel groups
  r1 = RStoolbox::unsuperClass(r.img1,nClasses = 3)
  # determine average group
  clust.center = apply(r1$model$centers,1,mean)
  # silver tape should have highest luminance across clusters -> select max lum.cluster (but not close to == 1 [pure white?])
  clust= which(clust.center ==max(clust.center[clust.center > tape.brightness*max(values(r.img1))]))
  # identify the end of tape by rowsum threshold[]
  rr1 = r1$map == clust


  ## determine 0cm offset from the middle
  rsums = rowSums(rr1)
  # bin into two parts & take the middle one (my images have a complete tape part and a partial tape part.
  # Partial tape indicates the tube surface. Here, we assume that the partial tape is well placed !! (?)
  # needs manual calibration -- work in progress)
  bin = dplyr::ntile(rsums,2)
  zero.rotation.center = median(which(bin == 2 ))
  return(zero.rotation.center)

}


### Rotation Censoring

# Session = sampling campaign with scans taken in short succession
# the CI-600 doesn't have a full 360 degree rotation and thus sequential images have an region with no overlap
# these regions needs to be censored (next implementation)
# keep in mind that inner and outer tube diameter differ and a structural underestimation of rootlength needs a resize coefficient!!!

## input rotation matters for th conclusion!

#' Rotation Correlation
#'
#' @param img1 Reference Image
#' @param img2 Subsequent Image
#' @param cor.type Two correlation types available: "ccf" cross correlation, and "phase" phase correlation in frequency domain. See ??imagefx
#'
#' @return numeric value corresponding to rotation in rows. Will be nonsens if image is rotated 90 degrees.
#' @export
#'
#' @examples y.lag = RotShiftDet(img,img2,"phase")
RotShiftDet = function(img1,img2,cor.type = "phase"){
  ## import image
  # check image type


  # grayscale conversion
  img11 <- img1[[1]]*0.21 + img1[[2]]*0.72 + img1[[3]] * 0.07
  img11 <- img2[[1]]*0.21 + img2[[2]]*0.72 + img2[[3]] * 0.07



  ## resize gray image dimensions between sessions; just to be sure
  dif.dim = dim(img11)-dim(img22)
  if(dif.dim[1] != 0 | dif.dim[2] != 0 ){
    print(paste0("Subsequent images differ in size. Difference is ",dif.dim," px."))
    img22=EBImage::resize(img22,output.dim = dim(img11)[1:2],w=dim(img11)[1],h= dim(img11)[2])
  }


  ########### shift detection part on grayscale images  ###
  if(cor.type == "phase"){
    ## phase correlation which should be faster than cross correlation
    a=imagefx::pcorr3d(img11,img22)
  }else{
    if(cor.type == "ccf"){
      ## phase correlation which should be faster than cross correlation
      a=imagefx::xcorr3d(img11,img22)
    }
  }
  y.lag = c(a$max.shifts[1])# ["row"]
  x.lag = c(a$max.shifts[2])# ["col"]


  if(x.lag > 10){
    print(paste0("Careful, also depth shift present!\n About ",x.lag," pixel." ))
  }

  return(y.lag)
}



### Cuts images along rotational center
#' Rotation Censor
#'
#' @param img The image which should be censored
#' @param center.offset Rotational shift in rows. Can be retrieved from 'RotShiftDet()', reference image is 0 offset
#' @param cut.buffer fixed ratio that will be cut in addition
#'
#' @return a raster
#' @export
#'
#' @examples censored.raster = RotCensor(img,center.offset = 120, cut.buffer = 0.02)
RotCensor = function(img, center.offset=0,  cut.buffer = 0.02){

  ###### uses the detected shift to whiten the region that doesnt appear in another image (y-dimension) or slide them up and down the tube (x-dimension)
  img.c=raster::brick(img)
  offset.ratio = (abs(center.offset))/  dim(img)[1]
  cut.ratio = offset.ratio + cut.buffer


  if(center.offset >0){
    # positive lag (upper rows cut)
    ex = raster::extent(img.c)
    ex[4] = 1-cut.ratio
    img.cc = raster::crop(img.c, ex)
  }

  if(center.offset <0){
    # negative lag (lower rows cut)
    ex = raster::extent(img.c)
    ex[3] = cut.ratio
    img.cc = raster::crop(img.c, ex)
  }

  if(center.offset == 0){
    ex = raster::extent(img.c)
    ex[3] = cut.ratio/2
    ex[4] = 1-cut.ratio/2
    img.cc = raster::crop(img.c, ex)
  }

  return(img.cc)
}








