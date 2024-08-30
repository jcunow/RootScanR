#' Estimates Rotation from Tape Coverage. Assumes that more tape is present on the tube upside.
#'
#' @param img takes raster, file name, or array input.
#' @param tape.brightness used for clustering. Silver Dukt tape appears bright e.g., > 0.66.
#' @param tape.quantile aligns extra.rows brightness with tape to ensure the same cluster class. The default uses Silver Tape as reference.
#' @param extra.rows In case no tape is present. Best leave unchanged - some extra.rows are recommended and will be subtracted from the output anyway.
#' @param search.area portion of image to perform the tape search on. Potential speed increase.
#' @param nclasses number of classes to group pixels with unsupervised Clustering from the RStoolbox package. Leave at default unless weird results are being displayed.
#' @importFrom "stats" "quantile"
#' @importFrom "stats" "median"
#'
#' @return numeric value corresponding to the center of extruding tape
#' @export
#'
#' @examples
#' img = terra::rast(seg_Oulanka2023_Session01_T067)
#' r0 = RotationE(img)
RotationE = function(img,tape.brightness = 0.66,extra.rows = 100,search.area = 0.45,tape.quantile = 0.98, nclasses = 3){

  if(is.character(img)){
    im = terra::rast(img)
    im = terra::as.array(im)
  }else{
    if(!is.array(img) ){
      im = terra::as.array(img)
    }else{
      if(is.array(img) ){
        im = img
      }else{
        print("not the right input format. Try RasterBrick or array type of format.")
      }

    }
  }


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

#' Rotation Correlation
#'
#' @param img1 Reference Image with 3 rgb channel
#' @param img2 Subsequent Image with 3 rgb channel
#' @param cor.type Two correlation types available: "ccf" cross correlation, and "phase" phase correlation in frequency domain. See ??imagefx
#'
#' @return numeric value corresponding to rotation in rows. Ensure correct image rotation i.e., rows == rotation
#' @export
#'
#' @examples
#' img1 = seg_Oulanka2023_Session01_T067
#' img2 = seg_Oulanka2023_Session03_T067
#' y.lag = RotShiftDet(img1,img2,"phase")
RotShiftDet = function(img1,img2,cor.type = "phase"){
  ## import image
  # check image type for image 1
  if(is.character(img1)){
    im1 = terra::rast(img1)
    im1 = terra::as.array(im1)
  }else{
    if(!is.array(img1) ){
      im1 = terra::as.array(img1)
    }else{
      im1 = img1
    }
  }


  # check image type for image 2
  if(is.character(img2)){
    im2 = terra::rast(img2)
    im2 = terra::as.array(im2)
  }else{
    if(!is.array(img2) ){
      im2 = terra::as.array(img2)
    }else{
      im2 = img2
    }
  }


  # grayscale conversion
  img11 <- im1[,,1]*0.21 + im1[,,2]*0.72 + im1[,,3] * 0.07
  img22 <- im2[,,1]*0.21 + im2[,,2]*0.72 + im2[,,3] * 0.07



  ## resize gray image dimensions between sessions; just to be sure
  dif.dim1 = nrow(img11) - nrow(img22)
  dif.dim2 = ncol(img11) - ncol(img22)
  if(dif.dim1 != 0 | dif.dim2 != 0 ){
    print(paste0("Subsequent images differ in size. If the difference is large, reconsider if this operation is valid! Difference is ",dif.dim1," px & ",dif.dim1," px."))
    #img22=resize(img22,output.dim = dim(img11)[1:2],w=dim(img11)[1],h= dim(img11)[2]) # Used to be EBImage function
    img22=OpenImageR::resizeImage(img22,width=dim(img11)[1],height= dim(img11)[2], method = "bilinear")
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
#' @param cut.buffer ratio of image dimensions that will be cut if fixed.rotation=FALSE
#' @param fixed.rotation specifies whether censoring is applied to fixed output dimensions (=TRUE) or proportional to input dimensions (=FALSE)
#' @param fixed.width if fixed.rotation is TRUE, fixed.width specifies the final amount of rows centered at center.offset
#'
#' @return a raster
#' @export
#'
#' @examples
#' data(seg_Oulanka2023_Session01_T067)
#' img = terra::rast(seg_Oulanka2023_Session01_T067)
#' censored.raster = RotCensor(img,center.offset = 120, cut.buffer = 0.02)
RotCensor = function(img, center.offset=0,  cut.buffer = 0.02, fixed.rotation = TRUE,fixed.width  =1000 ){

  ###### uses the detected shift to whiten the region that doesnt appear in another image (y-dimension) or slide them up and down the tube (x-dimension)
  img.c=img
  offset.ratio = (abs(center.offset))/  dim(img)[1]
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




#' Estimate the position of the soil surface by tape presence
#'
#' @param img raster,filename, or array input
#' @param search.area ratio of image which is used to look for tape cover. Speeds up computation.
#' @param tape.tresh ratio of how much of the tube rotation needs to covered in tape
#' @param dpi image resolution
#' @param nclasses number of clusters to discern tape, roots, background. 3 for rgb images. 2 can be appropriate for segmented 0-1 images.
#' @param tape.overlap assumes a safety margin on the tape. The soil surface will be shifted by this amount in cm
#' @param tape.brightness used for clustering. Tape appears bright e.g., 0.66
#' @param tape.quantile aligns extra.rows brightness with the tape. The default uses Silver Tape as reference.
#' @param extra.rows In case no tape is present. Best leave unchanged - some extra.rows are recommended and will be subtracted from the output anyway.
#' @importFrom "stats" "quantile"
#'
#' @return data.frame with tape end and soil surface estimation in rows
#' @export
#'
#' @examples
#' img = seg_Oulanka2023_Session01_T067
#' Soil0Estimates = SoilSurfE(img)
SoilSurfE = function(img,search.area = 0.45, tape.tresh = 0.33,dpi = 300, nclasses = 3, inverse = FALSE,
                     tape.overlap = 0.5,tape.brightness = 0.75,extra.rows = 100,tape.quantile = 0.98 ){

  if(is.character(img)){
    im = terra::rast(img)
    im = terra::as.array(im)
  }else{
    if(!is.array(img) ){
      im = terra::as.array(img)
    }else{
      im = img
    }
  }


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
    clust= which(clust.center ==min(clust.center[clust.center < tape.brightness*terra::global(r.img1,"min")[[1]]]))
  }else{
    # silver tape should have highest luminance across clusters -> select max lum.cluster (but not close to == 1 [pure white?])
    clust= which(clust.center ==max(clust.center[clust.center > tape.brightness*terra::global(r.img1,"max")[[1]]]))
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






