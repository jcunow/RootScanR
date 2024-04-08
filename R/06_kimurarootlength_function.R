#' Estimate Kimura Root Length
#'
#' @param im a skeletonized root image
#' @param unit output unit
#' @param dpi image resolution
#'
#' @return numeric value of root length
#' @export
#'
#' @examples value = RootLength(im)
RootLength = function(im,unit="cm",dpi=300){
  im2 = im / max(raster::values(im),na.rm=T)
  ## kimura image
  k0 = matrix(c(0,1,0,0,1,0,0,0,0), nrow = 3, ncol = 3)
  k1 = matrix(c(0,0,0,1,1,0,0,0,0), nrow = 3, ncol = 3)
  r0 <- terra::focal(im2,w = k0, fun = "sum")
  r1 <- terra::focal(im2,w = k1, fun = "sum")
  orth.img = sum((r0 == 2) | (r1 == 2))

  k0 = matrix(c(1,0,0,0,1,0,0,0,0), nrow = 3, ncol = 3)
  k1 = matrix(c(0,0,1,0,1,0,0,0,0), nrow = 3, ncol = 3)
  r0 <- terra::focal(im2 ,w = k0, fun = "sum")
  r1 <- terra::focal(im2 ,w = k1, fun = "sum")
  diag.img = sum((r0 == 2) | (r1 == 2))

  kimura.sum.diag <- raster::cellStats(diag.img,"sum")
  kimura.sum.orth <- raster::cellStats(orth.img,"sum")

  if(unit == "px"){
  rootlength = round(( kimura.sum.diag**2 + (kimura.sum.diag + kimura.sum.orth/2)**2 )**0.5   + kimura.sum.orth/2)
  }else{
    if(unit=="cm"){
      rootlength = round((( kimura.sum.diag**2 + (kimura.sum.diag + kimura.sum.orth/2)**2 )**0.5   + kimura.sum.orth/2) / dpi/2.54)
    }else{
      print("only 'px' or 'cm' acceptable as unit")
    }
  }


  return(rootlength)


}
