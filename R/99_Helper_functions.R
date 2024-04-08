## helper functions

#' Takes a continues depth map and bins it to a specified range
#'
#' @param depthmap image with depth information as illumination
#' @param nn bandwidth in cm
#'
#' @return depths in bins
#' @export
#'
#' @examples image = binning(depthmap,nn = 5)
binning = function(depthmap,nn){
  im = round(round(depthmap*(1/nn),0)*nn)
  #return(im)
}



#' Cuts out a zone of interest
#'
#' @param rootpic the image which should be cut
#' @param binned.map image which supplies the cut condition
#' @param indexD the condition
#' @param nn bin width in binned.map. Should be the same as used in 'binning()'
#' @param silent verbose
#'
#' @return a cutout of the root pic image
#' @export
#'
#' @examples raster_image = zone.fun(rootpic, binned.map, nn=5, silent = F)
zone.fun = function(rootpic,binned.map,indexD,nn = 5,silent =F){
  max.NAs = round((1-(1/(max(raster::values(binned.map),na.rm=T)/nn)))+0.02,4)
  # match dimensions
  ex.bm = c(0,dim(binned.map)[1],0,dim(binned.map)[2])
  ex.rp = c(0,dim(rootpic)[1],0,dim(rootpic)[2])
  if(any(ex.bm!=ex.rp)){
    rootpic = raster::t(rootpic)
  }
  #raster::extent(binned.map) <- ex.bm
  #raster::extent(rootpic) <- ex.rp

  r = (rootpic)
  if(raster::extent(rootpic)[2] != raster::extent(binned.map)[2] | raster::extent(rootpic)[4] != raster::extent(binned.map)[4]){
    r = raster::crop(r,ex.bm)
  }
  # set uninterested zones NA
  raster::values(r)[raster::values(binned.map) != indexD ] <- NA
  raster::values(r)[ is.na(raster::values(binned.map))] <- NA
  r = terra::rast(r)
  # non NA values need to cover at least ... % of average depth slice pixel number
  coverNA = round(sum(is.na(raster::values(r[[1]]))) / (raster::ncell(r)),2)
  if( coverNA  < max.NAs){
    r = r %>% terra::trim()
  }else{
    raster::values(r) = NA # set the remaining values NA ?
    if(silent == F){
      print( paste0("Depth: ",indexD,"cm. Not enough informative pixels.",
                    " In the whole Image, ", coverNA*100,"% are NAs after cutting this Depth Slice."," Expected NA% is: ~",(max.NAs-0.02)*100 ))
    }

  }
  r = raster::raster(r)
  return(r)
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
