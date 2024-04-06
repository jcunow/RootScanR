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
  max.NAs = round((1-(1/(max(values(bm),na.rm=T)/nn)))+0.02,4)
  # match dimensions
  ex.bm = c(0,dim(binned.map)[2],0,dim(binned.map)[1])
  ex.rp = c(0,dim(rootpic)[2],0,dim(rootpic)[1])
  extent(binned.map) <- ex.bm
  extent(rootpic) <- ex.rp
  r = rootpic
  if(extent(rootpic)[2] != extent(binned.map)[2] | extent(rootpic)[4] != extent(binned.map)[4]){
    r = crop(r,ex.bm)
  }
  # set uninterested zones NA
  cutout.map = binned.map
  r[binned.map != indexD | is.na(binned.map)] = NA
  r = terra::rast(r)
  # non NA values need to cover at least ... % of average depth slice pixel number
  coverNA = round(sum(is.na(values(r[[1]]))) / (ncell(r)),2)
  if( coverNA  < max.NAs){
    r = r %>% terra::trim()
  }else{
    values(r) = NA # set the remaining values NA ?
    if(silent == F){
      print( paste0("Depth: ",indexD,"cm. Not enough informative pixels.",
                    " In the whole Image, ", coverNA*100,"% are NAs after cutting this Depth Slice."," Expected NA% is: ~",(max.NAs-0.02)*100 ))
    }

  }
  return(r)
}
