## helper functions

#' Takes a continues depth map and bins it to a specified range
#'
#' @param depthmap image with depth information as illumination
#' @param nn bandwidth in cm
#'
#' @return depths in bins
#' @export
#'
#' @examples image = binning(depthmap,5)
binning = function(depthmap,nn){
  im = round(round(depthmap*(1/nn),0)*nn)
  #return(im)
}


#' Zoning
#'
#' nullifies all depth except the specified target depth in the root image
#' @param rootpic the segmented root image
#' @param binned.map depthmap with bands
#' @param indexD the target depth
#'
#' @return roots, but only in the specified depth band
#'
#' @examples zone.fun(root.pic, binned.map,IndexD) = root.pic2
zone.fun = function(rootpic,binned.map,indexD){
  nulled.depth.map = binned.map
  nulled.depth.map[nulled.depth.map != indexD] = 0
  r = terra::rast(rootpic) * terra::rast(nulled.depth.map)
  return(r)
}
