



#' counts all pixels in a binarized root image
#'
#' @param root.zone a one layer raster layer
#'
#' @return a numeric value
#' @export
#'
#' @examples value  = px.sum(root.zone)
px.sum = function(root.zone){
  rootpx = as.numeric(raster::values(root.zone)/ max(raster::values(root.zone),na.rm=T) ) ==1
  srpx = sum(rootpx,na.rm=T)
  return(srpx)
}
