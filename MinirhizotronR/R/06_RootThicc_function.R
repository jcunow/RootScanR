### Root Thickness


#' Approximate average Root Thickness
#'
#' @param kimuralength length of roots in image section
#' @param rootpx amount of rootpx in the image section
#' @param dpi image resolution
#'
#' @return a value in units cm
#' @export
#'
#' @examples root.ticc = root.tickness(300,9500,300)
root.tickness = function(kimuralength,rootpx,dpi){
  px.per.length = rootpx / kimuralength
  thiccness = (px.per.length/118) / (dpi/2.54)
  return(thiccness) # cm
}
