## Image Coloration


#' Color Coordinates and Brightness of an image segment
#'
#' @param img raster with 3 layers, first is red, second is green, third is blue
#' @param r grayscale factor red
#' @param g grayscale factor green
#' @param b grayscale factor blue
#'
#' @return vector with red chromatic coordinate, green chromatic coordinate, blue chromatic coordinate, mean brightness, mean luminosity
#' @export
#'
#' @examples c(rcc, gcc, bcc, mean.bright, mean.lum)  = Tube.coloration(img)
Tube.coloration = function(img,r=0.2126,g=0.7152,b=0.0722){
  vr = values(img[[1]])
  vg = values(img[[2]])
  vb = values(img[[3]])
  bright = vr + vg + vb
  lum.gray =  vr * r + vg * g + vb *b
  mean.bright = mean(bright,na.rm=T) %>% round(4)
  mean.lum = mean(lum.gray,na.rm=T) %>% round(4)
  rcc = mean(vr / bright,na.rm = T) %>% round(4)
  gcc = mean(vg / bright,na.rm = T) %>% round(4)
  bcc = mean(vb / bright,na.rm = T) %>% round(4)
  return(c(rcc,gcc,bcc,mean.bright,mean.lum))
}
