#' Calculate Root Length using Kimura's Method with optimizations
#'
#' @param img A skeletonized root image raster
#' @param unit Output unit ("px" or "cm")
#' @param dpi Image resolution (required when unit = "cm")
#' @param select.layer Integer. Specifies which layer to use if the input is a multi-band image. Default is `2`.
#' @return Numeric value representing root length in specified unit
#' @export
#'
#' @examples
#' data(skl_Oulanka2023_Session01_T067)
#' img = terra::rast(skl_Oulanka2023_Session01_T067)
#' RL = RootLength(img = img,unit = "cm", dpi = 300)
RootLength = function(img,unit="cm",dpi=300,select.layer = NULL){

  # flexible input
  img <- load_flexible_image(img, select.layer = select.layer, output_format = "spatrast", normalize = TRUE  )


  ## kimura image
  k0 = matrix(c(0,1,0,0,1,0,0,0,0), nrow = 3, ncol = 3)
  k1 = matrix(c(0,0,0,1,1,0,0,0,0), nrow = 3, ncol = 3)
  r0 <- terra::focal(img,w = k0, fun = "sum",na.rm=TRUE)
  r1 <- terra::focal(img,w = k1, fun = "sum",na.rm=TRUE)
  orth.img = sum((r0 == 2) | (r1 == 2))

  g0 = matrix(c(1,0,0,0,1,0,0,0,0), nrow = 3, ncol = 3)
  g1 = matrix(c(0,0,1,0,1,0,0,0,0), nrow = 3, ncol = 3)
  u0 <- terra::focal(img ,w = g0, fun = "sum",na.rm=TRUE)
  u1 <- terra::focal(img ,w = g1, fun = "sum",na.rm=TRUE)
  diag.img = sum((u0 == 2) | (u1 == 2))

  kimura.sum.diag <- terra::global(diag.img,"sum",na.rm=TRUE)
  kimura.sum.orth <- terra::global(orth.img,"sum",na.rm=TRUE)

  if(unit == "px"){
    rootlength = round(( kimura.sum.diag**2 + (kimura.sum.diag + kimura.sum.orth/2)**2 )**0.5   + kimura.sum.orth/2)
  }else{
    if(unit=="cm"){
      rootlength = round((( kimura.sum.diag**2 + (kimura.sum.diag + kimura.sum.orth/2)**2 )**0.5   + kimura.sum.orth/2) / dpi/2.54,3)
    }else{
      print("only 'px' or 'cm' acceptable as unit")
    }
  }

  return(rootlength[[1]])

}




## RootScape Metrics

# input image should be segmented raster

#' RootScapeMetric relies on Landscapemetrics to extract 'Root Scape' Features akin to landscape analysis.
#'
#' @param img segmented raster  (values = 0,1). Consider whether skeletonized raster is appropriate.
#' @param indexD please specify depth. Will only affect the output column = "depth". Useful when used in a loop.
#' @param metrics which ,metrics should be calculated from the available ones in 'landscapemetrics::calculate_lsm()'.
#' @param select.layer Integer. Specifies which layer to use if the input is a multi-band image. Default is `2`.
#' @import dplyr
#' @return a bunch of metric values
#' @export
#'
#' @examples
#' data(seg_Oulanka2023_Session01_T067)
#' img = terra::rast(seg_Oulanka2023_Session01_T067)[[2]]
#' RootScapeObject  = RootScapeMetrics(img,indexD = 80, metrics = c("lsm_c_ca"))
RootScapeMetrics = function(img,indexD=NA, select.layer = NULL, metrics = c( "lsm_c_ca","lsm_l_ent","lsm_c_pd","lsm_c_np","lsm_c_pland",
                                                    "lsm_c_area_mn","lsm_c_area_cv","lsm_c_enn_mn","lsm_c_enn_cv")){

  # flexible input
  img <- load_flexible_image(img, select.layer = select.layer, output_format = "spatrast", normalize = FALSE  )


  rsm = landscapemetrics::calculate_lsm(img, directions = 8, neighbourhood = 8,what = metrics)
  t.object = ifelse(rsm$class == 0,"deletable","root")
  t.object = ifelse(is.na(rsm$class),"root",t.object )
  rsm$object = t.object
  rsm$depth = indexD
  rsm = dplyr::distinct(rsm)
  rsm = dplyr::filter(rsm, rsm$object != "deletable")
  rsm$id = NULL
  rsm$class = NULL
  rsm$level = NULL
  rsm$layer = NULL

  return(rsm)
}






#' counts all pixels in a segmented image
#'
#' @param img one layer image
#'
#' @return a numeric value
#' @export
#'
#' @examples
#' data(seg_Oulanka2023_Session01_T067)
#' img = terra::rast(seg_Oulanka2023_Session01_T067)[[2]]
#' rootpixel  = px.sum(img)
px.sum = function(img){
  # flexible input
  img <- load_flexible_image(img, output_format = "spatrast", normalize = TRUE  )

  srpx = terra::global(img,"sum",na.rm = TRUE)[[1]]
  return(srpx)
}




#' Calculate Image Coloration Metrics
#'
#' @param img Three-band raster (RGB) or path to image
#' @param r Red channel weight
#' @param g Green channel weight
#' @param b Blue channel weight
#' @return Data frame of color metrics
#' @export
#'
#' @examples
#' data(rgb_Oulanka2023_Session03_T067)
#' img = terra::rast(rgb_Oulanka2023_Session03_T067)
#' colorvector = Tube.coloration(img)
Tube.coloration = function(img,r=0.2126,g=0.7152,b=0.0722){
  # flexible input
  img <- load_flexible_image(img, output_format = "spatrast", normalize = FALSE  )


  vr = terra::values(img[[1]])
  vg = terra::values(img[[2]])
  vb = terra::values(img[[3]])
  mean.r = mean(vr,na.rm=T )
  mean.g = mean(vg,na.rm=T )
  mean.b = mean(vb,na.rm=T )

  hsl = grDevices::rgb2hsv(r = mean.r, g = mean.g, b = mean.b)
  intensity = vr + vg + vb
  lum.gray =  vr * r + vg * g + vb *b
  mean.intensity = round(mean(intensity,na.rm=T),4)
  mean.lum = round(mean(lum.gray,na.rm=T) ,4)
  rcc = round(mean(vr / intensity,na.rm = T) ,4)
  gcc = round(mean(vg / intensity,na.rm = T) ,4)
  bcc = round(mean(vb / intensity,na.rm = T),4)
  colordf = data.frame(rcc = rcc,gcc = gcc,bcc = bcc,
                       hue = hsl[1], saturation = hsl[2], luminosity = hsl[3],
                       red = mean.r,green = mean.g,blue = mean.b)
  return(colordf)
}



#' Enhanced texture calculation
#'
#' @param img.color Three-band raster or path to image
#' @param grays Number of gray levels
#' @param window Window size for GLCM
#' @param metrics Texture metrics to calculate
#' @return Raster with texture metrics
#' @export
#'
#' @import raster
#'
#' @examples
#' data(rgb_Oulanka2023_Session03_T067)
#' img = raster::brick(rgb_Oulanka2023_Session03_T067)
#' texture(img, 7, c(9,9), metrics = "second_moment")
texture = function(img.color,grays = 7, window = c(9,9), metrics = c("variance","second_moment") ){

  # flexible input
  img <- load_flexible_image(img, output_format = "spatrast", normalize = FALSE  )

  mx = max(raster::values(img.color),na.rm=TRUE)
  if(mx > 1){
    mx = 255
  }else{
    mx = 1
  }
  img.gray = (img.color[[1]]*0.21 + img.color[[2]]*0.72 + img.color[[3]]*0.07) / 255
  tx.im  =  glcm::glcm(img.gray,
                       n_grey = grays,
                       statistics = metrics)
}



#' Approximate average Root Thickness
#'
#' @param kimuralength length of roots in image section, input unit is cm.
#' @param rootpx amount of rootpx in the image section
#' @param dpi image resolution
#'
#' @return a value in units cm
#' @export
#'
#' @examples root.ticc = root.thickness(kimuralength = 300,rootpx = 9500, dpi = 300)
root.thickness = function(kimuralength,rootpx,dpi=300){
  px.length =  kimuralength * (dpi/2.54) # from cm_RL -> px_RL
  thiccness = rootpx / px.length
  px.thicc = thiccness / (dpi/2.54) # from px_thiccnes/px_RL -> cm_thickness/px_RL

  return(px.thicc) # output is root thickness in cm per smallest unit i.e., one pixel
}
