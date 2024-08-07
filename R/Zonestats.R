




#' Estimate Kimura Root Length
#'
#' @param im a skeletonized root image raster
#' @param unit output unit
#' @param dpi image resolution
#'
#' @return root length
#' @export
#'
#' @examples
#' data(skl_Oulanka2023_Session01_T067)
#' img = terra::rast(skl_Oulanka2023_Session01_T067[[2]])
#' RL = RootLength(im = img,unit = "cm", dpi = 300)
RootLength = function(im,unit="cm",dpi=300){
  im2 = im / terra::global(im,"max",na.rm=TRUE)[[1]]
  ## kimura image
  k0 = matrix(c(0,1,0,0,1,0,0,0,0), nrow = 3, ncol = 3)
  k1 = matrix(c(0,0,0,1,1,0,0,0,0), nrow = 3, ncol = 3)
  r0 <- terra::focal(im2,w = k0, fun = "sum",na.rm=TRUE)
  r1 <- terra::focal(im2,w = k1, fun = "sum",na.rm=TRUE)
  orth.img = sum((r0 == 2) | (r1 == 2))

  g0 = matrix(c(1,0,0,0,1,0,0,0,0), nrow = 3, ncol = 3)
  g1 = matrix(c(0,0,1,0,1,0,0,0,0), nrow = 3, ncol = 3)
  u0 <- terra::focal(im2 ,w = g0, fun = "sum",na.rm=TRUE)
  u1 <- terra::focal(im2 ,w = g1, fun = "sum",na.rm=TRUE)
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

  return(rootlength)

}


#' Estimate the directional of a skeleton pixels
#'
#' @param im skeletonized image raster. Roots must be 1, background 0. Consider the input image rotation to interpret the output
#' @param rotate should the image be rotated 90 degrees before entering?
#' @param kernelsize number of connecting cells
#' @param leeway amount of gaps allowed in connecting lines
#'
#' @return percentage of pixels with given neighbor pixel position
#' @export
#'
#' @examples
#' library(terra);library(raster)
#'
#' data("skl_Oulanka2023_Session01_T067")
#' img = terra::rast(skl_Oulanka2023_Session01_T067)
#' direction.frame = Directionality(img,kernelsize = 3)
Directionality = function(im,rotate = TRUE, kernelsize = 3, leeway =0){
  if(terra::global(im[[2]],"max",na.rm = TRUE)[[1]] > 1){
    im2 = im / 255
  }else{
    im2 = im
  }


  # the rotation of the image matters !
  if(rotate == TRUE){
    im2 = terra::t(im2)
  }

  # background white or objects white matters - we want to count white objects (1's not 0's)

  # Kernel
  kerns = adaptive.Kernelsize.Directionality(n = kernelsize)
  D_horizontal = kerns[[2]]
  D_vertical = kerns[[1]]
  D_dia_topleft = kerns[[3]]
  D_dia_botleft = kerns[[4]]

  # orthogonal
  r_Dho <- terra::focal(im2,w = D_horizontal, fun = "sum",na.rm = TRUE)
  r_Dve <- terra::focal(im2,w = D_vertical, fun = "sum",na.rm = TRUE)
  rr_Dho = sum(r_Dho >= (kernelsize-leeway))
  rr_Dve = sum(r_Dve >= (kernelsize-leeway))
  #sum.Dho <- raster::cellStats(rr_Dho,"sum")
  #sum.Dve <- raster::cellStats(rr_Dve,"sum")
  sum.Dho <- terra::global(rr_Dho,"sum",na.rm = TRUE)[[1]]
  sum.Dve <- terra::global(rr_Dve,"sum",na.rm = TRUE)[[1]]


  # diagonal
  r_Dtopleft <- terra::focal(im2,w = D_dia_topleft, fun = "sum",na.rm = TRUE)
  r_Dbotleft <- terra::focal(im2,w = D_dia_botleft, fun = "sum",na.rm = TRUE)
  rr_Dtopleft = sum(r_Dtopleft >= (kernelsize-leeway))
  rr_Dbotleft = sum(r_Dbotleft >= (kernelsize-leeway))
  # sum.Dtl <- raster::cellStats(rr_Dtopleft,"sum")
  # sum.Dbl <- raster::cellStats(rr_Dbotleft,"sum")
  sum.Dtl <- terra::global(rr_Dtopleft,"sum",na.rm = TRUE)[[1]]
  sum.Dbl <- terra::global(rr_Dbotleft,"sum",na.rm = TRUE)[[1]]

# pixel neigbourhood sums
  all.px = sum.Dho + sum.Dve + sum.Dtl + sum.Dbl



Directions = data.frame(
                        horizontal = round(sum.Dho / all.px, 4),
                        vertical = round(sum.Dve / all.px, 4),
                        topleft.bottomright = round(sum.Dtl / all.px, 4),
                        topright.bottom.left = round(sum.Dbl / all.px, 4),
                        horizontalpx = sum.Dho,
                        verticalpx = sum.Dve,
                        topleftpx = sum.Dtl,
                        botleftpx = sum.Dbl,
                        avg.down.angle = circular_mean(angles = c(rep(225,sum.Dbl),rep(135,sum.Dtl),rep(180,sum.Dve)),
                                                  input_units = "degrees",output_units = "degrees"))
colnames(Directions) <- c("a90_270%","a180%","a135%","a225%","a90_270px","a180px","a135px","a225px","avg.down_angle")
  return(Directions)

}



## RootScape Metrics

# input image should be segmented raster

#' RootScapeMetric relies on Landscapemetrics to extract 'Root Scape' Features akin to landscape analysis.
#'
#' @param im segmented raster  (values = 0,1). Consider whether skeletonized raster is appropriate.
#' @param indexD please specify depth. Will only affect the output column = "depth". Useful when used in a loop.
#' @param metrics which ,metrics should be calculated from the available ones in 'landscapemetrics::calculate_lsm()'.
#' @import dplyr
#' @return a bunch of metric values
#' @export
#'
#' @examples
#' data(seg_Oulanka2023_Session01_T067)
#' img = terra::rast(seg_Oulanka2023_Session01_T067)[[2]]
#' RootScapeObject  = RootScapeMetrics(img,indexD = 80, metrics = c("lsm_c_ca"))
RootScapeMetrics = function(im,indexD=NA, metrics = c( "lsm_c_ca","lsm_l_ent","lsm_c_pd","lsm_c_np","lsm_c_pland",
                                                    "lsm_c_area_mn","lsm_c_area_cv","lsm_c_enn_mn","lsm_c_enn_cv")){

  rsm = landscapemetrics::calculate_lsm(im, directions = 8, neighbourhood = 8,what = metrics)
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
#' @param root.zone one layer raster
#'
#' @return a numeric value
#' @export
#'
#' @examples
#' data(seg_Oulanka2023_Session01_T067)
#' img = terra::rast(seg_Oulanka2023_Session01_T067[[2]])
#' rootpixel  = px.sum(img)
px.sum = function(root.zone){
  rootpx = as.numeric(terra::values(root.zone)/ terra::global(root.zone,"max",na.rm = TRUE)[[1]] ) ==1
  rootpx = rootpx * 1
  srpx = sum(rootpx,na.rm=T)
  return(srpx)
}




#' Coloration of the image
#'
#' @param img raster with 3 color bands
#' @param r weight for first channel - typically red
#' @param g weight for second channel - typically green
#' @param b weight for third channel - typically blue
#' @importFrom "grDevices" "rgb2hsv"
#' @import dplyr
#'
#' @return a vector with chromatic coordinates,luminosity, brightness, luminosity, color values, saturation
#' @export
#'
#' @examples
#' data(rgb_Oulanka2023_Session03_T067)
#' img = terra::rast(rgb_Oulanka2023_Session03_T067)
#' colorvector = Tube.coloration(img = img)
Tube.coloration = function(img,r=0.2126,g=0.7152,b=0.0722){
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



#### Image Texture

# uses the glcm package
#' Texture corresponds to glmc::glmc()
#'
#' @param img.color image with three color channels. Will be converted in gray scale. raster format required.
#' @param grays number of gray shades. Documentation is lacking, see: ?glmc::glmc()
#' @param window convolution window size e.g., c(3,3)
#' @param metrics texture metrics based on illumination differences. see:: ?glmc::glmc() for available methods
#'
#' @return raster with each layer correspond to on metric
#' @export
#'
#' @examples
#' data(rgb_Oulanka2023_Session03_T067)
#' img = rgb_Oulanka2023_Session03_T067
#' texture(img, 7, c(9,9), metrics = "second_moment")
texture = function(img.color,grays = 7, window = c(9,9), metrics = c("variance","second_moment") ){

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


### Root Thickness


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
