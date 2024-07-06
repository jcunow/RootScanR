#' Estimate Kimura Root Length
#'
#' @param im a skeletonized root image
#' @param unit output unit
#' @param dpi image resolution
#'
#' @return root length
#' @export
#'
#' @examples RL = RootLength(im = image,unit = "cm", dpi = 300)
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
      rootlength = round((( kimura.sum.diag**2 + (kimura.sum.diag + kimura.sum.orth/2)**2 )**0.5   + kimura.sum.orth/2) / dpi/2.54,3)
    }else{
      print("only 'px' or 'cm' acceptable as unit")
    }
  }


  return(rootlength)

}


#' Estimate the directional of a skeleton pixels
#'
#' @param im a skeletonized image. Roots must be 1, background 0. The rotation of the image determines the output. The default is
#'
#' @return percentage of pixels with given neighbour pixel position
#' @export
#'
#' @examples direction.frame = Directionality(im)
Directionality = function(im,rotate = TRUE){
  im2 = im / max(raster::values(im),na.rm=T)
  # the rotation of the image matters !
  if(rotate == T){
    im2 = raster::t(im2)
  }
  # background white or objects white matters - we want to count white objects (1's not 0's)

  ## kimura image
  D_horizontal = matrix(c(0,1,0,0,1,0,0,0,0), nrow = 3, ncol = 3)
  D_vertical = matrix(c(0,0,0,1,1,0,0,0,0), nrow = 3, ncol = 3)
  D_dia_topleft = matrix(c(1,0,0,0,1,0,0,0,0), nrow = 3, ncol = 3)
  D_dia_botleft = matrix(c(0,0,1,0,1,0,0,0,0), nrow = 3, ncol = 3)

  # orthogonal
  r_Dho <- terra::focal(im2,w = D_horizontal, fun = "sum")
  r_Dve <- terra::focal(im2,w = D_vertical, fun = "sum")

  rr_Dho = sum((r_Dho == 2))
  rr_Dve = sum((r_Dve == 2))

  sum.Dho <- raster::cellStats(rr_Dho,"sum")
  sum.Dve <- raster::cellStats(rr_Dve,"sum")
  # diagonal
  r_Dtopleft <- terra::focal(im2,w = D_dia_topleft, fun = "sum")
  r_Dbotleft <- terra::focal(im2,w = D_dia_botleft, fun = "sum")

  rr_Dtopleft = sum((r_Dtopleft == 2))
  rr_Dbotleft = sum((r_Dbotleft == 2))

  sum.Dtl <- raster::cellStats(rr_Dtopleft,"sum")
  sum.Dbl <- raster::cellStats(rr_Dbotleft,"sum")

  all.px = sum.Dho + sum.Dve + sum.Dtl + sum.Dbl
  diag.px = sum.Dtl + sum.Dbl
  orth.px = sum.Dho + sum.Dve


Directions = data.frame(vertical = sum.Dve / all.px,
                        horizontal = sum.Dho / all.px,
                        topleft.bottomright = sum.Dtl / all.px,
                        topright.bottom.left = sum.Dbl / all.px,
                        diag.px = diag.px,
                        orth.px = orth.px)

  return(Directions)

}



## RootScape Metrics

#¤ input image should be segmented, unskeletonized image

#' RootScape relies on Landscapemetrics as working horse
#'
#' @param im segmented, unskeletonized raster
#' @param indexD please specify depth
#' @param metrics which ,metrics should be calculated from the available ones in 'landscapemetrics::calculate_lsm()'. A selection is used as default
#' Interpretation here: ...
#'
#' @return a bunch of metric values
#' @export
#'
#' @examples RootScapeObject  = RootScapeMetrics(image)
RootScapeMetrics = function(im,indexD, metrics = c( "lsm_c_ca","lsm_l_ent","lsm_c_pd","lsm_c_np","lsm_c_pland",
                                                    "lsm_c_area_mn","lsm_c_area_cv","lsm_c_enn_mn","lsm_c_enn_cv")){

  rsm = landscapemetrics::calculate_lsm(im, directions = 8, neighbourhood = 8,what = metrics)
  rsm$object = ifelse(rsm$class == 0,"deletable","root")
  rsm$object = ifelse(is.na(rsm$class),"root",rsm$object )
  rsm$depth = indexD
  rsm = dplyr::distinct(rsm)
  rsm = dplyr::filter(rsm, object != "deletable")
  rsm$id = NULL
  rsm$class = NULL
  rsm$level = NULL
  rsm$layer = NULL

  return(rsm)
}






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




## Halo function

#' Halo creates a buffer around pixel bigger than 0
#'
#' @param im segmented raster
#' @param width buffer around roots in px, the rhizosphere extent (exudate diffusion distance) is cited as 2mm (1-12mm) (source), but higher values have been suggested (Finzi )
#' @param halo.only set TRUE if only the buffer around roots should be returned (the rhizosphere only)
#' @return raster output
#' @export
#'
#' @examples buffIMG = Halo(im = im, width = 10, halo.only = T)
Halo = function(im,width=1, halo.only = T){
  im = im / max(raster::values(im),na.rm=T)
  im2 = im
  ## circular kernel
  k0 = matrix(c(1,1,1,1,0,1,1,1,1), nrow = 3, ncol = 3)

  itr = 1
  while(itr <= width){
    im2 <- terra::focal(im2,w = k0, fun = "sum") #%>% suppressWarnings()
    itr = itr + 1
  }

  out.im = sum((im2 >= 1),silent = T)-1 #%>% suppressWarnings()

  if(halo.only  == TRUE){
    out.im = out.im - im
  }
  return(out.im)
}


#' Coloration of the image
#'
#' @param img raster with 3 color bands
#' @param r weight for first channel - typically red
#' @param g weight for second channel - typically green
#' @param b weight for third channel - typically blue
#'
#' @return a vector with chromatic coordinates,luminosity, brightness, luminosity, color values, saturation
#' @export
#'
#' @examples colorvector = Tube.coloration(img = imrgb)
Tube.coloration = function(img,r=0.2126,g=0.7152,b=0.0722){
  vr = raster::values(img[[1]])
  vg = raster::values(img[[2]])
  vb = raster::values(img[[3]])
  mean.r = mean(vr,na.rm=T )
  mean.g = mean(vg,na.rm=T )
  mean.b = mean(vb,na.rm=T )

  hsl =rgb2hsv(r = mean.r, g = mean.g, b = mean.b)
  intensity = vr + vg + vb
  lum.gray =  vr * r + vg * g + vb *b
  mean.intensity = mean(intensity,na.rm=T) %>% round(4)
  mean.lum = mean(lum.gray,na.rm=T) %>% round(4)
  rcc = mean(vr / intensity,na.rm = T) %>% round(4)
  gcc = mean(vg / intensity,na.rm = T) %>% round(4)
  bcc = mean(vb / intensity,na.rm = T) %>% round(4)
  colordf = data.frame(rcc = rcc,gcc = gcc,bcc = bcc,
                       hue = hsl[1], saturation = hsl[2], luminosity = hsl[3],
                       red = mean.r,green = mean.g,blue = mean.b)
  return(colordf)
}



#### Image Texture

# uses the glcm package
#' Texture corresponds to glmc::glmc()
#'
#' @param img.color image with three color chanels. Will be converted in grayscale.
#' @param grays number of gray shades. Documentation is lacking, see: ?glmc::glmc()
#' @param window convolution window size e.g., c(3,3)
#' @param metrics texture metrics based on illumination differences. see:: ?glmc::glmc() for available methods
#'
#' @return raster with each layer correspond to on metric
#' @export
#'
#' @examples texture(img.color, 7, c(9,9), metrics = "second_moment")
texture = function(img.color,grays = 7, window = c(9,9), metrics = c("variance","second_moment") ){
  img.gray = (img.color[[1]]*0.21 + img.color[[2]]*0.72 + img.color[[3]]*0.07)/255
  tx.im  =  glcm::glcm(img.gray,
                       n_grey = grays,
                       statistics = metrics)
}


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
#' @examples root.ticc = root.thickness(kimuralength = 300,rootpx = 9500,300)
root.thickness = function(kimuralength,rootpx,dpi){
  px.per.length = rootpx / kimuralength
  thiccness = round((px.per.length/118) / (dpi/2.54),5)
  return(thiccness) # cm
}
