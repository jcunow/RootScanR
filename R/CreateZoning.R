#' Create a buffer (halo) around non-zero pixels
#'
#' @param img SpatRaster/matrix/array - segmented image
#' @param width numeric - buffer width in pixels (default: 2)
#' @param halo.only logical - if TRUE, returns only the buffer zone (default: TRUE)
#' @param kernel character - shape of the thickening kernel: "circle" or "diamond"
#'
#' @return SpatRast - buffer zone around non-zero pixels
#' @export
#'
#'
#' @examples
#' data(seg_Oulanka2023_Session03_T067)
#' img <- terra::rast(seg_Oulanka2023_Session03_T067)
#' Halo(img)
Halo = function(img, width=2, halo.only = TRUE, kernel = "circle"){

  # flexible input
  im <- load_flexible_image(img, output_format = "spatrast", normalize = TRUE  )


  ## circular kernel
  if(kernel == "circle"){
    k0 = matrix(c(1,1,1,1,0,1,1,1,1), nrow = 3, ncol = 3)
  }else if(kernel == "diamond"){
    k0 = matrix(c(0,1,0,1,0,1,0,1,0), nrow = 3, ncol = 3)
  }else{
    print("kernel name does not match choices")
  }


  # Apply the focal operation iteratively
  for(itr in 1:width){
    im2 <- terra::focal(im, w = k0, fun = "sum", na.policy = "omit")
  }

  out.im = (im2 >= 1) *1

  if(halo.only  == TRUE){
    out.im = out.im - im
    out.im = terra::subst(out.im, from = -1, to = 0)

  }
  return(out.im)
}




#' Bin continuous depth values into discrete intervals
#'
#' @param depthmap SpatRaster/matrix/array - continuous depth values
#' @param nn numeric - bin width
#' @param round.option character - binning method: "rounding", "ceiling", or "floor"
#' @return SpatRaster - binned depth values
#' @export
#'
#' @examples
#' data(seg_Oulanka2023_Session01_T067)
#' img = terra::rast(seg_Oulanka2023_Session01_T067)
#' mask = img[[1]] - img[[2]]
#' mask[mask == 255] <- NA
#' img = img[[2]]
#' depthmap = create.depthmap(img,mask,start.soil = 290 )
#' binned.map = binning(depthmap,nn = 5)
binning = function(depthmap,nn,round.option = "rounding"){


  # flexible input
  img <- load_flexible_image(depthmap, output_format = "spatrast", normalize = FALSE  )

  if(round.option == "rounding"){
    im = nn * round(depthmap / nn,0)
  }
  if(round.option == "ceiling"){
    im = nn * ceiling(depthmap / nn)
  }
  if(round.option == "floor"){
    im = nn * floor(depthmap / nn)
  }
im
}



#' Extract a zone of interest based on depth values
#'
#' @param rootpic SpatRaster/matrix/array - source image
#' @param binned.map SpatRaster/matrix/array - binned depth map
#' @param indexD numeric - depth index to extract. Usually used to index in a loop to iterate over depths.
#' @param nn numeric - bin width that has been used in binning
#' @param select.layer Integer. Specifies which layer to use if the input is a multi-band image. Default is `2`.
#' @param silent logical - suppress messages
#' @return SpatRaster - extracted zone
#' @export
#'
#' @examples
#' data(seg_Oulanka2023_Session01_T067)
#' img = terra::rast(seg_Oulanka2023_Session01_T067)
#' mask = img[[1]] - img[[2]]
#' mask[mask == 255] <- NA

#' depthmap = create.depthmap(img,mask,start.soil = 0, tube.thicc = 3, select.layer = 2 )
#' binned.map = binning(depthmap,nn = 5)
#' image.zone = zone.fun(img, binned.map, indexD = 10, nn=5, silent = FALSE,select.layer = 2)
zone.fun = function(rootpic,binned.map,indexD= 0,nn = 5,silent =FALSE, select.layer = NULL){

  # flexible input
  rootpic <- load_flexible_image(rootpic, select.layer = select.layer, output_format = "spatrast", normalize = FALSE  )
  binned.map <- load_flexible_image(binned.map, select.layer = select.layer, output_format = "spatrast", normalize = FALSE  )

  # match dimensions
  ex.bm = c(0,dim(binned.map)[1],0,dim(binned.map)[2])
  ex.rp = c(0,dim(rootpic)[1],0,dim(rootpic)[2])
  if(any(ex.bm!=ex.rp)){
    rootpic = terra::t(rootpic)
  }
  ex.bm = c(0,dim(binned.map)[1],0,dim(binned.map)[2])
  ex.rp = c(0,dim(rootpic)[1],0,dim(rootpic)[2])
  terra::ext(binned.map) = ex.bm
  terra::ext(rootpic) = ex.rp


  r = (rootpic)
  if(terra::ext(rootpic)[2] != terra::ext(binned.map)[2] | terra::ext(rootpic)[4] != terra::ext(binned.map)[4]){
    r = terra::crop(r,ex.bm)
  }
  # set uninterested zones NA
  terra::values(r)[terra::values(binned.map != indexD) ] <- NA
  terra::values(r)[ is.na(terra::values(binned.map))] <- NA
  # non NA values need to cover at least ... % of average depth slice pixel number
  coverNA = round(sum(is.na(terra::values(r[[1]]))) / (terra::ncell(r)),2)
  if( coverNA  < 0.95){
    r = terra::trim(r)
  }else{
    terra::values(r) = NA # set the remaining values NA ?
    if(silent == FALSE){
      cat("Depth: ",indexD,"cm. Not enough informative pixels. ",
                    coverNA*100,"% are NAs in this Depth Slice.")
    }

  }
  return(r)
}



# Adressing rotational Bias

#' Extract zones based on rotation axis
#'
#' @param rootpic array/matrix/SpatRaster - source image
#' @param k numeric vector - which cuts to keep (length 2)
#' @param kk numeric - total number of cuts
#' @param mm numeric vector - region limits (length 2)
#' @param select.layer Integer. Specifies which layer to use if the input is a multi-band image. Default is `NULL`.
#' @return SpatRaster - extracted rotation zone
#' @export
#'
#' @examples
#' data(seg_Oulanka2023_Session01_T067)
#' img = terra::rast(seg_Oulanka2023_Session01_T067)
#' rotationZone = zone.rotation.fun(img, k = c(1,2), kk = 7, mm = c(1500,3000))
zone.rotation.fun = function(rootpic,k=c(3,4),kk = 5,mm = c(2000,5000),select.layer = NULL){

  img0 <- load_flexible_image(rootpic, select.layer = select.layer, output_format = "array", normalize = FALSE  )

if(exists("mm")){
  m = seq(mm[1],mm[2],length = mm[2])
  ## rotation wise test
  img0 = img0[,m,]
}

  # lower row for the kth bin
  q = floor((k[1]*(nrow(img0)/kk)))
  q = q +1
  p = floor((k[2]*(nrow(img0)/kk)))

  img00 = img0


  if(is.na(dim(img0)[3] > 1 )){
    img00[1:dim(img00)[1],1:dim(img00)[2]] <- NA
  }else{
    img00[1:dim(img00)[1],1:dim(img00)[2],] <- NA
  }

  if(is.na(dim(img0)[3] > 1 )){
    img00[c(p:q),] = img0[c(p:q),]
  }else{
    img00[c(p:q),,] = img0[c(p:q),,]
  }


  img00 = terra::rast(img00)
  img00 = terra::trim(img00)

  return(img00)
}
