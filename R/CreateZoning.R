#' Create a buffer (halo) around non-zero pixels
#'
#' @param im SpatRaster/matrix/array - segmented image
#' @param width numeric - buffer width in pixels (default: 2)
#' @param halo.only logical - if TRUE, returns only the buffer zone (default: TRUE)
#' @param kernel character - shape of the thickening kernel: "circle" or "diamond"
#'
#' @return SpatRast - buffer zone around non-zero pixels
#' @export
#'
#' @import raster
#'
#' @examples
#' data(skl_Oulanka2023_Session03_T067)
#' im <- terra::rast(skl_Oulanka2023_Session03_T067)
#' Halo(im = im)
Halo = function(im, width=2, halo.only = TRUE, kernel = "circle"){

  if(terra::global(im,"max",na.rm=TRUE)$max[1] > 1){
    im = ceiling(im / 255)
  }

  im2 = im
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
    im2 <- terra::focal(im2, w = k0, fun = "sum", na.policy = "omit")
  }

  out.im = sum(im2 >= 1)

  if(halo.only  == TRUE){
    out.im = out.im - im
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
#' img = terra::rast(seg_Oulanka2023_Session01_T067[[2]])
#' mask = seg_Oulanka2023_Session01_T067[[1]] - seg_Oulanka2023_Session01_T067[[2]]
#' mask[mask == 255] <- NA
#' depthmap = create.depthmap(img,mask,start.soil = 290 )
#' binned.map = binning(depthmap,nn = 5)
binning = function(depthmap,nn,round.option = "rounding"){

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
#' @param indexD numeric - depth index to extract
#' @param nn numeric - bin width used in binning
#' @param silent logical - suppress messages
#' @return SpatRaster - extracted zone
#' @export
#'
#' @examples
#' data(seg_Oulanka2023_Session01_T067)
#' img = terra::rast(seg_Oulanka2023_Session01_T067[[2]])
#' mask = seg_Oulanka2023_Session01_T067[[1]] - seg_Oulanka2023_Session01_T067[[2]]
#' mask[mask == 255] <- NA
#' depthmap = create.depthmap(img,mask,start.soil = 290 )
#' binned.map = binning(depthmap,nn = 5)
#' image.zone = zone.fun(img, binned.map, indexD = 0, nn=5, silent = FALSE)
zone.fun = function(rootpic,binned.map,indexD= 0,nn = 5,silent =FALSE){
  max.NAs = round((1-(1/(terra::global(binned.map,"max",na.rm=TRUE)[[1]]/nn)))+0.02,4)
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
  if( coverNA  < max.NAs){
    r = terra::trim(r)
  }else{
    terra::values(r) = NA # set the remaining values NA ?
    if(silent == FALSE){
      print( paste0("Depth: ",indexD,"cm. Not enough informative pixels.",
                    " In the whole Image, ", coverNA*100,"% are NAs after cutting this Depth Slice."," Expected NA% is: ~",(max.NAs-0.02)*100 ))
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
#' @return SpatRaster - extracted rotation zone
#' @export
#'
#' @examples
#' data(seg_Oulanka2023_Session01_T067)
#' img = terra::rast(seg_Oulanka2023_Session01_T067)
#' rotationZone = zone.rotation.fun(img, k = c(1,2), kk = 7, mm = c(1500,3000))
zone.rotation.fun = function(rootpic,k=c(3,4),kk = 5,mm = c(2000,5000)){

  if(!is.array(rootpic) ){
    img0 = terra::as.array(rootpic)
  }else{
    img0 = rootpic
  }

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
