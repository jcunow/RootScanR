## Halo function

#' Halo creates a buffer around pixel bigger than 0
#'
#' @param im segmented raster
#' @param width buffer around roots in px, the rhizosphere extent.
#' exudate diffusion distance is reported as 2mm (1-12mm) (Finzi et al. 2015, https://doi.org/10.1111/gcb.12816), but higher values have been suggested.
#' At 300 dpi, 1mm = 11.8px
#' @param halo.only set TRUE if only the buffer around roots should be returned (the rhizosphere only)
#' @return SpatRaster
#' @export
#'
#' @examples
#' library(terra)
#'
#' data(seg_Oulanka2023_Session01_T067)
#' img = terra::rast(seg_Oulanka2023_Session01_T067)
#' buffIMG = Halo(im = img, width = 2, halo.only = TRUE)
Halo = function(im,width=2, halo.only = TRUE){
  im = im / terra::global(im,"max",na.rm = TRUE)[[1]]
  im2 = im
  ## circular kernel
  k0 = matrix(c(1,1,1,1,0,1,1,1,1), nrow = 3, ncol = 3)

  itr = 1
  while(itr <= width){
    im2 <- terra::focal(im2,w = k0, fun = "sum") #%>% suppressWarnings()
    itr = itr + 1
  }

  out.im = sum(im2 >= 1) #%>% suppressWarnings()

  if(halo.only  == TRUE){
    out.im = out.im - im
  }
  return(out.im)
}




#' Takes a continues depth map and bins it to a specified range
#'
#' @param depthmap 1-layer raster, takes output from create.depthmap()
#' @param nn bin width
#' @param round.option choose the binning operation. Available are "rounding", "ceiling", and "floor".
#'
#' @return raster with input depths in bins
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



#' Cuts out a zone of interest
#'
#' @param rootpic the image which should be cut
#' @param binned.map image which supplies the cut condition
#' @param indexD the condition
#' @param nn bin width in binned.map. Should be the same as used in 'binning()'
#' @param silent verbose
#' @import dplyr
#'
#' @return a cutout of the root pic image
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

#' RotationZones
#'
#' @param rootpic the "to be cut" image
#' @param kk  number of total cuts along rotation axis
#' @param k specify which cuts to keep. Must be <= nn
#' @param mm limit the region along the tube = c(start,end). Adjust to your tube dimensions!
#' @import dplyr
#'
#' @return raster, cut along rotation axis
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
