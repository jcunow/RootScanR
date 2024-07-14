#' Takes a continues depth map and bins it to a specified range
#'
#' @param depthmap 1-layer raster, takes output from create.depthmap()
#' @param nn bin width
#' @param round.option choose the binning operation. Available are "rounding", "ceiling", and "floor".
#'
#' @return raster with input depths in bins
#' @export
#'
#' @examples image = binning(depthmap,nn = 5)
binning = function(depthmap,nn,round.option = "standard"){

  if(round.option == "rounding"){
    im = nn * round(depthmap / nn,0)
  }
  if(round.option == "ceiling"){
    im = nn * ceiling(x / nn)
  }
  if(round.option == "floor"){
    im = nn * floor(x / nn)
  }

  return(im)
}



#' Cuts out a zone of interest
#'
#' @param rootpic the image which should be cut
#' @param binned.map image which supplies the cut condition
#' @param indexD the condition
#' @param nn bin width in binned.map. Should be the same as used in 'binning()'
#' @param silent verbose
#'
#' @return a cutout of the root pic image
#' @export
#'
#' @examples raster_image = zone.fun(rootpic, binned.map, nn=5, silent = F)
zone.fun = function(rootpic,binned.map,indexD,nn = 5,silent =F){
  max.NAs = round((1-(1/(max(raster::values(binned.map),na.rm=T)/nn)))+0.02,4)
  # match dimensions
  ex.bm = c(0,dim(binned.map)[1],0,dim(binned.map)[2])
  ex.rp = c(0,dim(rootpic)[1],0,dim(rootpic)[2])
  if(any(ex.bm!=ex.rp)){
    rootpic = raster::t(rootpic)
  }


  r = (rootpic)
  if(raster::extent(rootpic)[2] != raster::extent(binned.map)[2] | raster::extent(rootpic)[4] != raster::extent(binned.map)[4]){
    r = raster::crop(r,ex.bm)
  }
  # set uninterested zones NA
  raster::values(r)[raster::values(binned.map) != indexD ] <- NA
  raster::values(r)[ is.na(raster::values(binned.map))] <- NA
  r = terra::rast(r)
  # non NA values need to cover at least ... % of average depth slice pixel number
  coverNA = round(sum(is.na(raster::values(r[[1]]))) / (raster::ncell(r)),2)
  if( coverNA  < max.NAs){
    r = r %>% terra::trim()
  }else{
    raster::values(r) = NA # set the remaining values NA ?
    if(silent == F){
      print( paste0("Depth: ",indexD,"cm. Not enough informative pixels.",
                    " In the whole Image, ", coverNA*100,"% are NAs after cutting this Depth Slice."," Expected NA% is: ~",(max.NAs-0.02)*100 ))
    }

  }
  r = raster::raster(r)
  return(r)
}



# Adressing rotational Bias

#' RotationZones
#'
#' @param rootpic the "to be cut" image
#' @param nn  number of total cuts along rotation axis
#' @param k specify which cuts to keep. Must be <= nn
#' @param mm limit the region along the tube = c(start,end). Adjust to your tube dimensions!
#'
#' @return raster, cut along rotation axis
#' @export
#'
#' @examples rotationZone1 = zone.rotation.fun(rootpic, k = c(1,2), nn = 7, mm = c(1500,3000))
zone.rotation.fun = function(rootpic,k=c(3,4),nn = 5,mm = c(2000,5000)){

  if(!is.array(rootpic) ){
    img0 = raster::as.array(rootpic)
  }else{
    img0 = rootpic
  }


  m = seq(mm[1],mm[2],length = mm[2])
  ## rotation wise test

  img0 = img0[,m,]


  # lower row for the kth bin
  q = (k[1]*(nrow(img0)/nn)) %>% floor()
  q = q +1
  p = (k[2]*(nrow(img0)/nn)) %>% floor()

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
  img00 = img00 %>% terra::trim()
  img00 = raster::raster(img00)

  return(img00)
}
