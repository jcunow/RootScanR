## Rotation Estimate

RotationE = function(path,Session,file.name,tape.brightness = 0.95){
  # print
  print(paste0("currently processing Tube: ",file.name))


  ### read image
  img1 <- readTIFF(paste(path,Session,file.name,sep=''))

  ## add one row of red tape pixel
  red.line = array(dim = c(dim(img1)[1],extra.rows,dim(img1)[3]))
  red.line[,,2:3] <- 0
  red.line[,,1] <- max(img1[,,1])
  img1 = abind::abind(red.line,img1,along = 2)

  r.img1 = brick(img1)
  ### make the search are smaller
  r.img1 = crop(r.img1,extent(0,search.area,0,1))

  # identify distinct pixel groups
  r1 = RStoolbox::unsuperClass(r.img1,nClasses = 3)
  # determine average group
  clust.center = apply(r1$model$centers,1,mean)
  # silver tape should have highest luminance across clusters -> select max lum.cluster (but not close to == 1 [pure white?])
  clust= which(clust.center ==max(clust.center[clust.center < tape.brightness]))
  # identify the end of tape by rowsum threshold[]
  rr1 = r1$map == clust


  ## determine 0cm offset from the middle
  rsums = rowSums(rr1)
  # bin into two parts & take the middle one (my images have a complete tape part and a partial tape part.
  # Partial tape indicates the tube surface. Here, we assume that the partial tape is well placed !! (?)
  # needs manual calibration -- work in progress)
  bin = dplyr::ntile(rsums,2)
  zero.rotation2[ii] = median(which(bin == 2 ))
  zero.rotation.prop.offset2[ii] =   round(zero.rotation2[ii]/ (dim(img1)[1]),3)
  # end
  if(ii == length(file.ls)){
    print("### Complete !!!")
  }
  print(paste0("columns cecked: ",i))
}

