## create Depthmaps

SoilSurfE = function(path,Session,file.name,search.area = 0.33, tape.tresh = 0.33,dpi = 300,
                     tape.overlap = 0.5,tape.brightness = 0.95,extra.rows = 100 ){
  # print
  print(paste0("currently processing Tube: ",file.name))


  ### read image
  img1 <- readTIFF(paste(path,Session,"/",file.name,sep=''))

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



  # iterate over depth and check proportion covered by tape (>0.2 is considered >0 soil depth)
  # checks also 0.1 cm (+12 rows) and 0.2 cm (+24rows) and 0.3 cm (+36) down if the tape reappears  (in case a row impurities or other reasons for missclassification)
  i=1
  while (
    #sum(rr1[,i])/length(rr1[,i]) >tresh &!
    sum(rr1[,i])/length(rr1[,i]) > tape.tresh |
    sum(rr1[,i+24])/length(rr1[,i]) > tape.tresh |
    sum(rr1[,i+48])/length(rr1[,i]) > tape.tresh) {
    i=i+1
    # when loop exceeds image limits, we assume there is no tape = no offset
    if(i > dim(rr1)[2]){
      i = 1
    }
  }

# store row index to determine soil 0cm offset
rw.ind = i
# remove a fixed amount of offset which corresponds to the tape that was installed deeper than 0 cm to prevent light intrusion
rw.ind = rw.ind - extra.rows - round(tape.overlap*dpi/2.54) # cm tape in the ground converted as pixel offset from soil 0
#start.soil = rw.ind # in case
tape.end = rw.ind + round(0.5*300/2.54)

out = data.fram(soil0 = rw.ind, tape.end = tape.end  )
return(out)
print(paste0("columns cecked: ",i))
}


