####### create the depthmap

### manual measurements
## cases where tape ends followed by a gap, followed by soil 0cm depth
# add extra rows from manual estimation
# soil.extra = data.frame(Tube =  paste0("T0",seq(37,72)),
#                         Plus.0 = c(rep(0,3),340,0,310,320,rep(0,10),590,300,530,0,920,920,0,270,1100,750,0,1320,0,220,410,0,500,290,0))
# soil.0plus = soil.0
# soil.0plus$tape.end = soil.0$tape.end + soil.extra$Plus.0



create.depthmap = function(path,Session,file.name,output.name,output.path,
                           tube.thicc = 7,tilt = 45,dpi = 300,
                           start.soil = 1,center.offset = 0){

  # used for sin()
  radiant = pi/180
  # law of sines
  tilt.factor = sin((180-90-tilt) * radiant) / sin(90 * radiant)
  # tilt correct
  tube.thicc.tilted = round(tube.thicc * tilt.factor,3 )

  # dpi-depth realtion
  px.to.cm.h = 1/(dpi/2.54)

  # get dims
  im = readTIFF(paste0(path,Session,file.name))
  target.col = dim(im)[1]
  target.row = dim(im)[2]

  # create mask to set tape px NA later
  mask = brick(im)
  mask = mask$layer.1 - mask$layer.2
  mask = t(mask)

  # simulate a sine wave function across one row
  df1 = seq(0*3.141,2*3.141,2*3.141/(target.col-1))
  # apply the function with the amplitude corresponding to the tilt
  # creates a cosine shaped curved shifted by the amount of rotation offset
  df11 = (cos(df1+(3.141*(1-center.offset))))*(tube.thicc.tilted/2)


  # stack up rows and adding flat depth to each row
  df = array(dim = c(target.row,target.col))
  for (ii in 1:target.row) {
    df[ii,] = df11+(ii*px.to.cm.h) # adds progressive depth to each row
  }
  # add soil surface offset estimated from tape cover
  df.depthmap = df - (soil.start* px.to.cm.h * tilt.factor)

  # masking tape
  masked.depthmap= raster(df.depthmap)
  # set mask NA
  values(masked.depthmap)[values(mask) == 1] <- NA

  ## write depth maps to disk (merge later)
  writeRaster(masked.depthmap,paste0(outputpath,"/",output.name),overwrite = T)

}

