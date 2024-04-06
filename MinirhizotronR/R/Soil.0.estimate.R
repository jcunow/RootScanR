##### create a unique depth map for each Tube & Session


## Create a base depth layer for minirhizotron tubes
# rm(list = ls()); gc()
# disclaimer:
# This only accounts for the curvature and the angle of the tube, however,
# the issue of where the soil starts == 0cm == offset still needs a solution.

# library(raster)
# library(ggplot2)
# library(scico)
# library(EBImage)
# library(metR)
# #library(EBImageExtra)
# library(tiff)
# library(png)
# library(jpeg)
# library(dplyr)
# library(RStoolbox)
# library(stringr)

## --> not all needed !

# Tube and Scan specs (300 dpi)
# column and row assignment can change depending on the object class


### for each file of each the sessions
# path = "C:/Users/jocu0013/Desktop/Oulanka/Scans_FullSkeleton/"
# path.depthmap = "C:/Users/jocu0013/Desktop/Oulanka/DepthMaps/"
# Session = "Oulanka2023_03/"



# one scanner window dims to estimate conversion rates

# its actually not thickness but the distance from the bottom to the top part of the tube;
# should be equal to tube thickness at 0degree insertion angle

# andthe rotation currently doesn't cover 360 degrees
# target.col = 2273 # scanned columns = rotation
# target.row = 2550 # scanned rows = depth
# image.width = 19.6 # cm (19.24)
# image.height = 21.6 # cm (21.59)

# px.to.cm.h = image.height/target.row # cm / px
# px.to.cm.w = image.width/target.col # cm / px


## calibration, distance from tape to true 0 in px
# soil.extra = data.frame(Tube =  paste0("T0",seq(37,72)),
#                         Plus.0 = c(rep(0,3),340,0,310,320,rep(0,10),590,300,530,0,920,920,0,270,1100,750,0,1320,0,220,410,0,500,290,0))

#' Soil surface estimation
#'
#'the function locates the soil surface by the presence of tape
#' @param path where segmented images are located
#' @param Session measurement campaign
#' @param thresh how much of the image needs to be covered in tape to be considered purposfully obscured %
#' @param dpi scan dpi to convert tape-soil overlap into px units
#' @param overlap how much of the soil is not visible due to overlaping tape extending into the soil cm
#' @param soil.extra enter a dataframe with first column Tube idientifier, and second column the calibration offset. This allows the user to specify how close the tape matches the measured soil surface in px.
#'
#' @return values indicating where the tape ends, also estimates the rotational position based on tape proportions
#' @export
#'
#' @examples soil.0.estimate(path,Session,
#' soil.extra = data.frame(
#' Tube =  paste0("T0",seq(37,72)),
#' soil.0plus = 0)) = values...

soil.0.estimate = function(path = NULL,
                           Session = NULL,
                           thresh = 0.33,
                           dpi = 300,
                           overlap = 0.5,
                           soil.extra = data.frame(Tube =  paste0("T0",seq(37,72)),
                                                   soil.0plus = 0)
                           ){

  file.ls = list.files(paste0(path,Session))
  Tube.ID.start = stringr::str_locate(file.ls, patter = "T0")[,1]
  Tube.ID.end = stringr::str_locate(file.ls, patter = ".tif")[,1]-1
  Tube.ID = stringr::str_sub(file.ls,start=Tube.ID.start,end = Tube.ID.end)

  # row index storage
  rw.ind = c()
  # store proportion of tape where the threshold occurred
  tape.prop = c()
  tape.end = c()
  zero.rotation = c()
  zero.rotation.prop.offset = c()
  # if median is better
  zero.rotation2 = c()
  zero.rotation.prop.offset2 = c()
  t1 = Sys.time()
  soil.start = c()
  ## revise if necessary!!!!
  ls.of.exceptions = c(21)#plot57 # tape does not reflect soil surface
  ls.of.exceptions = 0
  extra.rows = 100
  for(ii in 1:(length(file.ls))
  ){
    # print
    print(paste0("currently processing Tube: ",ii))


    ### Session 1
    img1 <- tiff::readTIFF(paste(path,Session,file.ls[ii],sep=''))
    ## add one row of red tape pixel
    red.line = array(dim = c(dim(img1)[1],extra.rows,dim(img1)[3]))
    red.line[,,2:3] <- 0
    red.line[,,1] <- max(img1[,,1])
    img1 = abind::abind(red.line,img1,along = 2)

    r.img1 = raster::brick(img1)
    ### make the search are smaller
    r.img1 = raster::crop(r.img1,extent(0,0.33,0,1))

    # identify distinct pixel groups
    r1 = RStoolbox::unsuperClass(r.img1,nClasses = 3)
    # determine average group
    clust.center = apply(r1$model$centers,1,mean)
    # silver tape should have highest luminance across clusters -> select max lum.cluster (but not close to == 1 [pure white?])
    clust= which(clust.center ==max(clust.center[clust.center < 0.95]))
    # identify the end of tape by rowsum threshold[]
    rr1 = r1$map == clust

    ######## exceptional cases ####
    if(ii == ls.of.exceptions){
      clust= which(clust.center ==min(clust.center))
      # identify the end of tape by rowsum threshold[]
      rr1 = r1$map != clust
      i=1
      # tresh = 0.33 # how much of the rotation needs to be covered in identified tape
      while (
        #sum(rr1[,i])/length(rr1[,i]) >tresh &!
        sum(rr1[,i])/length(rr1[,i]) > tresh |
        sum(rr1[,i+24])/length(rr1[,i]) > tresh |
        sum(rr1[,i+48])/length(rr1[,i]) > tresh)
        i=i+1
      # when loop exceeds image limits, we assume there is no tape = no offset
      if(i > dim(rr1)[2]){
        i = 1
      }
    }

    ##### regular cases
    if(ii != ls.of.exceptions){
      # iterate over depth and check proportion covered by tape (>0.2 is considered >0 soil depth)
      # checks also 0.1 cm (+12 rows) and 0.2 cm (+24rows) and 0.3 cm (+36) down if the tape reappears  (in case a row impurities or other reasons for missclassification)
      i=1
      #tresh = 0.33 # how much of the rotation needs to be covered in identified tape
      while (
        #sum(rr1[,i])/length(rr1[,i]) >tresh &!
        sum(rr1[,i])/length(rr1[,i]) > tresh |
        sum(rr1[,i+24])/length(rr1[,i]) > tresh |
        sum(rr1[,i+48])/length(rr1[,i]) > tresh) {
        i=i+1
        # when loop exceeds image limits, we assume there is no tape = no offset
        if(i > dim(rr1)[2]){
          i = 1
        }
      }
    }
    # store row index to determine soil 0cm offset
    rw.ind[ii] = i
    # remove a fixed amount of offset which corresponds to the tape that was installed deeper than 0 cm to prevent light intrusion
    rw.ind[ii] = rw.ind[ii] - round(overlap*dpi/2.54) # cm tape in the ground converted as pixel offset from soil 0
    #start.soil = rw.ind # in case
    tape.end[ii] = rw.ind[ii] + round(overlap*dpi/2.54)
    tape.prop[ii]= sum(rr1[,i])/length(rr1[,i])

    ## determine 0cm offset from the middle
    rsums = rowSums(rr1)
    # bin into two parts & take the middle one (my images have a complete tape part and a partial tape part.
    # Partial tape indicates the tube surface. Here, we assume that the partial tape is well placed !! (?)
    # needs manual calibration -- work in progress)
    bin = dplyr::ntile(rsums,2)
    zero.rotation[ii] = floor(max(which(bin == 2 )) - min(which(bin == 2 )))
    zero.rotation2[ii] = median(which(bin == 2 ))
    zero.rotation.prop.offset[ii] =   round(zero.rotation[ii]/ (dim(img1)[1]),3)
    zero.rotation.prop.offset2[ii] =   round(zero.rotation2[ii]/ (dim(img1)[1]),3)
    # end
    if(ii == length(file.ls)){
      print("### Complete !!!")
    }
    print(paste0("columns cecked: ",i))



  }



  # object which contains soil == 0cm offset in px
  soil.0 = data.frame(rw.ind = rw.ind-extra.rows,
                      tape.end = tape.end-extra.rows,
                      tape.prop = tape.prop,
                      file.ls = file.ls,
                      Tube = stringr::str_sub(file.ls,start = stringr::str_locate(file.ls,pattern = "T")[,1],end = stringr::str_locate(file.ls,pattern = "T")[,1]+3), #
                      zero.rotation = zero.rotation,
                      zero.rotation.prop.offset = zero.rotation.prop.offset,
                      zero.rotation2 = zero.rotation2,
                      zero.rotation.prop.offset2 = zero.rotation.prop.offset2)
  t2 = Sys.time();round(t2-t1)



  #### cases where tape ends followed by a gap, followed by soil 0cm depth
  # add extra rows from manual estimation/ calibration [px]
  soil.0plus = soil.0
  soil.0plus$tape.end = soil.0$tape.end + soil.extra$Plus.0


}



#' Create individual cosine shifted depth maps for each of the tubes
#'
#' @param path path
#' @param path.depth path to store depth maps
#' @param Session sampling campaign identifier
#' @param soil.0 comes from the soil.0.estimation function
#' @param tube.thicc diameter of the minirhizotron tube
#' @param tilt insertion angle of the minirhiztron
#'
#' @return a depthmap which accounts for varying tube insertion depths, insertion angle, and tube diamter
#' @export
#'
#' @examples create.depthmap(path,Session,soil.0) = depthmap
create.depthmap = function(path = NULL,
                           path.depth = NULL,
                           Session = NULL,
                           soil.0 = soil.0, # from soil.0.estimate
                           tube.thicc = 7,
                           tilt = 45
                           ){
  file.ls = list.files(paste0(path,Session))
  Tube.ID.start = stringr::str_locate(file.ls, patter = "T0")[,1]
  Tube.ID.end = stringr::str_locate(file.ls, patter = ".tif")[,1]-1
  Tube.ID = stringr::str_sub(file.ls,start=Tube.ID.start,end = Tube.ID.end)

  tube.thicc = tube.thicc # cm outer tube diameter (6.35 inner diameter,7 outer diameter)
  tilt = tilt # minirhizotron insertion angle 0-90; 90 == vertical, 0 == horizontally
  tube.thicc.tilted = round(tube.thicc * sin((180-90-tilt) * ( pi / 180 )) / sin(90 * ( pi / 180 )),3 )
  print(paste0("with tube diameter of: ",tube.thicc," and a tilt of: ",tilt,", the "))

  ts3 = Sys.time()
  ### for each file of each the sessions
  for(i in 1:length(file.ls)
  ){
    # get dims
    im = tiff::readTIFF(paste0(path,Session,file.ls[i]))
    target.col = dim(im)[1]
    target.row = dim(im)[2]

    # simulate a sine wave function across one row
    df1 = seq(0*3.141,2*3.141,2*3.141/(target.col-1))
    # apply the function with the amplitude corresponding to the tilt
    center.offset = soil.0$zero.rotation.prop.offset2[i]
    # creates a cosine shaped curved shifted by the amount of rotation offset
    df11 = (cos(df1+(3.141*(1-center.offset))))*tube.thicc.tilted/2


    # stack up rows and adding flat depth to each row
    df = array(dim = c(target.row,target.col))
    for (ii in 1:target.row) {
      df[ii,] = df11+(ii*px.to.cm.h) # adds progressive depth to each row
    }
    # add soil surface offset estimated from tape cover
    df.depthmap = df - (soil.0$tape.end[i]* px.to.cm.h)
    #norm.img= (df.depthmap+20)/120
    norm.img= raster::raster(df.depthmap)
    ## write depth maps to disk (merge later)
    raster::writeRaster(norm.img,paste0(path.depthmap,Session,Tube.ID[i]),overwrite = T)
  }
  ts4 = Sys.time();ts4-ts3

}







# ##### check and load pre-exisitng data
# tryCatch(load("C:/Users/jocu0013/Desktop/Oulanka/Data/soil0_Data.RData"),
#          error = function(e){
#            print("No dataset found. Define a new one.")
#            soil.0.exp = data.frame(rw.ind = NA,tape.end = NA,tape.prop=NA,
#                                    file.ls = NA,Tube = NA,zero.rotation = NA,
#                                    zero.rotation.prop.offset = NA,zero.rotation2 = NA,
#                                    zero.rotation.prop.offset2 = NA)
#          }) # export data
#
#
# ## merge new + old data
# soil.0.exp = rbind(soil.0.exp,soil.0plus)
#
# ## save data
# save(soil.0.exp,file = "C:/Users/jocu0013/Desktop/Oulanka/Data/soil0_Data.RData")

