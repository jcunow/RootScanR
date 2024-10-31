
#' Deep Drive Estimate
#'
#' @param DepthMap SpatRast with Depth information
#' @param AngleMap SpatRast with D8 root angles. See e.g., terra::terrain( v= "flowdir")
#' @param RootMap  If no Angle MAp is supplied, an AngleMap is calculated from SpatRast containing a segmented, presence-abscence root image.
#'
#' @return ratio of root px with an angle == to the steepest depth slope to the sum of all root px
#' @export
#'
#' @examples
#' data(skl_Oulanka2023_Session01_T067)
#' im = terra::rast(skl_Oulanka2023_Session01_T067)
#' DepthMap = terra::t(create.depthmap(im,center.offset=0,tube.thicc=3.5))
#'
#' deep.drive(DepthMap = DepthMap, RootMap = im)
deep.drive = function(DepthMap=NULL,AngleMap=NULL,RootMap = NULL){
  ## if not root angles are supplied than rootmap is used to create an AngleMap
  if(is.null(AngleMap)){
    # default chooses second layer
    if(dim(RootMap)[3] != 1){
      RootMap = RootMap[[2]]
    }

    if(terra::global(RootMap,"max",na.rm=TRUE)$max[1] > 1){
      RootMap = ceiling(RootMap / 255)
    }

    # align extents
    terra::ext(RootMap) <- terra::ext(DepthMap)
    # ensure positive Depth increments
    DepthMap = abs(DepthMap)
    # the D8 flowdir algorithm needs decreasing values
    dem = -DepthMap
    dem[RootMap != 1] <- NA
    dem = terra::t(terra::flip(dem))
    AngleMap = terra::terrain(dem,v = "flowdir")
    AngleMap = terra::subst(AngleMap, from = c(0, 1, 2, 4, 8, 16, 32, 64, 128), to = c( NA, 90,135,180,225,270,315,0,45))

  }

  # align orientation with AngleMap
  DepthMap = terra::t(terra::flip(DepthMap))
  # which pixel to go to reach the next deepest pixel in 8px neighbourhood
  mxslope <- terra::focal(DepthMap, w=c(3,3), fun=function(x)x[c(4,6,7:9)] - x[5])
  # correct for diagonal length
  mxslope[[c(3,5)]] = mxslope[[c(3,5)]] / sqrt(2)
  gg = terra::which.max(mxslope)
  # give meaningful labels
  gg = terra::subst(gg, from = c(1:5), to = c(270,90,225,180,135))


  a3 = AngleMap == 225 *1
  g<-gg
  terra::ext(g) <- terra::ext(a3)
  g[is.na(a3)]<-NA
  s3 = terra::zonal(a3,g,"sum")

  a1 = AngleMap == 270 *1
  g<-gg
  terra::ext(g) <- terra::ext(a1)
  g[is.na(a1)]<-NA
  s1 = terra::zonal(a1,g,"sum")

  a2 = AngleMap == 90 *1
  g<-gg
  terra::ext(g) <- terra::ext(a2)
  g[is.na(a2)]<-NA
  s2 = terra::zonal(a2,g,"sum")

  a4 = (AngleMap == 180)  *1
  g<-gg
  terra::ext(g) <- terra::ext(a4)
  g[is.na(a4)]<-NA
  s4 = terra::zonal(a4,g,"sum")

  a5 = (AngleMap == 135)  *1
  g<-gg
  terra::ext(g) <- terra::ext(a5)
  g[is.na(a5)]<-NA
  s5 = terra::zonal(a5,g,"sum")


  t.all = s1 + s2 + s3 + s4 +s5
  t.all[,1] = t.all[,1] / 5

  depthdrivepx = sum(c(s1$flowdir[s1$which.max == 270],s2$flowdir[s2$which.max == 90],
                       s3$flowdir[s3$which.max == 225],s4$flowdir[s4$which.max == 180],
                       s5$flowdir[s5$which.max == 135]),na.rm=TRUE)
  zonepx = sum(t.all$flowdir)


  deep.drive.v = depthdrivepx / zonepx
  return(deep.drive.v)

}
