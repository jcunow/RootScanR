
#' Calculate Deep Drive Estimate for Root Systems
#'
#' This function analyzes root growth directions in relation to depth gradients.
#' It calculates the ratio of root pixels growing in the direction of the steepest
#' depth slope compared to all root pixels.
#'
#' @param DepthMap SpatRast object containing depth information
#' @param AngleMap Optional SpatRast with D8 root angles (from terra::terrain(v="flowdir"))
#' @param RootMap Optional SpatRast containing segmented, presence-absence root image
#' @param select.layerRM Integer. Specifies which layer to use if the input is a multi-band image. Default is `2`.
#' @param select.layerDM Integer. Specifies which layer to use if the input is a multi-band image. Default is `NULL`.
#' @param select.layerAM Integer. Specifies which layer to use if the input is a multi-band image. Default is `NULL`.
#'
#' @return If return_diagnostics=FALSE, returns numeric value representing the ratio of aligned root pixels
#'         to total root pixels. If TRUE, returns a list containing the ratio and diagnostic information.
#' @export
#'
#'
#' @examples
#' data(skl_Oulanka2023_Session01_T067)
#' im = ceiling(terra::rast(skl_Oulanka2023_Session01_T067)/255)
#' DepthMap = terra::t(create_depthmap(im,center.offset=0,tube.thicc=3.5))
#'
#' deep_drive(DepthMap = DepthMap, RootMap = im, select.layerRM = 2, select.layerAM = 2)
deep_drive = function(DepthMap=NULL, AngleMap=NULL, RootMap=NULL,
                      select.layerRM=NULL, select.layerDM=NULL, select.layerAM=NULL) {

  # Validation module
  tryCatch({
    # Check required inputs
    if (is.null(DepthMap)) {
      stop("DepthMap is required")
    }

    # Input type validation
    if (!is.null(select.layerRM) && !is.numeric(select.layerRM)) {
      stop("select.layerRM must be numeric")
    }
    if (!is.null(select.layerDM) && !is.numeric(select.layerDM)) {
      stop("select.layerDM must be numeric")
    }
    if (!is.null(select.layerAM) && !is.numeric(select.layerAM)) {
      stop("select.layerAM must be numeric")
    }

    # Layer selection validation
    if (!is.null(select.layerRM) && select.layerRM < 1) {
      stop("select.layerRM must be positive")
    }
    if (!is.null(select.layerDM) && select.layerDM < 1) {
      stop("select.layerDM must be positive")
    }
    if (!is.null(select.layerAM) && select.layerAM < 1) {
      stop("select.layerAM must be positive")
    }

    # Check if RootMap is provided when AngleMap is NULL
    if (is.null(AngleMap) && is.null(RootMap)) {
      stop("Either AngleMap or RootMap must be provided")
    }

    # Load and validate DepthMap
    DepthMap <- load_flexible_image(DepthMap, select.layer=select.layerDM,
                                    output_format="spatrast", normalize=FALSE)
    if (is.null(DepthMap)) {
      stop("Failed to load DepthMap")
    }

    # Ensure DepthMap has valid values
    if (terra::global(DepthMap, "isNA", na.rm=TRUE)[1] == terra::ncell(DepthMap)) {
      stop("DepthMap contains only NA values")
    }

  }, error = function(e) {
    stop(paste("Validation error:", e$message))
  })

  # Main function logic with additional error handling
  tryCatch({
    if(is.null(AngleMap)){
      RootMap <- load_flexible_image(RootMap, select.layer=select.layerRM,
                                     output_format="spatrast", normalize=T)
      if (is.null(RootMap)) {
        stop("Failed to load RootMap")
      }

      # align extents
      terra::ext(RootMap) <- terra::ext(DepthMap)

      # ensure positive Depth increments
      DepthMap = abs(DepthMap)

      # the D8 flowdir algorithm needs decreasing values
      dem = -DepthMap
      dem[RootMap != 1] <- NA
      dem = terra::t(terra::flip(dem))
      AngleMap = terra::terrain(dem, v="flowdir")
      AngleMap = terra::subst(AngleMap,
                              from=c(0,1,2,4,8,16,32,64,128),
                              to=c(NA,90,135,180,225,270,315,0,45))
    } else {
      AngleMap <- load_flexible_image(AngleMap, select.layer=select.layerAM,
                                      output_format="spatrast", normalize=FALSE)
      if (is.null(AngleMap)) {
        stop("Failed to load AngleMap")
      }
    }

    # [Rest of the original function code remains the same...]
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




    # Add validation for final calculation
    if (zonepx == 0) {
      warning("No valid pixels found for calculation")
      return(NA)
    }

    deep.drive.v = depthdrivepx / zonepx
    return(deep.drive.v)

  }, error = function(e) {
    stop(paste("Error in main function:", e$message))
  }, warning = function(w) {
    warning(paste("Warning in main function:", w$message))
  })
}
