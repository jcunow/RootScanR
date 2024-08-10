

# D8 Flow Direction using normalized slope-based kernels
#' D8 Flow Direction using slope-based kernels
#'
#' @param skl skeletonized SpatRaster Image
#' @param depthmap created by create.depthmap()
#' @param kernelsize argument used for adaptive.Kernelsize.Directionality(). Indicates the amount of neighboring pixels evaluated i.e., kernelsize 3 = 5x5 kernel
#' Consider readjusting normalized.center when changing kernel size.
#' @param normalized.center kernel center value
#' @param diag.weighted diagonal pixels are weighted 1/sqrt(2) to reflect increased diagonal pixel distance.
#' @param depth.positive depth values positive sign
#'
#' @return SpatRaster with deepest D8 neighbor-pixel angles:  0, 45, 90, 135, 180, 225, 270, 315
#' @export
#'
#' @examples
#' library(terra)
#' library(raster)
#'
#' data(skl_Oulanka2023_Session03_T067)
#' skl = terra::rast(skl_Oulanka2023_Session03_T067)[[2]]
#' depthmap = create.depthmap(skl,skl)
#' D8_FlowPath(skl,depthmap,kernelsize =2,normalized.center=-1)
D8_FlowPath <- function(skl,depthmap, diag.weighted = TRUE, kernelsize = 2,normalized.center = -1, depth.positive = TRUE) {


  if(terra::global(skl,"max",na.rm=TRUE) > 1){
    skl = ceiling(skl / 255)
  }

  # thin depth map to only contain vlaues where the skeleton is
  terra::values(depthmap)[terra::values(skl) != 1] <- NA


  kernels = adaptive.Kernelsize.Directionality(n=kernelsize,fill.value = NA,normalized.center,diag.weighted = diag.weighted)

  # Initialize a list to store the slope results for each direction
  results <- list()


  # Convolve the DEM with each kernel to compute slope in each direction
  for (i in 1:8) {
    slope <- terra::focal(depthmap, w = kernels[[i]], fun = sum, pad = TRUE, padValue = NA)
    results[[i]] <- slope
  }

  # Stack the slope results into a SpatRaster
  stack_slope <- terra::rast(results)
  if(depth.positive == TRUE){
    flow_dir <- terra::which.max(stack_slope)
  }else{
    flow_dir <- terra::which.min(stack_slope)
  }




  # Encode the flow direction using powers of 2
  directions <- c(0, 45, 90, 135, 180, 225, 270, 315)
  flow_dir[] <- directions[flow_dir[]]

  # Return the flow direction SpatRaster
  return(flow_dir)
}




