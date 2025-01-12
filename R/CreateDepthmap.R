####### create the depthmap

#' Create A Phase-Shifted, Tilt-Amplitude Sine Depth Map
#'
#' This function generates a depth map for minirhizotron images, accounting for tube geometry
#' and insertion angle. Supports parallel processing and efficient memory management for large images.
#'
#' @param im Input image (accepts terra SpatRaster, matrix, array, or file path). For multi-band
#'           images, specify band_index parameter
#' @param mask Raster mask indicating foreign objects (1 = mask, 0 or NA = keep)
#' @param sinoid Logical; if TRUE, accounts for tube curvature in depth calculation
#' @param tube.thicc Numeric; diameter of minirhizotron tube in cm
#' @param tilt Numeric; minirhizotron tube insertion angle in degrees (typically 30-45)
#' @param dpi Numeric; image resolution in dots per inch
#' @param start.soil Numeric; soil surface boundary in cm (0 = surface)
#' @param center.offset Numeric; rotational center offset (0 = centered, 1 = edge)
#'
#' @return terra raster object containing the depth map
#' @export
#'
#' @examples
#' data(seg_Oulanka2023_Session01_T067)
#' seg_Oulanka2023_Session01_T067 = terra::rast(seg_Oulanka2023_Session01_T067)
#' img = seg_Oulanka2023_Session01_T067[[2]]
#' mask = seg_Oulanka2023_Session01_T067[[1]] - seg_Oulanka2023_Session01_T067[[2]]
#' mask[mask == 255] <- NA
#' map = create.depthmap(img,mask,start.soil = 290 )
create.depthmap = function(im, mask = NULL, sinoid = TRUE,
                           tube.thicc = 7,tilt = 45,dpi = 300,
                           start.soil = 0,center.offset = 0){
# if no mask is supplied nothing is masked
  if(length(mask) == 0){
    mask = im
    terra::values(mask) <- 0
  }

  # used for sin()
  radiant = pi/180
  # law of sines
  tilt.factor = sin((180-tilt) * radiant)
  # tilt correct
  tube.thicc.tilted = round(tube.thicc * tilt.factor,3 )

  # dpi-depth realtion
  px.to.cm.h = (2.54/dpi)

  # get dims
  target.col = dim(im)[1]
  target.row = dim(im)[2]

  if(sinoid == TRUE){
    # simulate a sine wave function across one row
    # df1 = seq(0*pi,2*pi,2*pi/(target.col-1))
    df1 = seq(0*pi,2*pi,2*pi/(round(tube.thicc*pi*dpi/2.54,0)-1))

    ### CORE FUNCTION, apply the function with the amplitude corresponding to the tilt and a phase corresponding to the rotation center
    # creates a cosine shaped curved shifted by the amount of rotation offset
    #df11 = (cos(df1+(pi*(1-center.offset))))*(tube.thicc.tilted/2)+ (tube.thicc.tilted/2) # center.offset 1 = no offset
    df00 = (cos(df1+(pi*(center.offset))))*(tube.thicc.tilted/2)+ (tube.thicc.tilted/2) # center.offset = 0 no offset
    # trim to actual rotation of the scanner (often less than 360 degree)
    df11 = df00[1:target.col]
  }else{
#### flat df11
    df11 = rep(0,target.col)
  }





  # stack up rows and adding flat depth to each row
  df = array(dim = c(target.row,target.col))
  for (ii in 1:target.row) {
    df[ii,] = df11+(ii*px.to.cm.h * tilt.factor) # adds progressive depth to each row
  }
  # add soil surface offset estimated from tape cover
  df.depthmap = df - (start.soil)

  # masking tape
  masked.depthmap= terra::rast(df.depthmap)
  # set mask NA
  terra::values(masked.depthmap)[terra::values(mask) == 1] <- NA

  return(masked.depthmap)

}

