####### create the depthmap

#' Create A Phase-Shifted, Tilt-Amplitude Sine Depth Map
#'
#' @param im input image
#' @param mask indicating which pixels are foreign objects like tape. The mask will be: mask = im[[1]] - im[[2]] if 'RootDetector' format is used
#' @param tube.thicc Diameter of minirhizotron tubes
#' @param tilt Minirhiztron tube insertion angle (typically 30-45 degrees)
#' @param dpi Image resolution
#' @param start.soil indicates soil boundary 0cm. Can be retrieved from 'SoilSurfE()' but in-situ calibration is recommended
#' @param center.offset rotational center. Set 0 if top facing tube side is perfectly in the middle
#'
#' @return raster image
#' @export
#'
#' @examples map = create.depthmap(im,mask,start.soil = 290 )
create.depthmap = function(im, mask,
                           tube.thicc = 7,tilt = 45,dpi = 300,
                           start.soil = 0,center.offset = 0){

  # used for sin()
  radiant = pi/180
  # law of sines
  tilt.factor = sin((180-90-tilt) * radiant) / sin(90 * radiant)
  # tilt correct
  tube.thicc.tilted = round(tube.thicc * tilt.factor,3 )

  # dpi-depth realtion
  px.to.cm.h = 1/(dpi/2.54)

  # get dims
  target.col = dim(im)[1]
  target.row = dim(im)[2]



  # simulate a sine wave function across one row
  df1 = seq(0*pi,2*pi,2*pi/(target.col-1))

  ### CORE FUNCTION, apply the function with the amplitude corresponding to the tilt and a phase corresponding to the rotation center
  # creates a cosine shaped curved shifted by the amount of rotation offset
  #df11 = (cos(df1+(pi*(1-center.offset))))*(tube.thicc.tilted/2)+ (tube.thicc.tilted/2) # center.offset 1 = no offset
  df11 = (cos(df1+(pi*(center.offset))))*(tube.thicc.tilted/2)+ (tube.thicc.tilted/2) # center.offset = 0 no offset


  # stack up rows and adding flat depth to each row
  df = array(dim = c(target.row,target.col))
  for (ii in 1:target.row) {
    df[ii,] = df11+(ii*px.to.cm.h * tilt.factor) # adds progressive depth to each row
  }
  # add soil surface offset estimated from tape cover
  df.depthmap = df + (start.soil* px.to.cm.h)

  # masking tape
  masked.depthmap= raster::raster(df.depthmap)
  # set mask NA
  raster::values(masked.depthmap)[raster::values(mask) == 1] <- NA

  return(masked.depthmap)

}

