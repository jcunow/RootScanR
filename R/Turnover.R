#' Calculate global root production and root turnover from temporal comparison
#'
#' @param im.t1 SpatRaster object for timepoint 1
#' @param im.t2 SpatRaster object for timepoint 2
#' @param method Analysis method: "kimura" or "rootpx"
#' @param unit Unit of root length measurement (only for method = "kimura"). Default: "cm"
#' @param dpi Image resolution (only for method = "kimura"). Default: 300
#'
#' @return data.frame containing:
#'   - standingroot_t1: Standing roots at first timepoint
#'   - standingroot_t2: Standing roots at second timepoint
#'   - production: Root production between timepoints
#'   - newroot%per_t1: Percentage of new roots compared to starting conditions
#'   - newroot%per_t2: Percentage of new roots at second timepoint
#' @export
#'
#' @import raster
#'
#' @examples
#'   data(skl_Oulanka2023_Session01_T067)
#'   data(skl_Oulanka2023_Session03_T067)
#'   time1 <- terra::rast(skl_Oulanka2023_Session01_T067)
#'   time2 <- terra::rast(skl_Oulanka2023_Session03_T067)
#'   turnover.values <- Turnover.TC(
#'     im.t1 = time1,
#'     im.t2 = time2,
#'     method = "kimura",
#'     verbose = TRUE)
Turnover.TC = function(im.t1, im.t2,method="kimura",unit = "cm",dpi = 300){

  if(method == "rootpx"){
    px1 = px.sum(im.t1)
    px2 = px.sum(im.t2)
    px.prod =  px2-px1
    px.turn1 =  round(px.prod / px1,4)
    px.turn2 =  round(px.prod / px2,4)
  }
  if(method == "kimura"){

    px1 = RootLength(im.t1,dpi = dpi,unit = unit)
    px2 = RootLength(im.t2,dpi = dpi,unit = unit)
    px.prod =  px2-px1
    px.turn1 =  round(px.prod / px1,4)
    px.turn2 =  round(px.prod / px2,4)
  }

  return(data.frame("standingroot_t1" = px1,"standingroot_t2" = px2,"production" = px.prod,"newroot%per_t1" = px.turn1,"newroot%per_t2" = px.turn2 ))
}



#' Extract Root Decay, New Root Production, and No-Change Roots (only 'RootDetector' images)
#'
#' @param img SpatRaster with three layers for production, decay, and stagnation
#' @param product.layer Integer indicating the production layer index (1-3)
#' @param decay.layer Integer indicating the decay & tape layer index (1-3)
#' @param blur.capture Threshold for pixel inclusion (0-1). Default: 0.95
#' @param im.return Logical: return images instead of values? Default: FALSE
#' @param include.virtualroots Logical: consider all roots present at any timepoint? Default: FALSE
#'
#' @return If im.return = FALSE: tibble with pixel sums and ratios
#'         If im.return = TRUE: list of SpatRaster layers for tape, constant, production, and decay
#' @export
#'
#' @examples
#' data(TurnoverDPC_data)
#' img = terra::rast(TurnoverDPC_data)
#' DPCs = Turnover.DPC(img = img, im.return = FALSE)
Turnover.DPC = function(img,product.layer = 2, decay.layer = 1, blur.capture = 0.95, im.return = FALSE, include.virtualroots = FALSE){
  l.indx = 1:3
  no.change.layer = which(!l.indx %in% c(product.layer,decay.layer))
  l.pr = img[[product.layer]]
  l.no = img[[no.change.layer]]
  l.dc = img[[decay.layer]]

  tape = l.no - l.dc
  tape = tape <= terra::global(tape,"min")[[1]]*blur.capture
  tape = tape * terra::global(l.pr,"max")[[1]]

  l.pr2 = l.pr - l.no - tape
  l.pr2 = l.pr2 >= terra::global(l.pr2,"max")[[1]]*blur.capture
  l.pr2 = l.pr2 * 1

  l.dc2 = l.dc - l.no - tape
  l.dc2 = l.dc2 >= terra::global(l.dc2)[[1]]*blur.capture
  l.dc2 = l.dc2 * 1

  l.no2 = l.no >= terra::global(l.no,"max")[[1]]*blur.capture
  l.no2 = l.no2 * 1

  tape = tape / terra::global(l.pr,"max")[[1]]

  if(im.return == TRUE){
    return(list("tape" = tape, "constant" = l.no2, "production" = l.pr2, "decay" = l.dc2))
  }else{
    tape.px = px.sum(tape)
    const.px = px.sum(l.no2)
    prodc.px = px.sum(l.pr2)
    decay.px = px.sum(l.dc2)
if(include.virtualroots == TRUE){
  newgrowth.ratio  = prodc.px / (prodc.px + const.px + decay.px)
  decay.ratio  = decay.px / (decay.px + const.px + prodc.px)
}else{
  newgrowth.ratio  = prodc.px / (prodc.px + const.px)
  decay.ratio  = decay.px / (decay.px + const.px)
}
  constant.ratio  = const.px / (prodc.px + decay.px + const.px)
  state.change = (decay.px + prodc.px) / (prodc.px + decay.px + const.px)
  return(dplyr::tibble("tape" = tape.px, "constant" = const.px,
                  "production" = prodc.px, "decay" = decay.px,
                  "newgrowth.ratio" = round(newgrowth.ratio,4),
                  "decay.ratio" = round(decay.ratio,4),
                  "constant.ratio" = round(constant.ratio,4)))
  }
}

