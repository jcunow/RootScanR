#' Global root production and root turnover from temporal comparison
#'
#' @param im.t1 brick raster timepoint 1
#' @param im.t2 brick raster timepoint 2
#' @param method choose between "kimura" & "rootpx". Adjust the input image accordingly
#' @param unit unit of root length. only applies if method = "kimura"
#' @param dpi image resolution. only applies if method = "kimura"
#'
#' @return standing roots at the start, standing roots at the end, root production, % of new roots compared to starting conditions, % of new roots at the second timepoint
#' @export
#' @details
#' Chosing kimura as method will return root length with default settings from RootLength()
#'
#'
#' @examples tunrover.values = TurnoverE(im.t1 = time1, im.t2 = time2, method= "kimura")
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



#' Estimates New Root Production, Root Decay, and Roots without change. Relies on 'RootDetector'
#'
#' @param img image in the 'RootDetector' format - one layer for production, one layer for decay, one layer for stagnation
#' @param product.layer layer indicating production
#' @param decay.layer layer indicating decay & tape
#' @param blur.capture pixel are included if:  value >= max value * blur.capture. Ensures that attenuated pixels (as a result of blurring or resizing) are also included
#' @param im.return return images instead of values?
#' @param include.virtualroots should all roots which were present at some point in one of the two time steps be considered? Consider decay in production ratio, and production in decay ratio?
#'
#' @return either pixel sums if im.return = F, or individual layers corresponding to tape, production, decay, and no change
#' @export
#'
#' @examples  PDCs = Turnover.PDC(img = img, im.return = F)
Turnover.PDC = function(img,product.layer = 2, decay.layer = 1, blur.capture = 0.95, im.return = FALSE, include.virtualroots = FALSE){
  l.indx = 1:3
  no.change.layer = which(!l.indx %in% c(product.layer,decay.layer))
  l.pr = img[[product.layer]]
  l.no = img[[no.change.layer]]
  l.dc = img[[decay.layer]]

  tape = l.no - l.dc
  tape = tape <= min(raster::values(tape))*blur.capture
  tape = tape * max(raster::values(l.pr),na.rm=T)

  l.pr2 = l.pr - l.no - tape
  l.pr2 = l.pr2 >= max(raster::values(l.pr2),na.rm=T)*blur.capture
  l.pr2 = l.pr2 * 1

  l.dc2 = l.dc - l.no - tape
  l.dc2 = l.dc2 >= max(raster::values(l.dc2),na.rm=T)*blur.capture
  l.dc2 = l.dc2 * 1

  l.no2 = l.no >= max(raster::values(l.no),na.rm=T)*blur.capture
  l.no2 = l.no2 * 1

  tape = tape / max(raster::values(l.pr),na.rm=T)

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

