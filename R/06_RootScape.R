## RootScape Metrics

#¤ input image should be segmented, unskeletonized image

#' RootScape relies on Landscapemetrics as working horse
#'
#' @param im segmented, unskeletonized raster
#' @param indexD please specify depth
#' @param metrics which ,metrics should be calculated from the available ones in 'landscapemetrics::calculate_lsm()'. A selection is used as default
#' Interpretation here: ...
#'
#' @return a bunch of metric values
#' @export
#'
#' @examples RootScapeObject  = RootScapeMetrics(image)
RootScapeMetrics = function(im,indexD, metrics = c( "lsm_c_ca","lsm_l_ent","lsm_c_pd","lsm_c_np","lsm_c_pland",
                                      "lsm_c_area_mn","lsm_c_area_cv","lsm_c_enn_mn","lsm_c_enn_cv")){

  rsm = landscapemetrics::calculate_lsm(im, directions = 8, neighbourhood = 8,what = metrics)
  rsm$object = ifelse(rsm$class == 0,"deletable","root")
  rsm$object = ifelse(is.na(rsm$class),"root",rsm$object )
  rsm$depth = indexD
  rsm = dplyr::distinct(rsm)
  rsm = dplyr::filter(rsm, object != "deletable")
   rsm$id = NULL
   rsm$class = NULL
   rsm$level = NULL
   rsm$layer = NULL

  return(rsm)
}
