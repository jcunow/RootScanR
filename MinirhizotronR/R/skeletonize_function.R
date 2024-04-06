#' Skeletonize a segmented image
#'
#' @param img binary iamge
#' @param method which approach should be used
#'
#' @return binary image with area turned into linear features
#' @export
#'
#' @examples skeletonize(img,"gonzales") = skeleton_image
skeletonize_function = function(img,method = "gonzales"){
  if(method == "gonzales"){
    skeleton = dipr::thinning(img)
  }
  if(method == "erode"){
    skeleton = dipr::skeletonize(img)
  }
  if(method == "lantuejoul"){
    skeleton = mmand::skeletonise(method = method )
  }
  if(method == "beucher"){
    skeleton = mmand::skeletonise(method = method )
  }
  if(method == "hitormiss"){
    skeleton = mmand::skeletonise(method = method )
  }
  if(method %in% any(c("gonzales","erode","lantuejoul","beucher","hitormiss"))  ){
    print("method not appropriate")
  }
}
