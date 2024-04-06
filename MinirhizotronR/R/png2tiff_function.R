#' PNG2TIFF file translation
#'
#' Reads .png images and writes them as .tiff files elsewhere
#' @param path path wHERE png files are expected
#' @param path.out path where tiff files going to be stored
#'
#' @return .tiff file
#' @export
#'
#' @examples png2tiff(image.png) = image.tiff
png2tiff = function(path,path.out){
  for (i in file.ls) {
    file.ls = list.files(path = path, pattern = ".png")
    im = png::readPNG(paste0(path,i))
    tiff::writeTIFF(im,where = paste0(dir.out,stringr::str_replace( i, pattern = ".png",".tiff")))
  }
}
#
