## Image Resolution Alignment ##

# import.path= "C:/Users/jocu0013/Desktop/Oulanka/ExtraScan/Oulanka2022_July_highRes/"
# output.path= paste0(getwd(),"/")


#' Resize the Image
#'
#' useful if some scanning campaigns have used different dpi
#' @param import.path Input Path
#' @param output.path Output Path
#' @param height target height
#' @param width target width
#'
#' @return a image with new dimensions and resolution
#' @export
#'
#' @examples resize.image(import.path, output.path) = resized.image
resize.image = function(import.path,output.path,height = 2550,width= 2273){
  img.list = list.files(import.path, pattern = ".tiff")

  # sequentially read images, lower resolution and restore them
  t1 = Sys.time()
  for (i in 1:length(img.list)) {
    # read
    temp.im = OpenImageR::readImage(paste0(import.path,img.list[i]))
    # resample bilinear
    temp.im = OpenImageR::resizeImage(temp.im,
                                      method = "bilinear",
                                      height = 2550, width = 2273)
    # transformation to 0-1 range
    temp.im = imagefx::range01(temp.im)
    # write as tiff
    OpenImageR::writeImage(temp.im,file_name = paste0(output.path,img.list[i]))

  }
  t2 = Sys.time();floor(t2-t1)

}






