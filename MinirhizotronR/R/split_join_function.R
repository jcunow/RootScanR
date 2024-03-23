# split and join images
#### Translate PNG to TIFF
library(png);library(tiff);library(stringr)
#path = "C:/Users/jocu0013/Desktop/blended_images/AffineStitcher/2023_025/"
path = "C:/Users/jocu0013/Desktop/Oulanka/Scans_Blended/Full/Oulanka2023_03/"
dir.out = "C:/Users/jocu0013/Desktop/Oulanka/Scans_Blended/Split/Oulanka2023_03/"

# split
#' Split Images
#'
#' Use if too many roots are present for the Root Detector AI image segmentation
#' @param path Input path
#' @param dir.out Output path
#' @param pattern only include images with this pattern
#'
#' @return two sets of images
#' @export
#'
#' @examples split_im(path,dir.out, ".tiff) = top.im; bot.im
split_im = function(path,dir.out,pattern = ".tiff"){
  file.ls = list.files(path = path, pattern = pattern,ratio = 0.5)
  ## split the image
  for (i in file.ls) {
    im = readTIFF(paste0(path,i))
    im.top = im[,1:(dim(im)[2]*ratio),]
    im.dwn = im[,(dim(im)[2]*ratio+1):dim(im)[2],]
    writeTIFF(im.top,where = paste0(dir.out,"Split_top_",i))
    writeTIFF(im.dwn,where = paste0(dir.out,"Split_dwn_",i))
  }
}


# join
## fuse the image
#' Fuse two Images
#'
#' this function is intended to to be used after using the root detector
#' @param path Input path
#' @param dir.out Output path
#' @param pattern only include images with this pattern
#'
#' @return one complete image
#' @export
#'
#' @examples join_im(path,out.path) = img
join_im = function(path,dir.out,pattern = "skeleton"){
  ### segmented
  dir.ls = paste0(list.dirs(path = path)[-1],"/")
  file.ls = list.files(path = dir.ls, pattern = pattern)
  top.files =  file.ls[stringr::str_detect(file.ls,pattern = "Split_top")]
  dwn.files =  file.ls[stringr::str_detect(file.ls,pattern = "Split_dwn")]

  for (i in 1:length(top.files)) {
    im.top = readPNG(paste0(dir.ls[i+36],top.files[i]))
    im.dwn = readPNG(paste0(dir.ls[i],dwn.files[i]))
    im.all = abind::abind(im.top,im.dwn, along = 2 )
    file.name= top.files[i] %>% stringr::str_remove(pattern="Split_top_") %>% stringr::str_remove(pattern=".tiff.segmentation.png")
    writeTIFF(im.all,where = paste0(dir.out,"FullSegmented_",file.name,".tiff"))

  }
}

