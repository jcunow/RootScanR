#### read, rename,  and rewrite
path= "C:/Users/jocu0013/Desktop/Oulanka/Scan_Raw/Oulanka2020_November/"
ls.file = list.files(path);ls.file

for (img in ls.file) {
  im = tiff::readTIFF(paste0(path,img))
  name = stringr::str_replace(img, pattern = "Ecfg",replacement = "Oulanka_2020")
  tiff::writeTIFF(im,paste0(path,name))
}

