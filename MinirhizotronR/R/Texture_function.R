#### Image Texture

function(img.color,grays = 7, window = c(9,9), metrics = c("variance","second_moment") ){
  img.gray = (img.color[[1]]*0.21 + img.color[[2]]*0.72 + img.color[[3]]*0.07)/255
  tx.im  =  glcm::glcm(img.gray,
                       n_grey = grays,
                       statistics = metrics)
  }
