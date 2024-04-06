### Root Thickness


root.tickness = function(kimuralength,rootpx,dpi){
  px.per.length = rootpx / kimuralength
  thiccness = (px.per.length/118) / (dpi/2.54)
  return(thiccness) # cm
}
