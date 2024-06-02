
library(MinirhizotronR)
library(dplyr)

path = "C:/Users/jocu0013/Desktop/Oulanka/Scans_FullSkeleton/Oulanka2023_03/"
im.ls = list.files(path)
s0 = read.csv("C:/Users/jocu0013/Desktop/Oulanka/OulankaTubeInsertionAngle.csv")
s0$soil0 = round(s0$Soil.Excess * (300/2.54))

k=36
nn = 5
root.frame = data.frame(depth = seq(-20,140, nn))
for (i in im.ls) {
# Depthmapping
 k=k+1
 im =  raster::brick(paste0(path,i))
 mask = (im[[1]] - im[[2]] ) / 255
 mask = raster::t(mask)
 r0 = RotationE(im)
 center.offset = r0  / dim(im)[1]
 DepthMap = create.depthmap(im = im, mask = mask,start.soil = s0$soil0[k], center.offset = center.offset, tilt = s0$InsertionAngle[k ])
 # Stats

 bm = binning(depthmap = DepthMap,nn = nn)
 bm = t(bm)
 raster::extent(bm)<- raster::extent(im)
 nb = unique(raster::values(bm)) %>% na.omit() %>% sort()

 new.col = data.frame( rep(NA,length(root.frame$depth)))
 Tube.name =  stringr::str_sub(i, start =  stringr::str_locate( i, pattern = "T0")[1],end = stringr::str_locate(i,".tiff")[1]-1)
 root.frame = data.frame(root.frame,"rootpx" = rootpixel,"vertical" = new.col,"horizontal" = new.col,"topright" = new.col,"botright" = new.col)
 ## rename the newly added 5 columns last columns
 colnames(root.frame)[dim(root.frame)[2]-4] = paste0(Tube.name,"_","rootpx")
 colnames(root.frame)[dim(root.frame)[2]-3] = paste0(Tube.name,"_","vertical")
 colnames(root.frame)[dim(root.frame)[2]-2] = paste0(Tube.name,"_","horizontal")
 colnames(root.frame)[dim(root.frame)[2]-1] = paste0(Tube.name,"_","topright")
 colnames(root.frame)[dim(root.frame)[2]] = paste0(Tube.name,"_","botright")
 for (j in nb) {
  rootzone = zone.fun(rootpic = im[[2]], binned.map = bm, indexD = j, nn = nn )
  rootpixel = px.sum(root.zone = rootzone)
  Directs = Directionality(rootzone)
  print(j)
  ## assign values to new columns
  root.frame[root.frame$depth == j,dim(root.frame)[2]-4] = rootpixel
  root.frame[root.frame$depth == j,dim(root.frame)[2]-3] = Directs$vertical
  root.frame[root.frame$depth == j,dim(root.frame)[2]-2] = Directs$horizontal
  root.frame[root.frame$depth == j,dim(root.frame)[2]-1] = Directs$topright.bottom.left
  root.frame[root.frame$depth == j,dim(root.frame)[2]] =   Directs$topleft.bottomright

}
}

ggplot(data = root.frame,aes(depth, T037_vertical)) + geom_point() + geom_smooth()

root.frame2 = reshape2::melt(root.frame,id = "depth", value.name = "rootpx")
root.frame2$Plot = stringr::str_sub(root.frame2$variable,start=3) %>% as.numeric()
Oulanka.treat <- read.table(paste0("C:/Users/jocu0013/Desktop/Oulanka/Data/trt_ls.csv"), header = TRUE, sep = ";")
root.frame3 = merge(root.frame2,Oulanka.treat,by.x = "Plot", by.y = "Plotid")
root.frame3$Block = factor(root.frame3$Block)
root.frame3$Grazing = factor(root.frame3$Grazing)
#plot
root.frame3 %>% na.omit() %>% filter(depth >= 0) %>% filter(rootpx != 0 | depth == 0) %>%
  ggplot(aes(depth,rootpx,color = Grazing)) + geom_point() + geom_smooth(span=0.5,se=T,size = 2)

root.frame.exp = root.frame3
save(root.frame.exp, file =  "C:/Users/jocu0013/Documents/rootdepthDATA.RData")

root.frame4 = root.frame3 %>% na.omit() %>% filter(depth >= 0) %>% filter(rootpx != 0 | depth == 0)

library(mgcv)
mod.gam = gam(data= root.frame4, rootpx ~ Grazing + s(depth, bs = "cr", by = Grazing, k=15) +  s(Block, bs = 're')  )
mod.gam %>% summary()
mod.gam %>% gam.check()
flexplot::visualize(mod.gam, formula = rootpx ~ depth)
performance::check_model(mod.gam)
performance(mod.gam)

