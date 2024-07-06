### DemoFile

### note
# A sneak peak into the functions I want to make availble in an R package (name suggestions r welcome)
# You can find all the stuff with a tini tiny bit more documentation on: https://github.com/jcunow/MinirhizotronR
#
# Call to action:
# 1) You might need to install some of these package dependencies:
# raster, imagefx,  terra, landscapemetrics, RStoolbox, dplyr, stars, OpenImageR,  glcm, mmand, png, tiff, sf, stringr, cowplot, nloptr, minpack.lm
# 2) Adjust the Diameter of the Minirhizotron Tube to your size and specifiy the insertion angle. Change scanning resolution (default is 300 dpi). See create.depthmap()
# 3) Change path to the directory containing the segmented images





library(MinirhizotronR)
library(dplyr)
library(ggplot2)

path = "C:/Users/UserName/Images/"
#path.skeleton = "C:/Users/jocu0013/Desktop/Oulanka/Scans_FullSkeleton/Oulanka2022_03/"
im.ls = list.files(path)
#im.ls.skeleton = list.files(path.skeleton)

#s0 = read.csv("C:/Users/jocu0013/Desktop/Oulanka/OulankaTubeInsertionAngle.csv")
#s0$soil0 = round(s0$Soil.Excess * (300/2.54))

### Action !!!
s0 = data.frame(soil0 = 0, InsertionAngle = 30)
r0 = 2160/2 # image center row/ upside of the tube
tube.thicc = 5 # cm diameter

#k=0
k=1 # in case you have real calibrations on insertion angles and soil beginning for each tube
nn = 1 # depth slice intervall in cm
root.frame = data.frame(depth = seq(-10,99, nn))
for (i in im.ls) {
  # Depthmapping
  #k=k+1
  im =  raster::brick(paste0(path,i))
  #im.skeleton =  raster::brick(paste0(path.skeleton,im.ls.skeleton[k-36]))
  # mask = (im[[1]] - im[[2]] ) / 255
  mask = im[[1]]
  mask = raster::t(mask)
  im = im[[2]]/255

  #r0 = RotationE(im)
  r0 = dim(im[1])/2
  center.offset = r0  / (dim(im)[1]-r0) # 1 means the tube upside in perfectly in the middle of the image
  DepthMap = create.depthmap(im = im, mask = mask,start.soil = s0$soil0[k], center.offset = center.offset, tilt = s0$InsertionAngle[k ],tube.thicc = tube.thicc)


  bm = binning(depthmap = DepthMap,nn = nn)
  bm = raster::t(bm)
  raster::extent(bm)<- raster::extent(im)
  nb = unique(raster::values(bm)) %>% na.omit() %>% sort()

  new.col = data.frame( rep(NA,length(root.frame$depth)))
  ### adjust to however your naming convension works
  #Tube.name =  stringr::str_sub(i, start =  stringr::str_locate( i, pattern = "T0")[1],end = stringr::str_locate(i,".tiff")[1]-1)
  Tube.name =  paste0("T", stringr::str_sub(i,start=-8,end=-5))
  #root.frame = data.frame(root.frame,"rootlength" = new.col,"voidpx" = new.col,"rootpx" = new.col,"vertical" = new.col,"horizontal" = new.col,"topright" = new.col,"botright" = new.col,diag.px = new.col)
  root.frame = data.frame(root.frame,"voidpx" = new.col,"rootpx" = new.col)
  print(Tube.name)
  ## rename the newly added columns
  # colnames(root.frame)[dim(root.frame)[2]-7] = paste0(Tube.name,"_","rootlength")
  # colnames(root.frame)[dim(root.frame)[2]-6] = paste0(Tube.name,"_","voidpx")
  # colnames(root.frame)[dim(root.frame)[2]-5] = paste0(Tube.name,"_","rootpx")
  # colnames(root.frame)[dim(root.frame)[2]-4] = paste0(Tube.name,"_","vertical")
  # colnames(root.frame)[dim(root.frame)[2]-3] = paste0(Tube.name,"_","horizontal")
  # colnames(root.frame)[dim(root.frame)[2]-2] = paste0(Tube.name,"_","topright")
  # colnames(root.frame)[dim(root.frame)[2]-1] = paste0(Tube.name,"_","botright")
  # colnames(root.frame)[dim(root.frame)[2]] = paste0(Tube.name,"_","diag.px")


   colnames(root.frame)[dim(root.frame)[2]-1] = paste0(Tube.name,"_","voidpx")
   colnames(root.frame)[dim(root.frame)[2]] = paste0(Tube.name,"_","rootpx")

  for (j in nb) {
    # rootpx
    rootzone = zone.fun(rootpic = im, binned.map = bm, indexD = j, nn = nn )
    rootpixel = px.sum(root.zone = rootzone)
    # void
    void = rootzone
    raster::values(void) = 1-raster::values(rootzone)
    voidpixel = px.sum(root.zone = void)
    # # root length
    # rootzone.skeleton = zone.fun(rootpic = im.skeleton[[2]], binned.map = bm, indexD = j, nn = nn )
    # rootlength = RootLength(rootzone.skeleton,"cm",300)
    # directions
    # Directs = Directionality(raster::t(rootzone.skeleton),rotate = T)



    print(j)
    ## assign values to new columns
    # root.frame[root.frame$depth == j,dim(root.frame)[2]-7] = rootlength
    # root.frame[root.frame$depth == j,dim(root.frame)[2]-6] = voidpixel
    # root.frame[root.frame$depth == j,dim(root.frame)[2]-5] = rootpixel
    # root.frame[root.frame$depth == j,dim(root.frame)[2]-4] = Directs$vertical
    # root.frame[root.frame$depth == j,dim(root.frame)[2]-3] = Directs$horizontal
    # root.frame[root.frame$depth == j,dim(root.frame)[2]-2] = Directs$topright.bottom.left
    # root.frame[root.frame$depth == j,dim(root.frame)[2]-1] =   Directs$topleft.bottomright
    # root.frame[root.frame$depth == j,dim(root.frame)[2]] =   Directs$diag.px

    root.frame[root.frame$depth == j,dim(root.frame)[2]-1] = voidpixel
    root.frame[root.frame$depth == j,dim(root.frame)[2]] = rootpixel


  }
}


ggplot(data = root.frame,aes(depth, root.frame$T2c_1_voidpx)) + geom_point() + geom_smooth(span = 0.85)


### voidpx
meas.variable = colnames(root.frame)[seq(2,450,2)] %>% na.omit() # normal
root.frame2.vo = reshape2::melt(root.frame,id = "depth", value.name = "rootpx",measure.vars = meas.variable)
#root.frame2.vo$Plot = stringr::str_sub(root.frame2.vo$variable,start=3,end=4)  %>% as.numeric()
root.frame2.vo$Plot =stringr::str_sub(root.frame2.vo$variable,end = 5)
colnames(root.frame2.vo)[3] <- "void_px"
root.frame2.vo = root.frame2.vo[,-2]
### rootpx
meas.variable = colnames(root.frame)[seq(3,450,2)] %>% na.omit() # normal
root.frame2.ro = reshape2::melt(root.frame,id = "depth", value.name = "rootpx",measure.vars = meas.variable)
#root.frame2.ro$Plot = stringr::str_sub(root.frame2.ro$variable,start=3,end=4)  %>% as.numeric()
root.frame2.ro$Plot =stringr::str_sub(root.frame2.ro$variable,end = 5)
colnames(root.frame2.ro)[3] <- "root_px"
root.frame2.ro = root.frame2.ro[,-2]


root.frame2 = merge(root.frame2.vo,root.frame2.ro,by = c("Plot","depth"))


## add treatment ID
# Oulanka.treat <- read.table(paste0("C:/Users/jocu0013/Desktop/Oulanka/Data/trt_ls.csv"), header = TRUE, sep = ";")
# root.frame3 = merge(root.frame2,Oulanka.treat,by.x = "Plot", by.y = "Plotid")
# root.frame3$Block = factor(root.frame3$Block)
# root.frame3$Grazing = factor(root.frame3$Grazing)



#plotting
root.frame2 %>% na.omit() %>% filter(depth >= -10)  %>% filter(root_px > 0 | depth == 0) %>%
  ggplot(aes(depth,root_px)) + geom_text(aes(label = Plot),size = 1.9) + geom_smooth(span=0.55,se=F,size = 1)+scale_x_continuous(limits = c(-10,110))+ theme_classic()
root.frame.exp = root.frame4

## import export
# save(root.frame.exp, file =  "C:/Users/jocu0013/Documents/rootdepthDATA.RData")
# load("C:/Users/jocu0013/Documents/rootdepthDATA.RData")

tr = root.frame2 %>% na.omit() %>%  filter(root_px > 1 | (depth == 0 & root_px == 0 ))

#### accumulative curves
crf = tr %>% group_by(Plot) %>%
  #filter(depth>30) %>%
  arrange(depth) %>%
  mutate(cs = cumsum(scale(root_px,center = F,scale = F)))
crf = crf %>% group_by(Plot) %>%
  #filter(depth>50) %>%
  arrange(depth) %>%
  mutate(cs_percentroots = cumsum(root_px/(void_px+root_px)))
crf %>%
  #filter(Snow == "Control") %>%
  ggplot(aes(y = cs,x = depth,group = factor(Plot))) +
  geom_text(aes(label = factor(Plot)),size= 1.8) +
  #scale_y_continuous(limits = c(400000,1300000))+
  #geom_point() +
  geom_line()+
  ylab("cummulative roots")+
  #geom_smooth( span = 1, aes(color = "loess"),se = F)+
  ## verse rootpx
  #  geom_smooth(method="nls",formula = y ~Vmax *(1-exp(-x/tau)) ,se=F,
  #                method.args = list(start=c(tau=30,Vmax=400000)),
  #                aes(group=factor(Plot),color="asymptotic growth"),size = 0.1)+
  #
  # geom_smooth(method="nls",formula = y ~ a * exp(-b * exp(-c * x)),se=F,
  #               method.args = list(start=c(c=0.2,a=400000,b=2),nls.control(maxiter =200 )),
  #               aes(group=factor(Plot),color="gompertz"),size = 0.1)+
  #
  # geom_smooth(method="nls",formula = y ~ a * (1-exp(-b *(x)))**2 + c,se=F,
  #             method.args = list(start=c(a=0.6,b=0.1,c=0.05),nls.control(maxiter =200 )),
  #             aes(group=factor(Plot),color="Bertalanffy"),size = 0.1)+
  #
  # geom_smooth(method="nls",formula = y ~ (a*x)/ (b+x) ,se=F,
  #             method.args = list(start=c(a=0.6,b=0.1),nls.control(maxiter =200 )),
  #             aes(group=factor(Plot),color="Monod"),size = 0.1)+

  scale_x_continuous(limits = c(0,80)) + theme_classic()







######## model individual non-linear growth curves (and compare shape coefficients) #####

## estimate tau
tdf = data.frame(Plot=unique(crf$Plot),a_asymp.growth = NA,b_asymp.growth = NA,mae_asymp.growth = NA ,a_gompertz = NA, b_gompertz = NA, c_gompertz = NA,mae_gompertz = NA)
tdf = merge(tdf,Oulanka.treat, by.x = "Plot", by.y = "Plotid")
for(i in unique(tdf$Plot)){

  tryCatch({
    print(i)
    m1= crf %>% filter(Plot == i) %>% nls( formula = cs ~a*(1-exp(-depth/b)) ,start = c(b=50,a=40000),algorithm = "port",control = nls.control(maxiter=500))
    tdf$a_asymp.growth[tdf$Plot == i] <- coef(m1)["a"]
    tdf$b_asymp.growth[tdf$Plot == i] <- coef(m1)["b"]
    tdf$mae_asymp.growth[tdf$Plot == i] <- mae(m1)

    m2= crf %>% filter(Plot == i) %>% nls( formula = cs ~ a * exp(-b * exp(-c * depth)),start = c(c=0.2,a=400000,b=2),algorithm = "port",control = nls.control(maxiter=500) )
    tdf$a_gompertz[tdf$Plot == i] <- coef(m2)["a"]
    tdf$b_gompertz[tdf$Plot == i] <- coef(m2)["b"]
    tdf$c_gompertz[tdf$Plot == i] <- coef(m2)["c"]
    tdf$mae_gompertz[tdf$Plot == i] <- mae(m2)

  },error = function(e){})

}


