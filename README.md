# MinirhizotronR - R package

<br />

## Overview
MinirhizotronR is a package which has been designed to make the analysis of minirhizotron root scans just a little less stupid. 
The package relies on prior Image segmentation. We recommend using *RootDetector* or *RootPainter*. This packages enables the user to map various root features to a continues depth distributions. 
The user has the choice to include the tube insertion angle and tube diameter - or ignore it (?). Rotational bias of root prevalence can be tested. Let me know if you have questions, improvements, or found bugs & errors.

````
# install the package
if (!require("devtools")) {
  install.packages("devtools")
}
devtools::install_github("jcunow/MinirhizotronR")
````


____________________________________


<br />

## Calibration 
More accurate results require easy to gather *in-situ* calibration of each tube. Traditionally, the assumption is that all tube have been installed perfectly, i.e., exact same insertion angle, insertion depth, same scanner rotation position, and tape overshoot for each tube. MinirhizotronR offers the possibility  to approximate the soil start and rotation center based an assumption of well-taped tubes.

1. Soil Surface Position Estimation 
````
start.soil = SoilSurfE(img)
````

2. Scanner Rotation Position, Rotation Differences Between Two Scans, and Removal of Scan Edges only present in one Scan and not the other.
````  
RotationE(img)
RotShiftDet(img1,img2)
RotCensor(img,center.offset = 0.3,fixed.rotation = FALSE)

# under construction: estimating rotation by root distribution (assuming top-tube bias)
LR_rhythmicity(rotation.cuts, root.px, period= number.rotation.slices)$phase * (dim(roation.cuts)[1]/number.rotations.slices) %>% round()
````

<br />

## The core function: Create a depthmap
requires: soil surface estimate, rotation center estimate, tube insertion angle, tube diameter, scan resolution

1. Phase shifted, trimmed sine depth mapping including tube diameter and tube insertion angle 
````
depthmap = create.depthmap(img,mask, tube.thicc = 7, tilt = 45, dpi = 300, start.soil = 0,center.offset = 0)
````

<br />

## Zoning - Create analytical units
Creates discrete sub-units of the whole image. This determines the resolution of the depth distributions. Currently, the package only supports analysis within discrete depth bands.
Future releases will tackle response functions on continuous depth.

1. Depth-wise zoning
````
binned.map = binning(depthmap, nn = 1,"rounding")
root.zone = zone.fun(img, binned.map, indexD = 5, nn = 1)
````

<br />

##  Feature Extractions
1. Parameter from segmented image cut
````
RootScapeMetrics(root.zone, metrics = c("lsm_c_ca","lsm_c_pland","lsm_c_enn_mn"))
Root.px = px.sum(root.zone)
````

2. Parameter from skeletonized image cut
````
root.zone.skl = skeletonize(root.zone)
RL = RootLength(root.zone.skeleton)
````

3. Combine both
````
root.thickness(Root.px, RL)
````

4. Parameter from RGB image cut
````
Rhizosphere = Halo(root.zone, width = 5, halo.only = FALSE)
soil = rgb.img[Rhizosphere == 1] < - NA  
Soil.texture(soil)
soil.color = Tube.coloration(soil)
````

5. Turnover - either two images from different time points, or the 'RootDetector' root tracking output 
````
Turnover.TC(root.zone1,root.zone2, method = "kimura",dpi = 300, unit = "cm")
Turnover.DPC(root.zone.dpc, product.layer = 2, decay.layer = 1, im.return = F)
````

<br />


## Rotational Bias
1. Rotational zoning,  Feature Extraction, 
````
rotation.rootpx = data.frame(kk = 1:12,px = NA)

for(i in rotation.rootpx){
rotation.zone = zone.rotation.fun(img,k = c(i-1,i),kk = 12, mm = c(1500,5000))
rotation.rootpx$px[i,] = px.sum(rotation.zone)
}


````

2. Test Amplitude Difference
Can be used to text e.g., differences in top-bottom-tube bias between varying insertion angle or tube diameter
To test differences between sinus fits, we make use of the diffCircadian package https://github.com/diffCircadian/diffCircadian
````
fitSinCurve(rotation.rootpx$px,period = 12)
LTTest_diff_amp(img1_rotation.rootpx$px,img2_rotation.rootpx$px)
````
