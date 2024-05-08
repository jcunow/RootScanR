# MinirhizotronR
a collection of R code to analyze iamges from minirhiztron tube

 ## The workflow could look like this:
 ### Step 1: Bringing seperate Images together
  * Option 1)  (Not yet available in R) Affine Keyfeature Stitching of all Scans corresponding to a particular Tube in Time; see Try AffineStitcher.py in StitchR repository. Works well for feature rich images and less good for homogneneous images.
  * Option 2)  (Not yet available in R) Use ImageJ's Addon Stitcher with predetermined Image position (ca. 100px overlap)
        
### Step 2: Image segmentation and root detection (+Turnover)
 * Option 1)  (Not in R) “RootDetector” applied to stitched Images;
   Peters et al. (2023)  *Sci Rep* **13**,
   Software at https://github.com/ExPlEcoGreifswald/RootDetector
   
### Step 3: Calibration 
  * Soil Surface Position Estimation; it's recommended to measure soil  position *in-situ*; otherwise, tape-cover is used
  * Option 1) Scanner Rotation Position Estimation; it's recommended to measure soil rotation *in-situ* instead; otherwise, tape-cover is used here
  * Option 2) Rotation shift detection; estimates rotation in relation to a reference image with known rotation offset
  * Rotation Censoring; remove unique image sections only present in at timepoint

### Step 4: Create a depthmap
requires: soil surface estimate, rotation center estimate, tube insertion angle, tube diameter, scan resolution
 * Option 1) Traditional Depthmap ignoring tube diameter depth but considering insertion angle
 * Option 2) Phase shifted sine depth mapping including tube diameter and tube insertion angle 

### Step 5: Zoning - Create analytical units
Determine the resolution of anaylysis. Currently, the package only supports analysis within discrete depthbands. 
Future releases will tackle response functions on continous depth.

### Step 6: Description and Stats
Estimate parameters:
* Root Cover
* Root Length
* Root Thickness
* Peat Coloration
* Peat Texture
* Rotation Bias
* Turnover




