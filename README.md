# MinirhizotronR
a collection of R code to analyze minirhiztron tubes

> ### The workflow typically goes like this:
> * (Not yet available in R) Affine Keyfeature Stitching of all Scans corresponding to a particular Tube in Time
>     1) Try AffineStitcher.py; works well for feature rich images like forest soils
>     2) Use ImageJ's Addon Stitcher with predetermined Image position (ca. 100px overlap)
> * (Not in R) RootDetector applied to stitched Images;  Peters et al. (2023) As good as human experts in detecting plant roots in minirhizotron images but efficient and reproducible: the convolutional neural network “RootDetector”. *Sci Rep* **13**, Software at https://github.com/ExPlEcoGreifswald/RootDetector
> * (optional) Soil Surface Position Estimation & Rotational Position; it's recommended to measure soil & rotation position *in-situ*; otherwise, tape-cover is used instead as proxy 
>    1) Function: SoilSurfPosE()
>    2) Function: RotPosE()
> *   Rotation Censoring keeps the image area which is present in all sequential images and removes stand-alone areas. This is based on the rotational center and removes edge areas.
>    1) Function: RotCensor()
> * Phase Shifted Cosine Depthmap
>    1) Function: PhShCosDepthmap()
>  * Zonal Statistics calcluates a bunch of metrics for each depth layer
>    1) Function: Binning() - creates depth bands with specified range, small diameter tubes need wider bands 
>    2) Function: KimuraLength() - estimates Root Length based on Kimura et al. 1999 [...]
>    3) Function: TurnTurn() - Estimates amount of new, constant, and disappearing root pixels & root length
>    4) Function: Texture() - Calculates Texture Metrics
>    5) Function: ColorCoords() - retrieves color coordinates, brightness, and luminosity
>       




