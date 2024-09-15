##### Simple Lines ######

#' @title Example Image with circa vertical lines
#'
#'@description
#' Use this dataset to test directionality
#'
#' @author Johannes Cunow \email{johannes.cunow@gmail.com}
#' @usage data(lines_vertical)
#' @docType data
#' @format a RasterBrick with 144 columns x 99 rows x 4 layer
#' @source Images by J.Cunow
"lines_vertical"


#' @title Example Image with circa horizontal lines
#'
#'@description
#' Use this dataset to test directionality
#'
#' @author Johannes Cunow \email{johannes.cunow@gmail.com}
#' @usage data(lines_horizontal)
#' @docType data
#' @format a RasterBrick with 184 columns x 90 rows x 4 layer
#' @source Images by J.Cunow
"lines_horizontal"


#' @title Example Image with circa left to right diagonal lines
#'
#'@description
#' Use this dataset to test directionality
#'
#' @author Johannes Cunow \email{johannes.cunow@gmail.com}
#' @usage data(lines_diagonal_topleft)
#' @docType data
#' @format a RasterBrick with 194 columns x 100 rows  x 4 layer
#' @source Images by J.Cunow
"lines_diagonal_topleft"

#' @title Example Image with circa right to left diagonal lines
#'
#' @description
#' Use this dataset to test directionality
#'
#' @author Johannes Cunow \email{johannes.cunow@gmail.com}
#' @usage data(lines_diagonal_botleft)
#' @docType data
#' @format a RasterBrick with 191 columns x 99 rows x 4 layer
#' @source Images by J.Cunow
"lines_diagonal_botleft"


########## Root Scans ######




#' @title An original Minirhizotron Root Scan from Timepoint 2
#'
#' @description
#' This is a blended image of multiple scans stitched together.
#' The scan originates from a sedge fen in northern Finland in late October.
#' Columns correspond to tube length and rows to tube rotation
#'
#' @author Johannes Cunow \email{johannes.cunow@gmail.com}
#' @usage data(rgb_Oulanka2023_Session03_T067)
#' @docType data
#' @format a RasterBrick with 4900 columns x 1161 rows x 3 layer
#' @source Images by J.Cunow
"rgb_Oulanka2023_Session03_T067"


#' @title A Minirhizotron Root Scan after Segmentation from Timepoint 1
#'
#'@description
#' The image is derived from a rgb root scan after using "RootDetector" to segment.
#' Roots are represented as = 1, Background = 0.
#' The first layer also shows foreign objects such as tape as = 1.
#' The scan originates from a sedge fen in northern Finland in early June.
#' Columns correspond to tube length and rows to tube rotation
#'
#' @author Johannes Cunow \email{johannes.cunow@gmail.com}
#' @usage data(seg_Oulanka2023_Session01_T067)
#' @docType data
#' @format a RasterBrick with 4900 columns x 1144 rows x 3 layer
#' @source Images by J.Cunow
"seg_Oulanka2023_Session01_T067"



#' @title A Minirhizotron Root Scan after Segmentation from Timepoint 2
#'
#'@description
#' The image is derived from a rgb root scan after using "RootDetector" to segment.
#' Roots are represented as = 1, Background = 0.
#' The first layer also shows foreign objects such as tape as = 1.
#' The scan originates from a sedge fen in northern Finland in early June.
#' Columns correspond to tube length and rows to tube rotation
#'
#' @author Johannes Cunow \email{johannes.cunow@gmail.com}
#' @usage data(seg_Oulanka2023_Session03_T067)
#' @docType data
#' @format a RasterBrick with 4900 columns x 1161 rows x 3 layer
#' @source Images by J.Cunow
"seg_Oulanka2023_Session03_T067"


#' @title A Minirhizotron Root Scan after Segmentation and Skeletonization from Timepoint 1
#'
#'@description
#' The image is derived from a rgb root scan after using "RootDetector" to segment and skeletonize.
#' Roots are represented as = 1, Background = 0.
#' The first layer also shows foreign objects such as tape as = 1.
#' The scan originates from a sedge fen in northern Finland in early June.
#' Columns correspond to tube length and rows to tube rotation
#'
#' @author Johannes Cunow \email{johannes.cunow@gmail.com}
#' @usage data(skl_Oulanka2023_Session01_T067)
#' @docType data
#' @format a RasterBrick with 4900 columns x 1144 rows x 3 layer
#' @source Images by J.Cunow
"skl_Oulanka2023_Session01_T067"


#' @title A Minirhizotron Root Scan after Segmentation and Skeletonization from Timepoint 2
#'
#'@description
#' The image is derived from a rgb root scan after using "RootDetector" to segment and skeletonize.
#' Roots are represented as = 1, Background = 0.
#' The first layer also shows foreign objects such as tape as = 1.
#' The scan originates from a sedge fen in northern Finland in early June.
#' Columns correspond to tube length and rows to tube rotation
#'
#' @author Johannes Cunow \email{johannes.cunow@gmail.com}
#' @usage data(skl_Oulanka2023_Session03_T067)
#' @format a RasterBrick with 4900 columns x 1161 rows x 3 layer
#' @source Images by J.Cunow
"skl_Oulanka2023_Session03_T067"



#' @title Turnover Estimation of Individual Roots from "RootDetector"
#'
#'@description
#' The image is derived from a rgb root scan after using "RootDetector" to segment and compare against a second time point using the root tracking feature.
#'
#' Roots are represented as max value, Background = 0.
#' First layer shows root decay (disappeared between time points), foreign objects such as tape, and no-change roots
#' Second Layer shows new root growth and no-change roots
#' Third Layer shows no-change-roots
#' The scan originates from a sedge fen in northern Finland.
#' Columns correspond to tube length and rows to tube rotation
#'
#'
#'
#' @author Johannes Cunow \email{johannes.cunow@gmail.com}
#' @usage data(TurnoverDPC_data)
#' @format a RasterBrick with 2550 columns x 2273 rows x 3 layer
#' @source Images by J.Cunow
"TurnoverDPC_data"








