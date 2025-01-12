#### open issue: get away from raster and save data as other type

#' @title Example Image with Vertical Lines
#' @description
#' Test dataset containing an image with approximately vertical lines for
#' directionality analysis.
#'
#' @details
#' The image consists of multiple vertical lines suitable for testing
#' directionality detection algorithms.
#'
#' @author Johannes Cunow \email{johannes.cunow@gmail.com}
#' @format A RasterBrick object with dimensions:
#'   \itemize{
#'     \item 144 columns (width)
#'     \item 99 rows (height)
#'     \item 4 layers (channels)
#'   }
#' @usage data(lines_vertical)
#' @source Images by J.Cunow
#' @examples
#' \donttest{
#'   data(lines_vertical)
#'   plot(lines_vertical)
#' }
"lines_vertical"

#' @title Example Image with Horizontal Lines
#' @description
#' Test dataset containing an image with approximately horizontal lines for
#' directionality analysis.
#'
#' @details
#' The image consists of multiple horizontal lines suitable for testing
#' directionality detection algorithms.
#'
#' @author Johannes Cunow \email{johannes.cunow@gmail.com}
#' @format A RasterBrick object with dimensions:
#'   \itemize{
#'     \item 184 columns (width)
#'     \item 90 rows (height)
#'     \item 4 layers (channels)
#'   }
#' @usage data(lines_horizontal)
#' @source Images by J.Cunow
#' @examples
#' \donttest{
#'   data(lines_horizontal)
#'   plot(lines_horizontal)
#' }
"lines_horizontal"

#' @title Example Image with Diagonal Lines (Top-Left to Bottom-Right)
#' @description
#' Test dataset containing an image with diagonal lines running from
#' top-left to bottom-right for directionality analysis.
#'
#' @details
#' The image consists of multiple diagonal lines suitable for testing
#' directionality detection algorithms.
#'
#' @author Johannes Cunow \email{johannes.cunow@gmail.com}
#' @format A RasterBrick object with dimensions:
#'   \itemize{
#'     \item 194 columns (width)
#'     \item 100 rows (height)
#'     \item 4 layers (channels)
#'   }
#' @usage data(lines_diagonal_topleft)
#' @source Images by J.Cunow
#' @examples
#' \donttest{
#'   data(lines_diagonal_topleft)
#'   plot(lines_diagonal_topleft)
#' }
"lines_diagonal_topleft"

#' @title Example Image with Diagonal Lines (Bottom-Left to Top-Right)
#' @description
#' Test dataset containing an image with diagonal lines running from
#' bottom-left to top-right for directionality analysis.
#'
#' @details
#' The image consists of multiple diagonal lines suitable for testing
#' directionality detection algorithms.
#'
#' @author Johannes Cunow \email{johannes.cunow@gmail.com}
#' @format A RasterBrick object with dimensions:
#'   \itemize{
#'     \item 191 columns (width)
#'     \item 99 rows (height)
#'     \item 4 layers (channels)
#'   }
#' @usage data(lines_diagonal_botleft)
#' @source Images by J.Cunow
#' @examples
#' \donttest{
#'   data(lines_diagonal_botleft)
#'   plot(lines_diagonal_botleft)
#' }
"lines_diagonal_botleft"

#' @title Original Minirhizotron Root Scan - Session 3, Tube 67
#' @description
#' Original RGB root scan image from a sedge fen in northern Finland (October 2023).
#' This image represents a composite of multiple stitched scans.
#'
#' @details
#' The image represents a complete tube scan where:
#' \itemize{
#'   \item Columns correspond to tube length
#'   \item Rows correspond to tube rotation
#'   \item RGB channels represent true color information
#' }
#'
#' @author Johannes Cunow \email{johannes.cunow@gmail.com}
#' @format A RasterBrick object with dimensions:
#'   \itemize{
#'     \item 4900 columns (width)
#'     \item 1161 rows (height)
#'     \item 3 layers (RGB channels)
#'   }
#' @usage data(rgb_Oulanka2023_Session03_T067)
#' @source Images by J.Cunow
#' @examples
#' \donttest{
#'   data(rgb_Oulanka2023_Session03_T067)
#'   plotRGB(rgb_Oulanka2023_Session03_T067)
#' }
"rgb_Oulanka2023_Session03_T067"

#' @title Segmented Minirhizotron Root Scan - Session 1, Tube 67
#' @description
#' Segmented root scan image from a sedge fen in northern Finland (June 2023).
#' The image was processed using RootDetector for root segmentation.
#'
#' @details
#' Binary mask representation where:
#' \itemize{
#'   \item Roots = 1
#'   \item Background = 0
#'   \item Layer 1 includes foreign objects (e.g., tape) marked as 1
#' }
#' Spatial dimensions correspond to physical tube measurements:
#' columns = tube length, rows = tube rotation.
#'
#' @author Johannes Cunow \email{johannes.cunow@gmail.com}
#' @format A RasterBrick object with dimensions:
#'   \itemize{
#'     \item 4900 columns (width)
#'     \item 1144 rows (height)
#'     \item 3 layers (channels)
#'   }
#' @usage data(seg_Oulanka2023_Session01_T067)
#' @source Images by J.Cunow
#' @examples
#' \donttest{
#'   data(seg_Oulanka2023_Session01_T067)
#'   plot(seg_Oulanka2023_Session01_T067)
#' }
"seg_Oulanka2023_Session01_T067"

#' @title Segmented Minirhizotron Root Scan - Session 3, Tube 67
#' @description
#' Segmented root scan image from a sedge fen in northern Finland (October 2023).
#' The image was processed using RootDetector for root segmentation.
#'
#' @details
#' Binary mask representation where:
#' \itemize{
#'   \item Roots = 1
#'   \item Background = 0
#'   \item Layer 1 includes foreign objects (e.g., tape) marked as 1
#' }
#' Spatial dimensions correspond to physical tube measurements:
#' columns = tube length, rows = tube rotation.
#'
#' @author Johannes Cunow \email{johannes.cunow@gmail.com}
#' @format A RasterBrick object with dimensions:
#'   \itemize{
#'     \item 4900 columns (width)
#'     \item 1161 rows (height)
#'     \item 3 layers (channels)
#'   }
#' @usage data(seg_Oulanka2023_Session03_T067)
#' @source Images by J.Cunow
#' @examples
#' \donttest{
#'   data(seg_Oulanka2023_Session03_T067)
#'   plot(seg_Oulanka2023_Session03_T067)
#' }
"seg_Oulanka2023_Session03_T067"

#' @title Skeletonized Root Scan - Session 1, Tube 67
#' @description
#' Skeletonized root scan image from a sedge fen in northern Finland (June 2023).
#' The image was processed using RootDetector for segmentation and skeletonization.
#'
#' @details
#' Binary mask representation where:
#' \itemize{
#'   \item Root skeletons = 1
#'   \item Background = 0
#'   \item Layer 1 includes foreign objects (e.g., tape) marked as 1
#' }
#' Skeletonization reduces root width to single-pixel lines while preserving
#' the root system topology.
#'
#' @author Johannes Cunow \email{johannes.cunow@gmail.com}
#' @format A RasterBrick object with dimensions:
#'   \itemize{
#'     \item 4900 columns (width)
#'     \item 1144 rows (height)
#'     \item 3 layers (channels)
#'   }
#' @usage data(skl_Oulanka2023_Session01_T067)
#' @source Images by J.Cunow
#' @examples
#' \donttest{
#'   data(skl_Oulanka2023_Session01_T067)
#'   plot(skl_Oulanka2023_Session01_T067)
#' }
"skl_Oulanka2023_Session01_T067"

#' @title Skeletonized Root Scan - Session 3, Tube 67
#' @description
#' Skeletonized root scan image from a sedge fen in northern Finland (October 2023).
#' The image was processed using RootDetector for segmentation and skeletonization.
#'
#' @details
#' Binary mask representation where:
#' \itemize{
#'   \item Root skeletons = 1
#'   \item Background = 0
#'   \item Layer 1 includes foreign objects (e.g., tape) marked as 1
#' }
#' Skeletonization reduces root width to single-pixel lines while preserving
#' the root system topology.
#'
#' @author Johannes Cunow \email{johannes.cunow@gmail.com}
#' @format A RasterBrick object with dimensions:
#'   \itemize{
#'     \item 4900 columns (width)
#'     \item 1161 rows (height)
#'     \item 3 layers (channels)
#'   }
#' @usage data(skl_Oulanka2023_Session03_T067)
#' @source Images by J.Cunow
#' @examples
#' \donttest{
#'   data(skl_Oulanka2023_Session03_T067)
#'   plot(skl_Oulanka2023_Session03_T067)
#' }
"skl_Oulanka2023_Session03_T067"

#' @title Root Turnover Analysis Data
#' @description
#' Root turnover analysis from a sedge fen in northern Finland, comparing
#' root presence between two time points (June 2023 and October 2023) using
#' RootDetector root tracking feature.
#'
#' @details
#' Multi-layer representation of root dynamics where:
#' \itemize{
#'   \item Layer 1: Root decay (value = 1)
#'     \itemize{
#'       \item Roots that disappeared between time points
#'       \item Foreign objects (e.g., tape)
#'       \item Persistent roots
#'     }
#'   \item Layer 2: Root growth (value = 1)
#'     \itemize{
#'       \item New roots that appeared
#'       \item Persistent roots
#'     }
#'   \item Layer 3: Persistent roots (value = 1)
#'     \itemize{
#'       \item Roots present in both time points
#'     }
#' }
#' Background is represented as 0 in all layers.
#'
#' @author Johannes Cunow \email{johannes.cunow@gmail.com}
#' @format A RasterBrick object with dimensions:
#'   \itemize{
#'     \item 2550 columns (width)
#'     \item 2273 rows (height)
#'     \item 3 layers (decay, growth, persistent)
#'   }
#' @usage data(TurnoverDPC_data)
#' @source Images by J.Cunow
#' @examples
#' \donttest{
#'   data(TurnoverDPC_data)
#'   img = terra::rast(TurnoverDPC_data)
#'   # Plot individual layers
#'   terra::plot(img)
#'   # Calculate turnover statistics
#'   decay <- sum(terra::values(img[[1]]) == 1)
#'   growth <- sum(terra::values(img[[2]]) == 1)
#'   persistent <- sum(terra::values(img[[3]]) == 1)
#' }
"TurnoverDPC_data"
