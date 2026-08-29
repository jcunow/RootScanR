# ============================================================
# VALIDATION AGAINST SYNTHETIC ROOT IMAGES WITH KNOWN PROPERTIES
# ------------------------------------------------------------
# The branching pipeline reports length, diameter and topology, but on a real
# scan none of those are known, so a wrong answer looks like a plausible one.
# The functions here draw root systems whose geometry is prescribed -- exact
# axis length, stroke width, tip count, branch-point count and branch order --
# and score the pipeline's output against that ground truth.
# ============================================================


#' Rasterise a line between two pixels (Bresenham)
#'
#' @param p0,p1 Integer length-2 (row, col) endpoints.
#' @return An Nx2 integer matrix of (row, col) pixels, 8-connected and 1 px wide.
#' @keywords internal
#' @noRd
.bresenham <- function(p0, p1) {
  r0 <- as.integer(p0[1]); c0 <- as.integer(p0[2])
  r1 <- as.integer(p1[1]); c1 <- as.integer(p1[2])
  dr <- abs(r1 - r0); dc <- abs(c1 - c0)
  sr <- sign(r1 - r0); sc <- sign(c1 - c0)
  err <- dr - dc
  out <- matrix(NA_integer_, dr + dc + 2L, 2L); k <- 0L
  repeat {
    k <- k + 1L; out[k, ] <- c(r0, c0)
    if (r0 == r1 && c0 == c1) break
    e2 <- 2L * err
    if (e2 > -dc) { err <- err - dc; r0 <- r0 + sr }
    if (e2 <  dr) { err <- err + dr; c0 <- c0 + sc }
  }
  out[seq_len(k), , drop = FALSE]
}

#' Pixels of a polyline
#' @keywords internal
#' @noRd
.polyline_px <- function(v) {
  if (nrow(v) < 2L) return(v)
  do.call(rbind, lapply(seq_len(nrow(v) - 1L),
                        function(i) .bresenham(v[i, ], v[i + 1L, ])))
}

#' Exact Euclidean length of a polyline
#' @keywords internal
#' @noRd
.polyline_length <- function(v) {
  if (nrow(v) < 2L) return(0)
  sum(sqrt(rowSums(diff(v)^2)))
}

#' Stamp a filled disk of radius r at every pixel of a polyline
#' @keywords internal
#' @noRd
.stamp_stroke <- function(m, px, radius) {
  rr <- max(0L, as.integer(round(radius)))
  g   <- expand.grid(dr = -rr:rr, dc = -rr:rr)
  off <- as.matrix(g[g$dr^2 + g$dc^2 <= radius^2, , drop = FALSE])
  if (nrow(off) == 0L) off <- matrix(0L, 1L, 2L)
  for (k in seq_len(nrow(px))) {
    R <- px[k, 1] + off[, 1]; C <- px[k, 2] + off[, 2]
    ok <- R >= 1L & R <= nrow(m) & C >= 1L & C <= ncol(m)
    m[cbind(R[ok], C[ok])] <- 1L
  }
  m
}


#' Stroke definitions for the built-in phantom designs
#'
#' Each stroke is \code{list(v, radius, order, parent)} where \code{v} is an Nx2
#' matrix of (row, col) vertices and \code{parent} indexes the stroke this one
#' attaches to (\code{NA} for a free-standing root). A child's first vertex is
#' the attachment point and sits exactly on the parent's centre line, so the
#' ground-truth length of the union is the sum of the stroke lengths.
#'
#' Stamping radius \code{r} yields a stroke \code{2r + 1} px wide.
#'
#' @keywords internal
#' @noRd
.phantom_strokes <- function(design, size) {
  n <- size; mid <- round(n / 2)
  lo <- round(0.05 * n); hi <- round(0.95 * n)
  switch(
    design,
    comb = {
      # one straight main axis, five perpendicular laterals attached at interior
      # points of it. Everything is axis-aligned, so the chain-code length of the
      # drawn skeleton equals its Euclidean length exactly.
      s <- list(list(v = rbind(c(lo, mid), c(hi, mid)), radius = 4, order = 1L, parent = NA_integer_))
      for (r in round(seq(0.2, 0.8, length.out = 5) * n))
        s[[length(s) + 1L]] <- list(v = rbind(c(r, mid), c(r, round(0.85 * n))),
                                    radius = 2, order = 2L, parent = 1L)
      s
    },
    herringbone = {
      # 45-degree laterals: still exact under the sqrt(2) chain code, but they
      # exercise diagonal tracing and the continuation rule at oblique junctions.
      s <- list(list(v = rbind(c(lo, mid), c(hi, mid)), radius = 4, order = 1L, parent = NA_integer_))
      at <- round(seq(0.25, 0.75, length.out = 4) * n); d <- round(0.18 * n)
      for (i in seq_along(at)) {
        side <- if (i %% 2L == 1L) 1L else -1L
        s[[length(s) + 1L]] <- list(v = rbind(c(at[i], mid), c(at[i] + d, mid + side * d)),
                                    radius = 2, order = 2L, parent = 1L)
      }
      s
    },
    hierarchical = {
      # three generations, every attachment a T on an interior point, so the
      # continuation rule has one unambiguous answer at every junction.
      s <- list(list(v = rbind(c(lo, mid), c(hi, mid)), radius = 4, order = 1L, parent = NA_integer_))
      at  <- round(seq(0.25, 0.7, length.out = 3) * n)
      len <- round(0.22 * n); sub <- round(0.10 * n)
      for (i in seq_along(at)) {
        s[[length(s) + 1L]] <- list(v = rbind(c(at[i], mid), c(at[i], mid + len)),
                                    radius = 2, order = 2L, parent = 1L)
        pid <- length(s)
        # distinct attachment points: two sub-laterals leaving one point in
        # opposite directions would be an X in outline, not a pair of T's
        for (k in seq_along(c(-1L, 1L))) {
          side <- c(-1L, 1L)[k]
          col  <- mid + round(len * c(0.4, 0.7)[k])
          s[[length(s) + 1L]] <- list(
            v = rbind(c(at[i], col), c(at[i] + side * sub, col)),
            radius = 1, order = 3L, parent = pid)
        }
      }
      s
    },
    cross = {
      # two independent roots that merely overlap in projection: the pipeline
      # must recover TWO roots and ZERO branch points, not an X-shaped fork.
      list(list(v = rbind(c(lo, mid), c(hi, mid)), radius = 3, order = 1L, parent = NA_integer_),
           list(v = rbind(c(mid, lo), c(mid, hi)), radius = 3, order = 1L, parent = NA_integer_))
    },
    fork = {
      # a symmetric dichotomous Y. Length and topology are well defined, but
      # "which arm continues the parent" is not, so root_phantom() marks this
      # design as having no ground truth for the per-root metrics.
      d <- round(0.30 * n)
      f <- c(round(0.40 * n), mid)
      s <- list(list(v = rbind(c(lo, mid), f), radius = 4, order = 1L, parent = NA_integer_))
      for (side in c(-1L, 1L))
        s[[length(s) + 1L]] <- list(v = rbind(f, c(f[1] + d, mid + side * d)),
                                    radius = 4, order = 2L, parent = 1L)
      s
    },
    stop("Unknown design: ", design)
  )
}


#' Ground-truth topology implied by a set of strokes
#'
#' A child attaches at its own first vertex. That point is a branch junction; the
#' parent end that coincides with it, if any, stops being a free end. Counting
#' this way is design-agnostic: it is right for a T-attachment and for a
#' symmetric fork alike.
#'
#' @keywords internal
#' @noRd
.phantom_topology <- function(strokes) {
  key <- function(p) paste(round(p[1]), round(p[2]), sep = "_")
  attach <- unique(vapply(strokes, function(s)
    if (is.na(s$parent)) NA_character_ else key(s$v[1L, ]), character(1)))
  attach <- attach[!is.na(attach)]
  tips <- 0L
  for (s in strokes) {
    ends <- if (is.na(s$parent)) list(s$v[1L, ], s$v[nrow(s$v), ]) else list(s$v[nrow(s$v), ])
    for (e in ends) if (!key(e) %in% attach) tips <- tips + 1L
  }
  list(n_tips = tips, n_branch_points = length(attach))
}


#' Synthetic root image with known geometry
#'
#' Draws a root system whose length, width and topology are prescribed, so the
#' branching pipeline can be scored against ground truth. Returns both the
#' filled mask (as a scanner would see it) and the exact one-pixel centre line,
#' which lets graph errors be told apart from thinning errors.
#'
#' @details
#' Every design is built from strokes -- polylines stamped with a disk of radius
#' \code{r}, giving a root \code{2r + 1} px wide -- and each lateral starts
#' exactly on its parent's centre line. The ground truth is therefore analytic:
#' \describe{
#'   \item{\code{total_length}}{Sum of the stroke centre-line lengths (px).}
#'   \item{\code{n_tips}, \code{n_branch_points}}{Counted from the stroke
#'     attachment graph, so they hold for T-attachments and forks alike.}
#'   \item{\code{n_roots}, \code{max_branch_order}, \code{by_order}}{Per-root
#'     truth, defined only where the continuation rule is unambiguous (see
#'     \code{continuation_defined}).}
#' }
#' Two conventions are worth stating, because the truth is written to match the
#' package and not the other way round:
#' \itemize{
#'   \item \strong{Diameter.} Rootopia reports \code{2 * EDT}, the inscribed
#'     diameter measured to the nearest \emph{background} pixel, exactly as
#'     \code{\link{root_diameter}} does. For a stroke \code{W} px wide that is
#'     \code{W + 1}, so \code{by_order$mean_diameter} carries the \code{+1}.
#'   \item \strong{Length.} The designs run along the axes or at 45 degrees, so
#'     the drawn skeleton's chain-code length equals its Euclidean length. At
#'     other angles the pipeline would \emph{correctly} report the well-known
#'     chain-code overestimate, which is a property of the measure rather than a
#'     defect.
#' }
#' \code{"fork"} is a symmetric dichotomous Y. Which arm continues the parent is
#' genuinely undefined there, so \code{truth$continuation_defined} is
#' \code{FALSE} and \code{\link{validate_branching}} scores only length and
#' topology for it.
#'
#' @param design One of \code{"comb"} (main axis + 5 perpendicular laterals),
#'   \code{"herringbone"} (45-degree laterals), \code{"hierarchical"} (three
#'   generations of T-attachments), \code{"cross"} (two roots overlapping in an
#'   X) or \code{"fork"} (symmetric dichotomous Y).
#' @param size Image side length in pixels (square image).
#' @param dpi Nominal scan resolution recorded in \code{$truth$dpi}; used to
#'   express the ground truth in cm as well as pixels.
#' @param as \code{"matrix"} (default) or \code{"spatraster"} for the returned
#'   images.
#' @return A list with \code{$mask} (filled binary root image), \code{$skeleton}
#'   (exact 1-px centre line), \code{$truth} (ground-truth list, including
#'   \code{$by_order}) and \code{$strokes} (the stroke definitions).
#' @seealso \code{\link{validate_branching}}, \code{\link{branch_order_map}}
#' @examples
#' ph <- root_phantom("comb", size = 200)
#' ph$truth$n_tips
#' ph$truth$total_length
#' @export
root_phantom <- function(design = c("comb", "herringbone", "hierarchical", "cross", "fork"),
                         size = 400, dpi = 300, as = c("matrix", "spatraster")) {
  design <- match.arg(design); as <- match.arg(as)
  if (size < 60) stop("'size' must be at least 60 px.")
  strokes <- .phantom_strokes(design, size)

  mask <- matrix(0L, size, size)
  skel <- matrix(0L, size, size)
  for (s in strokes) {
    px <- .polyline_px(s$v)
    ok <- px[, 1] >= 1L & px[, 1] <= size & px[, 2] >= 1L & px[, 2] <= size
    px <- px[ok, , drop = FALSE]
    skel[px] <- 1L
    mask <- .stamp_stroke(mask, px, s$radius)
  }

  ord <- vapply(strokes, function(s) s$order, integer(1))
  len <- vapply(strokes, function(s) .polyline_length(s$v), numeric(1))
  # Diameter truth in the package's own 2 * EDT convention, measured on the
  # phantom we just drew. For an axis-aligned stroke 2r + 1 px wide this is
  # exactly 2r + 2; for a 45-degree stroke the disk stamp gives slightly less,
  # and near a tip the rounded cap gives less again -- all of which a fixed
  # formula would get wrong while the phantom's own EDT gets right.
  DT0 <- .distance_transform_edt(mask, backend = "baseR")
  dia <- vapply(strokes, function(s) {
    px <- .polyline_px(s$v)
    ok <- px[, 1] >= 1L & px[, 1] <= size & px[, 2] >= 1L & px[, 2] <= size
    2 * mean(DT0[px[ok, , drop = FALSE]])
  }, numeric(1))
  topo <- .phantom_topology(strokes)
  lv <- sort(unique(ord)); ch <- as.character(lv)

  truth <- list(
    design                = design,
    dpi                   = dpi,
    continuation_defined  = design != "fork",
    total_length          = sum(len),
    total_length_cm       = sum(len) * 2.54 / dpi,
    n_tips                = as.integer(topo$n_tips),
    n_branch_points       = as.integer(topo$n_branch_points),
    n_roots               = length(strokes),
    max_branch_order      = max(ord),
    by_order = data.frame(
      order         = lv,
      n_roots       = as.integer(tapply(len, ord, length)[ch]),
      total_length  = as.numeric(tapply(len, ord, sum)[ch]),
      mean_diameter = as.numeric(tapply(dia * len, ord, sum)[ch] /
                                   tapply(len, ord, sum)[ch]),
      stringsAsFactors = FALSE)
  )

  if (as == "spatraster") {
    if (!requireNamespace("terra", quietly = TRUE))
      stop("Package 'terra' is required for as = 'spatraster'.")
    mask <- terra::rast(mask); skel <- terra::rast(skel)
  }
  list(mask = mask, skeleton = skel, truth = truth, strokes = strokes)
}


#' Default tolerances for validate_branching()
#' @keywords internal
#' @noRd
.branch_tolerance <- function(from) {
  # Counts must be exact. Lengths carry an unavoidable negative bias: junction
  # contraction dissolves a couple of pixels per branch point, and thinning
  # additionally erodes each root tip by about one radius, so the mask route
  # gets more room than the skeleton route. Diameter is scored against the
  # 2 * EDT convention the package uses, so only sampling effects remain --
  # a lateral is fatter where it enters its parent and thinner at its tip.
  if (from == "skeleton")
    list(count = 0, length = 0.03, diameter = 0.10)
  else
    list(count = 0, length = 0.06, diameter = 0.10)
}


#' Score the branching pipeline against a phantom with known properties
#'
#' Runs \code{\link{branch_order_map}} on a synthetic root image whose geometry
#' is known exactly and reports, metric by metric, what was expected, what was
#' measured and whether the difference is within tolerance. This is the check to
#' run after touching the tracing, ordering or length code, and the one to cite
#' when reporting what the package's numbers mean.
#'
#' @details
#' \code{from = "skeleton"} feeds the exact one-pixel centre line, so the score
#' isolates the graph: tracing, junction contraction, crossing resolution,
#' ordering and length integration. \code{from = "mask"} skeletonises the filled
#' image first and therefore also carries the thinning error -- chiefly the
#' erosion of about one root radius at every tip, which shortens the total by a
#' few percent. Both are legitimate; they answer different questions, and the
#' default tolerances differ accordingly.
#'
#' Counts (\code{n_tips}, \code{n_branch_points}, \code{n_roots},
#' \code{max_branch_order}, \code{n_unordered}) are required to be exact.
#' Lengths and diameters are scored as relative error.
#'
#' @param design Phantom design, passed to \code{\link{root_phantom}}; ignored
#'   when \code{phantom} is supplied. For \code{"fork"} the per-root metrics
#'   (\code{n_roots}, \code{max_branch_order}, per-order length and diameter)
#'   have no ground truth and are omitted from the report.
#' @param phantom A phantom from \code{\link{root_phantom}} (or any list with
#'   \code{$mask}, \code{$skeleton} and \code{$truth}).
#' @param from \code{"skeleton"} to score the graph alone, \code{"mask"} to
#'   score the whole pipeline including skeletonisation.
#' @param size,dpi Passed to \code{\link{root_phantom}} when building a phantom.
#' @param tolerance Named list with \code{count}, \code{length} and
#'   \code{diameter} relative tolerances; defaults depend on \code{from}.
#' @param overlay_png Optional path for the order-coloured validation image.
#' @param verbose Print the report.
#' @param ... Passed to \code{\link{branch_order_map}}.
#' @return A data.frame with one row per metric (\code{metric}, \code{expected},
#'   \code{observed}, \code{abs_error}, \code{rel_error}, \code{tolerance},
#'   \code{pass}). \code{attr(., "passed")} is \code{TRUE} when every row passes;
#'   \code{attr(., "result")} holds the \code{branchOrderMap} object.
#' @seealso \code{\link{root_phantom}}, \code{\link{branch_order_map}}
#' @examples
#' \donttest{
#' v <- validate_branching("comb", size = 200, verbose = FALSE)
#' v[, c("metric", "expected", "observed", "pass")]
#' attr(v, "passed")
#' }
#' @export
validate_branching <- function(design = c("comb", "herringbone", "hierarchical", "cross", "fork"),
                               phantom = NULL,
                               from = c("skeleton", "mask"),
                               size = 400, dpi = 300,
                               tolerance = NULL, overlay_png = NULL,
                               verbose = TRUE, ...) {
  from <- match.arg(from)
  if (is.null(phantom)) {
    design  <- match.arg(design)
    phantom <- root_phantom(design, size = size, dpi = dpi)
  }
  tol <- .branch_tolerance(from)
  if (!is.null(tolerance)) tol[names(tolerance)] <- tolerance
  truth <- phantom$truth

  skel <- if (from == "skeleton") phantom$skeleton else
    .to_binary_matrix(skeletonize_image(phantom$mask, verbose = FALSE))

  res <- branch_order_map(skel = skel, mask = phantom$mask, order = "branch_order",
                          unit = "px", dpi = truth$dpi, return_map = FALSE,
                          overlay_png = overlay_png, verbose = FALSE, ...)
  et <- res$edges

  obs_orders <- stats::na.omit(et$branch_order)
  rows <- list(
    .vrow("total_length",    truth$total_length,    sum(et$length),          tol$length),
    .vrow("n_tips",          truth$n_tips,          sum(et$n_tips),          tol$count),
    .vrow("n_branch_points", truth$n_branch_points, sum(et$n_branch_points), tol$count),
    .vrow("n_unordered",     0,                     sum(is.na(et$tip_order)), tol$count)
  )

  if (isTRUE(truth$continuation_defined)) {
    rows <- c(rows, list(
      .vrow("n_roots",          truth$n_roots,          length(unique(et$root_id)), tol$count),
      .vrow("max_branch_order", truth$max_branch_order,
            if (length(obs_orders)) max(obs_orders) else NA_real_,                  tol$count)))
    # per-order length and diameter, read off the pipeline's own summary
    smry <- res$summary
    for (k in seq_len(nrow(truth$by_order))) {
      o <- truth$by_order$order[k]
      hit <- smry[smry$order == o, , drop = FALSE]
      rows[[length(rows) + 1L]] <- .vrow(
        sprintf("order%d_total_length", o), truth$by_order$total_length[k],
        if (nrow(hit)) hit$total_length[1] else NA_real_, tol$length)
      rows[[length(rows) + 1L]] <- .vrow(
        sprintf("order%d_mean_diameter", o), truth$by_order$mean_diameter[k],
        if (nrow(hit)) hit$mean_diameter[1] else NA_real_, tol$diameter)
    }
  }

  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  attr(out, "passed") <- all(out$pass)
  attr(out, "design") <- truth$design
  attr(out, "from")   <- from
  attr(out, "result") <- res

  if (verbose) {
    cat(sprintf("\n== validate_branching: design '%s', from '%s' ==\n", truth$design, from))
    print(within(out, {
      expected  <- round(expected, 2); observed <- round(observed, 2)
      rel_error <- ifelse(is.na(rel_error), NA, sprintf("%+.1f%%", 100 * rel_error))
    }))
    cat(sprintf("-> %s (%d/%d metrics within tolerance)\n\n",
                if (all(out$pass)) "PASS" else "FAIL", sum(out$pass), nrow(out)))
  }
  out
}

#' One row of the validation report
#' @keywords internal
#' @noRd
.vrow <- function(metric, expected, observed, tol) {
  expected <- as.numeric(expected); observed <- as.numeric(observed)
  abs_err  <- observed - expected
  rel_err  <- if (isTRUE(expected != 0)) abs_err / expected else NA_real_
  pass <- if (is.na(observed)) FALSE
          else if (tol == 0) isTRUE(abs_err == 0)
          else isTRUE(abs(rel_err) <= tol)
  data.frame(metric = metric, expected = expected, observed = observed,
             abs_error = abs_err, rel_error = rel_err, tolerance = tol,
             pass = pass, stringsAsFactors = FALSE)
}
