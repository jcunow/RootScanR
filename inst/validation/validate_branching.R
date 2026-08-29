#!/usr/bin/env Rscript
# ============================================================================
# Rootopia -- public validation of the branching pipeline
# ----------------------------------------------------------------------------
# Scores branch_order_map() against synthetic root images whose length, width
# and topology are known exactly (see ?root_phantom). Nothing here depends on
# field data, so anyone can reproduce the numbers the package reports.
#
#   Rscript inst/validation/validate_branching.R [--size 400] [--out DIR] [--quick]
#
# Exits 0 when every metric of every design is within tolerance, 1 otherwise,
# so it can be used as a CI gate.
# ============================================================================

suppressPackageStartupMessages(library(Rootopia))

args <- commandArgs(trailingOnly = TRUE)
arg  <- function(flag, default) {
  i <- match(flag, args); if (is.na(i) || i == length(args)) default else args[i + 1L]
}
size  <- as.integer(arg("--size", "400"))
outdir <- arg("--out", NA_character_)
quick <- "--quick" %in% args

designs <- c("comb", "herringbone", "hierarchical", "cross", "fork")
routes  <- if (quick) "skeleton" else c("skeleton", "mask")

if (!is.na(outdir)) dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

cat(strrep("=", 78), "\n")
cat("Rootopia branching validation --  image size ", size, " px\n", sep = "")
cat("  route 'skeleton': exact 1-px centre line in  -> scores the graph alone\n")
cat("  route 'mask'    : filled roots in, thinned here -> scores the whole pipeline\n")
cat(strrep("=", 78), "\n")

all_rows <- list()
for (route in routes) {
  for (d in designs) {
    png <- if (is.na(outdir)) NULL else file.path(outdir, sprintf("%s_%s.png", d, route))
    v <- validate_branching(d, from = route, size = size,
                            overlay_png = png, verbose = TRUE)
    v$design <- d; v$route <- route
    all_rows[[length(all_rows) + 1L]] <- v[, c("design", "route", "metric", "expected",
                                               "observed", "rel_error", "tolerance", "pass")]
  }
}

report <- do.call(rbind, all_rows)
rownames(report) <- NULL

cat(strrep("=", 78), "\n")
cat("SUMMARY\n")
by_case <- aggregate(pass ~ design + route, report, function(x) sprintf("%d/%d", sum(x), length(x)))
names(by_case)[3] <- "within_tolerance"
print(by_case, row.names = FALSE)

failed <- report[!report$pass, , drop = FALSE]
if (nrow(failed)) {
  cat("\nFAILING METRICS\n")
  print(failed, row.names = FALSE)
}
cat(sprintf("\n%s -- %d of %d metrics within tolerance\n",
            if (nrow(failed) == 0L) "ALL DESIGNS PASS" else "VALIDATION FAILED",
            sum(report$pass), nrow(report)))

if (!is.na(outdir)) {
  csv <- file.path(outdir, "validation_report.csv")
  utils::write.csv(report, csv, row.names = FALSE)
  cat("Report written to ", csv, "\n", sep = "")
}

quit(status = if (nrow(failed) == 0L) 0L else 1L)
