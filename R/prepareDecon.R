
#' Prepare expression and background matrices for SpatialDecon
#'
#' @param spe SpatialExperiment object.
#' @param assay2use The assay to use when `recomputeCPM = FALSE`. Defaults to
#'   `"logcounts"`, preserving any normalization or batch correction in that assay.
#' @param negProbeName The name of the negative probe gene. By default is NegProbe-WTX.
#' @param pool A vector of probe pools aligned to the rows of `spe`. Required
#'   when different probe pools need separate background estimates. By default,
#'   all negative probes are averaged within each ROI.
#' @param recomputeCPM Logical. If TRUE, recompute linear CPM from the `"counts"`
#'   assay, using the same library sizes for expression and negative probes.
#'   This replaces any normalization or batch correction in `assay2use`.
#'   Defaults to FALSE.
#' @param isLog Logical. Whether the selected assay contains log2 values, which
#'   are converted to the linear scale using `2^x`. Defaults to TRUE for
#'   `assay2use = "logcounts"` and FALSE otherwise. Set to FALSE if the assay
#'   already contains linear values, including output from `geomxNorm(log = FALSE)`.
#'   Ignored when `recomputeCPM = TRUE`.
#'
#' @return A list containing two linear-scale matrices: `normCount` (expression)
#'   and `backGround` (expected background), with matching dimensions.
#'
#' @details SpatialDecon requires expression and background on the linear scale.
#'   By default, negative probes must be retained in the selected assay so they
#'   receive the same transformation as the expression data. Raw negative probes
#'   stored in `metadata(spe)$NegProbes` can be used with `recomputeCPM = TRUE`,
#'   or with `assay2use = "counts", isLog = FALSE`. They cannot reconstruct the
#'   background of an arbitrary normalized or batch-corrected assay.
#'
#'   Converting log2 values with `2^x` retains any prior-count contribution from
#'   log-CPM normalization; it does not recover the original raw counts.
#'   Preparing batch-corrected input does not establish that correction preserves
#'   cell-type signals or that corrected negative probes adequately model background.
#' @export
#'
#' @examples
#' countFile <- data.frame(
#'   TargetName = c("GeneA", "GeneB", "NegProbe-WTX"),
#'   ROI1 = c(10, 20, 1),
#'   ROI2 = c(30, 40, 2)
#' )
#' sampleAnnoFile <- data.frame(
#'   SegmentDisplayName = c("ROI1", "ROI2"),
#'   ROICoordinateX = c(1, 2),
#'   ROICoordinateY = c(3, 4)
#' )
#' 
#' spe <- readGeoMx(countFile, sampleAnnoFile, rmNegProbe = FALSE)
#' 
#' # Use the existing logcounts assay, converting expression and background to linear values.
#' out <- prepareSpatialDecon(spe)
#' # Alternatively, recompute linear CPM from raw counts.
#' out_cpm <- prepareSpatialDecon(spe, recomputeCPM = TRUE)
#' 
prepareSpatialDecon <- function(spe, assay2use = "logcounts", negProbeName = "NegProbe-WTX",
                                pool = NA, recomputeCPM = FALSE,
                                isLog = identical(assay2use, "logcounts")) {

  stopifnot(is.logical(recomputeCPM), length(recomputeCPM) == 1L, !is.na(recomputeCPM))
  if (recomputeCPM) {
    assay2use <- "counts"
    isLog <- FALSE
  }
  stopifnot(is.logical(isLog), length(isLog) == 1L, !is.na(isLog))

  norm <- as.matrix(assay(spe, assay2use))
  neg_norm <- .get_neg_probe_values(spe, norm, assay2use, negProbeName, isLog)

  if (is.null(neg_norm)){
    stop(
      paste0(
        paste(negProbeName, collapse = ", "),
        " must be included in the assay rows or stored in metadata(spe)$NegProbes. ",
        "If using readGeoMx(), set rmNegProbe=FALSE to keep negative probes in assays."
      )
    )
  }

  if (recomputeCPM) {
    # Include metadata-only negative probes once when calculating library sizes.
    extra_neg <- neg_norm[!rownames(neg_norm) %in% rownames(norm), , drop = FALSE]
    lib_size <- colSums(rbind(norm, extra_neg))
    norm <- edgeR::cpm(norm, lib.size = lib_size, log = FALSE)
    neg_norm <- edgeR::cpm(neg_norm, lib.size = lib_size, log = FALSE)
  } else if (isLog) {
    # SpatialDecon fits linear mixtures, even when the input assay is log2 scaled.
    norm <- 2^norm
    neg_norm <- 2^neg_norm
  }

  if (all(is.na(pool))){
    mean_bg <- colMeans(neg_norm, na.rm = TRUE)
    bg <- sweep(norm*0, 2, mean_bg, "+")
  } else {
    if (nrow(norm) != length(pool)){
      stop("length(pool) should be equal to nrow(spe).")
    }

    bg <- norm*0

    for (p in unique(pool)){
      genes_p <- rownames(norm)[pool == p]
      negp_p <- base::intersect(rownames(neg_norm), genes_p)

      stopifnot("Negative probes not found in one of the pool." = length(negp_p) > 0)

      mean_bg_p <- colMeans(neg_norm[negp_p, , drop = FALSE])

      bg[pool == p, ] <- sweep(bg[pool == p, ], 2, mean_bg_p, "+")
    }
  }

  out <- list(norm, bg)

  names(out) <- c("normCount","backGround")

  return(out)
}

.get_neg_probe_values <- function(spe, norm, assay2use, negProbeName, isLog) {
  neg_rows <- .match_neg_probe_rows(rownames(norm), negProbeName)
  if (length(neg_rows) > 0) {
    return(norm[neg_rows, , drop = FALSE])
  }

  neg_counts <- S4Vectors::metadata(spe)$NegProbes
  if (is.null(neg_counts)) {
    return(NULL)
  }

  neg_rows <- .match_neg_probe_rows(rownames(neg_counts), negProbeName)
  if (length(neg_rows) == 0) {
    return(NULL)
  }
  neg_counts <- as.matrix(neg_counts[neg_rows, , drop = FALSE])
  neg_counts <- neg_counts[, colnames(norm), drop = FALSE]

  if (identical(assay2use, "counts") && !isLog) {
    return(neg_counts)
  }

  stop(
    "Negative probes stored in metadata(spe)$NegProbes are raw counts; ",
    "keep negative probes in the selected assay to preserve its normalization or batch correction, ",
    "or set recomputeCPM=TRUE to recompute expression and background on the same linear CPM scale."
  )
}

.match_neg_probe_rows <- function(row_names, negProbeName) {
  exact <- which(row_names %in% negProbeName)
  if (length(exact) > 0) {
    return(exact)
  }

  if (length(negProbeName) != 1) {
    return(integer())
  }

  which(startsWith(row_names, paste0(negProbeName, "_")))
}
