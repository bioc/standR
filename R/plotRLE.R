#' Compute and plot relative log expression (RLE) values of gene expression data
#'
#' @param ordannots variables or computations to sort samples by (tidy style).
#'
#' @inheritParams drawPCA
#' @param sce_thresh Integer value. The threshold of sample size for using 
#' dot plot instead of box plot.
#' @return a ggplot2 object, containing the RLE plot.
#' @export
#'
#' @examples
#' data("dkd_spe_subset")
#' plotRLExpr(dkd_spe_subset)
#'
setGeneric(
  "plotRLExpr",
  function(object,
           ordannots = c(),
           ...) {
    standardGeneric("plotRLExpr")
  }
)

#' @rdname plotRLExpr
setMethod(
  "plotRLExpr",
  signature("DGEList", "ANY"),
  function(object, ordannots, ...) {
    # extract sample data
    sdata <- object$samples
    # extract expression data (and transform)
    object <- edgeR::cpm(object, log = TRUE)
    # create data structure
    samporder <- orderSamples(sdata, ordannots)
    rledf <- pdataRLE_intl(object, samporder)
    p1 <- plotRLExpr_intl(rledf, sdata, isSCE = FALSE, ...)

    return(p1)
  }
)

#' @rdname plotRLExpr
setMethod(
  "plotRLExpr",
  signature("ExpressionSet", "ANY"),
  function(object, ordannots, ...) {
    # extract sample data
    sdata <- Biobase::pData(object)
    # extract expression data (and transform)
    object <- Biobase::exprs(object)
    # create data structure
    samporder <- orderSamples(sdata, ordannots)
    rledf <- pdataRLE_intl(object, samporder)
    p1 <- plotRLExpr_intl(rledf, sdata, isSCE = FALSE, ...)

    return(p1)
  }
)

#' @rdname plotRLExpr
setMethod(
  "plotRLExpr",
  signature("SummarizedExperiment", "ANY"),
  function(object, ordannots, assay = 1, sce_thresh = 1000, ...) {
    isSCE <- is(object, "SingleCellExperiment")

    # extract sample data
    sdata <- BiocGenerics::as.data.frame(colData(object), optional = TRUE)

    # extract expression data (and transform)
    object <- assay(object, i = assay)
    # create data structure
    samporder <- orderSamples(sdata, ordannots)
    rledf <- pdataRLE_intl(object, samporder)
    p1 <- plotRLExpr_intl(rledf, sdata, isSCE = isSCE, 
                          sce_thresh = sce_thresh, ...)

    return(p1)
  }
)

# plot data preparation using MDS results
pdataRLE_intl <- function(emat, sampord) {
  # compute RLE
  rle <- emat - Biobase::rowMedians(as.matrix(emat))
  # order samples
  rle <- rle[, sampord]

  # compute boxplot
  rledf <- t(apply(rle, 2, function(x) grDevices::boxplot.stats(x)$stats))
  rledf <- as.data.frame(rledf)
  colnames(rledf) <- c("ymin", "lower", "middle", "upper", "ymax")
  rledf$x <- seq(nrow(rledf))
  rledf$RestoolsMtchID <- rownames(rledf)

  return(rledf)
}

plotRLExpr_intl <- function(plotdf, sdata, isSCE = FALSE, textScale = 1, 
                            sce_thresh = 1000, ...) {

  # constant - sample size at which standard plot becomes dense
  dense_thresh <- 50
  sce_thresh <- sce_thresh

  # extract aes
  aesmap <- rlang::enquos(...)

  # annotate samples
  plotdf <- addSampleAnnot(plotdf, sdata)

  # compute plot
  aesmap <- aesmap[!names(aesmap) %in% c("x", "ymin", "ymax", 
                                         "upper", "middle", "lower")] # remove fixed mappings if present

  # split aes params into those that are not aes i.e. static parametrisation
  if (length(aesmap) > 0) {
    is_aes <- vapply(aesmap, rlang::quo_is_symbolic, FUN.VALUE = logical(1))
    defaultmap <- lapply(aesmap[!is_aes], rlang::eval_tidy)
    aesmap <- aesmap[is_aes]
  } else {
    defaultmap <- list()
  }

  # build plot
  if (isSCE & nrow(plotdf) > sce_thresh) {
    p1 <- ggplot2::ggplot(plotdf, aes(x = upper - lower, y = middle, !!!aesmap)) +
      ggplot2::geom_point() +
      ggplot2::geom_hline(yintercept = 0, colour = 2, lty = 2) +
      ggplot2::labs(y = "Relative log expression (median)", x = "Relative log expression (IQR)") +
      do.call(ggplot2::geom_point, defaultmap) +
      bhuvad_theme(textScale)
  } else {
    p1 <- ggplot2::ggplot(plotdf, aes(x = x, y = middle, group = x, !!!aesmap)) +
      ggplot2::geom_boxplot(
        aes(ymin = ymin, ymax = ymax, upper = upper, middle = middle, lower = lower),
        stat = "identity"
      ) +
      ggplot2::geom_hline(yintercept = 0, colour = 2, lty = 2) +
      ggplot2::ylab("Relative log expression") +
      do.call(ggplot2::geom_boxplot, defaultmap) +
      bhuvad_theme(textScale) +
      ggplot2::theme(axis.text.x = element_blank())

    # update plot if too many samples are plot
    if (nrow(plotdf) > dense_thresh) {
      ## geom_point will inherit relevant aesthetics from top `aes`, include y=middle
      p1 <- p1 + ggplot2::geom_point()
    }
  }

  return(p1)
}


utils::globalVariables(c("upper", "lower", "middle", "x", "ymin", "ymax"))
