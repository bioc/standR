% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotDR.R
\name{plotMDS}
\alias{plotMDS}
\alias{plotMDS,DGEList-method}
\alias{plotMDS,ExpressionSet-method}
\alias{plotMDS,SummarizedExperiment-method}
\alias{plotMDS,SingleCellExperiment-method}
\alias{plotMDS,SpatialExperiment-method}
\title{Compute and plot the results of a MDS analysis on gene expression data}
\usage{
plotMDS(
  object,
  dims = c(1, 2),
  precomputed = NULL,
  textScale = 1,
  assay = 1,
  ...
)

\S4method{plotMDS}{DGEList}(
  object,
  dims = c(1, 2),
  precomputed = NULL,
  textScale = 1,
  assay = 1,
  ...
)

\S4method{plotMDS}{ExpressionSet}(
  object,
  dims = c(1, 2),
  precomputed = NULL,
  textScale = 1,
  assay = 1,
  ...
)

\S4method{plotMDS}{SummarizedExperiment}(
  object,
  dims = c(1, 2),
  precomputed = NULL,
  textScale = 1,
  assay = 1,
  ...
)

\S4method{plotMDS}{SingleCellExperiment}(
  object,
  dims = c(1, 2),
  precomputed = NULL,
  textScale = 1,
  assay = 1,
  ...
)

\S4method{plotMDS}{SpatialExperiment}(
  object,
  dims = c(1, 2),
  precomputed = NULL,
  textScale = 1,
  assay = 1,
  ...
)
}
\arguments{
\item{object}{a DGEList, SummarizedExperiment or ExpressionSet object
containing gene expression data.}

\item{dims}{a numeric, containing 2 values specifying the dimensions to plot.}

\item{precomputed}{a dimensional reduction results from either
\code{limma::plotMDS}.}

\item{textScale}{a numeric, specifying the relative scale factor to apply to text on
the plot.}

\item{assay}{a numeric or character, specifying the assay to use (for
\code{SummarizedExperiment} and its derivative classes).}

\item{...}{aesthetic mappings to pass to \code{ggplot2::aes_string()}.}
}
\value{
a ggplot2 object
}
\description{
Compute and plot the results of a MDS analysis on gene expression data
}
\examples{
data("dkd_spe_subset")
standR::plotMDS(dkd_spe_subset)

}
