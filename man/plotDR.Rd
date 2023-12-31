% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotDR.R
\name{plotDR}
\alias{plotDR}
\alias{plotDR,SingleCellExperiment-method}
\alias{plotDR,SpatialExperiment-method}
\title{Compute and plot the results of any dimension reduction methods on gene expression data}
\usage{
plotDR(object, dims = c(1, 2), ...)

\S4method{plotDR}{SingleCellExperiment}(object, dims, dimred = "PCA", textScale = 1, ...)

\S4method{plotDR}{SpatialExperiment}(object, dims, dimred = "PCA", textScale = 1, ...)
}
\arguments{
\item{object}{a DGEList, SummarizedExperiment or ExpressionSet object
containing gene expression data.}

\item{dims}{a numeric, containing 2 values specifying the dimensions to plot.}

\item{...}{aesthetic mappings to pass to \code{ggplot2::aes_string()}.}

\item{dimred}{a string or integer scalar indicating the reduced dimension
result in \code{reducedDims(object)} to plot.}

\item{textScale}{a numeric, specifying the relative scale factor to apply to text on
the plot.}
}
\value{
a ggplot2 object
}
\description{
Compute and plot the results of any dimension reduction methods on gene expression data
}
\examples{
library(scater)
data("dkd_spe_subset")
spe <- scater::runPCA(dkd_spe_subset)
plotDR(spe, dimred = "PCA")

}
