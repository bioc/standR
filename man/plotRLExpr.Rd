% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotRLE.R
\name{plotRLExpr}
\alias{plotRLExpr}
\alias{plotRLExpr,DGEList-method}
\alias{plotRLExpr,ExpressionSet-method}
\alias{plotRLExpr,SummarizedExperiment-method}
\title{Compute and plot relative log expression (RLE) values of gene expression data}
\usage{
plotRLExpr(object, ordannots = c(), ...)

\S4method{plotRLExpr}{DGEList}(object, ordannots = c(), ...)

\S4method{plotRLExpr}{ExpressionSet}(object, ordannots = c(), ...)

\S4method{plotRLExpr}{SummarizedExperiment}(object, ordannots, assay = 1, sce_thresh = 1000, ...)
}
\arguments{
\item{object}{a DGEList, SummarizedExperiment or ExpressionSet object
containing gene expression data.}

\item{ordannots}{variables or computations to sort samples by (tidy style).}

\item{...}{aesthetic mappings to pass to \code{ggplot2::aes_string()}.}

\item{assay}{a numeric or character, specifying the assay to use (for
\code{SummarizedExperiment} and its derivative classes).}

\item{sce_thresh}{Integer value. The threshold of sample size for using
dot plot instead of box plot.}
}
\value{
a ggplot2 object, containing the RLE plot.
}
\description{
Compute and plot relative log expression (RLE) values of gene expression data
}
\examples{
data("dkd_spe_subset")
plotRLExpr(dkd_spe_subset)

}
