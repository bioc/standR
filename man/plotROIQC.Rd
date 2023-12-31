% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotLibQC.R
\name{plotROIQC}
\alias{plotROIQC}
\title{Plot Sample-wise QC plot}
\usage{
plotROIQC(
  spe_object,
  x_axis = "AOINucleiCount",
  y_axis = "lib_size",
  x_lab = "AOINucleiCount",
  y_lab = "Library size",
  x_threshold = NULL,
  y_threshold = NULL,
  regression_col = "purple",
  hist_col = "black",
  hist_fill = "white",
  bin_num = 50,
  threshold_col = "red",
  threshold_linetype = "dashed",
  layout_ncol = 2,
  layout_nrow = 2,
  leyout_height = c(0.8, 2.5),
  layout_width = c(2.5, 0.8),
  ...
)
}
\arguments{
\item{spe_object}{A SpatialExperiment object.}

\item{x_axis}{Numeric feature to plot as x axis.}

\item{y_axis}{Numeric feature to plot as y axis.}

\item{x_lab}{Label name for x axis.}

\item{y_lab}{Label name for y axis.}

\item{x_threshold}{Threshold to draw.}

\item{y_threshold}{Threshold to draw.}

\item{regression_col}{Color for the regression line.}

\item{hist_col}{Color for the histograms.}

\item{hist_fill}{Fill for the histograms.}

\item{bin_num}{Bin numbers for the histograms.}

\item{threshold_col}{Threshold line color.}

\item{threshold_linetype}{Threshold line type.}

\item{layout_ncol}{Column number layout.}

\item{layout_nrow}{Row number layout.}

\item{leyout_height}{Height layout.}

\item{layout_width}{Width layout.}

\item{...}{aesthetic mappings to pass to \code{ggplot2::aes()} of the dot plots.}
}
\value{
A ggplot object.
}
\description{
Plot Sample-wise QC plot
}
\examples{
library(ggplot2)
library(patchwork)
data("dkd_spe_subset")
spe <- addPerROIQC(dkd_spe_subset)

plotROIQC(spe)

}
