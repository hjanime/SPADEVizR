% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/viewers.R
\name{streamgraphViewer}
\alias{streamgraphViewer}
\title{Visualization of cluster abundance dynamics}
\usage{
streamgraphViewer(Results, samples = NULL, clusters = NULL,
  use.relative = FALSE, show.on_device = TRUE)
}
\arguments{
\item{Results}{a SPADEResults or Results object}

\item{samples}{a character vector providing the sample names to used (all samples by default)}

\item{clusters}{a character vector containing the clusters names to be visualized}

\item{use.relative}{a logical specifying if the visualization should be performed on relative abundance}

\item{show.on_device}{a logical specifying if the representation will be displayed on device}
}
\value{
a 'ggplot' object
}
\description{
Generate a streamgraph representation showing the dynamic evolution of the number of cells in clusters across samples.
The 'clusters' parameter is required.
}
\details{
The order of samples in the 'samples' vector correspond to the order where the sample will be displayed
}

