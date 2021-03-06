% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/stats.R
\name{computeEigenCellClusters}
\alias{computeEigenCellClusters}
\title{Internal - Eigen vector classification}
\usage{
computeEigenCellClusters(data, eigencell.correlation.th = 0.8)
}
\arguments{
\item{data}{a numeric matrix with all clusters in rownames}

\item{eigencell.correlation.th}{a numeric value indicating the correlation coefficient threshold}
}
\value{
a dataframe containing for each cluster, its name and class
}
\description{
This function is used internally to classify clusters abundance profiles or phenotype profiles using eigen vector decomposition.
}
\details{
This method compute the performs a eigen vector decomposition and then calculate the correlations between the matrix rows and these vectors.
Clusters which correlate above a specific threshold with the same eigen vector are classified together.
This correlation threshold (a numeric double between 0 and 1 included, default is 0.8) provided using the 'eigencell.correlation.th' parameter.
}

