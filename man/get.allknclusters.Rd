% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/knutils.R
\name{get.allknclusters}
\alias{get.allknclusters}
\title{Extract indices of the areas in the clusters detected}
\usage{
get.allknclusters(spdf, knresults)
}
\arguments{
\item{spdf}{Spatial object with data used in the detection of clusters.}

\item{knresults}{Table with the clusters detected.}
}
\value{
A categorical vector with value the cluster to which area belongs.
It is the empty string for regions not in a cluster.
}
\description{
This function returns a categorical vector that identifies
to which cluster a given areas belongs. It is the empty string for 
areas not in a cluster.
}
