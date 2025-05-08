\name{Navarre}
\encoding{latin1}
\docType{data}
\alias{Navarre}
\alias{brainnav}
\alias{brainnavnb}
\title{Brain cancer in males in Navarre, Spain, 1988-1994.}
\description{This data set contains the male mortality due to brain cancer in the 40 basic health zones (BHZ) in Navarre over the period 1988-1994, and the neighborhood structure of the BHZ in Navarre.
In addition, the data set also includes information about
the location of the BHZ,
the expected cases, the Standardized Mortality Ratio (SMR),
relative risk estimates and 95\% confidence intervals.
}
\usage{data(Navarre)}
\format{brainnav: A SpatialPolygonsDataFrame with 40 polygons representing the basic health zones (BHZ) in Navarre, and the following information about each BHZ:
\tabular{rl}{
ZBS        \tab  \cr Basic Health Zone Code
NAME       \tab  Name \cr
OBSERVED   \tab  Number of observed brain cancer cases in males \cr
EXPECTED   \tab  Number of expected brain cancer cases in males. They are computed using indirect age-standardization using Navarre as a standard population. \cr
RISK       \tab   Relative Risk Estimates \cr
RISKLL     \tab   Relative 95\% confidence interval, lower limit \cr
RISKUL     \tab   Relative 95\% confidence interval, upper limit \cr
SMR        \tab Standardized Mortality Ratio (OBSERVED/EXPECTED) \cr
x          \tab x coordinate \cr
y          \tab y coordinate \cr
}
brainnavnb: A neighbor (nb) object which contains the index numbers of the neighbors of each BHZ.
}
\source{ 
Data set obtained from Ugarte et al. (2004). Boundaries downloaded 
in shapefile format from \url{https://geoportal.navarra.es/es/idena}. These have been thinned to reduce space use.
}
\references{
Ugarte, M. D., B. Ibáñez, and A. F. Militino (2004). Testing for poisson zero a
inflation in disease mapping. Biometrical Journal 46 (5), 526-539.

Ugarte, M. D., B. Ibáñez, and A. F. Militino (2006). Modelling risks in a
disease mapping. Statistical Methods in Medical Research 15, 21-35.
}
\keyword{datasets}
