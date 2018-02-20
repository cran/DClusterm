\name{NY8}
\encoding{latin1}
\docType{data}
\alias{NY8}
\alias{NY8_clusters}
\alias{ny.m0}
\alias{ny.m1}
\alias{ny.cl0}
\alias{ny.cl1}
\alias{TCE}
\title{Leukemia in an eight-county region of upstate New York, 1978-1982.}
\description{This data set provides the number of incident leukemia cases per census tract in an eight-county
region of upstate New York in the period 1978-1982. In addition, the data set also includes information about
the location of the census tracts,
the population in 1980,
the inverse of the distance to the nearest Trichloroethene (TCE) site,
the percentage of people aged 65 or more, and
the percentage of people who own their home.

The dataset also provides the locations of the TCE sites.

File NY8_clusters contains the results of running DetectClustersModel on a null
model ('ny.m0') and another one with covariates ('ny.m1'). The results are in
'ny.cl0' and 'ny.cl1', respectively.

}
\usage{data(NY8)}
\format{A SpatialPolygonsDataFrame with 281 polygons representing the census tracts, and the following information about each census tract:
\tabular{rl}{
AREANAME    \tab  Name\cr
AREAKEY     \tab  Identifier \cr
X           \tab  x coordinate \cr
Y           \tab  y coordinate \cr
POP8        \tab  Population in 1980\cr
TRACTCAS    \tab  Number of leukemia cases rounded to 2 decimals \cr
PROPCAS     \tab  Ratio of the number of leukemia cases to the population in 1980 \cr
PCTOWNHOME  \tab  Proportion of people who own their home \cr
PCTAGE65P   \tab  Proportion of people aged 65 or more \cr
Z           \tab  \cr
AVGIDIST    \tab  \cr
PEXPOSURE   \tab  Inverse of the distance to the nearest TCE site \cr
Cases       \tab  Number of leukemia cases \cr
Xm          \tab  x coordinate (in meters)\cr
Ym          \tab  y coordinates(in meters)\cr
Xshift      \tab  Shifted Xm coordinate\cr
Yshift      \tab  Shifted Ym coordinate\cr
}}
\source{Waller and Gotway (2004) and Bivand et al. (2008)}

\references{
Bivand, R.S., E. J. Pebesma and V. Gómez-Rubio (2008). Applied Spatial Data
Analysis with R. Springer.

Waller, L., B. Turnbull, L. Clark, and P. Nasca (1992). Chronic disease
surveillance and testing of clustering of disease and exposure: application
to leukemia incidence in tce-contamined dumpsites in upstate New York.
Environmetrics 3, 281-300

Waller, L. A. and C. A. Gotway (2004). Applied Spatial Statistics for Public
Health Data. John Wiley & Sons, Hoboken, New Jersey. \url{http://web1.sph.emory.edu/users/lwaller/WGindex.htm}


}

\keyword{datasets}
