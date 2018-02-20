\name{brainNM}
\docType{data}
\alias{brainNM}
\alias{brainst}
\alias{losalamos}
\alias{brainNM_clusters}
\alias{nm.m0}
\alias{nm.m1}
\alias{nm.cl0}
\alias{nm.cl1}
\title{Brain cancer in New Mexico, USA, 1973-1991.}
\description{This data set contains the number of incident brain cancer cases in the 32 counties of New Mexico, USA, and each year of the period 1973-1991,
and the location of Los Alamos National Laboratory.
In addition, the data set also includes for each county and year information about
the expected cases, the Standardized Morbidity Ratio (SMR),
the FIPS, ...

File brainNM_clusters contains the results of running DetectClustersModel on a
null model ('nm.m0') and another one with covariates ('nm.m1'). The results are
in 'nm.cl0' and 'nm.cl1', respectively.


}
\usage{data(brainNM)}
\format{brainst: A STFDF object containing the following information for each county and year:
\tabular{rl}{
Observed  \tab  Number of observed brain cancer cases \cr
Expected  \tab  Number of expected brain cancer cases. Standardisation is done taking the whole time-period and not year-ly to keep any temporal trend. \cr
SMR        \tab Standardized Morbidity Ratio (observed/expected) \cr
Year       \tab Year \cr
FIPS      \tab  FIPS Code \cr
ID         \tab  ID (from 1 to 32)\cr
IDLANL    \tab  Inverse distance to Los Alamos National Laboratory\cr
IDLANLre   \tab  Re-scaled Inverse distance to Los Alamos National Laboratory 
(i.e., IDLANL/mean(IDLANL))\cr
}
losalamos: A SpatialPoints object which contains the location (in long/lat) of Los Alamos National Laboratory obtained from the Wikipedia: -106.298333, 35.881667.
}
\source{ 
Data have been downlodad from the SatScan  website. Boundaries have been
obtained from the U.S. Census Bureau. Cibola and Valencia counties has been
merged together.
}
\references{
SatScan (c). \url{http://www.spatstat.org}

Kulldorff, M., W. F. Athas, E. J. Feurer, B. A. Miller, and C. R. Key (1998).
Evaluating cluster alarms: a space-time scan statistic and brain cancer in
los alamos, new mexico. American Journal of Public Health 88, 1377-1380.

}
\keyword{datasets}









 








