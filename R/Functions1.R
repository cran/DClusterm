##' @title Detects clusters and computes their significance.
##' 
##' @description Searches all possible clusters with start and end dates within minDateUser
##' and maxDateUser, so that the maximum fraction of the total population inside
##' the cluster is less than fractpop, and the maximum distance to the center is
##' less than radius.
##' The search can be done for spatial or spatio-temporal clusters.
##' The significance of the clusters is obtained with a Monte Carlo procedure
##' or based on the chi-square distribution (\link{glm}, \link{glmer} or \link{zeroinfl} models)
##' or DIC (\link[INLA]{inla} models).
##'
##' @param stfdf object containing the data.
##' If data is spatial, stfdf is a \link{SpatialPolygonsDataFrame} object from sp.
##' If data is spatio-temporal, stfdf is a \link{STFDF} object from spacetime.
##' The data contain a \link{SpatialPolygons} object with the coordinates,
##' and if applicable, a time object holding time information,
##' an endTime vector of class \link{POSIXct} holding end points of time intervals.
##' It also contain a data.frame with the Observed, Expected and potential covariates
##' in each location and time (if applicable). Note that the function DetectClustersModel
##' does not use the endTime vector. We can define endTime, for example,
##' as the vector of class \link{POSIXct} which contains the same dates as the ones
##' contained in the time object. 
##' @param thegrid two-columns matrix containing the points of the grid to be
##' used. If it is null, a rectangular grid is built.
##' @param radius maximum radius of the clusters.
##' @param step step of the thegrid built.
##' @param fractpop maximum fraction of the total population inside the cluster.
##' @param alpha significance level used to determine the existence of clusters.
##' @param typeCluster type of clusters to be detected. "ST" for spatio-temporal
##' or "S" spatial clusters.
##' @param minDateUser start date of the clusters.
##' @param maxDateUser end date of the clusters.
##' @param R If the cluster's significance is calculated based on the chi-square
##' distribution or DIC, R is NULL. If the cluster's significance is calculated using a
##' Monte Carlo procedure, R represents the number replicates under the null hypothesis.
##' @param model0 Initial model (including covariates).
##' @param ClusterSizeContribution Indicates the variable to be used as the population at risk in the cluster. This is the variable name to be used by 'fractpop' when checking the fraction of the population inside the cluster. The default column name is 'Population'.
##' This can be "glm" for generalized linear models (\link{glm} {stats}),
##' "glmer" for generalized linear mixed model (\link{glmer} {lme4}),
##' "zeroinfl" for zero-inflated models (\link{zeroinfl} {pscl}), or
##' "inla" for generalized linear, generalized linear mixed or zero-inflated models fitted with \link[INLA]{inla}.
##'
##' @return data frame with information of the detected clusters ordered by its
##' log-likelihood ratio value or DIC. Each row represents the information of
##' one of the clusters. It contains the coordinates of the center, the size,
##' the start and end dates, the log-likelihood ratio or DIC, the p-value,
##' the risk of the cluster, and a boolean indicating if it is a
##' cluster (TRUE in all cases).
##' It also returns alpha_bonferroni which is the level of significance
##' adjusted for multiple testing using Bonferroni correction.
##' Thus, rows that should be considered clusters are the ones with
##' p-value less than alpha_bonferroni.
##'
##' @references
##'
##' Bilancia M, Demarinis G (2014) Bayesian scanning of spatial disease rates
##' with the Integrated Nested Laplace Approximation (INLA). Statistical
##' Methods & Applications 23(1): 71 - 94. \url{http://dx.doi.org/10.1007/s10260-013-0241-8}
##'
##' Jung I (2009) A generalized linear models approach to spatial scan 
##' statistics for covariate adjustment. Statistics in Medicine 28(7): 1131 - 1143.
##'
##' Fast Bayesian classification for disease mapping and the detection of 
##' disease clusters (2017) Gomez-Rubio V, Molitor J, Moraga P. Submitted.
##'
##' @examples
##' library("DClusterm")
##' data("NY8")
##'
##' NY8$Observed <- round(NY8$Cases)
##' NY8$Expected  <- NY8$POP8 * sum(NY8$Observed) / sum(NY8$POP8)
##'
##' NY8$x <- coordinates(NY8)[, 1]
##' NY8$y <- coordinates(NY8)[, 2]
##'
##'
##' #Model to account for covariates
##' ny.m1 <- glm(Observed ~ offset(log(Expected)) + PCTOWNHOME + PCTAGE65P +
##' PEXPOSURE, family = "poisson", data = NY8)
##'
##' #Indices of areas that are possible cluster centres
##' idxcl <- c(120, 12, 89, 139, 146)
##' 
##' #Cluster detection adjusting for covariates
##' ny.cl1 <- DetectClustersModel(NY8,
##' thegrid = as.data.frame(NY8)[idxcl, c("x", "y")],
##' fractpop = 0.15, alpha = 0.05,
##' typeCluster = "S", R = NULL, model0 = ny.m1,
##' ClusterSizeContribution = "POP8")
##'
##' #Display results
##' ny.cl1
##'
##' @export
##'
##' @import DCluster
##' @import parallel
##' @import sp
##' @import spacetime
##' 
##' @importFrom methods as 
##' @importFrom xts xts 
##' @importFrom stats time pchisq rpois formula fitted
##' @importFrom stats coef deviance family glm logLik model.matrix
##' @importFrom grDevices dev.new
##' @importFrom pscl zeroinfl
DetectClustersModel <- function(stfdf, thegrid = NULL, radius = Inf,
  step = NULL, fractpop, alpha, typeCluster = "S",
  minDateUser = NULL, maxDateUser = NULL,
  R = NULL, model0, ClusterSizeContribution = "Population") {
  
  
  #############################
  # If data is spatial, stfdf is a \link{SpatialPolygonsDataFrame} object. We need to convert it to \link{STFDF}. We add date as.Date("1970-01-01").
  # If data is spatio-temporal, stfdf is a \link{STFDF} object.

  #Spatial object?
  is.sp <- FALSE
  if(class(stfdf) == "SpatialPolygonsDataFrame"){
    is.sp <- TRUE
    stfdf <- STFDF(as(stfdf, "SpatialPolygons"), xts(1, as.Date("1970-01-01")),
                   stfdf@data, endTime = as.POSIXct(strptime(c("1970-01-01"), "%Y-%m-%d"), tz = "GMT"))
  }
  if(is.null(minDateUser)){
    minDateUser <- min(time(stfdf@time))
  } 
  if(is.null(maxDateUser)){
    maxDateUser <- max(time(stfdf@time))
  }
  
  #############################

  # Create column with ID. Unique identifier
  stfdf[['ID']] <- 1:nrow(stfdf@data)

  sortDates <- sort(unique(time(stfdf@time)))

  # Check minDateUser and maxDateUser make sense
  minDateData <- min(time(stfdf@time))
  maxDateData <- max(time(stfdf@time))

  if(minDateUser > maxDateUser) {
    stop('Cluster minimum date is greater than cluster maximum date')
  }

  if(minDateUser > maxDateData) {
    stop('Cluster minimum date is greater than data set maximum date')
  }

  if(maxDateUser < minDateData) {
    stop('Cluster maximum date is smaller than data set minimum date')
  }

  # Closest dates to minDateUser and maxDateUser
  idMinDateCluster <- min(which(sortDates >= minDateUser))
  idMaxDateCluster <- max(which(sortDates <= maxDateUser))


  # If grid is null, create a new grid
  if(is.null(thegrid)) {
    CreateGridDClusterm(stfdf, radius, step)
  }

  # Radius
  rr <- radius * radius

  ## Init snowfall
  #if(!is.null(numCPUS)){
  #sfInit( parallel=TRUE, cpus=numCPUS )
  #sfLibrary(spdep)
  #sfLibrary(splancs)
  #sfLibrary(spacetime)
  #sfLibrary(DCluster)
  #sfLibrary(pscl)
  #sfLibrary(INLA)
  #sfLibrary(DClusterm)
  ##sfSource("R/Functions1PAU.R")
  ##sfSource("R/Functions2PAU.R")
  ##sfSource("R/glm.isclusterPAU.R")
  ##sfSource("R/knutils.R")
  #}

  #Init parallel
  numCPUS <- getOption("mc.cores")
  if(!is.null(numCPUS)) {
    cl <- makeCluster(numCPUS)
    setDefaultCluster(cl)
  }


  # Statistic of each cluster
  statsAllClusters <- CalcStatsAllClusters(thegrid, CalcStatClusterGivenCenter,
    stfdf, rr,  typeCluster, sortDates, idMinDateCluster, idMaxDateCluster,
    fractpop, model0, ClusterSizeContribution, numCPUS)
  
  #######################
  # Bonferroni correction to deal with multiple testing problem
  # reject null hypothesis if the p-value is less than alpha/num_tests
  alphaoriginal <- alpha
  #alpha <- alpha/nrow(statsAllClusters)
  alpha_bonferroni <- alpha/nrow(statsAllClusters)
  #######################
  


  # Remove rows where sizeCluster == -1
  idRemove <- which(statsAllClusters$sizeCluster == -1)
  if(length(idRemove) > 0) {
    statsAllClusters <- statsAllClusters[-idRemove, ]
  }

  # If there are no clusters return "No clusters found"
  #FIXME: Retiurn something different
  if(dim(statsAllClusters)[1] == 0) {
    print("No clusters found")
    return("No clusters found")
  }

  if(class(model0)[1] == "inla") {
    vecpvalue <- statsAllClusters$pvalue
    veccluster <- vecpvalue < alpha
  }
  else {

    # p-value of each cluster
    vecpvalue <- matrix(NA, dim(statsAllClusters)[1], 1)
    veccluster <- matrix(NA, dim(statsAllClusters)[1], 1)

    ##############################################################################################

  # 1. p-value without Monte Carlo
  if(is.null(R)) {
    for(i in 1:nrow(statsAllClusters)) {
      vecpvalue[i] <- 1 - pchisq(2 * statsAllClusters$statistic[i], 1)
      veccluster[i] <- vecpvalue[i] < alpha
    }
  } else {

    # 2. p-value with Monte Carlo
    maxStatisticRReplicas <- matrix(NA, R, 1)
    stfdfMC <- stfdf

    for(i in 1:R) {
      # Generate data set under H_0
      #stfdfMC$Observed <- 
        rpois(length(stfdf$Observed), lambda = stfdf$Expected)

      obslab <- as.character(formula(model0))[2]

      stfdfMC[[obslab]] <- rpois(nrow(stfdf@data), lambda = fitted(model0))

      # Statistic of each cluster
      statsAllClustersMC <- CalcStatsAllClusters(thegrid, 
        CalcStatClusterGivenCenter, stfdfMC, rr, typeCluster, sortDates, 
        idMinDateCluster, idMaxDateCluster, fractpop, model0,
        ClusterSizeContribution,numCPUS)

      maxStatisticRReplicas[i] <- max(statsAllClustersMC$statistic)

      print(paste("replica",i))
    }

    # p-value according to rank
    for(i in 1:(dim(statsAllClusters)[1])) {
      vecpvalue[i] <- 
        (sum(maxStatisticRReplicas > statsAllClusters$statistic[i]) + 1)/(R + 1)
      veccluster[i]<-vecpvalue[i]<alpha
    }
  }

  ##############################################################################################

  }

  ## End snowfall
  #if(!is.null(numCPUS)){
  #sfStop()
  #}

  #Stop parallel cluster
  if(!is.null(numCPUS)){
    stopCluster(cl)
  }

  statsAllClusters$pvalue <- vecpvalue
  statsAllClusters$cluster <- veccluster

  #statsAllClusters <- cbind(statsAllClusters, veccluster, vecpvalue)
  names(statsAllClusters) <- c("x", "y", "size", "minDateCluster",
    "maxDateCluster", "statistic", "pvalue", "risk", "cluster")


  # Selection of significant clusters
  indpvalueNA <- which(is.na(statsAllClusters$pvalue))
  if(length(indpvalueNA) > 0) {
    statsAllClusters <- statsAllClusters[-indpvalueNA, ]
  }
  statsAllClusters <- statsAllClusters[statsAllClusters$pvalue < alpha, ]

  # If there are no clusters return "No clusters found"
  #FIXME: Return a different value, i.e., empty data.frame
  if(dim(statsAllClusters)[1] == 0) {
    print(paste("No significant clusters found with alpha =", alphaoriginal))
    return("No clusters found")
  }


  # Ordered results by statistic value
  if(class(model0)[1] == "inla") {
    statsAllClusters <- statsAllClusters[order(statsAllClusters$statistic), ]
  } else {
    statsAllClusters <-
      statsAllClusters[rev(order(statsAllClusters$statistic)), ]
  }

  # Return
  statsAllClusters$minDateCluster <- as.POSIXct(statsAllClusters$minDateCluster,
    origin = "1970-01-01", tz = "GMT")
  statsAllClusters$maxDateCluster <- as.POSIXct(statsAllClusters$maxDateCluster,
    origin = "1970-01-01", tz = "GMT")
  
  
  #######################
  # Bonferroni correction to deal with multiple testing problem
  # Add column alpha_bonferroni = alpha/num_tests
  # reject null hypothesis if the p-value is less than alpha_bonferroni
  if(nrow(statsAllClusters) > 0){
    statsAllClusters$alpha_bonferroni <- alpha_bonferroni
  }
  #######################
  
  #If stdf is a spatial object, then remove dates
  if(is.sp) {
    statsAllClusters <- statsAllClusters[, 
      -which(names(statsAllClusters) %in% c("minDateCluster", "maxDateCluster")
   )]
  }
  
  return(statsAllClusters)
}





##' @title Removes the overlapping clusters.
##' 
##' @description Function DetectClustersModel() detects duplicated clusters.
##' This function reduces the number of clusters by removing the overlapping
##' clusters.
##'
##' @param stfdf spatio-temporal class object containing the data.
##' @param statsAllClusters data frame with information of the detected
##' clusters obtained with DetectClustersModel().
##'
##' @return data frame with the same information than statsAllClusters but only
##' for clusters that do not overlap.
##'
##' @examples
##' library("DClusterm")
##' data("brainNM")
##' data("brainNM_clusters")
##' 
##' SelectStatsAllClustersNoOverlap(brainst, nm.cl1)
##' @export

SelectStatsAllClustersNoOverlap <- function(stfdf, statsAllClusters) {
  # statsAllClusters is ordered by statistic value
  coordx <- coordinates(stfdf@sp)[, 1]
  coordy <- coordinates(stfdf@sp)[, 2]

  idSpaceAllClustersNoOverlap <- NULL
  idTimeAllClustersNoOverlap <- NULL
  statsAllClustersNoOverlap <- NULL

  for(i in 1:(dim(statsAllClusters)[1])) {
    xd <- (coordx-statsAllClusters$x[i])
    yd <- (coordy-statsAllClusters$y[i])

    dist <- xd * xd + yd * yd

    idSpaceOneCluster <- order(dist)[1:statsAllClusters$size[i]]
    idTimeOneCluster <- 
      which(as.Date(time(stfdf@time)) >= as.Date(statsAllClusters$minDateCluster[i]) & as.Date(time(stfdf@time)) <= as.Date(statsAllClusters$maxDateCluster[i]))

    if(sum(idSpaceOneCluster %in% idSpaceAllClustersNoOverlap) ==0 || 
      sum(idTimeOneCluster %in% idTimeAllClustersNoOverlap) == 0) {

      statsAllClustersNoOverlap <-
        rbind(statsAllClustersNoOverlap, statsAllClusters[i, ])

      idSpaceAllClustersNoOverlap <-
        c(idSpaceAllClustersNoOverlap, idSpaceOneCluster)

      idTimeAllClustersNoOverlap <-
        c(idTimeAllClustersNoOverlap, idTimeOneCluster)
    }
  }

  return(statsAllClustersNoOverlap)
}

