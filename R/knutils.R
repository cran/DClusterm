##' @title Computes the probability that a model parameter is <=k from inla marginals
##'
##' @description This function will be used to calculate the P(coeficient variable cluster <=0)
##'
##' @param func is the inla marginals of the model parameter
##' @param k is the cutoff
##' 
##' @return probability model coefficient <=k
##'
computeprob <- function(func, k) {
  idx <- which(func[, 1] <= k)
  if(length(idx) > 0) {
    weights <- func[idx + 1, 1] - func[idx, 1]

    prob <- sum(func[idx, 2] * weights)
  }
  else
  {
    prob<-0
  }
  return(prob)
}



##' @title Gets areas in a spatio-temporal cluster
##
##' @description This function is similar to get.knclusters but it also allows
##' for spatio-temporal clusters.
##
##' @param stfdf A sp or spacetime object with the information about the data.
##' @param results Results from a call to DetectClusterModel
##' 
##' @return A list with as many elements as clusters in 'results'
##'
get.stclusters <- function(stfdf, results) {
  if(inherits(stfdf, "Spatial")) {
    d <- as.data.frame(coordinates(stfdf))
    names(d) <- c("x", "y")
    return(get.knclusters(d, results))
  } else {
    d <- as.data.frame(coordinates(stfdf@sp))
    names(d) <- c("x", "y")
    knres <- get.knclusters(d, results)

    res <- as.list(rep(NA, nrow(results)))

    tms <- stfdf@time

    nsp <- nrow(d)
    ntms <- length(tms)

    for(i in 1:length(res)) {
      tidx <- which(as.Date(time(tms)) >= as.Date(results$minDateCluster[i]) &
       as.Date(time(tms)) <= as.Date(results$maxDateCluster[i]) )

      res[[i]] <- as.vector(sapply(tidx, function(X){
       (X-1) * nsp + knres[[i]]}
      ))

    }
  }

  return(res)
}



##' @title Constructs data frame with clusters in binary format.
##'
##' @description This function constructs a data frame with number of columns equal to the
##' number of clusters. Each column is a binary representation of one of the
##' clusters. The position i of the column is equal to 1 if the polygon i is
##' in the cluster or 0 if it is not in the cluster.
##'
##' @param datamap data of the SpatialPolygonsDataFrame with the polygons
##' of the map.
##' @param knresults data frame with information of the detected clusters.
##' Each row represents the information of one of the clusters.
##' It contains the coordinates of the center, the size, the start
##' and end dates, the log-likelihood ratio, a boolean indicating if it is a
##' cluster (TRUE in all cases), and the p-value of the cluster.
##'
##' @return data frame where the columns represent the clusters in binary format.
##' The position i of the column is equal to 1 if the polygon i is in the cluster
##' or 0 if it is not in the cluster.
##'
knbinary <- function(datamap, knresults) {
  clusters <- get.stclusters(datamap, knresults)
  res <- lapply(clusters, function(X, n) {
   v <- rep(0,n)
   v[X] <- 1
   return(v)
   }, n = dim(datamap@data)[1])

  res <- data.frame(matrix(unlist(res), nrow = dim(datamap@data)[1]))
  names(res) <- paste("CL", 1:length(clusters), sep = "")
  return(res)
}


##' @title Merges clusters so that they are identifed as levels of a factor.
##'
##' @description Given a data frame with clusters that do not overlap 
##' this function merges the clusters and construct a factor.
##' The levels of the factor are "NCL" if the polygon of the map is not
##' in any cluster, and "CL" if the polygon i is in cluster i.
##'
##' @param datamap data of the SpatialPolygonsDataFrame with the polygons
##' of the map.
##' @param knresults Data frame with information of the detected clusters.
##' Each row represents the information of one of the clusters.
##' It contains the coordinates of the center, the size, the start
##' and end dates, the log-likelihood ratio, a boolean indicating if it is a
##' cluster (TRUE in all cases), and the p-value of the cluster.
##' @param indClustersPlot rows of knresults that denote the clusters to be plotted.
##'
##' @return factor with levels that represent the clusters.
##'
mergeknclusters <- function(datamap, knresults, indClustersPlot) {
  n <- nrow(knresults)
  knbin <- as.matrix(knbinary(datamap, knresults))	
  res <- as.factor(knbin %*% matrix(1:n))
  levels(res) <- c("NCL", paste("CL", indClustersPlot, sep = "") )
  return(res)
}



##' @title Plots the clusters that do not overlap.
##' 
##' @description This function plots the detected clusters that do not overlap.
##' There are as many windows as different start dates. All clusters
##' with the same start date are represented in the same window.
##'
##' @param statsAllClustersNoOverlap data frame with information of the detected clusters
##' that no overlap. Each row represents the information of one of the clusters.
##' It contains the coordinates of the center, the size, the start
##' and end dates, the log-likelihood ratio, a boolean indicating if it is a
##' cluster (TRUE in all cases), and the p-value of the cluster.
##' @param colors vector with the colors of the clusters.
##' @param map SpatialPolygonsDataFrame with the polygons of the map.
##'
##' @return plots of the detected clusters for each start date.
##'
PlotClustersNoOverlap <- function(statsAllClustersNoOverlap, colors, map) {

  # Name of the clusters. Cluster's order by their significance
  nameClusters <- 1:nrow(statsAllClustersNoOverlap)

  # First maps are the ones with clusters with lower minDateCluster
  sortMinDateCluster <- sort(unique(statsAllClustersNoOverlap$minDateCluster))

  lsort <- length(sortMinDateCluster)
  for(i in 1:lsort) {

    indClustersPlot <-
     which(statsAllClustersNoOverlap$minDateCluster == sortMinDateCluster[i])
    knslim <- statsAllClustersNoOverlap[indClustersPlot, ]
    map$clusters <-
     mergeknclusters(as(map, "data.frame"), knslim, indClustersPlot)

    textLegend<-c("       End date:",
     paste(nameClusters, ". ", 
     statsAllClustersNoOverlap$maxDateCluster)[indClustersPlot])
    colLegend<-c("white", colors[indClustersPlot])

    p1 = spplot(map, "clusters",  col = "gray", col.regions = colLegend,
     colorkey = FALSE, 
     key = list(space = "right", 
      points = list(pch = 19, cex = 1.8, col = colLegend), 
      text = list(textLegend)), 
     sp.layout = list(list("sp.text", as.matrix(knslim[, 1:2]), 
      nameClusters[indClustersPlot], cex=1, font=2)),
     main = paste("Start date", sortMinDateCluster[i]))

    dev.new()

    plot(p1)
    #print(p1, position = c(0+(i-1)/lsort,0,i/lsort,1), more=T)
  }
}

#' @title Remove overlapping clusters
#
#' @description This function slims the number of clusters down.
#' The spatial scan statistic is known to detect duplicated
#' clusters. This function aims to reduce the number of clusters
#' by removing duplicated and overlapping clusters.
#'
#' @param d Data.frame with data used in the detection of clusters.
#' @param knresults Object returned by function opgam() with the clusters detected.
#' @param minsize Minimum size of cluster (default to 1).
#'
#' @return A subset of knresults with non-overlaping clusters of at least
#' minsize size.
#'

slimknclusters<-function(d, knresults, minsize = 1)
{
        #Filter by minsize
        knresults <- knresults[which(knresults$size >= minsize), ]

        #Ordering according to the test statistic
        idxcl <- rev(order(knresults$statistic))

        knbin <- knbinary(d, knresults)

        clusters <- c()
        while(length(idxcl) > 0)
        {
                cl <- idxcl[1]

                if(length(idxcl) > 0)
                {
                  res <- apply(as.matrix(knbin[, idxcl]), 2,
                        function(X, clbin){sum(X * clbin)},
                        clbin = knbin[, cl])
                   #Here is where we decide what clusters to remove
                   idxrem <- which(res > 0)
                   idxcl <- idxcl[-c(1, idxrem)]
                }
                else
                {
                  idxcl <- c() #idxcl[-c(1)]
                }

                clusters <- c(clusters, cl)
        }

        return(knresults[clusters, ])
}

