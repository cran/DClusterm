### R code from vignette source 'jss_paper.Rnw'

###################################################
### code chunk number 1: jss_paper.Rnw:5-7
###################################################
#JSS style
options(prompt = "R> ", continue = "+  ", width = 70, useFancyQuotes = FALSE)


###################################################
### code chunk number 2: jss_paper.Rnw:291-295
###################################################
library("DClusterm")
library("xts")

data("NY8")


###################################################
### code chunk number 3: jss_paper.Rnw:306-311
###################################################
NY8$Observed <- round(NY8$Cases)

NY8$Expected  <- NY8$POP8 * sum(NY8$Observed) / sum(NY8$POP8)

NY8$SMR <- NY8$Observed / NY8$Expected


###################################################
### code chunk number 4: jss_paper.Rnw:320-322
###################################################
NY8$x <- coordinates(NY8)[, 1]
NY8$y <- coordinates(NY8)[, 2]


###################################################
### code chunk number 5: jss_paper.Rnw:334-337
###################################################
NY8st <- STFDF(as(NY8, "SpatialPolygons"), xts(1, as.Date("1972-01-01")),
  NY8@data, endTime = as.POSIXct(strptime(c("1972-01-01"), "%Y-%m-%d"),
  tz = "GMT"))


###################################################
### code chunk number 6: jss_paper.Rnw:343-362
###################################################
library("RColorBrewer")

#Save plots
p1 <- spplot(NY8, "SMR", cuts = 8, col.regions = brewer.pal(9, "Oranges"),
  main = "Standardised Mortality Ratio", 
  sp.layout = list(list("sp.points", TCE, col = "red")))
p2 <- spplot(NY8, "PCTOWNHOME", cuts = 8, col.regions = brewer.pal(9, "Blues"),
  main = "PCTOWNHOME")
p3 <- spplot(NY8, "PCTAGE65P", cuts = 8, col.regions = brewer.pal(9, "Greens"),
  main = "PCTAGE65P")
p4 <- spplot(NY8, "PEXPOSURE", cuts = 8, col.regions = brewer.pal(9, "Reds"),
  main = "PEXPOSURE")

#Display plots
print(p1, position = c(0, 0.5, 0.5, 1), more = TRUE)
print(p2, position = c(0.5, 0.5, 1, 1), more = TRUE)
print(p3, position = c(0, 0, 0.5, 0.5), more = TRUE)
print(p4, position = c(0.5, 0, 1, 0.5))



###################################################
### code chunk number 7: jss_paper.Rnw:380-382
###################################################
ny.m0 <- glm(Observed ~ offset(log(Expected)) + 1, family = "poisson",
  data = NY8)


###################################################
### code chunk number 8: jss_paper.Rnw:412-413
###################################################
options(mc.cores = 1)


###################################################
### code chunk number 9: jss_paper.Rnw:423-429
###################################################
idxcl <- c(120, 12, 89, 139, 146)
ny.cl0 <- DetectClustersModel(NY8st,
  thegrid = as.data.frame(NY8)[idxcl, c("x", "y")], 
  fractpop = 0.15, alpha = 0.05, radius = Inf, step = NULL,
  typeCluster = "S", R = NULL, model0 = ny.m0,
  ClusterSizeContribution = "POP8")


###################################################
### code chunk number 10: jss_paper.Rnw:445-446
###################################################
ny.cl0


###################################################
### code chunk number 11: jss_paper.Rnw:454-467
###################################################
get.allknclusters <- function (spdf, knresults) {
  clusters <- rep("", nrow(spdf))

  knclusters <- get.knclusters(spdf, knresults)
  if(length(knclusters) >0 ) {
    clusters[unique(unlist(knclusters))] <- "CLUSTER"
    clusters <- as.factor(clusters)
  }

  return(clusters)
}
  
NY8$CLUSTER0 <- get.allknclusters(NY8, ny.cl0)


###################################################
### code chunk number 12: jss_paper.Rnw:477-480 (eval = FALSE)
###################################################
## cl0.centres <- SpatialPoints(ny.cl0[, 1:2])
## spplot(NY8, "CLUSTER0", col.regions = c("white", "gray"), col = "#4D4D4D",
##   sp.layout = list(list("sp.points", cl0.centres, pch = 19, col = "black")))


###################################################
### code chunk number 13: jss_paper.Rnw:492-495
###################################################
ny.m1 <- glm(Observed ~ offset(log(Expected)) + PCTOWNHOME + PCTAGE65P + 
  PEXPOSURE, family = "poisson", data = NY8)
summary(ny.m1)


###################################################
### code chunk number 14: jss_paper.Rnw:503-508
###################################################
ny.cl1 <- DetectClustersModel(NY8st, 
  thegrid = as.data.frame(NY8)[idxcl, c("x", "y")], 
  fractpop = 0.15, alpha = 0.05,
  typeCluster = "S", R = NULL, model0 = ny.m1,
  ClusterSizeContribution = "POP8")


###################################################
### code chunk number 15: jss_paper.Rnw:511-512
###################################################
ny.cl1


###################################################
### code chunk number 16: jss_paper.Rnw:519-520
###################################################
NY8$CLUSTER1 <- get.allknclusters(NY8, ny.cl1)


###################################################
### code chunk number 17: jss_paper.Rnw:535-549
###################################################
library("gridExtra")
library("latticeExtra")

ny.map0 <- spplot(NY8, c("CLUSTER0"), col.regions = c("white", "gray"),
  col = "#4D4D4D", colorkey = FALSE) +
  layer(lpoints(ny.cl0$x, ny.cl0$y, pch = 19))

ny.map1 <- spplot(NY8, c("CLUSTER1"), col.regions = c("white", "gray"),
  col = "#4D4D4D", colorkey = FALSE) +
  layer(lpoints(ny.cl1$x, ny.cl1$y, pch = 19))

grid.arrange(ny.map0, ny.map1, ncol = 2)




###################################################
### code chunk number 18: jss_paper.Rnw:596-598
###################################################
DeanB(ny.m0)
DeanB2(ny.m0)


###################################################
### code chunk number 19: jss_paper.Rnw:607-609
###################################################
DeanB(ny.m1)
DeanB2(ny.m1)


###################################################
### code chunk number 20: jss_paper.Rnw:626-630
###################################################
library("lme4")
ny.mm0 <- glmer(Observed ~ offset(log(Expected)) + (1 | AREANAME), 
  data = as(NY8, "data.frame"), family = "poisson")
summary(ny.mm0)


###################################################
### code chunk number 21: jss_paper.Rnw:633-637
###################################################
ny.clmm0 <- DetectClustersModel(NY8st,
  thegrid = as.data.frame(NY8)[idxcl, c("x", "y")], fractpop = 0.15,
  alpha = 0.05, typeCluster = "S", R = NULL, model0 = ny.mm0,
  ClusterSizeContribution = "POP8")


###################################################
### code chunk number 22: jss_paper.Rnw:640-641
###################################################
ny.clmm0


###################################################
### code chunk number 23: jss_paper.Rnw:653-660 (eval = FALSE)
###################################################
## #Compute SMR in clusters
## cls <- get.knclusters(NY8, ny.cl0)
## lapply(cls, function(X) {
##   res <- c(sum(NY8$Observed[X]), sum(NY8$Expected[X]))
##   c(res, res[1]/res[2])
## })
## 


###################################################
### code chunk number 24: jss_paper.Rnw:665-669
###################################################
ny.mm1 <- glmer(Observed ~ offset(log(Expected)) + PCTOWNHOME + 
  PCTAGE65P + PEXPOSURE + (1 | AREANAME),
  data = as(NY8, "data.frame"), family = "poisson")
summary(ny.mm1)


###################################################
### code chunk number 25: jss_paper.Rnw:673-677
###################################################
ny.clmm1 <- DetectClustersModel(NY8st,
  thegrid = as.data.frame(NY8)[idxcl, c("x", "y")], fractpop = 0.15,
  alpha = 0.05, typeCluster = "S", R = NULL, model0 = ny.mm1,
  ClusterSizeContribution = "POP8")


###################################################
### code chunk number 26: jss_paper.Rnw:680-681
###################################################
ny.clmm1


###################################################
### code chunk number 27: jss_paper.Rnw:694-696
###################################################
NY8$CLUSTERMM0 <- get.allknclusters(NY8, ny.clmm0)
NY8$CLUSTERMM1 <- get.allknclusters(NY8, ny.clmm1)


###################################################
### code chunk number 28: jss_paper.Rnw:702-713
###################################################

ny.mapmm0 <- spplot(NY8, c("CLUSTERMM0"), col.regions = c("white", "gray"),
  col = "#4D4D4D", colorkey = FALSE) + 
  layer(lpoints(ny.clmm0$x, ny.clmm0$y, pch = 19))

ny.mapmm1 <- spplot(NY8, c("CLUSTERMM1"), col.regions = c("white", "gray"),
  col = "#4D4D4D", colorkey = FALSE) +
  layer(lpoints(ny.clmm1$x, ny.clmm1$y, pch = 19))

grid.arrange(ny.mapmm0, ny.mapmm1, ncol = 2)



###################################################
### code chunk number 29: jss_paper.Rnw:770-775
###################################################
library("DClusterm")
#library(snowfall)
#library(pscl)

data("Navarre")


###################################################
### code chunk number 30: jss_paper.Rnw:780-782
###################################################
print(spplot(brainnav, "SMR", cuts = 8, 
  col.regions = brewer.pal(9, "Oranges")))


###################################################
### code chunk number 31: jss_paper.Rnw:797-800
###################################################
nav.m0 <- glm(OBSERVED ~ offset(log(EXPECTED)) + 1, family = "poisson",
  data = brainnav)
summary(nav.m0)


###################################################
### code chunk number 32: jss_paper.Rnw:806-809
###################################################
nav.m0q <- glm(OBSERVED ~ offset(log(EXPECTED)) + 1, data = brainnav,
  family = "quasipoisson")
summary(nav.m0q)


###################################################
### code chunk number 33: jss_paper.Rnw:821-825
###################################################
library("pscl")
nav.m0zip <- zeroinfl(OBSERVED ~ offset(log(EXPECTED)) + 1 | 1,
  data = brainnav, dist = "poisson", x = TRUE)
summary(nav.m0zip)


###################################################
### code chunk number 34: jss_paper.Rnw:835-840
###################################################
brainnav$Expected <- brainnav$EXPECTED

brainnavst <- STFDF(as(brainnav, "SpatialPolygons"),
  xts(1,as.Date("1990-01-01")), as(brainnav, "data.frame"),
  endTime = as.POSIXct(strptime(c("1990-01-01"), "%Y-%m-%d"), tz = "GMT"))


###################################################
### code chunk number 35: jss_paper.Rnw:848-851
###################################################
nav.cl0 <- DetectClustersModel(brainnavst, coordinates(brainnav),
  fractpop = 0.25, alpha = 0.05, typeCluster = "S", R = NULL,
  model0 = nav.m0zip, ClusterSizeContribution = "EXPECTED")


###################################################
### code chunk number 36: jss_paper.Rnw:856-857
###################################################
nav.cl0


###################################################
### code chunk number 37: jss_paper.Rnw:875-879
###################################################
nav.clusters <- get.knclusters(brainnav, nav.cl0)
brainnav$CLUSTER <- ""
brainnav$CLUSTER [ nav.clusters[[1]] ] <- "CLUSTER"
brainnav$CLUSTER <- as.factor(brainnav$CLUSTER)


###################################################
### code chunk number 38: jss_paper.Rnw:884-886
###################################################
print(spplot(brainnav, "CLUSTER",  col = "#4D4D4D", 
  col.regions = c("white", "grey")) )


###################################################
### code chunk number 39: jss_paper.Rnw:937-942 (eval = FALSE)
###################################################
## library("DClusterm")
## #debug(DetectClustersModel)
## #debug(glmAndZIP.iscluster)
## #debug(CalcStatsAllClusters)
## #library(snowfall)


###################################################
### code chunk number 40: jss_paper.Rnw:945-946
###################################################
data("brainNM")


###################################################
### code chunk number 41: jss_paper.Rnw:956-958
###################################################
print(stplot(brainst[, , "SMR"], cuts = 8, 
  col.regions = brewer.pal(9, "Oranges")))


###################################################
### code chunk number 42: jss_paper.Rnw:973-976
###################################################
nm.m0 <- glm(Observed ~ offset(log(Expected)) + 1, family = "poisson", 
  data = brainst)
summary(nm.m0)


###################################################
### code chunk number 43: jss_paper.Rnw:984-985
###################################################
NM.coords <- coordinates(brainst@sp)


###################################################
### code chunk number 44: jss_paper.Rnw:995-999
###################################################
nm.cl0 <- DetectClustersModel(brainst, NM.coords,
  minDateUser = "1985-01-01", maxDateUser = "1989-01-01",
  fractpop = 0.15, alpha = 0.05, typeCluster = "ST", R = NULL, 
  model0 = nm.m0, ClusterSizeContribution = "Expected")


###################################################
### code chunk number 45: jss_paper.Rnw:1007-1008
###################################################
nm.cl0[1:5,]


###################################################
### code chunk number 46: jss_paper.Rnw:1019-1020
###################################################
dst <- spDistsN1(pts = NM.coords, pt = losalamos, longlat = TRUE)


###################################################
### code chunk number 47: jss_paper.Rnw:1028-1030
###################################################
nyears <- length(unique(brainst$Year))
brainst$IDLANL <- rep(1 / dst, nyears)


###################################################
### code chunk number 48: jss_paper.Rnw:1035-1038
###################################################
nm.m1 <- glm(Observed ~ offset(log(Expected)) + IDLANL,
  family = "poisson", data = brainst)
summary(nm.m1)


###################################################
### code chunk number 49: jss_paper.Rnw:1047-1051
###################################################
nm.cl1 <- DetectClustersModel(brainst, NM.coords, fractpop = 0.15, 
  alpha = 0.05, minDateUser = "1985-01-01", maxDateUser = "1989-01-01",
  typeCluster = "ST", R = NULL, model0 = nm.m1,
  ClusterSizeContribution = "Expected")


###################################################
### code chunk number 50: jss_paper.Rnw:1060-1061
###################################################
nm.cl1[1:5,]


###################################################
### code chunk number 51: jss_paper.Rnw:1072-1075
###################################################
stcl <- get.stclusters(brainst, nm.cl0)
brainst$CLUSTER <- ""
brainst$CLUSTER[ stcl[[1]] ] <- "CLUSTER"


###################################################
### code chunk number 52: jss_paper.Rnw:1080-1082
###################################################
print(stplot(brainst[, , "CLUSTER"], at = c(0, 0.5, 1.5), col = "#4D4D4D",
  col.regions = c("white", "gray")))


