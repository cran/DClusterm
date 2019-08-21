### R code from vignette source 'DClusterm.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: DClusterm.Rnw:5-7
###################################################
#JSS style
options(prompt = "R> ", continue = "+", width = 70, useFancyQuotes = FALSE)


###################################################
### code chunk number 2: DClusterm.Rnw:315-318
###################################################
library("DClusterm")

data("NY8")


###################################################
### code chunk number 3: DClusterm.Rnw:329-334
###################################################
NY8$Observed <- round(NY8$Cases)

NY8$Expected  <- NY8$POP8 * sum(NY8$Observed) / sum(NY8$POP8)

NY8$SMR <- NY8$Observed / NY8$Expected


###################################################
### code chunk number 4: DClusterm.Rnw:343-345
###################################################
NY8$x <- coordinates(NY8)[, 1]
NY8$y <- coordinates(NY8)[, 2]


###################################################
### code chunk number 5: DClusterm.Rnw:353-356 (eval = FALSE)
###################################################
## NY8st <- STFDF(as(NY8, "SpatialPolygons"), xts(1, as.Date("1972-01-01")),
##   NY8@data, endTime = as.POSIXct(strptime(c("1972-01-01"), "%Y-%m-%d"),
##   tz = "GMT"))


###################################################
### code chunk number 6: DClusterm.Rnw:364-383
###################################################
library("RColorBrewer")

#Save plots
p1 <- spplot(NY8, "SMR", cuts = 8, col.regions = brewer.pal(9, "Oranges"),
  main = "Standardized mortality ratio", 
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
### code chunk number 7: DClusterm.Rnw:401-403
###################################################
ny.m0 <- glm(Observed ~ offset(log(Expected)) + 1, family = "poisson",
  data = NY8)


###################################################
### code chunk number 8: DClusterm.Rnw:457-458
###################################################
options(mc.cores = 2)


###################################################
### code chunk number 9: DClusterm.Rnw:466-472
###################################################
idxcl <- c(120, 12, 89, 139, 146)
ny.cl0 <- DetectClustersModel(NY8,
  thegrid = as.data.frame(NY8)[idxcl, c("x", "y")], 
  fractpop = 0.15, alpha = 0.05, radius = Inf, step = NULL,
  typeCluster = "S", R = NULL, model0 = ny.m0,
  ClusterSizeContribution = "POP8")


###################################################
### code chunk number 10: DClusterm.Rnw:488-489
###################################################
ny.cl0


###################################################
### code chunk number 11: DClusterm.Rnw:499-500
###################################################
NY8$CLUSTER0 <- get.allknclusters(NY8, ny.cl0)


###################################################
### code chunk number 12: DClusterm.Rnw:510-513 (eval = FALSE)
###################################################
## cl0.centres <- SpatialPoints(ny.cl0[, 1:2])
## spplot(NY8, "CLUSTER0", col.regions = c("white", "gray"), col = "#4D4D4D",
##   sp.layout = list(list("sp.points", cl0.centres, pch = 19, col = "black")))


###################################################
### code chunk number 13: DClusterm.Rnw:526-529
###################################################
ny.m1 <- glm(Observed ~ offset(log(Expected)) + PCTOWNHOME + PCTAGE65P + 
  PEXPOSURE, family = "poisson", data = NY8)
summary(ny.m1)


###################################################
### code chunk number 14: DClusterm.Rnw:537-541
###################################################
ny.cl1 <- DetectClustersModel(NY8, 
  thegrid = as.data.frame(NY8)[idxcl, c("x", "y")], fractpop = 0.15, 
  alpha = 0.05, typeCluster = "S", R = NULL, model0 = ny.m1,
  ClusterSizeContribution = "POP8")


###################################################
### code chunk number 15: DClusterm.Rnw:544-545
###################################################
ny.cl1


###################################################
### code chunk number 16: DClusterm.Rnw:552-553
###################################################
NY8$CLUSTER1 <- get.allknclusters(NY8, ny.cl1)


###################################################
### code chunk number 17: DClusterm.Rnw:568-593
###################################################
library("gridExtra")
library("latticeExtra")

ny.map0 <- spplot(NY8, c("CLUSTER0"), col.regions = c("white", "gray"),
  col = "#4D4D4D", colorkey = FALSE) +
  layer(lpoints(ny.cl0$x, ny.cl0$y, pch = 19))

ny.map1 <- spplot(NY8, c("CLUSTER1"), col.regions = c("white", "gray"),
  col = "#4D4D4D", colorkey = FALSE) +
  layer(lpoints(ny.cl1$x, ny.cl1$y, pch = 19))

#Syracuse City
syracuse <- which(NY8$AREANAME == "Syracuse city")

ny.map0.s <- spplot(NY8[syracuse, ], c("CLUSTER0"), col.regions = c("white", "gray"),
  col = "#4D4D4D", colorkey = FALSE) +
  layer(lpoints(ny.cl0$x, ny.cl0$y, pch = 19))

ny.map1.s <- spplot(NY8[syracuse, ], c("CLUSTER1"), col.regions = c("white", "gray"),
  col = "#4D4D4D", colorkey = FALSE) +
  layer(lpoints(ny.cl1$x, ny.cl1$y, pch = 19))

grid.arrange(ny.map0, ny.map1, ny.map0.s, ny.map1.s, ncol = 2)




###################################################
### code chunk number 18: DClusterm.Rnw:622-623
###################################################
slimknclusters(NY8, ny.cl1, 3)


###################################################
### code chunk number 19: DClusterm.Rnw:678-680
###################################################
DeanB(ny.m0)
DeanB2(ny.m0)


###################################################
### code chunk number 20: DClusterm.Rnw:689-691
###################################################
DeanB(ny.m1)
DeanB2(ny.m1)


###################################################
### code chunk number 21: DClusterm.Rnw:708-711
###################################################
library("lme4")
ny.mm0 <- glmer(Observed ~ offset(log(Expected)) + (1 | AREANAME), 
  data = as(NY8, "data.frame"), family = "poisson")


###################################################
### code chunk number 22: DClusterm.Rnw:714-718
###################################################
ny.clmm0 <- DetectClustersModel(NY8,
  thegrid = as.data.frame(NY8)[idxcl, c("x", "y")], fractpop = 0.15,
  alpha = 0.05, typeCluster = "S", R = NULL, model0 = ny.mm0,
  ClusterSizeContribution = "POP8")


###################################################
### code chunk number 23: DClusterm.Rnw:721-722
###################################################
ny.clmm0


###################################################
### code chunk number 24: DClusterm.Rnw:734-741 (eval = FALSE)
###################################################
## #Compute SMR in clusters
## cls <- get.knclusters(NY8, ny.cl0)
## lapply(cls, function(X) {
##   res <- c(sum(NY8$Observed[X]), sum(NY8$Expected[X]))
##   c(res, res[1]/res[2])
## })
## 


###################################################
### code chunk number 25: DClusterm.Rnw:748-752
###################################################
ny.mm1 <- glmer(Observed ~ offset(log(Expected)) + PCTOWNHOME + 
  PCTAGE65P + PEXPOSURE + (1 | AREANAME),
  data = as(NY8, "data.frame"), family = "poisson")
#summary(ny.mm1)


###################################################
### code chunk number 26: DClusterm.Rnw:760-764
###################################################
ny.clmm1 <- DetectClustersModel(NY8,
  thegrid = as.data.frame(NY8)[idxcl, c("x", "y")], fractpop = 0.15,
  alpha = 0.05, typeCluster = "S", R = NULL, model0 = ny.mm1,
  ClusterSizeContribution = "POP8")


###################################################
### code chunk number 27: DClusterm.Rnw:767-768
###################################################
ny.clmm1


###################################################
### code chunk number 28: DClusterm.Rnw:781-783
###################################################
NY8$CLUSTERMM0 <- get.allknclusters(NY8, ny.clmm0)
NY8$CLUSTERMM1 <- get.allknclusters(NY8, ny.clmm1)


###################################################
### code chunk number 29: DClusterm.Rnw:789-811
###################################################

ny.mapmm0 <- spplot(NY8, c("CLUSTERMM0"), col.regions = c("white", "gray"),
  col = "#4D4D4D", colorkey = FALSE) + 
  layer(lpoints(ny.clmm0$x, ny.clmm0$y, pch = 19))

ny.mapmm1 <- spplot(NY8, c("CLUSTERMM1"), col.regions = c("white", "gray"),
  col = "#4D4D4D", colorkey = FALSE) +
  layer(lpoints(ny.clmm1$x, ny.clmm1$y, pch = 19))

#Syracuse
ny.mapmm0.s <- spplot(NY8[syracuse, ], c("CLUSTERMM0"), col.regions = c("white", "gray"),
  col = "#4D4D4D", colorkey = FALSE) +
  layer(lpoints(ny.clmm0$x, ny.clmm0$y, pch = 19))

ny.mapmm1.s <- spplot(NY8[syracuse, ], c("CLUSTERMM1"), col.regions = c("white", "gray"),
  col = "#4D4D4D", colorkey = FALSE) +
  layer(lpoints(ny.clmm1$x, ny.clmm1$y, pch = 19))



grid.arrange(ny.mapmm0, ny.mapmm1, ny.mapmm0.s, ny.mapmm1.s, ncol = 2)



###################################################
### code chunk number 30: DClusterm.Rnw:870-873
###################################################
library("DClusterm")

data("Navarre")


###################################################
### code chunk number 31: DClusterm.Rnw:878-880
###################################################
print(spplot(brainnav, "SMR", cuts = 8, 
  col.regions = brewer.pal(9, "Oranges")))


###################################################
### code chunk number 32: DClusterm.Rnw:895-897
###################################################
nav.m0 <- glm(OBSERVED ~ offset(log(EXPECTED)) + 1, family = "poisson",
  data = brainnav)


###################################################
### code chunk number 33: DClusterm.Rnw:903-905
###################################################
nav.m0q <- glm(OBSERVED ~ offset(log(EXPECTED)) + 1, data = brainnav,
  family = "quasipoisson")


###################################################
### code chunk number 34: DClusterm.Rnw:919-923
###################################################
library("pscl")
nav.m0zip <- zeroinfl(OBSERVED ~ offset(log(EXPECTED)) + 1 | 1,
  data = brainnav, dist = "poisson", x = TRUE)
summary(nav.m0zip)


###################################################
### code chunk number 35: DClusterm.Rnw:939-942
###################################################
nav.cl0 <- DetectClustersModel(brainnav, coordinates(brainnav),
  fractpop = 0.25, alpha = 0.05, typeCluster = "S", R = NULL,
  model0 = nav.m0zip, ClusterSizeContribution = "EXPECTED")


###################################################
### code chunk number 36: DClusterm.Rnw:947-948
###################################################
nav.cl0


###################################################
### code chunk number 37: DClusterm.Rnw:966-970
###################################################
nav.clusters <- get.knclusters(brainnav, nav.cl0)
brainnav$CLUSTER <- ""
brainnav$CLUSTER [ nav.clusters[[1]] ] <- "CLUSTER"
brainnav$CLUSTER <- as.factor(brainnav$CLUSTER)


###################################################
### code chunk number 38: DClusterm.Rnw:975-977
###################################################
print(spplot(brainnav, "CLUSTER",  col = "#4D4D4D", 
  col.regions = c("white", "grey")) )


###################################################
### code chunk number 39: DClusterm.Rnw:1027-1028 (eval = FALSE)
###################################################
## library("DClusterm")


###################################################
### code chunk number 40: DClusterm.Rnw:1031-1032
###################################################
data("brainNM")


###################################################
### code chunk number 41: DClusterm.Rnw:1042-1044
###################################################
print(stplot(brainst[, , "SMR"], cuts = 8, 
  col.regions = brewer.pal(9, "Oranges")))


###################################################
### code chunk number 42: DClusterm.Rnw:1059-1062
###################################################
nm.m0 <- glm(Observed ~ offset(log(Expected)) + 1, family = "poisson", 
  data = brainst)
summary(nm.m0)


###################################################
### code chunk number 43: DClusterm.Rnw:1070-1071
###################################################
NM.coords <- coordinates(brainst@sp)


###################################################
### code chunk number 44: DClusterm.Rnw:1083-1087
###################################################
nm.cl0 <- DetectClustersModel(brainst, NM.coords,
  minDateUser = "1985-01-01", maxDateUser = "1989-01-01",
  fractpop = 0.15, alpha = 0.05, typeCluster = "ST", R = NULL, 
  model0 = nm.m0, ClusterSizeContribution = "Expected")


###################################################
### code chunk number 45: DClusterm.Rnw:1095-1097
###################################################
nm.cl0.s <- slimknclusters(brainst, nm.cl0)
nm.cl0.s


###################################################
### code chunk number 46: DClusterm.Rnw:1112-1113
###################################################
dst <- spDistsN1(pts = NM.coords, pt = losalamos, longlat = TRUE)


###################################################
### code chunk number 47: DClusterm.Rnw:1121-1123
###################################################
nyears <- length(unique(brainst$Year))
brainst$IDLANL <- rep(1 / dst, nyears)


###################################################
### code chunk number 48: DClusterm.Rnw:1128-1131
###################################################
nm.m1 <- glm(Observed ~ offset(log(Expected)) + IDLANL,
  family = "poisson", data = brainst)
summary(nm.m1)


###################################################
### code chunk number 49: DClusterm.Rnw:1140-1144
###################################################
nm.cl1 <- DetectClustersModel(brainst, NM.coords, fractpop = 0.15, 
  alpha = 0.05, minDateUser = "1985-01-01", maxDateUser = "1989-01-01",
  typeCluster = "ST", R = NULL, model0 = nm.m1,
  ClusterSizeContribution = "Expected")


###################################################
### code chunk number 50: DClusterm.Rnw:1153-1155
###################################################
nm.cl1.s <- slimknclusters(brainst, nm.cl1)
nm.cl1.s


###################################################
### code chunk number 51: DClusterm.Rnw:1167-1170
###################################################
stcl <- get.stclusters(brainst, nm.cl0.s)
brainst$CLUSTER <- ""
brainst$CLUSTER[ stcl[[1]] ] <- "CLUSTER"


###################################################
### code chunk number 52: DClusterm.Rnw:1175-1177
###################################################
print(stplot(brainst[, , "CLUSTER"], at = c(0, 0.5, 1.5), col = "#4D4D4D",
  col.regions = c("white", "gray")))


