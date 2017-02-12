setwd("C:/GSOC11/dclusterm/pkg/")

library(spdep)
library(RColorBrewer)
library(splancs)
library(spacetime)
library(DCluster)
library(pscl)
library(lme4)
library(snowfall)

# Functions
source("R/Functions1PAU.R")
source("R/Functions2PAU.R")
source("R/glm.isclusterPAU.R")
source("R/knutils.R")

AddCoordinatesToDataframe<-function(dataframe){
xy<-coordinates(dataframe)
dataframe$x<-xy[,1]
dataframe$y<-xy[,2]
names(xy)<-c("x", "y")
return(dataframe)
}


# Datasets

load("data/NY8.RData")
# Brain cancer in New Mexico, 1973-1991
load("data/nmf.RData")
# Brain cancer in Navarra, Spain
load("data/Navarra.RData")



# Example 1 (NY8)
NY8$Exp<- NY8$POP8 * sum(NY8$Cases)/sum(NY8$POP8)
NY8$SMR<-NY8$Cases/NY8$Exp

# STFDF object: sp, time and mydata
NY8<-AddCoordinatesToDataframe(NY8)
sp<-SpatialPoints(cbind(x = NY8$x, y = NY8$y))
# time must be ordered
time<-as.POSIXct(strptime(c("1972-01-01"), "%Y-%m-%d"), tz = "GMT")
mydata<-data.frame(
Observed = c(NY8$Cases),
Expected = c(NY8$Exp),
SMR      = c(NY8$SMR))
stfdf = STFDF(sp, time, mydata)


# Example 2 (Navarra)
# STFDF object: sp, time and mydata
brain<-AddCoordinatesToDataframe(brain)
sp<-SpatialPoints(cbind(x = brain$x, y = brain$y))
# time must be ordered
time<-as.POSIXct(strptime(c("1990-01-01"), "%Y-%m-%d"), tz = "GMT")
mydata<-data.frame(
Observed = c(brain$OBSERVED),
Expected = c(brain$EXPECTED),
SMR      = c(brain$SMR))
stfdf = STFDF(sp, time, mydata)


# Example 3 (New Mexico)
# STFDF object: sp, time and mydata
nmf<-AddCoordinatesToDataframe(nmf)
sp<-SpatialPoints(cbind(x = nmf$x, y = nmf$y))
# time must be ordered
time<-as.POSIXct(strptime(c("1973-01-01", "1974-01-01", "1975-01-01", "1976-01-01", "1977-01-01", "1978-01-01",
"1979-01-01", "1980-01-01", "1981-01-01", "1982-01-01", "1983-01-01", "1984-01-01", "1985-01-01", "1986-01-01",
"1987-01-01", "1988-01-01", "1989-01-01", "1990-01-01", "1991-01-01"), "%Y-%m-%d"), tz = "GMT")
mydata<-data.frame(
Observed = c(nmf$Observed73, nmf$Observed74, nmf$Observed75, nmf$Observed76, nmf$Observed77, nmf$Observed78,
nmf$Observed79, nmf$Observed80, nmf$Observed81, nmf$Observed82, nmf$Observed83, nmf$Observed84, nmf$Observed85,
nmf$Observed86, nmf$Observed87, nmf$Observed88, nmf$Observed89, nmf$Observed90, nmf$Observed91),
Expected = c(nmf$Expected73, nmf$Expected74, nmf$Expected75, nmf$Expected76, nmf$Expected77, nmf$Expected78,
nmf$Expected79, nmf$Expected80, nmf$Expected81, nmf$Expected82, nmf$Expected83, nmf$Expected84, nmf$Expected85,
nmf$Expected86, nmf$Expected87, nmf$Expected88, nmf$Expected89, nmf$Expected90, nmf$Expected91),
SMR      = c(nmf$SMR73, nmf$SMR74, nmf$SMR75, nmf$SMR76, nmf$SMR77, nmf$SMR78, nmf$SMR79,nmf$SMR80, nmf$SMR81,
nmf$SMR82, nmf$SMR83, nmf$SMR84, nmf$SMR85, nmf$SMR86, nmf$SMR87, nmf$SMR88, nmf$SMR89, nmf$SMR90, nmf$SMR91))
stfdf = STFDF(sp, time, mydata)



# Example 4 (SIDS)
nc_file <- system.file("etc/shapes/sids.shp", package = "spdep")[1]
nc <- readShapeSpatial(nc_file, ID = "FIPSNO", proj4string = CRS("+proj=longlat +datum=NAD27"))
nc$EXPECTED<-nc$BIR74*sum(nc$SID74)/sum(nc$BIR74)
nc$SIR<-nc$SID74/nc$EXPECTED

# STFDF object: sp, time and mydata
nc<-AddCoordinatesToDataframe(nc)
sp<-SpatialPoints(cbind(x = nc$x, y = nc$y))
# time must be ordered
time<-as.POSIXct(strptime(c("1990-01-01"), "%Y-%m-%d"), tz = "GMT")
mydata<-data.frame(
Observed = c(nc$SID74),
Expected = c(nc$EXPECTED),
SMR      = c(nc$SIR))
stfdf = STFDF(sp, time, mydata)











# typeCluster="ST" (Spatio-temporal) or "S" (Spatial)
# modelType="glm", "zeroinfl" or "glmer"
# R=NULL (p-value calculated with 1-pchisq(2*statsAllClusters$statistic[i], 1)) or
# R=number (p-value calculated with Monte Carlo)
# if modelType="zip", stfdf$Observed<-round(stfdf$Observed)


# Call method to detect clusters

statsAllClusters<-DetectClustersModel(stfdf=stfdf, thegrid=as.data.frame(stfdf@sp), radius=Inf, step=NULL, fractpop=0.15, alpha=0.05,
typeCluster="S", minDateUser=time(stfdf@time)[1], maxDateUser=time(stfdf@time)[1], R=NULL, numCPUS=NULL,
modelType="glm",
modelFormula="Observed ~ offset(log(Expected)) ",
modelFamilyGlmGlmer=poisson(link = "log"),
modelDistZeroinfl="poisson", modelLinkZeroinfl="logit")

statsAllClusters

# Select clusters that do not overlap
statsAllClustersNoOverlap<-SelectStatsAllClustersNoOverlap(stfdf, statsAllClusters)
statsAllClustersNoOverlap

# Plot detected clusters
colors<-brewer.pal(12, "Set3")
map<-nc
PlotClustersNoOverlap(statsAllClustersNoOverlap, colors, map)


