##' @title Obtains the cluster with the maximum log-likelihood ratio or minimum DIC
##' of all the clusters with the same center and start and end dates.
##' 
##' @description This function constructs all the clusters with start date equal to
##' minDateCluster, end date equal to maxDateCluster, and with center specified
##' by the first element of idxorder, so that the maximum fraction of the total
##' population inside the cluster is less than fractpop, and the maximum
##' distance to the center is less than radius.
##' For each one of these clusters, the log-likelihood ratio test statistic
##' for comparing the alternative model with the cluster versus the null model
##' of no clusters (if model is glm, glmer or zeroinfl),
##' or the DIC (if model is inla) is calculated.
##' The cluster with maximum value of the log-likelihood ratio or
##' minimum DIC is returned.
##'
##' @param stfdf a spatio-temporal class object containing the data.
##' @param idxorder a permutation of the regions according to their distance to
##' the current center.
##' @param minDateCluster start date of the cluster.
##' @param maxDateCluster end date of the cluster.
##' @param fractpop maximum fraction of the total population inside the cluster.
##' @param model0 Initial model (including covariates).
##' @param ClusterSizeContribution Variable used to check the fraction of the 
##' population at risk in the cluster
##' This can be "glm" for generalized linear models (glm {stats}),
##' "glmer" for generalized linear mixed model (glmer {lme4}),
##' "zeroinfl" for zero-inflated models (zeroinfl {pscl}), or
##' "inla" for generalized linear, generalized linear mixed or zero-inflated models.
##'
##' @return vector containing the size, the start and end dates,
##' the log-likelihood ratio or DIC, the p-value and the risk
##' of the cluster with the maximum log-likelihood ratio or minimum DIC.
##'
glmAndZIP.iscluster <- function(stfdf, idxorder, minDateCluster,
  maxDateCluster, fractpop, model0, ClusterSizeContribution) {

  d0 <- stfdf@data

  mclass <- class(model0)[1]

  if(inherits(model0, "glm") | inherits(model0, "mer") | 
   inherits(model0, "glmerMod")) {
    
    modelFormula <- paste(formula(model0)[c(2, 3)], collapse = "~")
    modelType <- ifelse(mclass == "glm", "glm", "glmer")
    modelFamilyGlmGlmer <- family(model0)
    modelDistZeroinfl <- NULL
    modelLinkZeroinfl <- NULL
  }

  if(inherits(model0, "zeroinfl")) {
    modelFormula <- paste(formula(model0)[c(2, 3)], collapse = "~")
    modelType<-mclass
    modelFamilyGlmGlmer<-NULL
    modelDistZeroinfl<-model0$dist
    modelLinkZeroinfl<-model0$link
  }
  if(inherits(model0, "inla")) {
    modelFormula <- paste(model0$.args$formula[c(2, 3)], collapse = "~")
    modelType <- mclass
    modelFamilyINLA <- model0$.args$family
  }


  # If model is inla, statistic= DIC. If model is not inla statistic=difL
  sizeCluster <- NA
  estadistico <- NA
  pvalue <- NA
  risk <- NA
  idTime <- which( (time(stfdf@time) >= minDateCluster) &
   (time(stfdf@time) <= maxDateCluster))

  # idTime, idSpace: indexes corresponding to the time and locations inside
  # the cluster

  if(length(idxorder) == 0) {
    print('length(idxorder) = 0')
    return(data.frame(size = sizeCluster, minDateCluster = minDateCluster,
     maxDateCluster = maxDateCluster, statistic = estadistico, pvalue = pvalue,
     risk=risk))
  }


  for(i in 1:length(idxorder)) {
    idSpace <- idxorder[1:i]
    d0$CLUSTER <- SetVbleCluster(stfdf, idTime, idSpace)

    # cluster size must be smaller than fractpop of the total
    if((sum(d0$CLUSTER * stfdf[[ClusterSizeContribution]]) - 
     fractpop * sum(stfdf[[ClusterSizeContribution]])) < 0) {


      # Formula to refit the cluster coefficient only
      if(modelType == "glm") {
        newformula <- formula(paste( strsplit(modelFormula, "~")[[1]][1], 
         "~ -1+CLUSTER "))
      }
      if(modelType == "glmer"){
        #newformula O ~ -1 + CLUSTER + random effects
        #need to keep random effects from model0 formula
        formularight<-strsplit(modelFormula, "~")[[1]][2]
        terms<-strsplit(formularight, "[+]")[[1]]
        indexrandomeffects<-grep("[|]",terms)
        reModelFormula<- paste(terms[indexrandomeffects],collapse="+")
        newformula <- formula(paste( strsplit(modelFormula, "~")[[1]][1], "~ -1+CLUSTER +",reModelFormula))
        
        
        #When I fit m1 I need to put offset equal to fitted fixedeffects model0 + offset model0
        #(fitted values without random effects)
        fixedeffectsmodel0<- model.matrix(model0) %*% model0@beta
        #if offset is in parameter offset (eg. offset=Expected)
        if("(offset)"%in%names(model0@frame)){
          offsetmodel0<-log(model0@frame[,"(offset)"])
        }
        #if offset is in the formula (eg. ~ offset(log(Expected)))
        if("offset"%in%substr(names(model0@frame),1,6)){
          offsetmodel0<-model0@frame[which(substr(names(model0@frame),1,6)=="offset")]
        }
        #This will be the offset in m1
        fittedwithoutre <- fixedeffectsmodel0+unlist(offsetmodel0)
      }
      
      if(modelType == "zeroinfl") {
        offzero <- as.vector(model0$x$zero %*% 
         matrix(model0$coefficients$zero, ncol=1))

        if(!is.null(model0$offset$zero)) {
          offzero <- offzero + model0$offset$zero
        }

        offcount <- 
         model0$x$count %*% matrix(model0$coefficients$count, ncol = 1)
        offcount <- as.vector(offcount)

        if(!is.null(model0$offset$count)) {
          offcount <- offcount + model0$offset$count
        }

        newformula <- formula(paste( strsplit(modelFormula, "~")[[1]][1],
         "~ -1 + CLUSTER + offset(offcount) | 1 + offset(offzero)"))
      }

      if(modelType=="inla") {
        #random effects terms in formula
        reModelFormula <- 
         paste(regmatches(modelFormula, gregexpr("(?=f\\().*?(?<=\\))", 
         modelFormula, perl = TRUE))[[1]], collapse = "+")

        if(reModelFormula != ""){
          reModelFormula <- paste("+", reModelFormula)
        }

        newformula <- formula(paste( strsplit(modelFormula, "~")[[1]][1],
         "~ -1 + CLUSTER", reModelFormula))
      }


      switch(modelType,
       glm = {
         m1 <- glm(newformula, data = d0, family = modelFamilyGlmGlmer,
          offset = log(fitted(model0)))
         riskAux <- coef(m1)[1]
         estadisticoAux <- ifelse(riskAux > 0,
          (deviance(model0) - deviance(m1))/2, 0)
       }, 
       glmer = {
         m1 <- lme4::glmer(newformula, data = d0, family = modelFamilyGlmGlmer,
                     offset=fittedwithoutre)
         riskAux <- ((coef(m1))[[1]][, 1])[1]
         estadisticoAux <- ifelse(riskAux > 0, 
          (deviance(model0) - deviance(m1))/2, 0)
       },
       zeroinfl = {
         m1 <- zeroinfl(newformula, data = d0, dist = modelDistZeroinfl,
          link = modelLinkZeroinfl)
         riskAux <- as.numeric(coef(m1)[1])
         estadisticoAux <- ifelse(riskAux > 0, 
          (-2*logLik(model0) + 2*logLik(m1))/2, 0)
       },
       inla = {
         if(is.null(model0$.args$E)) {
           esperados <- model0$summary.fitted.values[, "mean"]
         } else {
           esperados <- exp(model0$summary.linear.predictor[, "mean"] +
            log(model0$.args$E))
         }
         m1 <- INLA::inla(newformula, family = modelFamilyINLA, data = d0,
          E = esperados, control.inla = list(h = 0.01), 
          control.compute = list(dic = TRUE, mlik=TRUE), verbose = FALSE)

         riskAux <- m1$summary.fixed[1]
         estadisticoAux <- ifelse(riskAux > 0, as.numeric(m1$dic$dic), Inf)
         pvalueAux <- computeprob(m1$marginals.fixed$CLUSTER, 0)
       }
      )

      if(modelType != "inla") {
        if(estadisticoAux > estadistico | i==1) {
          sizeCluster <- i
          estadistico <- as.numeric(estadisticoAux)
          risk <- as.numeric(riskAux)
        }
      } else {
        if(estadisticoAux < estadistico | i==1) {
          sizeCluster <- i
          estadistico <- as.numeric(estadisticoAux)
          risk <- as.numeric(riskAux)
          pvalue <- as.numeric(pvalueAux)
        }
      }
    } else {
      return(data.frame(size = sizeCluster, minDateCluster = minDateCluster,
       maxDateCluster = maxDateCluster, statistic = estadistico,
       pvalue = pvalue, risk = risk))
    }
  }

  return(data.frame(size = sizeCluster, minDateCluster = minDateCluster,
   maxDateCluster = maxDateCluster, statistic = estadistico, pvalue = pvalue,
   risk=risk))

}





##' @title Constructs a variable that indicates the locations and times that pertain
##' to a cluster.
##' 
##' @description This function constructs a variable that indicates the locations and times
##' that pertain to a cluster. Each position of the variable is equal to 1 if
##' it corresponds to a location and time inside the cluster, and 0 otherwise.
##' This is one of the explanatory variables used in the glmAndZIP.iscluster
##' function to model the observed cases.
##'
##' @param stfdf spatio-temporal class object containing the data.
##' @param idTime vector with the indexes of the stfdf object corresponding to
##' the time inside the cluster.
##' @param idSpace vector with the indexes of the stfdf object corresponding to
##' the locations inside the cluster.
##'
##' @return vector with 1's or 0's that indicates the locations and times that
##' pertain to a cluster.
##'
SetVbleCluster <- function(stfdf,idTime,idSpace) {
  # vbleCluster, vbleCluster[i] = 1 if position i corresponds to space and time inside the cluster, 0 otherwise

  vbleCluster <- rep(0, nrow(stfdf@data))

  idSubsetTime <- stfdf[, idTime, drop = FALSE]@data$ID
  idSubsetSpace <- stfdf[idSpace, drop = FALSE]@data$ID

  vbleCluster[(stfdf$ID %in% idSubsetTime) & (stfdf$ID %in% idSubsetSpace)] <- 1

  return(vbleCluster)
}
