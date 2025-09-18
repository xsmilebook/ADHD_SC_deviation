library(tidyr)
library(mgcv)
library(gratia)
library(tidyverse)
library(ecostats)

#### Interaction analysis for linear effects (cognition) by  a categorical variable####
##Function to predict fitted values of a region for a each level of a categorical variable, using a varying coefficients linear-by-category interaction
gam.smooth.predict.interaction <- function(region, dataname, smooth_var, int_var, covariates, knots, set_fx = FALSE, stats_only=FALSE, increments=1000){
  #Fit the gam
  gam.data <- get(dataname)
  gam.data[[region]] <- scale(gam.data[[region]])
  parcel <- region
  
  #Fit the gam
  if (is.na(covariates)){
    modelformula <- as.formula(sprintf("%1$s ~ s(%2$s, by=%5$s, k=%3$s, fx=%4$s)+ %5$s", region, smooth_var, knots, set_fx, int_var))
    modelformula.null <- as.formula(sprintf("%1$s ~ %5$s+s(%2$s, k=%3$s, fx=%4$s)", region, smooth_var, knots, set_fx, int_var))
    modelformula.null.null <- as.formula(sprintf("%1$s ~ s(%2$s, k=%3$s, fx=%4$s)", region, smooth_var, knots, set_fx))
    
  }else{
    modelformula <- as.formula(sprintf("%1$s ~ s(%2$s, by=%5$s, k=%3$s, fx=%4$s) + %5$s + %6$s", region, smooth_var, knots, set_fx, int_var, covariates))
    modelformula.null <- as.formula(sprintf("%1$s ~ %5$s+s(%2$s, k=%3$s, fx=%4$s)+ %6$s", region, smooth_var, knots, set_fx, int_var, covariates))
    modelformula.null.null <- as.formula(sprintf("%1$s ~ s(%2$s, k=%3$s, fx=%4$s)+ %5$s", region, smooth_var, knots, set_fx, covariates))
  }

  gam.model <- gam(modelformula, method = "REML", data = gam.data)
  gam.results <- summary(gam.model)
  
  gam.model.null <- gam(modelformula.null, method = "REML", data = gam.data)
  gam.null.results <- summary(gam.model.null)
  
  gam.model.null.null <- gam(modelformula.null.null, method = "REML", data = gam.data)
  gam.null.null.results <- summary(gam.model.null.null)
  
  ##Full versus reduced model anova p-value
  #anova.int.pvalue <- anova(gam.model.null, gam.model,test='Chisq')$`Pr(>Chi)`[2]
  intresults <- anovaPB(gam.model.null, gam.model, n.sim = 10000,test='Chisq', ncpus=1)
  anova.int.chisq <- intresults$Deviance[2]
  anova.int.pvalue <- intresults$`Pr(>Chi)`[2]
  gam.int.pvalue <- gam.results$s.table[grep(x=rownames(gam.results$s.table),pattern = ":"),"p-value"]
  # interaction effect size
  sse.model <- sum((gam.model$y - gam.model$fitted.values)^2)
  sse.nullmodel <- sum((gam.model.null$y - gam.model.null$fitted.values)^2)
  IntpartialRsq <- (sse.nullmodel - sse.model)/sse.nullmodel
  
  # Disease effects
  T.disease <- gam.null.results$p.table[2:nlevels(gam.data[[int_var]]),3]
  P.disease <- gam.null.results$p.table[2:nlevels(gam.data[[int_var]]),4]
  diseaseresults <- anovaPB(gam.model.null.null, gam.model.null, n.sim = 10000,test='Chisq', ncpus=1)
  anova.disease.chisq <- diseaseresults$Deviance[2]
  anova.disease.pvalue <- diseaseresults$`Pr(>Chi)`[2]
  
  stats.reults<-cbind(parcel, int_var, anova.int.pvalue, anova.int.chisq, gam.int.pvalue, IntpartialRsq, 
                      T.disease, P.disease, anova.disease.pvalue, anova.disease.chisq)

  return(stats.reults)
  
}