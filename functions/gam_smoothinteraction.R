library(tidyr)
library(mgcv)
library(gratia)
library(tidyverse)
library(mgcv)
#library(car)


#### Fit GAMM containing a smooth by 
## The interest variable is continuous, and the interaction variable is factor.

gam.predict.smoothinteraction <- function(region, dataname, smooth_var, VOI, int_var, covariates=NA, knots, set_fx = FALSE){
  gam.data <- get(dataname)
  # tmp<-gam.data[,region]
  # outlierindx<-which(tmp<mean(tmp)-3*sd(tmp) | tmp>mean(tmp)+3*sd(tmp))
  # if (length(outlierindx)>0){
  #   gam.data<-gam.data[-outlierindx, ]
  # }
  parcel <- region
  gam.data <- gam.data %>% drop_na(region)
  
  #Fit the gam
  if (is.na(covariates)){
    modelformula <- as.formula(sprintf("%1$s ~ s(%2$s, by=%6$s, k=%4$s, fx=%5$s) + s(%3$s, k=%4$s, fx=%5$s) + %6$s", region, VOI, smooth_var, knots, set_fx, int_var))
    modelformula.null <- as.formula(sprintf("%1$s ~ s(%2$s, k=%4$s, fx=%5$s) + s(%3$s, k=%4$s, fx=%5$s) + %6$s", region, VOI, smooth_var, knots, set_fx, int_var))
    modelformula.null2 <- as.formula(sprintf("%1$s ~ s(%2$s, k=%3$s, fx=%4$s) + %5$s", region, smooth_var, knots, set_fx, int_var))
  }else{
    modelformula <- as.formula(sprintf("%1$s ~ s(%2$s, by=%6$s, k=%4$s, fx=%5$s) + s(%3$s, k=%4$s, fx=%5$s) + %6$s + %7$s", region, VOI, smooth_var, knots, set_fx, int_var, covariates))
    modelformula.null <- as.formula(sprintf("%1$s ~ s(%2$s, k=%4$s, fx=%5$s) + s(%3$s, k=%4$s, fx=%5$s) + %6$s + %7$s", region, VOI, smooth_var, knots, set_fx, int_var, covariates))
    modelformula.null2 <- as.formula(sprintf("%1$s ~ s(%2$s, k=%3$s, fx=%4$s) + %5$s + %6$s", region, smooth_var, knots, set_fx, int_var, covariates))
  }

  gam.model <- gam(modelformula, method="REML", data = gam.data)
  gam.results <- summary(gam.model)
  
  gam.model.null <- gam(modelformula.null, method="REML", data = gam.data)
  gam.null.results <- summary(gam.model.null)
  
  gam.model.null2 <- gam(modelformula.null2, method="REML", data = gam.data)
  gam.null2.results <- summary(gam.model.null2)
  
  # The full model contains the smooth by factor interaction term, while the null model contains 
  # the smooth term of VOI and main effect of factor variable separately.
  
  
  # stats
  # Interaction effect
  anova.result.int <- anova(gam.model.null, gam.model, test="Chisq")
  anova_pvalue.int <- anova.result.int$`Pr(>Chi)`[2]
  anova_chisq.int <- anova.result.int$Deviance[2]
  gam.VOI.pvalue.level1 <- gam.results$s.table[grep(x=rownames(gam.results$s.table),pattern = ":")[1],4]
  gam.VOI.Fvalue.level1 <- gam.results$s.table[grep(x=rownames(gam.results$s.table),pattern = ":")[1],3]
  gam.VOI.pvalue.level2 <- gam.results$s.table[grep(x=rownames(gam.results$s.table),pattern = ":")[2],4]
  gam.VOI.Fvalue.level2 <- gam.results$s.table[grep(x=rownames(gam.results$s.table),pattern = ":")[2],3]
  
  
  # VOI effect
  anova.result.VOI <- anova(gam.model.null2, gam.model.null, test="Chisq")
  anova_pvalue.VOI <- anova.result.VOI$`Pr(>Chi)`[2]
  anova_chisq.VOI <- anova.result.VOI$Deviance[2]
  gam.VOI.pvalue <- gam.null.results$s.table[grep(x=rownames(gam.null.results$s.table),pattern = VOI),4]
  gam.VOI.Fvalue <- gam.null.results$s.table[grep(x=rownames(gam.null.results$s.table),pattern = VOI),3]
  
  
  stats.reults<-cbind(parcel, int_var, VOI, anova_pvalue.int, anova_chisq.int, gam.VOI.pvalue.level1, gam.VOI.Fvalue.level1,
                      gam.VOI.pvalue.level2, gam.VOI.Fvalue.level2, anova_pvalue.VOI, anova_chisq.VOI, gam.VOI.pvalue, gam.VOI.Fvalue)
  
  return(stats.reults)
  
}
