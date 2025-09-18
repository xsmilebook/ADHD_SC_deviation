library(lme4)
library(psych)
library(tidyverse)
library(mediation)
library(Formula)
library(stats)


#### FIT GAM SMOOTH ####
##Function to fit a GAM (measure ~ s(X_var, k = knots, fx = set_fx) + covariates)) per each region in atlas and save out statistics and derivative-based characteristics
linearmediation <- function(mediator, dataname, out_var, X_var, covariates){
  
  #Fit the gam
  gam.data <- get(dataname)
  NonNANIndex <- which(!is.na(gam.data[ ,out_var]))
  gam.data <- gam.data[NonNANIndex,]
  # tmp<-gam.data[,mediator]
  # outlierindx<-which(tmp<mean(tmp)-3*sd(tmp), tmp>mean(tmp)+3*sd(tmp))
  # if (length(outlierindx)>0){
  #   gam.data<-gam.data[-outlierindx, ]
  # }
  # parcel <- as.character(mediator)
  # gam.data$tmp<-gam.data[ ,mediator]
  #fit models
  #X=X_var, M=mediator, Y=out_var
  
  # med.lmer <- lmer(as.formula(sprintf("%s ~ 1 + (1 | subID)", mediator)), data = gam.data)
  # out.lmer <- lmer(as.formula(sprintf("%s ~ 1 + (1 | subID)", out_var)), data = gam.data)
  # 
  # gam.data$mediator.res <- residuals(med.lmer)
  # gam.data$out_var.res <- residuals(out.lmer)
  # gam.data$agefactor <- (gam.data$age < 12)
  
  #modelformula_med <- as.formula(sprintf("%s ~s(%s, k = %s, fx = %s) + %s",mediator, X_var, knots, set_fx, covariates))
  #med.fit <- gam(mediator.res ~ s(age, k = 3, fx = TRUE) + sex + agefactor, method="REML", data = gam.data)
  med.fit <- lme4::lmer(as.formula(sprintf("%s ~ age + sex + (1|subID)", mediator)), data = gam.data)

  #modelformula_out <- as.formula(sprintf("%s ~ %s +s(%s, k = %s, fx = %s) + %s",out_var, mediator,X_var, knots, set_fx, covariates))
  #out.fit <- gam(out_var.res ~ mediator.res + s(age, k = 3, fx = F) + sex +agefactor, method="REML", data = gam.data)
  out.fit <- lme4::lmer(as.formula(sprintf("%s ~ %s + age + sex + (1|subID)", out_var, mediator)), data = gam.data)
  
  med.out <- mediation::mediate(med.fit, out.fit, treat = X_var, 
                                mediator = mediator, sims = 1000, boot=F)
  
  
  medparameter<-summary(med.out)
  indirect.Beta <- medparameter$d.avg
  indirect.lwth95CI <- medparameter$d.avg.ci[1]
  indirect.upth95CI <- medparameter$d.avg.ci[2]
  indirect.P <- medparameter$d.avg.p
  direct.Beta <- medparameter$z.avg
  direct.lwth95CI <- medparameter$z.avg.ci[1]
  direct.upth95CI <- medparameter$z.avg.ci[2]
  direct.P <- medparameter$z.avg.p
  total.Beta <- medparameter$tau.coef
  total.lwth95CI <- medparameter$tau.ci[1]
  total.upth95CI <- medparameter$tau.ci[2]
  total.P <- medparameter$tau.p
  prop.Med <- medparameter$n.avg
  prop.lwth95CI <- medparameter$n.avg.ci[1]
  prop.upth95CI <- medparameter$n.avg.ci[2]
  prop.P <- medparameter$n.avg.p
  
  stats.results <- cbind(mediator, out_var, indirect.Beta, indirect.lwth95CI, indirect.upth95CI,
                         indirect.P, direct.Beta, direct.lwth95CI, direct.upth95CI, direct.P, total.Beta,
                         total.lwth95CI, total.upth95CI, total.P, prop.Med, prop.lwth95CI, prop.upth95CI,
                         prop.P)
  return(stats.results)
}





