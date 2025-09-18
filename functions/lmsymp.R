library(tidyr)
library(psych)
library(dplyr)

#### Fit the linear effects of a continuous variable with linear covariates ####
## Function to fit evaluate the linear effects of a continuous variable on dependent variables per edge.
lm.fit.symp <- function(region, dataname, symp_var, covariates=NA, if_residual=F, IDvar="subID"){
  
  #Fit the gam
  lm.data <- get(dataname)
  symp<-lm.data[[symp_var]]
  if (is.numeric(symp)){
    outlierindx1<-which(symp<mean(symp,na.rm = T)-3*sd(symp,na.rm = T) | symp>mean(symp,na.rm = T)+3*sd(symp,na.rm = T))
    lm.data[outlierindx1 ,symp_var]<-NA

  }
  NonNANIndex <- which(!is.na(lm.data[ ,symp_var]) & !is.na(lm.data[ ,region]))

  lm.data <- lm.data[NonNANIndex,]
  tmp<-lm.data[[region]]
  outlierindx<-which(tmp<mean(tmp)-3*sd(tmp) | tmp>mean(tmp)+3*sd(tmp))
  if (length(outlierindx)>0){
    lm.data<-lm.data[-outlierindx, ]
  }

  parcel <- as.character(region)
  if (! is.na(covariates)){
    covariates.var <- unlist(str_split(covariates, "\\+"))
    lm.data <- lm.data %>% drop_na(covariates.var)
    modelformula <- as.formula(sprintf("%s ~ %s + %s",symp_var, region, covariates))
    modelformula.null<-as.formula(sprintf("%s ~ %s",symp_var, covariates))
    lm.model <- lm(modelformula, data = lm.data)
    lm.model.null <- lm(modelformula.null, data = lm.data)
    anova.cov.pvalue<-anova(lm.model.null, lm.model, test="Chisq")$`Pr(>Chi)`[2]
    #LM statistics
    lm.results <- summary(lm.model)
    adj.Rsq <- NA
  }else{
    modelformula <- as.formula(sprintf("%s ~ %s", symp_var,region))
    lm.model <- lm(modelformula, data = lm.data)
    anova.cov.pvalue<-anova(lm.model)$`Pr(>F)`[1]
    #LM statistics
    lm.results <- summary(lm.model)
    ##Full versus reduced model direction-dependent partial R squared
    ### effect size
    adj.Rsq <- lm.results$adj.r.squared
  }
  
  
  
  #F value for the smooth term and LM-based significance of the smooth term
  lm.symp.t <- lm.results$coefficients[2,3]
  lm.symp.pvalue <- lm.results$coefficients[2,4]

  ### effect direction
  if(lm.symp.t < 0){ #if the gam t-value for covariate of interest is less than 0, make the partialRsq negative
    adj.Rsq <- adj.Rsq*-1}
  
  df.new <- data.frame(int=rep(1, times=nrow(lm.data)))
  for (i in 1:length(covariates.var)){
    if (is.numeric(lm.data[[covariates.var[i]]])){
      df.new[[covariates.var[i]]] <- mean(lm.data[[covariates.var[i]]])
    }else if (is.factor(lm.data[[covariates.var[i]]])){
      df.new[[covariates.var[i]]] <- levels(lm.data[[covariates.var[i]]])[2]
    }
  }
  
  #residual correlation
  if (is.na(covariates)){
    res1 <- res2 <- NA
    res1.normtest <- ks.test(lm.data[[region]], "pnorm", mean = mean(lm.data[[region]]), sd = sd(lm.data[[region]]))
    res2.normtest <- ks.test(lm.data[[symp_var]], "pnorm", mean = mean(lm.data[[symp_var]]), sd = sd(lm.data[[symp_var]]))
    if (res1.normtest$p.value > 0.01 & res2.normtest$p.value > 0.01){
      corrmethod <- "pearson"
    }else{
      corrmethod <- "spearman"
    }
    PCorr_Test <- corr.test(lm.data[[region]], lm.data[[symp_var]], method=corrmethod)
    correstimate<-as.numeric(PCorr_Test$r)
    corrp <- as.numeric(PCorr_Test$p)
  }else{
    res1 <- residuals(lm.model.null)+predict(lm.model.null, newdata = df.new) # symptom
    modelformula.region<-as.formula(sprintf("%s ~ %s",region, covariates))
    lm.model.region <- lm(modelformula.region, data = lm.data)
    res2 <- residuals(lm.model.region)+predict(lm.model.region, newdata = df.new) # region
    
    res1.normtest <- ks.test(res1, "pnorm", mean = mean(res1), sd = sd(res1))
    res2.normtest <- ks.test(res2, "pnorm", mean = mean(res2), sd = sd(res2))
    if (res1.normtest$p.value > 0.01 & res2.normtest$p.value > 0.01){
      corrmethod <- "pearson"
    }else{
      corrmethod <- "spearman"
    }
    PCorr_Test <- corr.test(res1, res2, method=corrmethod)
    correstimate<-as.numeric(PCorr_Test$r)
    corrp <- as.numeric(PCorr_Test$p)
  }
  
  
  stats.results <- cbind(parcel, symp_var, lm.symp.t, lm.symp.pvalue,anova.cov.pvalue, adj.Rsq, correstimate, corrp)
  
  residual.sum <- data.frame(SCDelta.res=res2, SymptomDelta.res=res1, ID=lm.data[[IDvar]])
  
  if (if_residual == F){
    return(stats.results)
  }else{
    return(residual.sum)
  }
}





