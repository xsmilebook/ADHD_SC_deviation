library(lavaan)
library(lme4)

parallelmediation <- function(mediator, dataname, out_var, X_var, covariates, boottime=1000, longitudinal=T){
  
  gam.data <- get(dataname)
  covariates <- gsub(" ", "", covariates)
  covariates1 <- unlist(strsplit(covariates, "+", fixed=T))
  
  Data_Y_New <- data.frame(X=gam.data[[X_var]])
  
  for (c in 1:length(covariates1)){
    Data_Y_New[[covariates1[c]]] <- gam.data[[covariates1[c]]]
  }
  
  if (longitudinal == T){
    out.lmer <- lmer(as.formula(sprintf("%s ~ 1 + (1 | subID)", out_var)), data = gam.data)
    Data_Y_New$Y <- residuals(out.lmer)
    
    for (n in 1:length(mediator)){
      med.lmer <- lmer(as.formula(sprintf("%s ~ 1 + (1 | subID)", mediator[n])), data = gam.data)
      Data_Y_New[[paste0("M", n)]] <- residuals(med.lmer)
    }
  }else{
    Data_Y_New$Y <- gam.data[[out_var]]
    for (n in 1:length(mediator)){
      Data_Y_New[[paste0("M", n)]] <- gam.data[[mediator[n]]]
    }
  }
  
  
  directeffect <- '# direct effect
  Y ~ c*X '
  
  mediatoreffect <- '# mediator
  '
  
  indirecteffect <- '# indirect effect (IDE)
  '
  sumindirecteffect <- 'sumIDE := '
  
  totaleffect <- '# total effect
  total := c '
  
  for (n in 1:length(mediator)){
    directeffect <- paste0(directeffect, "+ b", n, "*M", n)
    mediatoreffect <- paste0(mediatoreffect, "\n M", n, " ~ a", n, "*X + ", covariates)
    indirecteffect <- paste0(indirecteffect, "\n M", n, "IDE := a", n, "*b", n)
    if (n == 1){
      sumindirecteffect <- paste0(sumindirecteffect, "(a", n , "*b", n, ")")
    }else{
      sumindirecteffect <- paste0(sumindirecteffect, " + (a", n , "*b", n, ")")
    }
    
    totaleffect <- paste0(totaleffect, "+ (a", n, "*", "b", n, ")")
  }
  
  model <- paste0(directeffect, " + ", covariates, "\n ", mediatoreffect, "\n", indirecteffect, "\n", sumindirecteffect, "\n ", totaleffect)
  
  
  # model <- ' # direct effect
  # Y ~ c*X + b1*M1 + sex + meanFD
  # #+ b2*M2 + sex
  # # mediator
  # M1 ~ a1*X + sex
  # M2 ~ a2*X + sex
  # 
  # # indirect effect (IDE)
  # M1IDE  := a1*b1
  # M2IDE  := a2*b2
  # sumIDE := (a1*b1) + (a2*b2)
  # 
  # # total effect
  # total := c + (a1*b1) + (a2*b2)
  # '
  
  set.seed(1234)
  LVSV_fit <- sem(model, data = Data_Y_New, test="bootstrap", bootstrap=boottime)
  summary(LVSV_fit)
  assign("LVSV_fit", LVSV_fit, envir = .GlobalEnv)
  
  LVSV_boot <- parameterEstimates(LVSV_fit, boot.ci.type = "perc")
  
  ### Select features
  LVSV_boot.IDE <- LVSV_boot %>% filter(str_detect(lhs, "IDE"), str_detect(lhs, "M"))
  
  selectM <- LVSV_boot.IDE$lhs[LVSV_boot.IDE$pvalue < 0.05]
  if (length(selectM) > 0){
    selectM <- gsub("IDE", "", selectM)
    selectM <- gsub("M", "", selectM)
    mediator2 <- mediator[as.numeric(selectM)]
    
    Data_Y_New <- data.frame(X=gam.data[[X_var]])
    
    for (c in 1:length(covariates1)){
      Data_Y_New[[covariates1[c]]] <- gam.data[[covariates1[c]]]
    }
    
    if (longitudinal == T){
      out.lmer <- lmer(as.formula(sprintf("%s ~ 1 + (1 | subID)", out_var)), data = gam.data)
      Data_Y_New$Y <- residuals(out.lmer)
      
      for (n in 1:length(mediator2)){
        med.lmer <- lmer(as.formula(sprintf("%s ~ 1 + (1 | subID)", mediator2[n])), data = gam.data)
        Data_Y_New[[paste0("M", n)]] <- residuals(med.lmer)
      }
    }else{
      Data_Y_New$Y <- gam.data[[out_var]]
      for (n in 1:length(mediator2)){
        Data_Y_New[[paste0("M", n)]] <- gam.data[[mediator2[n]]]
      }
    }
    
    
    directeffect <- '# direct effect
  Y ~ c*X '
    
    mediatoreffect <- '# mediator
  '
    
    indirecteffect <- '# indirect effect (IDE)
  '
    sumindirecteffect <- 'sumIDE := '
    
    totaleffect <- '# total effect
  total := c '
    
    for (n in 1:length(mediator2)){
      directeffect <- paste0(directeffect, "+ b", n, "*M", n)
      mediatoreffect <- paste0(mediatoreffect, "\n M", n, " ~ a", n, "*X + ", covariates)
      indirecteffect <- paste0(indirecteffect, "\n M", n, "IDE := a", n, "*b", n)
      if (n == 1){
        sumindirecteffect <- paste0(sumindirecteffect, "(a", n , "*b", n, ")")
      }else{
        sumindirecteffect <- paste0(sumindirecteffect, " + (a", n , "*b", n, ")")
      }
      
      totaleffect <- paste0(totaleffect, "+ (a", n, "*", "b", n, ")")
    }
    
    model <- paste0(directeffect, " + ", covariates, "\n ", mediatoreffect, "\n", indirecteffect, "\n", sumindirecteffect, "\n ", totaleffect)
    
    set.seed(1234)
    if (length(mediator) == length(mediator2)){
      final_mod <- "LVSV_fit"
      LVSV_fit.final <- LVSV_fit
      LVSV_boot.final <- LVSV_boot
      
    }else{
      LVSV_fit2 <- sem(model, data = Data_Y_New, test="bootstrap", bootstrap=boottime)
      summary(LVSV_fit2)
      assign("LVSV_fit2", LVSV_fit2, envir = .GlobalEnv)
      
      LVSV_boot2 <- parameterEstimates(LVSV_fit2, boot.ci.type = "perc")
      
      ### Full model VS Reduced model
      anova.result <- anova(LVSV_fit, LVSV_fit2)
      anova.P <- anova(LVSV_fit, LVSV_fit2)$`Pr(>Chisq)`[2]
      if (anova.P < 0.05){
        final_mod <- rownames(anova.result)[which.min(anova.result$BIC)]
        
      }
      
      LVSV_fit.final <- get(final_mod)
      LVSV_boot.final <- parameterEstimates(LVSV_fit.final, boot.ci.type = "perc")
    }
    
    
    
    if (str_detect(final_mod, "2")){
      mediator.f <- mediator2
    }else{
      mediator.f <- mediator
    }
    
    result.sum <- list(mediation=LVSV_boot.final, mediator=mediator.f)
    
    return(result.sum)
  }else{
    
    print("IDE for no mediator is significant!")
    
    result.sum <- list(mediation=LVSV_boot, mediator=mediator)
    
    return(result.sum)
  }
  
}

