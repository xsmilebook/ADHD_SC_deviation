#### Fit GAM containing interaction term of continuous by categorical variables, which does not contain the smooth var ####
## discrete interaction covariate
gam.linearvar.predict.interaction <- function(region, dataname, smooth_var, VOI, int_var, covariates=NA, knots, set_fx = FALSE){
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
    modelformula <- as.formula(sprintf("%1$s ~ %2$s * %6$s + s(%3$s, k=%4$s, fx=%5$s)", region, VOI, smooth_var, knots, set_fx, int_var))
    modelformula.null <- as.formula(sprintf("%1$s ~ %2$s + s(%3$s, k=%4$s, fx=%5$s) + %6$s", region, VOI, smooth_var, knots, set_fx, int_var))
    modelformula.null2 <- as.formula(sprintf("%1$s ~ s(%2$s, k=%3$s, fx=%4$s) + %5$s", region, smooth_var, knots, set_fx, int_var))
  }else{
    modelformula <- as.formula(sprintf("%1$s ~ %2$s * %6$s + s(%3$s, k=%4$s, fx=%5$s) + %7$s", region, VOI, smooth_var, knots, set_fx, int_var, covariates))
    modelformula.null <- as.formula(sprintf("%1$s ~ %2$s + s(%3$s, k=%4$s, fx=%5$s) + %6$s + %7$s", region, VOI, smooth_var, knots, set_fx, int_var, covariates))
    modelformula.null2 <- as.formula(sprintf("%1$s ~ s(%2$s, k=%3$s, fx=%4$s) + %5$s + %6$s", region, smooth_var, knots, set_fx, int_var, covariates))
  }
  
  gam.model <- gam(modelformula, method="REML", data = gam.data)
  gam.results <- summary(gam.model)
  
  gam.model.null <- gam(modelformula.null, method="REML", data = gam.data)
  gam.null.results <- summary(gam.model.null)
  
  gam.model.null2 <- gam(modelformula.null2, method="REML", data = gam.data)
  gam.null2.results <- summary(gam.model.null2)
  
  # stats
  # bootstrap.result.int <- anovaPB(gam.model.null,gam.model, n.sim = 1000,test='Chisq', ncpus=1)
  # bootstrap_pvalue <- bootstrap.result.int[1]
  # bootstrap_chisq <- bootstrap.result.int[2]
  bootstrap.result.int <- anova(gam.model.null,gam.model,test='Chisq')
  bootstrap_pvalue <- bootstrap.result.int$`Pr(>Chi)`[2]
  bootstrap_chisq <- bootstrap.result.int$Deviance[2]

  gam.int.pvalue <- gam.results$p.table[grep(x=rownames(gam.results$p.table),pattern = ":"),4]
  gam.int.T <- gam.results$p.table[grep(x=rownames(gam.results$p.table),pattern = ":"),3]
  # interaction effect size
  sse.model <- sum((gam.model$y - gam.model$fitted.values)^2)
  sse.nullmodel <- sum((gam.model.null$y - gam.model.null$fitted.values)^2)
  IntpartialRsq <- (sse.nullmodel - sse.model)/sse.nullmodel
  
  
  # VOI effect
  bootstrap.result.VOI <- anova(gam.model.null2,gam.model.null,test='Chisq')
  bootstrap_pvalue.VOI <- bootstrap.result.VOI$`Pr(>Chi)`[2]
  bootstrap_chisq.VOI <- bootstrap.result.VOI$Deviance[2]
  gam.VOI.pvalue <- gam.null.results$p.table[grep(x=rownames(gam.null.results$p.table),pattern = VOI),4]
  gam.VOI.T <- gam.null.results$p.table[grep(x=rownames(gam.null.results$p.table),pattern = VOI),3]
  
  
  stats.reults<-cbind(parcel, VOI, int_var, bootstrap_pvalue, bootstrap_chisq, gam.int.pvalue, gam.int.T, IntpartialRsq,
                      bootstrap_pvalue.VOI, bootstrap_chisq.VOI, gam.VOI.pvalue, gam.VOI.T)
  
  return(stats.reults)
  
}