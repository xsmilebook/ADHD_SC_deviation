library(tidyr)
library(mgcv)
library(gratia)
library(tidyverse)
library(lme4)
library(gamm4)
library(pbkrtest)
#library(car)

pbootint <- function(modelobj){
  numsims <- 100
  set.seed(925)
  df <- modelobj$gam$model
  thisResp <- as.character(modelobj$gam$terms[[2]])
  f1 <- modelobj$gam$formula
  theseVars <- attr(terms(f1),"term.labels")
  
  if (sum(str_detect(theseVars, "^ti"))==1){
    intindx <- which(str_detect(theseVars, "^ti"))
    f2 <- reformulate(theseVars[-intindx],response = thisResp)
  }else{
    f2 <- reformulate(theseVars[2:(length(theseVars))],response = thisResp)
  }
  
  g1 <- gam(f1,data = df)
  g2 <- gam(f2,data = df)
  
  mat1 <- model.matrix(g1)
  mat2 <- model.matrix(g2)
  
  subID<-df$subID
  y <- df[,thisResp]
  
  m1 <- lmer(y ~ -1 + mat1 + (1|subID))
  m2 <- lmer(y ~ -1 + mat2 + (1|subID))
  refdist <- PBrefdist(m1, m2, nsim=numsims)
  pb <- PBmodcomp(m1, m2, ref = refdist, details = 1)
  int_pval <- pb$test["PBtest","p.value"]
  int_chisq <- pb$test["PBtest","stat"]
  
  return(c(int_pval, int_chisq))
}

#### Fit GAMM containing a smooth by 
## Use tensor product for interaction term.

gamm.predict.tensorinteraction <- function(region, dataname, smooth_var, VOI, int_var, covariates=NA, knots, set_fx = FALSE){
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
    modelformula <- as.formula(sprintf("%1$s ~ ti(%2$s, %3$s, k=%4$s, fx=%5$s) + s(%2$s, k=%4$s, fx=%5$s) + s(%3$s, k=%4$s, fx=%5$s) + %6$s", region, VOI, smooth_var, knots, set_fx, int_var))
    modelformula.null <- as.formula(sprintf("%1$s ~ s(%2$s, k=%4$s, fx=%5$s) + s(%3$s, k=%4$s, fx=%5$s) + %6$s", region, VOI, smooth_var, knots, set_fx, int_var))
  }else{
    modelformula <- as.formula(sprintf("%1$s ~ ti(%2$s, %3$s, k=%4$s, fx=%5$s) + s(%2$s, k=%4$s, fx=%5$s) + s(%3$s, k=%4$s, fx=%5$s) + %6$s + %7$s", region, VOI, smooth_var, knots, set_fx, int_var, covariates))
    modelformula.null <- as.formula(sprintf("%1$s ~ s(%2$s, k=%4$s, fx=%5$s) + s(%3$s, k=%4$s, fx=%5$s) + %6$s + %7$s", region, VOI, smooth_var, knots, set_fx, int_var, covariates))
  }
  
  gamm.model <- gamm(modelformula, random=list(subID = ~1), REML=TRUE, data = gam.data)
  gamm.results <- summary(gamm.model$gam)
  
  gamm.model.null <- gamm(modelformula.null, random=list(subID = ~1), REML=TRUE, data = gam.data)
  gamm.null.results <- summary(gamm.model.null$gam)
  
  # The full model contains the smooth by factor interaction term, while the null model contains 
  # the smooth term of VOI and main effect of factor variable separately.
  
  
  # stats
  # Interaction effect
  bootstrap.result.int <- pbootint(gamm.model)
  bootstrap_pvalue.int <- bootstrap.result.int[1]
  bootstrap_chisq.int <- bootstrap.result.int[2]
  gam.Int.pvalue <- gamm.results$s.table[grep(x=rownames(gamm.results$s.table),pattern = "^ti")[1],4]
  gam.Int.Fvalue <- gamm.results$s.table[grep(x=rownames(gamm.results$s.table),pattern = "^ti")[1],3]
  
  
  stats.reults<-cbind(parcel, int_var, VOI, bootstrap_pvalue.int, bootstrap_chisq.int, gam.Int.pvalue, gam.Int.Fvalue)
  
  return(stats.reults)
  
}
