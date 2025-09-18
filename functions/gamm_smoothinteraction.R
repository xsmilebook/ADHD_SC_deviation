library(tidyr)
library(mgcv)
library(gratia)
library(tidyverse)
library(lme4)
library(gamm4)
library(pbkrtest)
#library(car)

pbootint <- function(modelobj){
  numsims <- 10000
  set.seed(925)
  df <- modelobj$gam$model
  thisResp <- as.character(modelobj$gam$terms[[2]])
  f1 <- modelobj$gam$formula
  theseVars <- attr(terms(f1),"term.labels")
  
  if (sum(str_detect(theseVars, "by"))==1){
    intindx <- which(str_detect(theseVars, "by"))
    intvar <- theseVars[intindx]
    addvar <- paste0(str_split(intvar, ",")[[1]][1], ",", str_split(intvar, ",")[[1]][3], ",",str_split(intvar, ",")[[1]][4])
    
    f2 <- reformulate(c(theseVars[-intindx], addvar),response = thisResp)
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
## The interest variable is continuous, and the interaction variable is factor.

gamm.predict.smoothinteraction <- function(region, dataname, smooth_var, VOI, int_var, covariates=NA, knots, set_fx = FALSE){
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
  }else{
    modelformula <- as.formula(sprintf("%1$s ~ s(%2$s, by=%6$s, k=%4$s, fx=%5$s) + s(%3$s, k=%4$s, fx=%5$s) + %6$s + %7$s", region, VOI, smooth_var, knots, set_fx, int_var, covariates))
    modelformula.null <- as.formula(sprintf("%1$s ~ s(%2$s, k=%4$s, fx=%5$s) + s(%3$s, k=%4$s, fx=%5$s) + %6$s + %7$s", region, VOI, smooth_var, knots, set_fx, int_var, covariates))
  }
  
  gamm.model <- gamm4(modelformula, random=~(1|subID), REML=TRUE, data = gam.data)
  gamm.results <- summary(gamm.model$gam)
  
  gamm.model.null <- gamm4(modelformula.null, random=~(1|subID), REML=TRUE, data = gam.data)
  gamm.null.results <- summary(gamm.model.null$gam)
  
  # The full model contains the smooth by factor interaction term, while the null model contains 
  # the smooth term of VOI and main effect of factor variable separately.
  
  
  # stats
  # Interaction effect
  bootstrap.result.int <- pbootint(gamm.model)
  bootstrap_pvalue.int <- bootstrap.result.int[1]
  bootstrap_chisq.int <- bootstrap.result.int[2]
  gam.VOI.pvalue.level1 <- gamm.results$s.table[grep(x=rownames(gamm.results$s.table),pattern = ":")[1],4]
  gam.VOI.Fvalue.level1 <- gamm.results$s.table[grep(x=rownames(gamm.results$s.table),pattern = ":")[1],3]
  gam.VOI.pvalue.level2 <- gamm.results$s.table[grep(x=rownames(gamm.results$s.table),pattern = ":")[2],4]
  gam.VOI.Fvalue.level2 <- gamm.results$s.table[grep(x=rownames(gamm.results$s.table),pattern = ":")[2],3]

  
  # VOI effect
  bootstrap.result.VOI <- pbootint(gamm.model.null)
  bootstrap_pvalue.VOI <- bootstrap.result.VOI[1]
  bootstrap_chisq.VOI <- bootstrap.result.VOI[2]
  gam.VOI.pvalue <- gamm.null.results$s.table[grep(x=rownames(gamm.null.results$s.table),pattern = VOI),4]
  gam.VOI.Fvalue <- gamm.null.results$s.table[grep(x=rownames(gamm.null.results$s.table),pattern = VOI),3]
  
  
  stats.reults<-cbind(parcel, int_var, VOI, bootstrap_pvalue.int, bootstrap_chisq.int, gam.VOI.pvalue.level1, gam.VOI.Fvalue.level1,
                      gam.VOI.pvalue.level2, gam.VOI.Fvalue.level2, bootstrap_pvalue.VOI, bootstrap_chisq.VOI, gam.VOI.pvalue, gam.VOI.Fvalue)
  
  return(stats.reults)
  
}
