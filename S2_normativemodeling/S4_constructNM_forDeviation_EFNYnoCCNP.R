rm(list=ls())
library(ggplot2)
library(tidyverse)
library(openxlsx)
library(parallel)
library(gamlss)
library(scales)
library(tableone)
library(psych)
# set resolution
Yeoresolution <- 17
if (Yeoresolution == 7){
  Yeoresolution.delLM = 6
}else if (Yeoresolution == 17){
  Yeoresolution.delLM = 15
}
element_num <- Yeoresolution.delLM*(Yeoresolution.delLM+1)/2
# input directory
wd <- getwd()
homepath <- str_split_i(wd, "Normative_model", 1)
demopath <- paste0(homepath, '/Normative_model/demography')
interfileFolder <- paste0(homepath, '/Normative_model/interfileFolder_EFNYnoCCNP')
interfileFolder_ABCD <- paste0(homepath, '/Normative_model/interfileFolder_ABCD')
functionFolder <- paste0(homepath, "/Normative_model/functions")
resultFolder <- paste0(homepath, "/Normative_model/results_EFNYnoCCNP")
functionFolder.SCDev <- paste0(homepath, "/SC_development/Rcode_SCdevelopment/gamfunction")
FigureFolder <- paste0(homepath, '/Normative_model/Figures_EFNYnoCCNP/Yeo', Yeoresolution,'/CV75')

# Load data
SCdata <- readRDS(paste0(interfileFolder, '/SCdata_Yeo',Yeoresolution, '_CV75_sumSCinvnode.sum.msmtcsd.merge.rds'))
SCdata$sex <- as.factor(SCdata$sex)
SCdata$diagnosis <- factor(SCdata$ADHD, levels=c(0,1), labels = c("TD", "ADHD"))
SCdata$visit <- "0"

SCdatapost <- readRDS(paste0(interfileFolder, '/SCdata_Yeo',Yeoresolution, '_CV75_sumSCinvnode.sum.msmtcsd.merge_postmedication.rds'))
SCdatapost$sex <- as.factor(SCdatapost$sex)
SCdatapost$diagnosis <- factor(SCdatapost$ADHD, levels=c(0,1), labels = c("TD", "ADHD"))
SCdatapost$visit <- "1"
SCdatapost$X <- NULL
SCdatapost$age <- SCdatapost$age + 12*7/365


SCdata <- rbind(SCdata, SCdatapost)

# source function
source(paste0(functionFolder, "/Construct_gamlss_set.R"))

SCdata.TD.trainset <- readRDS(paste0(interfileFolder, "/SCdata.allTD_EFNYnoCCNP.SCYeo", element_num, ".rds"))
SCdata.ADHD <- SCdata %>% filter(diagnosis=="ADHD")
SCdata.ADHD <- SCdata.ADHD %>% mutate(WB_SCmean = rowMeans(dplyr::select(.,starts_with("SC."))))
SCdata.ADHD.base <- SCdata.ADHD %>% filter(visit=="0")
SCdata.ADHD.base <- SCdata.ADHD.base %>% filter(WB_SCmean>=mean(WB_SCmean)-3*sd(WB_SCmean), WB_SCmean<=mean(WB_SCmean)+3*sd(WB_SCmean))
# n=482

SCdata.ADHD <- SCdata.ADHD %>% filter(ID %in% SCdata.ADHD.base$ID)

## 1. For ADHD dataset, all TD will be assigned as training set.
########################################

## Fit normative models
SCdata.TD.trainset <- as.data.frame(SCdata.TD.trainset)
SCdata.TD.trainset[,c("sex", "site")] <- lapply(SCdata.TD.trainset[,c("sex", "site")], as.factor)
dataname <- "SCdata.TD.trainset"
smoothterm <- "age"
covariates <- "sex+mean_fd"
randomvar <- "site"
mu.df <- degree <- 2
sigma.df <- 2
distribution.fam <- "BCPEo"
IDvar <- "ID"
quantile.vec <- c(0.025, 0.5, 0.975)
stratify <- c("sex", "site")
if (! file.exists(paste0(interfileFolder, "/GAMLSS_Yeo", Yeoresolution,".TDtraining.All.rds"))){
  mod120.training.sum <- mclapply(1:element_num, function(i){
    dependentvar <- paste0("SC.", i)
    sumlist <- construct_gamlss(dataname, dependentvar, smoothterm, covariates,randomvar, 
                                mu.df, sigma.df, degree, distribution.fam,IDvar, quantile.vec, stratify)
    
    return(sumlist)
  }, mc.cores = 10)
  saveRDS(mod120.training.sum, paste0(interfileFolder, "/GAMLSS_Yeo", Yeoresolution,".TDtraining.All.rds"))
}else{
  mod120.training.sum <- readRDS(paste0(interfileFolder, "/GAMLSS_Yeo", Yeoresolution,".TDtraining.All.rds"))
}

# compute deviations
summary(SCdata.ADHD[,c("sex", "age", "site", "mean_fd", "diagnosis", "visit")])

# compute deviation
## ADHD
SCdata.ADHD$sex <- as.factor(SCdata.ADHD$sex)
SCdata.ADHD_fixed <- SCdata.ADHD
SCdata.ADHD_fixed$mean_fd <- mean(SCdata.ADHD_fixed$mean_fd)
SCdata.ADHD_fixed.F <- SCdata.ADHD_fixed.M <- SCdata.ADHD_fixed
SCdata.ADHD_fixed.F$sex <- factor("F")
SCdata.ADHD_fixed.M$sex <- factor("M")
SCdata.TD.trainset$site <- droplevels(SCdata.TD.trainset$site)
sitelist <- unique(SCdata.TD.trainset$site)

gam.data2 <- SCdata.TD.trainset
gam.data2$visit <- "0"
if (! file.exists(paste0(interfileFolder, "/SCdata.sum75_Yeo", Yeoresolution,".ADHD_deviation_All.rds"))){
  deviations.sum <- list()
  for (n in 1:element_num){
    #i = parcel.num[n]
    i=n
    print(i)
    mod.tmp <- mod120.training.sum[[n]]$mod.tmp
    mu_pred <- predict(mod.tmp, newdata = SCdata.ADHD, what = "mu", type = "response")
    sigma_pred <- predict(mod.tmp, newdata = SCdata.ADHD, what = "sigma", type = "response")
    nu_pred <- predict(mod.tmp, newdata = SCdata.ADHD, what = "nu", type = "response")
    tau_pred <- predict(mod.tmp, newdata = SCdata.ADHD, what = "tau", type = "response")
    
    dependentvar <- paste0("SC.", i)
    deviation.df <- data.frame(ID=SCdata.ADHD$ID)
    observation <- SCdata.ADHD[[dependentvar]]
    centile <- pBCPEo(observation, mu = mu_pred, sigma = sigma_pred, nu = nu_pred, tau = tau_pred)
    # if (sum(centile==1) > 0){
    #   centile[which(centile==1)] <- 0.99999
    # }
    if (sum(centile==0) > 0){
      centile[which(centile==0)] <- 1-0.99999
    }else if (sum(centile==1) > 0){
      centile[which(centile==1)] <- 0.99999
    }
    deviation.df[[paste0("SC.", i, "_centile")]] <- centile
    deviation.df[[paste0("SC.", i, "_deviationZ")]] <- qnorm(centile)
    
    # compute the standard values while controlling for mean_fd and site
    standardvalue.mat <- matrix(NA, nrow(SCdata.ADHD_fixed.F), length(sitelist))
    for (j in 1:length(sitelist)){
      site = sitelist[j]
      SCdata.ADHD_fixed.F$site <- SCdata.ADHD_fixed.M$site <- site
      
      # Female
      mu_pred.fixed.F <- predict(mod.tmp, newdata = SCdata.ADHD_fixed.F, what = "mu", type = "response")
      sigma_pred.fixed.F <- predict(mod.tmp, newdata = SCdata.ADHD_fixed.F, what = "sigma", type = "response")
      nu_pred.fixed.F <- predict(mod.tmp, newdata = SCdata.ADHD_fixed.F, what = "nu", type = "response")
      tau_pred.fixed.F <- predict(mod.tmp, newdata = SCdata.ADHD_fixed.F, what = "tau", type = "response")
      standardvalue.F <- qBCPEo(centile, mu = mu_pred.fixed.F, sigma = sigma_pred.fixed.F, nu = nu_pred.fixed.F, tau = tau_pred.fixed.F)
      
      # Male
      mu_pred.fixed.M <- predict(mod.tmp, newdata = SCdata.ADHD_fixed.M, what = "mu", type = "response")
      sigma_pred.fixed.M <- predict(mod.tmp, newdata = SCdata.ADHD_fixed.M, what = "sigma", type = "response")
      nu_pred.fixed.M <- predict(mod.tmp, newdata = SCdata.ADHD_fixed.M, what = "nu", type = "response")
      tau_pred.fixed.M <- predict(mod.tmp, newdata = SCdata.ADHD_fixed.M, what = "tau", type = "response")
      standardvalue.M <- qBCPEo(centile, mu = mu_pred.fixed.M, sigma = sigma_pred.fixed.M, nu = nu_pred.fixed.M, tau = tau_pred.fixed.M)
      
      standardvalue <- (standardvalue.F + standardvalue.M) /2
      standardvalue.mat[,j] <- standardvalue
    }
    
    standardvalue <- rowMeans(standardvalue.mat)
    
    deviation.df <- deviation.df %>% drop_na()
    deviation.df[[paste0("SC.", i, "_standard")]] <- standardvalue
    deviation.df$visit <- SCdata.ADHD$visit
    
    deviations.sum[[n]] <- deviation.df
  }
  saveRDS(deviations.sum, paste0(interfileFolder, "/SCdata.sum75_Yeo", Yeoresolution,".ADHD_deviation_All.rds"))
}else{
  deviations.sum <- readRDS(paste0(interfileFolder, "/SCdata.sum75_Yeo", Yeoresolution,".ADHD_deviation_All.rds"))
  
  deviation.ADHD.df2 <- Reduce(function(x, y) merge(x, y, by = c("ID", "visit"), all = TRUE), deviations.sum)
}

deviation.ADHD.df2$scanID <- paste0(deviation.ADHD.df2$ID, "_", deviation.ADHD.df2$visit)
deviation.ADHD.df2 <- deviation.ADHD.df2 %>% distinct(scanID, .keep_all = T)
deviation.ADHD.df3 <- deviation.ADHD.df2 %>% left_join(SCdata.ADHD, join_by("ID", "visit"))

deviation.ADHD.df3$diagnosis <- as.factor(deviation.ADHD.df3$diagnosis)

describeBy(deviation.ADHD.df3$SC.109_deviationZ, group=deviation.ADHD.df3$diagnosis)
deviation.ADHD.df3 <- deviation.ADHD.df3 %>% drop_na(WB_SCmean)

saveRDS(deviation.ADHD.df3, paste0(interfileFolder, '/SCdata_Yeo', Yeoresolution,'_CV75.deviations_ADHD_All.rds'))
############################


## 2. For TD dataset, 10 fold CV will be conducted. 
## For each time, 9 folds will be set as training set, the rest set as test set. ref: X. Sun, 2023, BP
########################################
SCdata.TD <- SCdata.TD.trainset
SCdata.TD$site_sex <- paste0(SCdata.TD$site, "_", SCdata.TD$sex)
set.seed(12)

SCdata.TD.10fold <- SCdata.TD %>% group_by(site_sex) %>% 
  mutate(fold = sample(rep(1:10, length.out=n()))) %>% ungroup()
table(SCdata.TD.10fold$site_sex, SCdata.TD.10fold$fold)

## Fit normative models
dataname <- "SCdata.TD.trainset"
smoothterm <- "age"
covariates <- "sex+mean_fd"
randomvar <- "site"
mu.df <- degree <- 2
sigma.df <- 2
distribution.fam <- "BCPEo"
IDvar <- "ID"
quantile.vec <- c(0.025, 0.5, 0.975)
stratify <- c("sex", "site")
cl <- makeCluster(4)
clusterExport(cl, varlist = ls(), envir = environment())
invisible(clusterEvalQ(cl, {
  library(ggplot2)
  library(tidyverse)
  library(openxlsx)
  library(parallel)
  library(gamlss)
  library(scales)
  source(paste0(functionFolder, "/Construct_gamlss_set.R"))
}))

Normative_Model_fit <- function(f, cl){
  SCdata.TD.trainset <- as.data.frame(SCdata.TD.10fold %>% filter(fold!=f))
  SCdata.TD.trainset[,c("sex", "site")] <- lapply(SCdata.TD.trainset[,c("sex", "site")], as.factor)
  SCdata.test <- as.data.frame(SCdata.TD.10fold %>% filter(fold==f))
  assign("SCdata.TD.trainset", SCdata.TD.trainset, envir=.GlobalEnv)

  clusterExport(cl, varlist = ls(), envir = environment())
  mod120.training.sum <- parLapply(cl, 1:element_num, function(i){
    dependentvar <- paste0("SC.", i)
    sumlist <- construct_gamlss(dataname, dependentvar, smoothterm, covariates,randomvar, 
                                mu.df, sigma.df, degree, distribution.fam,IDvar, quantile.vec, stratify)
    
    return(sumlist)
  })
  saveRDS(mod120.training.sum, paste0(interfileFolder, "/GAMLSS_Yeo", Yeoresolution,".TDtraining_fold", f,".rds"))
  
  
  # compute deviation & draw plots
  ## test
  SCdata.test$sex <- as.factor(SCdata.test$sex)
  SCdata.test_fixed <- SCdata.test
  SCdata.test_fixed$mean_fd <- mean(SCdata.test_fixed$mean_fd)
  SCdata.test_fixed.F <- SCdata.test_fixed.M <- SCdata.test_fixed
  SCdata.test_fixed.F$sex <- factor("F")
  SCdata.test_fixed.M$sex <- factor("M")
  
  gam.data2 <- SCdata.TD.trainset
  assign("gam.data2", gam.data2, envir=.GlobalEnv)
  
  sitelist <- unique(SCdata.TD.trainset$site)
  deviations.sum <- list()
  for (n in 1:element_num){
    #i = parcel.num[n]
    i=n
    print(i)
    mod.tmp <- mod120.training.sum[[n]]$mod.tmp
    mu_pred <- predict(mod.tmp, newdata = SCdata.test, what = "mu", type = "response")
    sigma_pred <- predict(mod.tmp, newdata = SCdata.test, what = "sigma", type = "response")
    nu_pred <- predict(mod.tmp, newdata = SCdata.test, what = "nu", type = "response")
    tau_pred <- predict(mod.tmp, newdata = SCdata.test, what = "tau", type = "response")
    
    dependentvar <- paste0("SC.", i)
    deviation.df <- data.frame(ID=SCdata.test$ID)
    observation <- SCdata.test[[dependentvar]]
    centile <- pBCPEo(observation, mu = mu_pred, sigma = sigma_pred, nu = nu_pred, tau = tau_pred)
    # if (sum(centile==1) > 0){
    #   centile[which(centile==1)] <- 0.99999
    # }
    if (sum(centile==0) > 0){
      centile[which(centile==0)] <- 1-0.99999
    }else if (sum(centile==1) > 0){
      centile[which(centile==1)] <- 0.99999
    }
    deviation.df[[paste0("SC.", i, "_centile")]] <- centile
    deviation.df[[paste0("SC.", i, "_deviationZ")]] <- qnorm(centile)
    
    # compute the standard values while controlling for mean_fd and site
    standardvalue.mat <- matrix(NA, nrow(SCdata.test_fixed.F), length(sitelist))
    for (j in 1:length(sitelist)){
      site = sitelist[j]
      SCdata.test_fixed.F$site <- SCdata.test_fixed.M$site <- site
      
      # Female
      mu_pred.fixed.F <- predict(mod.tmp, newdata = SCdata.test_fixed.F, what = "mu", type = "response")
      sigma_pred.fixed.F <- predict(mod.tmp, newdata = SCdata.test_fixed.F, what = "sigma", type = "response")
      nu_pred.fixed.F <- predict(mod.tmp, newdata = SCdata.test_fixed.F, what = "nu", type = "response")
      tau_pred.fixed.F <- predict(mod.tmp, newdata = SCdata.test_fixed.F, what = "tau", type = "response")
      standardvalue.F <- qBCPEo(centile, mu = mu_pred.fixed.F, sigma = sigma_pred.fixed.F, nu = nu_pred.fixed.F, tau=tau_pred.fixed.F)
      
      # Male
      mu_pred.fixed.M <- predict(mod.tmp, newdata = SCdata.test_fixed.M, what = "mu", type = "response")
      sigma_pred.fixed.M <- predict(mod.tmp, newdata = SCdata.test_fixed.M, what = "sigma", type = "response")
      nu_pred.fixed.M <- predict(mod.tmp, newdata = SCdata.test_fixed.M, what = "nu", type = "response")
      tau_pred.fixed.M <- predict(mod.tmp, newdata = SCdata.test_fixed.M, what = "tau", type = "response")
      standardvalue.M <- qBCPEo(centile, mu = mu_pred.fixed.M, sigma = sigma_pred.fixed.M, nu = nu_pred.fixed.M, tau=tau_pred.fixed.M)
      
      standardvalue <- (standardvalue.F + standardvalue.M) /2
      standardvalue.mat[,j] <- standardvalue
    }
    
    standardvalue <- rowMeans(standardvalue.mat)
    
    deviation.df <- deviation.df %>% drop_na()
    deviation.df[[paste0("SC.", i, "_standard")]] <- standardvalue
    
    deviations.sum[[n]] <- deviation.df
  }
  
  deviation.test.unique <- lapply(deviations.sum, function(x) x %>% distinct(ID, .keep_all = T))
  
  deviation.test.df2 <- Reduce(function(x, y) merge(x, y, by = "ID", all = TRUE), deviation.test.unique)
  
  return(deviation.test.df2)
}


for (f in 1:10){
  print(f)
  deviation.test.df2 <- Normative_Model_fit(f, cl)
  saveRDS(deviation.test.df2, paste0(interfileFolder, "/SCdata.sum75_Yeo", Yeoresolution,".TD_deviation_fold", f,".rds"))

}

deviation.test.list <- list()
for (f in 1:10){
  deviation.test.list[[f]] <- readRDS(paste0(interfileFolder, "/SCdata.sum75_Yeo", Yeoresolution,".TD_deviation_fold", f,".rds"))
}

deviation.test <- do.call(rbind, deviation.test.list)

deviation.TD.df3 <- merge(deviation.test, SCdata.TD, by="ID")
deviation.TD.df3$diagnosis <- as.factor(deviation.TD.df3$diagnosis)
deviation.TD.df3$site_sex <- NULL

deviation.TD.df3$visit <- "0"
deviation.TD.df3$scanID <- paste0(deviation.TD.df3$ID, "_", deviation.TD.df3$visit)

deviation.all <- rbind(deviation.TD.df3, deviation.ADHD.df3)
describeBy(deviation.all$SC.109_deviationZ, group=deviation.all$diagnosis)

saveRDS(deviation.all, paste0(interfileFolder, '/SCdata_Yeo',Yeoresolution, '_CV75_deviations_all.rds'))


############################







