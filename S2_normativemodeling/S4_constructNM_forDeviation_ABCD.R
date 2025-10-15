rm(list=ls())
library(ggplot2)
library(tidyverse)
library(openxlsx)
library(parallel)
library(gamlss)
library(scales)
library(tableone)
#library(MatchIt)
# set resolution
Yeoresolution <- 17
if (Yeoresolution == 7){
  Yeoresolution.delLM = 6
}else if (Yeoresolution == 17){
  Yeoresolution.delLM = 15
}
element_num <- Yeoresolution.delLM*(Yeoresolution.delLM+1)/2
# input directory
homepath <- "/ibmgpfs/cuizaixu_lab/xuhaoshu/ADHD_SC_deviation"
demopath <- file.path(homepath, "data", 'demography')
interfileFolder <- file.path(homepath, "data", 'interfileFolder', "ABCD")
functionFolder <- file.path(homepath, "src", "functions")
resultFolder <- file.path(homepath, "reports", "results", "ABCD")
functionFolder.SCDev <- file.path(homepath, "src", "gamfunction")
FigureFolder <- paste0(homepath, '/reports/figures/ABCD/Yeo', Yeoresolution,'/CV75')

# Load data
SCdata <- readRDS(paste0(interfileFolder, '/SCdata_Yeo',Yeoresolution, '_CV75_sumSCinvnode.sum.msmtcsd.merge.rds'))
SCdata$sex <- as.factor(SCdata$demo_sex_v2)

# data clean
SCdata <- SCdata %>% mutate(WB_SCmean = rowMeans(dplyr::select(.,starts_with("SC."))))
#SCdata <- SCdata %>% filter(WB_SCmean>=mean(WB_SCmean)-3*sd(WB_SCmean), WB_SCmean<=mean(WB_SCmean)+3*sd(WB_SCmean))

SCdata.TD <- SCdata %>% filter(if_TD==1)
SCdata.ADHD <- SCdata %>% filter(if_TD==0)

# source function
source(paste0(functionFolder, "/Construct_gamlss_set.R"))

SCdata.TD.trainset <- readRDS(paste0(interfileFolder, "/SCdata.TD.trainset_SCYeo", element_num, ".rds"))
SCdata.TD.testset <- readRDS(paste0(interfileFolder, "/SCdata.TD.testset_SCYeo", element_num, ".rds"))

# Describe demographic information.
SCdata <- SCdata %>% mutate(
  group = case_when(
    scanID %in% SCdata.TD.trainset$scanID ~ "TDtrain",
    scanID %in% SCdata.TD.testset$scanID ~ "TDtest",
    if_ADHD == 1 ~ "ADHD",
    .default = NA
  ))

table(SCdata$group)
demovar <- c("sex", "age", "meanFD", "siteID", "ehi_y_ss_scoreb", "GENERAL", "cbcl_scr_syn_attention_t", "cbcl_scr_dsm5_adhd_t",
             "income.adj", "nihtbx_fluidcomp_uncorrected", "nihtbx_totalcomp_uncorrected",
             "nihtbx_cryst_uncorrected")
tableone.df <- CreateTableOne(demovar, strata="group", data=SCdata, 
                              factorVars = c("sex", "siteID", "ehi_y_ss_scoreb"), includeNA = T,
                              test=T)
tableone.df <- print(tableone.df)
write.csv(tableone.df, paste0(resultFolder, "/demoinfo_testtrainsets_Yeo17.csv"), row.names = T)

tableone.df <- CreateTableOne(demovar, strata="group", data=SCdata[SCdata$if_TD==1, ], 
                              factorVars = c("sex", "siteID", "ehi_y_ss_scoreb"), includeNA = T,
                              test=T)
tableone.df <- print(tableone.df)
write.csv(tableone.df, paste0(resultFolder, "/demoinfo_testtrainsetsTD_Yeo17.csv"), row.names = T)

## Fit normative models
SCdata.TD.trainset <- as.data.frame(SCdata.TD.trainset)
SCdata.TD.trainset[,c("sex", "siteID")] <- lapply(SCdata.TD.trainset[,c("sex", "siteID")], as.factor)
dataname <- "SCdata.TD.trainset"
smoothterm <- "age"
covariates <- "sex+meanFD"
randomvar <- "siteID"
mu.df <- sigma.df <- degree <- 2
distribution.fam <- "GG"
IDvar <- "scanID"
quantile.vec <- c(0.025, 0.5, 0.975)
stratify <- c("sex", "siteID")
if (! file.exists(paste0(interfileFolder, "/GAMLSS_Yeo", Yeoresolution,".TDtraining.sum.rds"))){
  mod.training.sum <- mclapply(1:element_num, function(i){
    dependentvar <- paste0("SC.", i)
    sumlist <- construct_gamlss(dataname, dependentvar, smoothterm, covariates,randomvar, 
                                mu.df, sigma.df, degree, distribution.fam,IDvar, quantile.vec, stratify)
    
    return(sumlist)
  }, mc.cores = 64)
  saveRDS(mod.training.sum, paste0(interfileFolder, "/GAMLSS_Yeo", Yeoresolution,".TDtraining.sum.rds"))
}else{
  mod.training.sum <- readRDS(paste0(interfileFolder, "/GAMLSS_Yeo", Yeoresolution,".TDtraining.sum.rds"))
}

mod.training.WB <- construct_gamlss(dataname, dependentvar="WB_SCmean", smoothterm, covariates,randomvar, 
                                    mu.df, sigma.df, degree, distribution.fam,IDvar, quantile.vec, stratify)

mod.training.sum[[element_num+1]] <- mod.training.WB
saveRDS(mod.training.WB, paste0(interfileFolder, "/GAMLSS_Yeo", Yeoresolution,"_WBSCmean.TDtraining.sum.rds"))

# compute deviations
SCdata.TD.trainset.subset <- SCdata.TD.trainset %>% dplyr::select(c(paste0("SC.", 1:element_num),"WB_SCmean", "age", "sex", "siteID","scanID",
                                                                    "meanFD", "if_TD", "eventname2")) %>% drop_na()
SCdata.ADHD <- SCdata.ADHD %>% filter(siteID %in% SCdata.TD.trainset$siteID)
SCdata.ADHD <- SCdata.ADHD %>% filter( WB_SCmean> mean(WB_SCmean)-3*sd(WB_SCmean), WB_SCmean< mean(WB_SCmean)+3*sd(WB_SCmean))
SCdata.ADHD.subset <- SCdata.ADHD %>% dplyr::select(c(paste0("SC.", 1:element_num),"WB_SCmean", "age", "sex", "siteID","scanID",
                                                      "meanFD","if_TD", "eventname2")) %>% drop_na()

SCdata.TD.testset.subset <- SCdata.TD.testset %>% dplyr::select(c(paste0("SC.", 1:element_num),"WB_SCmean", "age", "sex", "siteID","scanID",
                                                                  "meanFD", "if_TD", "eventname2")) %>% drop_na()
SCdata.requireZ <- rbind(SCdata.ADHD.subset, SCdata.TD.testset.subset)

# compute deviation & draw plots
SCdata.requireZ$sex <- as.factor(SCdata.requireZ$sex)
SCdata.requireZ_fixed <- SCdata.requireZ
SCdata.requireZ_fixed$meanFD <- mean(SCdata.requireZ_fixed$meanFD)
SCdata.requireZ_fixed.F <- SCdata.requireZ_fixed.M <- SCdata.requireZ_fixed
SCdata.requireZ_fixed.F$sex <- factor(2)
SCdata.requireZ_fixed.M$sex <- factor(1)
SCdata.TD.trainset$siteID <- droplevels(SCdata.TD.trainset$siteID)
sitelist <- unique(SCdata.TD.trainset$siteID)

gam.data2 <- SCdata.TD.trainset.subset
if (! file.exists(paste0(interfileFolder, "/SCdata.sum75_Yeo", Yeoresolution,".testset_ADHD_deviation.rds"))){
  
  deviations.sum <- mclapply(1:(element_num+1), function(i){
    mod.tmp <- mod.training.sum[[i]]$mod.tmp
    mu_pred <- predict(mod.tmp, newdata = SCdata.requireZ, what = "mu", type = "response")
    sigma_pred <- predict(mod.tmp, newdata = SCdata.requireZ, what = "sigma", type = "response")
    nu_pred <- predict(mod.tmp, newdata = SCdata.requireZ, what = "nu", type = "response")
    
    if (i == (element_num+1)){
      dependentvar <- "WB_SCmean"
      deviation.df <- data.frame(scanID=SCdata.requireZ$scanID)
      observation <- SCdata.requireZ[[dependentvar]]
      centile <- pGG(observation, mu = mu_pred, sigma = sigma_pred, nu = nu_pred)
      deviation.df[["WB_SCmean_centile"]] <- centile
      deviation.df[["WB_SCmean_deviationZ"]] <- qnorm(centile)
    }else{
      dependentvar <- paste0("SC.", i)
      deviation.df <- data.frame(scanID=SCdata.requireZ$scanID)
      observation <- SCdata.requireZ[[dependentvar]]
      centile <- pGG(observation, mu = mu_pred, sigma = sigma_pred, nu = nu_pred)
      deviation.df[[paste0("SC.", i, "_centile")]] <- centile
      deviation.df[[paste0("SC.", i, "_deviationZ")]] <- qnorm(centile)
    }
    
    
    # compute the standard values while controlling for meanFD and siteID
    standardvalue.mat <- matrix(NA, nrow(SCdata.requireZ_fixed.F), length(sitelist))
    for (j in 1:length(sitelist)){
      siteID = sitelist[j]
      SCdata.requireZ_fixed.F$siteID <- SCdata.requireZ_fixed.M$siteID <- siteID
      
      # Female
      mu_pred.fixed.F <- predict(mod.tmp, newdata = SCdata.requireZ_fixed.F, what = "mu", type = "response")
      sigma_pred.fixed.F <- predict(mod.tmp, newdata = SCdata.requireZ_fixed.F, what = "sigma", type = "response")
      nu_pred.fixed.F <- predict(mod.tmp, newdata = SCdata.requireZ_fixed.F, what = "nu", type = "response")
      standardvalue.F <- qGG(centile, mu = mu_pred.fixed.F, sigma = sigma_pred.fixed.F, nu = nu_pred.fixed.F)
      
      # Male
      mu_pred.fixed.M <- predict(mod.tmp, newdata = SCdata.requireZ_fixed.M, what = "mu", type = "response")
      sigma_pred.fixed.M <- predict(mod.tmp, newdata = SCdata.requireZ_fixed.M, what = "sigma", type = "response")
      nu_pred.fixed.M <- predict(mod.tmp, newdata = SCdata.requireZ_fixed.M, what = "nu", type = "response")
      standardvalue.M <- qGG(centile, mu = mu_pred.fixed.M, sigma = sigma_pred.fixed.M, nu = nu_pred.fixed.M)
      
      standardvalue <- (standardvalue.F + standardvalue.M) /2
      standardvalue.mat[,j] <- standardvalue
    }
    
    standardvalue <- rowMeans(standardvalue.mat)
    
    deviation.df <- deviation.df %>% drop_na()
    if (i == (element_num+1)){
      deviation.df[["WB_SCmean_standard"]] <- standardvalue
    }else{
      deviation.df[[paste0("SC.", i, "_standard")]] <- standardvalue
    }
    
    
    return(deviation.df)
  }, mc.cores = 64)
  saveRDS(deviations.sum, paste0(interfileFolder, "/SCdata.sum75_Yeo", Yeoresolution,".testset_ADHD_deviation.rds"))
}else{
  deviations.sum <- readRDS(paste0(interfileFolder, "/SCdata.sum75_Yeo", Yeoresolution,".testset_ADHD_deviation.rds"))
}

deviation.requireZ.unique <- lapply(deviations.sum, function(x) x %>% distinct(scanID, .keep_all = T))

deviation.requireZ.df2 <- Reduce(function(x, y) merge(x, y, by = "scanID", all = TRUE), deviation.requireZ.unique)

deviation.requireZ.df3 <- merge(deviation.requireZ.df2, SCdata, by="scanID") %>% 
  filter(scanID %in% c(SCdata.ADHD$scanID, SCdata.TD.testset$scanID))
deviation.requireZ.df3$if_TD <- factor(deviation.requireZ.df3$if_TD, levels = c(0, 1), labels = c("ADHD", "TD"))
describeBy(deviation.requireZ.df3$WB_SCmean_standard, group=deviation.requireZ.df3$if_TD)

saveRDS(deviation.requireZ.df3, paste0(interfileFolder, '/SCdata_Yeo',Yeoresolution, '_CV75_sumSCinvnode.deviations_TDtest_ADHD.rds'))


## ADHD
Behavior <- read.csv(paste0(demopath, "/demo_sublist7.csv"))
deviation.requireZ.ADHD <- deviation.requireZ.df3 %>% filter(if_TD=="ADHD")
deviation.requireZ.ADHD <- deviation.requireZ.ADHD %>% left_join(select(Behavior, scanID, GENERAL.norm, ATT.norm, HI.norm), by="scanID")
deviation.requireZ.ADHD$sex <- factor(deviation.requireZ.ADHD$sex, levels=c(1, 2), labels = c("M", "F"))
Interest.vars <- c(paste0("SC.", 1:element_num, "_deviationZ"), "cbcl_scr_syn_external_t", 
                   "cbcl_scr_syn_internal_t", "cbcl_scr_dsm5_adhd_t", "age", "meanFD")
tableone.df <- CreateTableOne(vars=Interest.vars, strata = "sex", data = deviation.requireZ.ADHD, test = TRUE)
tableone.df <- print(tableone.df)
write.csv(tableone.df, paste0(resultFolder, "/demoinfo_sexADHD.csv"), row.names = T)
## TD
Behavior <- read.csv(paste0(demopath, "/demo_sublist7.csv"))
deviation.requireZ.TD <- deviation.requireZ.df3 %>% filter(if_TD=="TD")
deviation.requireZ.TD <- deviation.requireZ.TD %>% left_join(select(Behavior, scanID, GENERAL.norm, ATT.norm, HI.norm), by="scanID")
deviation.requireZ.TD$sex <- factor(deviation.requireZ.TD$sex, levels=c(1, 2), labels = c("M", "F"))
Interest.vars <- c(paste0("SC.", 1:element_num, "_deviationZ"), "cbcl_scr_syn_external_t", 
                   "cbcl_scr_syn_internal_t", "cbcl_scr_dsm5_adhd_t", "age", "meanFD")
tableone.df <- CreateTableOne(vars=Interest.vars, strata = "sex", data = deviation.requireZ.TD, test = TRUE)
tableone.df <- print(tableone.df)
write.csv(tableone.df, paste0(resultFolder, "/demoinfo_sexTD.csv"), row.names = T)

