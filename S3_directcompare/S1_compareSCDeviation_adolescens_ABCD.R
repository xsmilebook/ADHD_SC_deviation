
rm(list=ls())
library(ggplot2)
library(corrplot)
library(RColorBrewer)
library(tidyverse)
library(mgcv)
library(lme4)
library(openxlsx)
library(parallel)
library(scales)
library(psych)
library(reshape)
library(patchwork)
library(ecostats)

# set resolution
Yeoresolution <- 17
if (Yeoresolution == 7){
  Yeoresolution.delLM = 6
}else if (Yeoresolution == 17){
  Yeoresolution.delLM = 15
}
element_num <- Yeoresolution.delLM*(Yeoresolution.delLM+1)/2

# input directory
homepath <- "D:/code/ADHD_SC_deviation"
demopath <- file.path(homepath, "data", 'demography')
interfileFolder <- file.path(homepath, "data", 'interfileFolder', "ABCD")
functionFolder <- file.path(homepath, "src", "functions")
resultFolder <- file.path(homepath, "reports", "results", "ABCD")
functionFolder.SCDev <- file.path(homepath, "src", "gamfunction")
FigureFolder <- paste0(homepath, '/reports/figures/ABCD/Yeo', Yeoresolution,'/CV75')


source(paste0(functionFolder, "/plotmatrix.R"))
source(paste0(functionFolder.SCDev, "/SCrankcorr.R"))
source(paste0(functionFolder.SCDev, "/gammsmooth.R"))
source(paste0(functionFolder.SCDev, "/plotdata_generate.R"))
source(paste0(functionFolder, "/gamm_varyingfactor.R"))
source(paste0(functionFolder, "/lmsymp.R"))
source(paste0(functionFolder.SCDev, "/gamm_factor_interaction.R"))
# load data
SCdata.sum75.merge.ADHDTD <- readRDS(paste0(interfileFolder, '/SCdata_Yeo',Yeoresolution, '_CV75_sumSCinvnode.deviations_TDtest_ADHD.rds'))
SCdata.sum75.merge.ADHDTD$sex <- as.factor(SCdata.sum75.merge.ADHDTD$sex)
SCdata.sum75.merge.ADHDTD$meanFD <- as.numeric(SCdata.sum75.merge.ADHDTD$meanFD)
SCdata.sum75.merge.ADHDTD$if_TD <- factor(SCdata.sum75.merge.ADHDTD$if_TD, levels = c("TD", "ADHD"))


SCdata.child <- SCdata.sum75.merge.ADHDTD %>% filter(tannerstage ==1)
SCdata.adoles <- SCdata.sum75.merge.ADHDTD %>% filter(tannerstage >1)
rbind(table(SCdata.child$if_TD, SCdata.child$sex), table(SCdata.adoles$if_TD, SCdata.adoles$sex))


dataname <- "SCdata.adoles"
data.tmp <- get(dataname)

# input
x <- commandArgs(trailingOnly = TRUE)

region <- paste0("SC.", x, "_deviationZ")
print(region)
if (dataname == "SCdata.adoles"){
  modelformula <- as.formula(sprintf("%1$s ~ if_TD+s(age, k=3, fx=T)+ sex+meanFD", region))
  gam.model <- gam(modelformula, REML=TRUE, data = data.tmp)
  modelformula.null <- as.formula(sprintf("%1$s ~ s(age, k=3, fx=T)+ sex+meanFD", region))
  gam.model.null <- gam(modelformula.null, REML=TRUE, data = data.tmp)
  
  anova.result <- anovaPB(gam.model.null, gam.model, n.sim = 10000, ncpus=1, test="Chisq")
  anova.cov.pvalue<-anova.result[2, 5]
  
  lm.model <- gam.model
  mod.result <- summary(gam.model)
  group_diff <- mod.result$p.table[2,1]
  tvalue <- mod.result$p.table[2,3]
}else{
  covariates <- "sex+age+meanFD"
  modelformula <- as.formula(sprintf("%s ~ %s + %s",region, "if_TD", covariates))
  modelformula.null<-as.formula(sprintf("%s ~ %s",region, covariates))
  lm.model <- lm(modelformula, data = data.tmp)
  lm.model.null <- lm(modelformula.null, data = data.tmp)
  anova.result <- anovaPB(lm.model.null, lm.model, n.sim = 10000, ncpus=1, test="Chisq")
  anova.cov.pvalue<-anova.result$`Pr(>F)`[2]
  mod.result <- summary(lm.model)
  tvalue=mod.result$coefficients[2,3]
  group_diff <- mod.result$coefficients[2,1]
}


resid_sd <- sigma(lm.model)
cohen_d <- group_diff / resid_sd

result.df <- data.frame(parcel=region, group="if_TD", dataname=dataname, tvalue=tvalue,
                        pvalue=anova.cov.pvalue, cohen_d=cohen_d)


if (! dir.exists(paste0(interfileFolder, "/GroupcompareDeviationAdolescens"))){
  dir.create(paste0(interfileFolder, "/GroupcompareDeviationAdolescens"))
}


saveRDS(result.df, paste0(interfileFolder, "/GroupcompareDeviationAdolescens/Yeo",
                                     element_num, "_", str_split_i(dataname, "SCdata.", 2), "_SC", x,".rds"))


