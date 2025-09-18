
rm(list=ls())
library(ggplot2)
library(corrplot)
library(RColorBrewer)
library(tidyverse)
library(mgcv)
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
wd <- getwd()
homepath <- str_split_i(wd, "Normative_model", 1)
demopath <- paste0(homepath, '/Normative_model/demography')
interfileFolder <- paste0(homepath, '/Normative_model/interfileFolder_EFNYnoCCNP')
functionFolder <- paste0(homepath, "/Normative_model/functions")
resultFolder <- paste0(homepath, "/Normative_model/results_EFNYnoCCNP")
functionFolder.SCDev <- paste0(homepath, "/SC_development/Rcode_SCdevelopment/gamfunction")
FigureFolder <- paste0(homepath, '/Normative_model/Figures_EFNYnoCCNP/Yeo', Yeoresolution,'/CV75')

source(paste0(functionFolder, "/plotmatrix.R"))
source(paste0(functionFolder.SCDev, "/SCrankcorr.R"))

# load data
SCdata.sum75.merge.ADHDTD <- readRDS(paste0(interfileFolder, '/SCdata_Yeo',Yeoresolution, '_CV75_deviations_all.rds'))
SCdata.sum75.merge.ADHDTD$sex <- as.factor(SCdata.sum75.merge.ADHDTD$sex)
SCdata.sum75.merge.ADHDTD$FD <- as.numeric(SCdata.sum75.merge.ADHDTD$FD)
SCdata.sum75.merge.ADHDTD$diagnosis <- factor(SCdata.sum75.merge.ADHDTD$diagnosis, levels = c("TD", "ADHD"))


SCdata.child <- SCdata.sum75.merge.ADHDTD %>% filter(age <=10)
SCdata.adoles <- SCdata.sum75.merge.ADHDTD %>% filter(age >10)
rbind(table(SCdata.child$diagnosis, SCdata.child$sex), table(SCdata.adoles$diagnosis, SCdata.adoles$sex))


dataname <- "SCdata.adoles"
data.tmp <- get(dataname)

# input
x <- commandArgs(trailingOnly = TRUE)

region <- paste0("SC.", x, "_deviationZ")
print(region)
if (dataname == "SCdata.adoles"){
  modelformula <- as.formula(sprintf("%1$s ~ diagnosis+s(age, k=3, fx=T)+ sex+FD", region))
  gam.model <- gam(modelformula, REML=TRUE, data = data.tmp)
  modelformula.null <- as.formula(sprintf("%1$s ~ s(age, k=3, fx=T)+ sex+FD", region))
  gam.model.null <- gam(modelformula.null, REML=TRUE, data = data.tmp)
  
  anova.result <- anovaPB(gam.model.null, gam.model, n.sim = 10000, ncpus=1, test="Chisq")
  anova.cov.pvalue<-anova.result[2, 5]
  
  lm.model <- gam.model
  mod.result <- summary(gam.model)
  group_diff <- mod.result$p.table[2,1]
  tvalue <- mod.result$p.table[2,3]
}else{
  covariates <- "sex+age+FD"
  modelformula <- as.formula(sprintf("%s ~ %s + %s",region, "diagnosis", covariates))
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

result.df <- data.frame(parcel=region, group="diagnosis", dataname=dataname, tvalue=tvalue,
                        pvalue=anova.cov.pvalue, cohen_d=cohen_d)


if (! dir.exists(paste0(interfileFolder, "/GroupcompareDeviationAdolescens"))){
  dir.create(paste0(interfileFolder, "/GroupcompareDeviationAdolescens"))
}


saveRDS(result.df, paste0(interfileFolder, "/GroupcompareDeviationAdolescens/Yeo",
                                     element_num, "_", str_split_i(dataname, "SCdata.", 2), "_SC", x,".rds"))


