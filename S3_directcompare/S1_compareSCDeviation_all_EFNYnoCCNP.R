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
library(MatchIt)
library(tableone)

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
resultFolder_ABCD <- paste0(homepath, "/Normative_model/results_ABCD")
functionFolder.SCDev <- paste0(homepath, "/SC_development/Rcode_SCdevelopment/gamfunction")
FigureFolder <- paste0(homepath, '/Normative_model/Figures_EFNYnoCCNP/Yeo', Yeoresolution,'/CV75')

source(paste0(functionFolder.SCDev, "/gam_factor_interaction_disease.R"))
# load data
SCdata.sum75.merge.ADHDTD <- readRDS(paste0(interfileFolder, '/SCdata_Yeo',Yeoresolution, '_CV75_deviations_all.rds'))
SCdata.sum75.merge.ADHDTD$sex <- as.factor(SCdata.sum75.merge.ADHDTD$sex)
SCdata.sum75.merge.ADHDTD$diagnosis <- factor(SCdata.sum75.merge.ADHDTD$diagnosis, levels = c("TD", "ADHD"))

set.seed(925)
knots = 3
int_var <- "diagnosis"
x <- commandArgs(trailingOnly = TRUE)
dependentvar <- paste0("SC.", x, "_deviationZ")
region <- dependentvar
dataname <- "SCdata.sum75.merge.ADHDTD"
smooth_var <- "age"
covariates <- "sex+mean_fd"
print(paste("Dependent var is", region, ", independent var is", smooth_var, ", interaction var is", int_var))
gam.res <- gam.smooth.predict.interaction(region, dataname, smooth_var, int_var, covariates, knots=3, set_fx = T, stats_only=FALSE, increments=100)
gam.res <- as.data.frame(gam.res)
gam.res$SClabel <- dependentvar

if (! dir.exists(paste0(interfileFolder, "/GroupcompareDeviationAll"))){
  dir.create(paste0(interfileFolder, "/GroupcompareDeviationAll"))
}


saveRDS(gam.res, paste0(interfileFolder, '/GroupcompareDeviationAll/gamresult_', x,'.rds'))




