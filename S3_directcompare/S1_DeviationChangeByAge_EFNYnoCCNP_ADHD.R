rm(list=ls())
library(mgcv)
library(parallel)
library(tidyverse)
library(patchwork)
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

source(paste0(functionFolder.SCDev, "/gamsmooth.R"))
source(paste0(functionFolder.SCDev, "/plotdata_generate.R"))
source(paste0(functionFolder.SCDev, "/SCrankcorr.R"))
source(paste0(functionFolder.SCDev, '/gamderivatives.R'))
# Load data
SCdata <- readRDS(paste0(interfileFolder, '/SCdata_Yeo',Yeoresolution, '_CV75_deviations_all.rds'))

Behavior <- read.csv(paste0(demopath, "/basic_demo_PKU6_addSymptom.csv"))
Behavior$ID <- paste0("sub-", Behavior$ID)
SCdata <- SCdata %>% left_join(dplyr::select(Behavior, ID, medication, IA_B, HI_B, TO_B, IA_F, HI_F, TO_F, IA_B.norm, HI_B.norm, TO_B.norm, RVPA, RVPML, RVPTFA, RVPTM, RVPTH, RVPPH, SSTSSRT, SSTDEG, SSTDES, SWMBE, ICV), by="ID")
SCdata$sex <- as.factor(SCdata$sex)
SCdata$diagnosis <- factor(SCdata$diagnosis, levels = c("TD", "ADHD"))
table(SCdata$diagnosis, SCdata$site)

SCdata.ADHD <- SCdata %>% filter(diagnosis == "ADHD", visit=="0")
SCdata.TD <- SCdata %>% filter(diagnosis == "TD")

covariates<-"sex+mean_fd"
dataname <- "SCdata.ADHD"
print(dataname)

x <- commandArgs(trailingOnly = TRUE)

region <- paste0("SC.", x, "_deviationZ")

gamresult<-gam.fit.smooth(region, dataname, smooth_var="age", covariates, knots=3, set_fx = T, stats_only = T, mod_only = FALSE)
gamresult<-as.data.frame(gamresult)

if (! dir.exists(paste0(interfileFolder, "/deviationdevelopmentADHD"))){
  dir.create(paste0(interfileFolder, "/deviationdevelopmentADHD"))
}


saveRDS(gamresult, paste0(interfileFolder, '/deviationdevelopmentADHD/gamresult_', x,'.rds'))



