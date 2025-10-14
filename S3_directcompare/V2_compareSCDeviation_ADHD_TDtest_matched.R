rm(list=ls())
library(tidyverse)
library(mgcv)
library(openxlsx)
library(parallel)
library(scales)
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
interfileFolder <- paste0(homepath, '/Normative_model/interfileFolder_ABCD')
functionFolder <- paste0(homepath, "/Normative_model/functions")
resultFolder <- paste0(homepath, "/Normative_model/results_ABCD")
functionFolder.SCDev <- paste0(homepath, "/SC_development/Rcode_SCdevelopment/gamfunction")
FigureFolder <- paste0(homepath, '/Normative_model/Figures_ABCD/Yeo', Yeoresolution,'/CV75')

source(paste0(functionFolder, "/plotmatrix.R"))
source(paste0(functionFolder.SCDev, "/SCrankcorr.R"))
source(paste0(functionFolder.SCDev, "/gammsmooth.R"))
source(paste0(functionFolder.SCDev, "/plotdata_generate.R"))
source(paste0(functionFolder, "/gamm_varyingfactor.R"))
source(paste0(functionFolder, "/lmsymp.R"))
source(paste0(functionFolder.SCDev, "/gamm_factor_interaction.R"))

# load data
seed <- commandArgs(trailingOnly = TRUE)

SCdata.sum75.merge.ADHDTD <- readRDS(paste0(interfileFolder, "/SCdata_Yeo", Yeoresolution, 
                                            "_CV75_sumSCinvnode.deviations_TDtest_ADHD_match_merge.rds"))
sigedge.all <- read.csv(paste0(resultFolder, "/ABCD_sigAgeChangedEdges.csv"))
sigedge.all <- sigedge.all$x
SCdata.sum75.merge.ADHDTD$sex <- as.factor(SCdata.sum75.merge.ADHDTD$sex)
SCdata.sum75.merge.ADHDTD$if_TD <- factor(SCdata.sum75.merge.ADHDTD$if_TD, levels = c("TD", "ADHD"))

if (str_detect(wd, "cuizaixu_lab")){
  set.seed(925)
  knots = 3
  int_var <- "if_TD"
  Int.results <- mclapply(sigedge.all, function(dependentvar){
    #dependentvar <- paste0("SC.", i, "_deviationZ")
    region <- dependentvar
    dataname <- "SCdata.sum75.merge.ADHDTD"
    smooth_var <- "age"
    covariates <- "sex+meanFD"
    print(paste("Dependent var is", region, ", independent var is", smooth_var, ", interaction var is", int_var))
    gam.res <- gamm.smooth.predict.interaction(region, dataname, smooth_var, int_var, covariates=covariates, knots, set_fx = T)
    gam.res <- as.data.frame(gam.res)
    gam.res$SClabel <- dependentvar
    
    return(gam.res)
  }, mc.cores=20)
  
  Int.results.df <- do.call(rbind, Int.results)
  saveRDS(Int.results.df, paste0(interfileFolder, "/Deviation_TDtest-ADHDMatched_SC", element_num,"_merge.rds"))
}else{
  Int.results.df <- readRDS(paste0(interfileFolder, "/Deviation_TDtest-ADHDMatched_SC", element_num,"_merge.rds"))
}



