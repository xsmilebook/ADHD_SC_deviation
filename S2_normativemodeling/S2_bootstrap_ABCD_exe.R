rm(list=ls())
library(tidyverse)

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

# load data
SCdataTD <- readRDS(paste0(interfileFolder, "/SCdata.TD.trainset_SCYeo", element_num, ".rds"))
source(paste0(wd, "/final_script/S2_normativemodeling/S2_bootstrap_ABCD.R"))

# get input
n <- commandArgs(trailingOnly = TRUE)
n <- as.numeric(n)

for (i in 1:element_num){
  SClabel <- paste0("SC.", i)
  SCdata.sum75.merge.TD <- SCdataTD %>% select(all_of(c(SClabel, "age", "sex", "meanFD", "siteID", "scanID")))
  SCdata.sum75.merge.TD$sex <- factor(SCdata.sum75.merge.TD$sex, levels=c(1,2), labels=c("M", "F"))
  SCdata.sum75.merge.TD$siteID <- as.factor(SCdata.sum75.merge.TD$siteID)
  
  execute_boot(n, SClabel)
}



