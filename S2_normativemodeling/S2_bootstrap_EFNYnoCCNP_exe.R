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
interfileFolder <- paste0(homepath, '/Normative_model/interfileFolder_EFNYnoCCNP')
functionFolder <- paste0(homepath, "/Normative_model/functions")
resultFolder <- paste0(homepath, "/Normative_model/results_EFNYnoCCNP")
functionFolder.SCDev <- paste0(homepath, "/SC_development/Rcode_SCdevelopment/gamfunction")
FigureFolder <- paste0(homepath, '/Normative_model/Figures_EFNYnoCCNP/Yeo', Yeoresolution,'/CV75')

# load data
SCdataTD <- readRDS(paste0(interfileFolder, "/SCdata.allTD_EFNYnoCCNP.SCYeo", element_num, ".rds"))
source(paste0(wd, "/final_script/S2_normativemodeling/S2_bootstrap_EFNYnoCCNP.R"))

# get input
n <- commandArgs(trailingOnly = TRUE)
n <- as.numeric(n)

for (i in 1:element_num){
  SClabel <- paste0("SC.", i)
  SCdata.sum75.merge.TD <- SCdataTD %>% select(all_of(c(SClabel, "age", "sex", "mean_fd", "site", "ID")))
  SCdata.sum75.merge.TD$sex <- as.factor(SCdata.sum75.merge.TD$sex)
  SCdata.sum75.merge.TD$site <- as.factor(SCdata.sum75.merge.TD$site)
  
  execute_boot(n, SClabel)
}



