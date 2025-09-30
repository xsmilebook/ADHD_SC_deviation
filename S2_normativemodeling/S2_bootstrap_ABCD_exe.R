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
homepath <- "/ibmgpfs/cuizaixu_lab/xuhaoshu/ADHD_SC_deviation"
demopath <- file.path(homepath, "data", "demography")
interfileFolder <- file.path(homepath, "data", "interfileFolder", "ABCD")
functionFolder <- file.path(homepath, "src", "functions")
resultFolder <- file.path(homepath, "reports", "results", "ABCD")
functionFolder.SCDev <- file.path(homepath, "src", "gamfunction")
FigureFolder <- paste0(homepath, '/reports/figures/ABCD/Yeo', Yeoresolution,'/CV75')

# load data
SCdataTD <- readRDS(paste0(interfileFolder, "/SCdata.TD.trainset_SCYeo", element_num, ".rds"))
source(file.path(homepath, "src", "S2_normativemodeling", "S2_bootstrap_ABCD.R"))

n <- commandArgs(trailingOnly = TRUE)
n <- as.numeric(n)

SClabels <- paste0("SC.", 1:element_num)

for (SClabel in SClabels) {
  cat("Processing", SClabel, "\n")
  SCdata.sum75.merge.TD <- SCdataTD %>% select(all_of(c(SClabel, "age", "sex", "meanFD", "siteID", "scanID")))
  SCdata.sum75.merge.TD$sex <- factor(SCdata.sum75.merge.TD$sex, levels=c(1,2), labels=c("M", "F"))
  SCdata.sum75.merge.TD$siteID <- as.factor(SCdata.sum75.merge.TD$siteID)
  execute_boot(n, SClabel)
}

# WB_SCmean
SCdata.sum75.merge.TD <- SCdataTD %>% select(all_of(c("WB_SCmean", "age", "sex", "mean_fd", "siteID", "ID")))
SCdata.sum75.merge.TD$sex <- as.factor(SCdata.sum75.merge.TD$sex)
SCdata.sum75.merge.TD$siteID <- as.factor(SCdata.sum75.merge.TD$siteID)
execute_boot(n, "WB_SCmean")