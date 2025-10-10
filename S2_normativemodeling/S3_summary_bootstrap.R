library(ggplot2)
library(RColorBrewer)
library(tidyverse)
library(mgcv)
library(openxlsx)
library(parallel)
library(gamlss)
library(scales)
library(psych)
library(reshape)

Yeoresolution <- 17
if (Yeoresolution == 7){
  Yeoresolution.delLM <- 6
}else if (Yeoresolution == 17){
  Yeoresolution.delLM <- 15
}
element_num <- Yeoresolution.delLM*(Yeoresolution.delLM+1)/2

homepath        <- "/ibmgpfs/cuizaixu_lab/xuhaoshu/ADHD_SC_deviation"
demopath        <- file.path(homepath, "data", "demography")
interfileFolder <- file.path(homepath, "data", "interfileFolder", "ABCD")
functionFolder  <- file.path(homepath, "src", "functions")
resultFolder    <- file.path(homepath, "reports", "results", "ABCD")
functionFolder.SCDev <- file.path(homepath, "src", "gamfunction")
FigureFolder    <- paste0(homepath, "/reports/figures/ABCD/Yeo", Yeoresolution, "/CV75")

SCdata.sum75.merge        <- readRDS(paste0(interfileFolder, "/SCdata_Yeo", Yeoresolution, "_CV75_sumSCinvnode.sum.msmtcsd.merge.rds"))
SCdata.sum75.merge.TD.trainset <- readRDS(paste0(interfileFolder, "/SCdata.TD.trainset_SCYeo", element_num, ".rds"))
SCdata.sum75.merge.TD.testset  <- readRDS(paste0(interfileFolder, "/SCdata.TD.testset_SCYeo", element_num, ".rds"))

source(paste0(functionFolder, "/Construct_gamlss_set.R"))

bootstrapdir <- paste0(interfileFolder, "/bootstrap/WB_SCmean")
boot_files   <- list.files(bootstrapdir,
                           pattern = "^centile_bootstrap_(\\d+)\\.rds$",
                           full.names = TRUE)
boot_ids     <- as.integer(gsub("^centile_bootstrap_(\\d+)\\.rds$", "\\1", basename(boot_files)))

centile_female_boot <- vector("list", 1000)
centile_male_boot   <- vector("list", 1000)
centile_boot        <- vector("list", 1000)

for (ii in seq_along(boot_files)){
  idx      <- boot_ids[ii]
  centile.tmp <- readRDS(boot_files[ii])
  
  if (centile.tmp$converged == TRUE){
    cent_F <- apply(centile.tmp$Centiles_female, c(2,3), mean)
    cent_M <- apply(centile.tmp$Centiles_male,   c(2,3), mean)
    centile_female_boot[[idx]] <- cent_F
    centile_male_boot[[idx]]   <- cent_M
    centile_boot[[idx]]        <- (cent_F + cent_M)/2
  }else{
    centile_female_boot[[idx]] <- NA
    centile_male_boot[[idx]]   <- NA
    centile_boot[[idx]]        <- NA
  }
}

n_missing <- sum(sapply(centile_female_boot, function(x) all(is.na(x))))
print(paste0(n_missing, " iterations are missing / not converged."))

centile_female_boot <- centile_female_boot[!sapply(centile_female_boot, function(x) all(is.na(x)))]
centile_male_boot   <- centile_male_boot[!sapply(centile_male_boot,   function(x) all(is.na(x)))]
centile_boot        <- centile_boot[!sapply(centile_boot,        function(x) all(is.na(x)))]

centile_female_boot_3d <- array(unlist(centile_female_boot),
                                dim = c(7, ncol(centile_female_boot[[1]]), length(centile_female_boot)))
centile_male_boot_3d   <- array(unlist(centile_male_boot),
                                dim = c(7, ncol(centile_male_boot[[1]]),   length(centile_male_boot)))
centile_boot_3d        <- array(unlist(centile_boot),
                                dim = c(7, ncol(centile_boot[[1]]),        length(centile_boot)))

quantiles <- c(0.01, 0.05, 0.25, 0.5, 0.75, 0.95, 0.99)
n_age     <- ncol(centile_female_boot_3d[1,,1])

quantile.median <- as.data.frame(matrix(NA, n_age, 7*3))
quantile.P2.5   <- as.data.frame(matrix(NA, n_age, 7*3))
quantile.P97.5  <- as.data.frame(matrix(NA, n_age, 7*3))

colnames(quantile.median) <- colnames(quantile.P2.5) <- colnames(quantile.P97.5) <-
  c(paste0("female_quantiles",   quantiles),
    paste0("male_quantiles",     quantiles),
    paste0("all_quantiles",      quantiles))

for (i in 1:7){
  # female
  qf.vec <- centile_female_boot_3d[i,,]
  quantile.median[,i]   <- apply(qf.vec, 1, median)
  quantile.P2.5[,i]     <- apply(qf.vec, 1, quantile, probs = 0.025)
  quantile.P97.5[,i]    <- apply(qf.vec, 1, quantile, probs = 0.975)
  
  # male
  qm.vec <- centile_male_boot_3d[i,,]
  quantile.median[,i+7]   <- apply(qm.vec, 1, median)
  quantile.P2.5[,i+7]     <- apply(qm.vec, 1, quantile, probs = 0.025)
  quantile.P97.5[,i+7]    <- apply(qm.vec, 1, quantile, probs = 0.975)
  
  # all
  qa.vec <- centile_boot_3d[i,,]
  quantile.median[,i+14]   <- apply(qa.vec, 1, median)
  quantile.P2.5[,i+14]     <- apply(qa.vec, 1, quantile, probs = 0.025)
  quantile.P97.5[,i+14]    <- apply(qa.vec, 1, quantile, probs = 0.975)
}

age_seq <- seq(min(SCdata.sum75.merge.TD.trainset$age),
               max(SCdata.sum75.merge.TD.trainset$age),
               length.out = n_age)
quantile.median$age <- age_seq
quantile.P2.5$age   <- age_seq
quantile.P97.5$age  <- age_seq

p_all <- ggplot() +
  geom_line(data = quantile.median, aes(x = age, y = all_quantiles0.5), linetype = "solid", linewidth = 1.5) +
  geom_line(data = quantile.P2.5,   aes(x = age, y = all_quantiles0.5), linetype = 2, linewidth = 1) +
  geom_line(data = quantile.P97.5,  aes(x = age, y = all_quantiles0.5), linetype = 2, linewidth = 1) +
  labs(x = "Age (years)", y = "SC strength",
       title = "Bootstrap 95% CI of Whole-brain averaged SC strength") +
  theme_classic() +
  theme(plot.title = element_text(size = 20, hjust = 0.5),
        axis.text  = element_text(size = 20),
        axis.title = element_text(size = 20))

ggsave(paste0(FigureFolder, "/ModelEvaluation/bootstrap_P50_95CI_WBmeanSC_all.tiff"),
       plot = p_all, width = 16, height = 14, units = "cm")

p_sex <- ggplot() +
  geom_line(data = quantile.median, aes(x = age, y = female_quantiles0.5), color = "red",  linetype = "solid", linewidth = 1.5) +
  geom_line(data = quantile.P2.5,   aes(x = age, y = female_quantiles0.5), color = "red",  linetype = 2, linewidth = 1) +
  geom_line(data = quantile.P97.5,  aes(x = age, y = female_quantiles0.5), color = "red",  linetype = 2, linewidth = 1) +
  geom_line(data = quantile.median, aes(x = age, y = male_quantiles0.5),   color = "blue", linetype = "solid", linewidth = 1.5) +
  geom_line(data = quantile.P2.5,   aes(x = age, y = male_quantiles0.5),   color = "blue", linetype = 2, linewidth = 1) +
  geom_line(data = quantile.P97.5,  aes(x = age, y = male_quantiles0.5),   color = "blue", linetype = 2, linewidth = 1) +
  labs(x = "Age (years)", y = "SC strength",
       title = "Bootstrap 95% CI of Whole-brain averaged SC strength") +
  theme_classic() +
  theme(plot.title = element_text(size = 20, hjust = 0.5),
        axis.text  = element_text(size = 20),
        axis.title = element_text(size = 20))

ggsave(paste0(FigureFolder, "/ModelEvaluation/bootstrap_P50_95CI_WBmeanSC_bysex.tiff"),
       plot = p_sex, width = 16, height = 14, units = "cm")