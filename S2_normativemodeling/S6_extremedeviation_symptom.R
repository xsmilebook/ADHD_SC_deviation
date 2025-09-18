rm(list=ls())
library(ggplot2)
library(tidyverse)
library(mgcv)
library(openxlsx)
library(parallel)
library(gamlss)
library(scales)
library(tableone)
library(bootES)
library(reshape)
# input directory
wd <- getwd()
if (str_detect(wd, "cuizaixu_lab")){
  demopath <- '/ibmgpfs/cuizaixu_lab/xuxiaoyu/Normative_model/demography'
  interfileFolder <- "/ibmgpfs/cuizaixu_lab/xuxiaoyu/Normative_model/interfileFolder"
  functionFolder <- "/ibmgpfs/cuizaixu_lab/xuxiaoyu/Normative_model/functions"
  resultFolder <- "/ibmgpfs/cuizaixu_lab/xuxiaoyu/Normative_model/results"
  functionFolder.SCDev <- "/ibmgpfs/cuizaixu_lab/xuxiaoyu/SCdevelopment/Rcode_SCdevelopment/gamfunction"
}else{
  demopath <- 'D:/xuxiaoyu/DMRI_network_development/Normative_model/demography'
  FigureFolder <- 'D:/xuxiaoyu/DMRI_network_development/Normative_model/Figures/SA12/CV75'
  interfileFolder <- "D:/xuxiaoyu/DMRI_network_development/Normative_model/interfileFolder"
  functionFolder <- "D:/xuxiaoyu/DMRI_network_development/Normative_model/functions"
  functionFolder.SCDev <- "D:/xuxiaoyu/DMRI_network_development/SC_development/Rcode_SCdevelopment/gamfunction"
  resultFolder <- "D:/xuxiaoyu/DMRI_network_development/Normative_model/results"
}
#D:/xuxiaoyu/DMRI_network_development
#/Users/xuxiaoyu_work/Cuilab/DMRI_network_development/Normative_model
# load data
# set resolution
ds.resolution <- 12
element_num <- ds.resolution*(ds.resolution+1)/2

# load data
SCdata.sum75.merge <- readRDS(paste0(interfileFolder, '/SCdata_SA',ds.resolution, '_CV75_sumSCinvnode.sum.msmtcsd.merge.rds'))
SCdata.sum25.merge <- readRDS(paste0(interfileFolder, '/SCdata_SA',ds.resolution, '_CV25_sumSCinvnode.sum.msmtcsd.merge.rds'))
source(paste0(functionFolder.SCDev, "/SCrankcorr.R"))
source(paste0(functionFolder, "/Construct_gamlss_set.R"))
sum_deviation <- readRDS(paste0(interfileFolder, '/SCdata_SA',ds.resolution, '_CV75_sumSCinvnode.deviations.rds'))
sum_deviation <- sum_deviation %>% filter(SC.78 > 1)

# Add extreme deviation number
sum_deviation <- sum_deviation %>% mutate(extreme_p_num = rowSums(across(ends_with("_deviationZ"), ~ . > 1.96)),
                                          extreme_n_num = rowSums(across(ends_with("_deviationZ"), ~ . < -1.96)))
sum_deviation$extreme_num <- sum_deviation$extreme_p_num + sum_deviation$extreme_n_num
describeBy(sum_deviation$extreme_num, group=sum_deviation$ADHD_condition, mat=T)
# item group1 vars    n     mean       sd median  trimmed    mad min max range     skew  kurtosis         se
# X11    1      0    1 6690 3.838714 4.891268      2 2.895366 2.9652   0  60    60 3.072828 15.715292 0.05980099
# X12    2      1    1  235 4.085106 5.116714      2 3.063492 2.9652   0  26    26 1.969529  4.002625 0.33377761
# X13    3      2    1  865 3.937572 4.703799      2 3.031746 2.9652   0  37    37 2.410250  7.845804 0.15993397
t.test(extreme_num ~ if_ADHD, data=sum_deviation)
# T = -0.834, P = 0.405
saveRDS(sum_deviation, paste0(interfileFolder, '/SCdata_SA',ds.resolution, '_CV75_sumSCinvnode.deviations.rds'))

extreme_proportion.ADHD <- sum(sum_deviation$extreme_num[sum_deviation$if_ADHD==1] > 0) / sum(sum_deviation$if_ADHD==1)
print(paste(round(extreme_proportion.ADHD,2), "ADHD all patients have at least one SC edge with extreme deviation."))
# 0.82 ADHD all patients have at least one SC edge with extreme deviation.
extreme_proportion.TD <- sum(sum_deviation$extreme_num[sum_deviation$if_TD==1] > 0) / sum(sum_deviation$if_TD==1)
print(paste(round(extreme_proportion.TD,2), "TD participants have at least one SC edge with extreme deviation."))
# 0.8 TD participants have at least one SC edge with extreme deviation.

# 1. Subjects with no extreme SC have less symptoms?
sum_deviation$ADHD_condition <- factor(sum_deviation$ADHD_condition, levels=c(0,1,2), labels=c("no", "partial", "present"))
sum_deviation$adhd_medication <- as.factor(sum_deviation$adhd_medication)
sum_deviation$sex <- factor(sum_deviation$sex, levels=c(1,2), labels=c("M", "F"))
sum_deviation_noextremeSC <- sum_deviation %>% filter(extreme_num == 0)
sum_deviation_noextremeSC_ADHD <- sum_deviation_noextremeSC %>% filter(if_ADHD == 1)
sum_deviation$if_extreme <- as.factor(sum_deviation$extreme_num > 0)
sum_deviation$if_extreme.p <- as.factor(sum_deviation$extreme_p_num > 0)
sum_deviation$if_extreme.n <- as.factor(sum_deviation$extreme_n_num > 0)
sum_deviation[,c("ADHD_condition_next1year", "ADHD_condition_next2year", "ADHD_condition_next3year")] <- lapply(sum_deviation[,c("ADHD_condition_next1year", "ADHD_condition_next2year", "ADHD_condition_next3year")], as.factor)

sum_deviation.TD <- sum_deviation %>% filter(if_TD == 1)
sum_deviation.ADHD <- sum_deviation %>% filter(if_ADHD == 1)

# CBCL vars
Symptom.var <- names(sum_deviation %>% select(starts_with("cbcl") & ends_with("_t")))
# trajectories, ADHD conditions
Interest.vars <- c(Symptom.var,"GENERAL", "ADHD_condition", "ADHD_trajectory", "adhd_medication", 
                   "age", "sex", "diff_value", "ADHD_condition_next1year", "ADHD_condition_next2year", "ADHD_condition_next3year",
                   "meanFD", "handedness", "income.adj")
factorVars <- c("ADHD_condition", "ADHD_trajectory", "adhd_medication", "sex", "handedness", "ADHD_condition_next1year", "ADHD_condition_next2year", "ADHD_condition_next3year")
### ADHD
## 1) ADHD with extreme (p / n) VS ADHD without extreme (p / n)
CBCL_ADHD.tab1.extremeB <- CreateTableOne(Interest.vars, strata = "if_extreme", data=sum_deviation.ADHD, 
                                          factorVars = factorVars, test = TRUE)
CBCL_ADHD.tab1.extremeBout <- print(CBCL_ADHD.tab1.extremeB, exact = c("status", "stage"), smd=T, quote = T)
write.csv(CBCL_ADHD.tab1.extremeBout, paste0(resultFolder, "/CBCL_ADHD.tab1.extremeBout_SC", element_num, "_CV75.csv"))

## 2) ADHD with extreme (p) VS ADHD without extreme (p)
CBCL_ADHD.tab1.extremeP <- CreateTableOne(Interest.vars, strata = "if_extreme.p", data=sum_deviation.ADHD, 
                                          factorVars = factorVars, test = TRUE)
CBCL_ADHD.tab1.extremePout <- print(CBCL_ADHD.tab1.extremeP, exact = c("status", "stage"), smd=T, quote = T)
write.csv(CBCL_ADHD.tab1.extremePout, paste0(resultFolder, "/CBCL_ADHD.tab1.extremePout_SC", element_num, "_CV75.csv"))

## 3) ADHD with extreme (n) VS ADHD without extreme (n)
CBCL_ADHD.tab1.extremeN <- CreateTableOne(Interest.vars, strata = "if_extreme.n", data=sum_deviation.ADHD, 
                                          factorVars = factorVars, test = TRUE)
CBCL_ADHD.tab1.extremeNout <- print(CBCL_ADHD.tab1.extremeN, exact = c("status", "stage"), smd=T, quote = T)
write.csv(CBCL_ADHD.tab1.extremeNout, paste0(resultFolder, "/CBCL_ADHD.tab1.extremeNout_SC", element_num, "_CV75.csv"))

### TD
## 1) TD with extreme (p / n) VS TD without extreme (p / n)
CBCL_TD.tab1.extremeB <- CreateTableOne(Interest.vars, strata = "if_extreme", data=sum_deviation.TD, 
                                          factorVars = factorVars, test = TRUE)
CBCL_TD.tab1.extremeBout <- print(CBCL_TD.tab1.extremeB, exact = c("status", "stage"), smd=T, quote = T)
write.csv(CBCL_TD.tab1.extremeBout, paste0(resultFolder, "/CBCL_TD.tab1.extremeBout_SC", element_num, "_CV75.csv"))

## 2) TD with extreme (p) VS TD without extreme (p)
CBCL_TD.tab1.extremeP <- CreateTableOne(Interest.vars, strata = "if_extreme.p", data=sum_deviation.TD, 
                                          factorVars = factorVars, test = TRUE)
CBCL_TD.tab1.extremePout <- print(CBCL_TD.tab1.extremeP, exact = c("status", "stage"), smd=T, quote = T)
write.csv(CBCL_TD.tab1.extremePout, paste0(resultFolder, "/CBCL_TD.tab1.extremePout_SC", element_num, "_CV75.csv"))

## 3) TD with extreme (n) VS TD without extreme (n)
CBCL_TD.tab1.extremeN <- CreateTableOne(Interest.vars, strata = "if_extreme.n", data=sum_deviation.TD, 
                                          factorVars = factorVars, test = TRUE)
CBCL_TD.tab1.extremeNout <- print(CBCL_TD.tab1.extremeN, exact = c("status", "stage"), smd=T, quote = T)
write.csv(CBCL_TD.tab1.extremeNout, paste0(resultFolder, "/CBCL_TD.tab1.extremeNout_SC", element_num, "_CV75.csv"))


# 2. Extreme SC effects are region-specific?
# Logistic regression will be used.
symptom_ROI <- c("cbcl_scr_syn_attention_t", "cbcl_scr_syn_external_t", "cbcl_scr_syn_totprob_t", "cbcl_scr_dsm5_adhd_t")

### 1) ADHD group
num_cores <- detectCores() - 8
cl <- makeCluster(num_cores)
result <- parLapply(cl, 1:10, function(x) x^2)
clusterEvalQ(cl, {
  library(tidyverse)
  library(bootES) # Load the package
})
permresult.df.list <- list()
for (j in 1:length(symptom_ROI)){
  Symptom.var.tmp <- symptom_ROI[j]
  clusterExport(cl, varlist = c("sum_deviation.ADHD", "Symptom.var.tmp"), envir = .GlobalEnv)
  bootresult <- parLapply(cl, 1:element_num, function(i){
    # positive
    sum_deviation.ADHD$independentvar.P <- as.factor(sum_deviation.ADHD[[paste0("SC.", i, "_deviationZ")]] > 1.96)
    data.tmp <- sum_deviation.ADHD %>% select(c(Symptom.var.tmp, "independentvar.P"))
    bootres.tmp <- bootES(data = data.tmp, R=10000, data.col = Symptom.var.tmp, ci.type = "norm",
                         group.col = "independentvar.P", contrast = c("FALSE", "TRUE"))
    p.value.tmp <- with(bootres.tmp, pnorm(abs((2*t0 - mean(t) - 1) / sqrt(var(t)[1,1])), lower.tail=F)*2)
    #bootres.df <- NA
    bootres.df <- data.frame(SClabel = paste0("SC.", i), boot.T = bootres.tmp$t0, boot.CI2.5 = bootres.tmp$bounds[[1]],
                             boot.CI97.5 = bootres.tmp$bounds[[2]], P.value = p.value.tmp, N_True = sum(data.tmp$independentvar.P==T),
                             N_False = sum(data.tmp$independentvar.P==F), CBCLvar = Symptom.var.tmp)
    
    return(bootres.df)
  })
  
  bootresult.df <- do.call(rbind, bootresult)
  bootresult.df$P.value.fdr <- p.adjust(bootresult.df$P.value, method="fdr")
  bootresult.df.list[[j]] <- bootresult.df
}


## plot
Matrix.tmp <- matrix(NA, ds.resolution, ds.resolution)
linerange_frame<-data.frame(x=c(0.5,12+0.5), ymin =rep(-12-0.5, times=2), ymax =rep(-0.5, times=2),
                            y=c(-0.5, -12-0.5), xmin=rep(0.5, times=2), xmax=rep(12+0.5, times=2))

for (j in 1:length(symptom_ROI)){
  bootresult.df <- bootresult.df.list[[j]]
  Symptom.var.tmp <- symptom_ROI[j]
  lmthr <- max(abs(bootresult.df$boot.T))
  
  Matrix.tmp[lower.tri(Matrix.tmp, diag = T)] <- bootresult.df$boot.T
  Matrix.tmp[upper.tri(Matrix.tmp)] <- t(Matrix.tmp)[upper.tri(Matrix.tmp)]
  colnames(Matrix.tmp) <-seq(1, 12)
  rownames(Matrix.tmp) <-seq(1, 12)
  matrixtmp.df <- as.data.frame(Matrix.tmp)
  matrixtmp.df$nodeid <- seq(1, 12)
  matrixtmp.df.melt <- melt(matrixtmp.df,id.vars=c("nodeid"))
  matrixtmp.df.melt$variable<-as.numeric(matrixtmp.df.melt$variable)
  matrixtmp.df.melt$nodeid<-0-matrixtmp.df.melt$nodeid
  matrixtmp.df.melt$value<-as.numeric(matrixtmp.df.melt$value)
  
  Matrix.tmp.sig <- matrix(NA, nrow = 12, ncol=12)
  Matrix.tmp.sig[lower.tri(Matrix.tmp.sig, diag = T)] <- (bootresult.df$P.value<0.05)
  Matrix.tmp.sig[upper.tri(Matrix.tmp.sig)] <- t(Matrix.tmp.sig)[upper.tri(Matrix.tmp.sig)]
  colnames(Matrix.tmp.sig) <-seq(1, 12)
  rownames(Matrix.tmp.sig) <-seq(1, 12)
  matrixtmp.df.sig <- as.data.frame(Matrix.tmp.sig)
  matrixtmp.df.sig$nodeid <- seq(1, 12)
  matrixtmp.df.sig.melt <- melt(matrixtmp.df.sig,id.vars=c("nodeid"))
  matrixtmp.df.sig.melt$variable<-as.numeric(matrixtmp.df.sig.melt$variable)
  matrixtmp.df.sig.melt$nodeid<-0-matrixtmp.df.sig.melt$nodeid
  matrixtmp.df.sig.melt$value<-as.numeric(matrixtmp.df.sig.melt$value)
  matrixtmp.df.sig.melt <- matrixtmp.df.sig.melt[-which(matrixtmp.df.sig.melt$value==0),]
  
  Symptom.var.core <- str_split_i(Symptom.var.tmp, "_", 4)
  plottitle <- paste("Difference of CBCL", Symptom.var.core, "between\n with and without extreme positive deviation.")
  
  ggplot(data=matrixtmp.df.melt)+
    geom_tile(aes(x=variable, y=nodeid, fill=value, color=value))+
    scale_fill_distiller(type="seq", palette = "RdBu",limits=c(-lmthr, lmthr), direction=-1,na.value = "grey")+
    scale_color_distiller(type="seq", palette = "RdBu",limits=c(-lmthr, lmthr),direction=-1, na.value = "grey")+
    geom_text(data =matrixtmp.df.sig.melt, aes(x=variable, y=nodeid, label = "*"), vjust = 0.7, hjust = 0.5, size=8)+
    geom_linerange(data=linerange_frame, aes(y=y, xmin =xmin, xmax =xmax), color="black", linewidth=0.5)+
    geom_linerange(data=linerange_frame, aes(x=x, ymin =ymin, ymax =ymax), color="black", linewidth=0.5)+
    geom_segment(aes(x = 0.5 , y = -0.5 , xend = 12+0.5 ,yend = -12-0.5), color="black", linewidth=0.5)+
    ggtitle(label = plottitle)+labs(x=NULL, y=NULL, color="T value", fill="T value")+
    scale_y_continuous(breaks=NULL, labels = NULL)+
    scale_x_continuous(breaks=NULL, labels = NULL)+
    theme(axis.line = element_blank(), 
          #axis.ticks=element_line(linewidth = 0),
          axis.text.x=element_text(size=12, angle=45, hjust=1), 
          axis.text.y=element_text(size=12, angle=315, hjust=1,vjust=1),
          axis.title =element_text(size=18),
          plot.title = element_text(size=18, hjust = 0.5),
          legend.title=element_text(size=18),
          legend.text=element_text(size=18), 
          panel.background=element_rect(fill=NA),
          panel.grid.major=element_line(linewidth = 0), 
          panel.grid.minor=element_line(linewidth = 1))
  
  ggsave(paste0(FigureFolder, "/deviationPlots/", Symptom.var.tmp, "_diff_ADHD_extreme.p.tiff"), height=18, width=20, units = "cm")

}








