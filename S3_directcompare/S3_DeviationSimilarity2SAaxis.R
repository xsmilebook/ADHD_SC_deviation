# Test whether the similarity to S-A axis differed by diagnosis
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
source(paste0(functionFolder.SCDev, "/gamminteraction.R"))
# load data
SCdata.sum75.merge.ADHDTD <- readRDS(paste0(interfileFolder, '/SCdata_Yeo',Yeoresolution, '_CV75_sumSCinvnode.deviations_TDtest_ADHD.rds'))
SCdata.sum75.merge.ADHDTD$sex <- as.factor(SCdata.sum75.merge.ADHDTD$sex)
SCdata.sum75.merge.ADHDTD$meanFD <- as.numeric(SCdata.sum75.merge.ADHDTD$meanFD)
SCdata.sum75.merge.ADHDTD$if_TD <- factor(SCdata.sum75.merge.ADHDTD$if_TD, levels = c("TD", "ADHD"))
SCdata.child <- SCdata.sum75.merge.ADHDTD %>% filter(age <= 10)
SCdata.adoles <- SCdata.sum75.merge.ADHDTD %>% filter(age > 10)

SCdeviationmat.child.ADHD <- SCdata.child %>% filter(if_TD=="ADHD") %>%
  dplyr::select(ends_with("_deviationZ") & starts_with("SC."))
sum(names(SCdeviationmat.child.ADHD)==paste0("SC.", 1:element_num, "_deviationZ"))
df <- data.frame(deviation=colMeans(SCdeviationmat.child.ADHD), SClabel=paste0("SC.", 1:element_num, "_deviationZ"))
axeslabels=c("VS.P", "SM.A", "VS.C", "SM.B", "DA.A", "DA.B", "FP.C", "DM.C", "VA.A", "TP", "FP.A", "VA.B", "DM.A", "FP.B", "DM.B")
matfig.ADHDdeviation <- plotmatrix("df", "deviation", Yeoresolution.delLM, Pvar=NA, NAcol="white", 
                      lmthr=NA, axeslabels=axeslabels, axeslabelsGap=F, linerange_frame=NA, PaletteSet=NA, Pvar.noFDR=NA)
SCrankcorr(df, "deviation", Yeoresolution.delLM)
# ds.resolution Interest.var r.spearman   p.spearman
# 1            15    deviation  0.4458701 3.334558e-07


# Compute the similarity with S-A axis
YeoSArank <- SCrankcorr(data.frame(rand=rnorm(element_num, 0, 1)), "rand", Yeoresolution.delLM, T)

SCdeviationmat <- SCdata.sum75.merge.ADHDTD %>% dplyr::select(ends_with("_deviationZ"), scanID)

SCdeviationmat.singleedge <- SCdeviationmat %>% dplyr::select(starts_with("SC."))
sum(names(SCdeviationmat.singleedge)==paste0("SC.", 1:element_num, "_deviationZ"))

similarity <- lapply(1:nrow(SCdeviationmat.singleedge), function(x){
  tmp <- SCdeviationmat.singleedge[x, ]
  corr <- cor(as.numeric(tmp), YeoSArank$SCrank, method = "spearman")
  return(corr)
})

similarity <- unlist(similarity)

df <- data.frame(similarity=similarity, scanID=SCdeviationmat$scanID)


SCdata.sum75.merge.ADHDTD <- SCdata.sum75.merge.ADHDTD %>% left_join(df, by="scanID")
describeBy(SCdata.sum75.merge.ADHDTD$similarity, group = SCdata.sum75.merge.ADHDTD$if_TD, mat = T)
# item group1 vars    n       mean        sd     median    trimmed       mad        min       max    range        skew
# X11    1     TD    1 2088 0.01124341 0.2652030 0.01379262 0.01309878 0.2947747 -0.6672894 0.6497326 1.317022 -0.04867316
# X12    2   ADHD    1 1114 0.04388723 0.2620671 0.04779499 0.04584808 0.2983116 -0.6608167 0.6468921 1.307709 -0.06415883


## 1. diagnostic effect
#################################
# All 8~15.5
region <- "similarity"
gamm.smooth.predict.interaction(region, dataname="SCdata.sum75.merge.ADHDTD", 
                                smooth_var="age", int_var="if_TD", covariates="sex+meanFD", knots=3, set_fx = T)
# T=2.765, P.disease=0.005724, bootstrap.P=0.004499550, bootstrap.chisq.disease=7.6484,  bootstrap_pvalue=0.287

# T=2.1655, bootstrap.P=0.03169683, bootstrap.chisq.disease=4.68884932078436,  bootstrap_pvalue=0.12708


# Childhood 8~10
dataname <- "SCdata.child"
smooth_var <- "age"
int_var <- "if_TD"

data.tmp <- get(dataname)
covariates <- "sex+age+meanFD"
modelformula <- as.formula(sprintf("%s ~ %s + %s",region, "if_TD", covariates))
modelformula.null<-as.formula(sprintf("%s ~ %s",region, covariates))
lm.model <- lm(modelformula, data = data.tmp)
lm.model.null <- lm(modelformula.null, data = data.tmp)
anova.result <- anovaPB(lm.model.null, lm.model, n.sim = 10000, ncpus=1)
anova.cov.pvalue<-anova.result$`Pr(>F)`[2]
mod.result <- summary(lm.model)
tvalue=mod.result$coefficients[2,3]
group_diff <- mod.result$coefficients[2,1]

resid_sd <- sigma(lm.model)
cohen_d <- group_diff / resid_sd

result.df <- data.frame(parcel=region, group="if_TD", dataname=dataname, tvalue=tvalue,
                        pvalue=anova.cov.pvalue, cohen_d=cohen_d)
#       parcel group     dataname   tvalue     pvalue   cohen_d
# 1 similarity if_TD SCdata.child 2.216511 0.02829717 0.1563991

# Adolescence 10~15.5
dataname <- "SCdata.adoles"
data.tmp <- get(dataname)
modelformula <- as.formula(sprintf("%1$s ~ if_TD+s(age, k=3, fx=T)+ sex+meanFD", region))
gam.model <- gam(modelformula, REML=TRUE, data = data.tmp)
modelformula.null <- as.formula(sprintf("%1$s ~ s(age, k=3, fx=T)+ sex+meanFD", region))
gam.model.null <- gam(modelformula.null, REML=TRUE, data = data.tmp)

anova.result <- anovaPB(gam.model.null, gam.model, n.sim = 10000, ncpus=1)
anova.cov.pvalue<-anova.result[2, 5]

lm.model <- gam.model
mod.result <- summary(gam.model)
group_diff <- mod.result$p.table[2,1]
tvalue <- mod.result$p.table[2,3]
resid_sd <- sigma(lm.model)
cohen_d <- group_diff / resid_sd

result.df <- data.frame(parcel=region, group="if_TD", dataname=dataname, tvalue=tvalue,
                        pvalue=anova.cov.pvalue, cohen_d=cohen_d)
#    parcel     group   dataname   tvalue    pvalue    cohen_d
# 1 similarity if_TD SCdata.adoles 1.549126 0.1245875 0.07039331

######################################


## 2. Age effect
######################################
SCdata.ADHD <- SCdata.sum75.merge.ADHDTD %>% filter(if_TD == "ADHD")
SCdata.TD <- SCdata.sum75.merge.ADHDTD %>% filter(if_TD == "TD")

result.age.ADHD <- gamm.fit.smooth(region, dataname="SCdata.ADHD", smooth_var="age", covariates = "sex+meanFD", knots=3, 
                                   set_fx = T, stats_only = FALSE, mod_only = FALSE)
# parcel       gamm.smooth.F      gamm.smooth.pvalue     partialRsq            bootstrap.zvalue    bootstrap_pvalue    
# cor "similarity" "8.77651473370244" "0.000165294375367105" "-0.0132340634886196" "-3.54011018851328" "0.0003999600039996"
# correstimate         corrp                  gamm.edf meanderv2              change.onset       peak.change       
# cor "-0.125696417893829" "2.58427485400727e-05" "2"      "-0.00233509915733154" "10.9624624624625" "13.4086586586587"

result.age.TD <- gamm.fit.smooth(region, dataname="SCdata.TD", smooth_var="age", covariates = "sex+meanFD", knots=3, 
                                   set_fx = T, stats_only = FALSE, mod_only = FALSE)
# parcel       gamm.smooth.F      gamm.smooth.pvalue    partialRsq             bootstrap.zvalue    bootstrap_pvalue    
# cor "similarity" "5.51543451879229" "0.00408316878902565" "-0.00261089297162481" "-2.90269817495577" "0.0036996300369963"
# correstimate          corrp                 gamm.edf meanderv2             change.onset       peak.change       
# cor "-0.0680838491373098" "0.00185311366692299" "2"      "0.00234756008917154" "8.91666666666667" "8.99904070737404"



################################

## 3. Associate with symptoms
######################################

symptomvar <- c("cbcl_scr_syn_external_t", "cbcl_scr_dsm5_adhd_t","cbcl_scr_syn_attention_t")
symptomresult <- list(); i=0
for (VOI in symptomvar){
  print(VOI)
  # resultsum <- mclapply(1:2, function(x){
  #   dataname <- c("SCdata.TD", "SCdata.ADHD")[x]
  #   gamresult<-gamm.smooth.predict.covariateinteraction(region=VOI, dataname, smooth_var="age", int_var="similarity", 
  #                                                       int_var.predict.percentile=0.1, covariates="sex+meanFD", knots=3, 
  #                                                       set_fx = T, increments=100, stats_only=F, if_residual=FALSE)
  #   gamresult<-as.data.frame(gamresult)
  #   gamresult$dataset <- dataname
  #   
  #   return(gamresult)
  # }, mc.cores = 5)
  
  resultsum <- list()
  for (x in 1:2){
    dataname <- c("SCdata.TD", "SCdata.ADHD")[x]
    gamresult<-gamm.smooth.predict.covariateinteraction(region=VOI, dataname, smooth_var="age", int_var="similarity",
                                                        int_var.predict.percentile=0.1, covariates="sex+meanFD", knots=3,
                                                        set_fx = T, increments=100, stats_only=F, if_residual=FALSE)
    gamresult<-as.data.frame(gamresult[[1]])
    gamresult$dataset <- dataname
    resultsum[[x]] <- gamresult
  }
  
  
  resultsum.df <- do.call(rbind, lapply(resultsum, function(x) data.frame(x)))
  resultsum.df[,c(3:18)]<-lapply(resultsum.df[,c(3:18)], as.numeric)
  i=i+1
  symptomresult[[i]] <- resultsum.df
}

symptomresult.df <- do.call(rbind, symptomresult)



saveRDS(symptomresult, paste0(interfileFolder, '/Deviation_similarity', VOI,'_Age_CV75_Yeo', Yeoresolution,'.rds'))











