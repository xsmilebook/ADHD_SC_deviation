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
homepath <- "D:/code/ADHD_SC_deviation"
demopath <- file.path(homepath, "data", 'demography')
interfileFolder <- file.path(homepath, "data", 'interfileFolder', "ABCD")
functionFolder <- file.path(homepath, "src", "functions")
resultFolder <- file.path(homepath, "reports", "results", "ABCD")
functionFolder.SCDev <- file.path(homepath, "src", "gamfunction")
FigureFolder <- paste0(homepath, '/reports/figures/ABCD/Yeo', Yeoresolution,'/CV75')

source(paste0(functionFolder, "/plotmatrix.R"))
source(paste0(functionFolder.SCDev, "/SCrankcorr.R"))
source(paste0(functionFolder.SCDev, "/gammsmooth.R"))
source(paste0(functionFolder, "/index2network.R"))

# load data
SCdata.sum75.merge.ADHDTD <- readRDS(paste0(interfileFolder, "/SCdata_Yeo", Yeoresolution, 
                                            "_CV75_sumSCinvnode.deviations_TDtest_ADHD_match_merge.rds"))
SCdata.sum75.merge.ADHDTD$sex <- as.factor(SCdata.sum75.merge.ADHDTD$sex)
SCdata.sum75.merge.ADHDTD$if_TD <- factor(SCdata.sum75.merge.ADHDTD$if_TD, levels = c("TD", "ADHD"))

SCdata.TD <- SCdata.sum75.merge.ADHDTD %>% filter(if_TD=="TD")
smooth_var<-"age"
#sigedge.8_10 <- diagnosis.rs.8_10$parcel[diagnosis.rs.8_10$p.fdr < 0.05]
#sigedge.all <- diffbyage$parcel
#sigedge.all <- gsub("_h", "_deviationZ", sigedge.all)
sigedge.all <- paste0("SC.", 1:element_num, "_deviationZ")
symptomvar <- c("cbcl_scr_syn_external_t", "cbcl_scr_dsm5_adhd_t",  "cbcl_scr_syn_attention_t")
dependent.var <- c(sigedge.all, symptomvar)
dataname <- "SCdata.TD"
if (str_detect(wd, "cuizaixu_lab")){
  resultsum <- mclapply(1:length(dependent.var), function(x){
    region<-dependent.var[x]
    
    if (str_detect(region, "cbcl_")){
      covariates<-"sex"
    }else{
      covariates<-"sex+meanFD"
    }
    
    gamresult<-gamm.fit.smooth(region, dataname, smooth_var, covariates, knots=3, set_fx=TRUE, stats_only = T, mod_only=FALSE)
    gamresult<-as.data.frame(gamresult)
    
    return(gamresult)
  }, mc.cores = 40)
  
  gamresultsum.df <- do.call(rbind, lapply(resultsum, function(x) data.frame(x)))
  gamresultsum.df[,c(2:10)]<-lapply(gamresultsum.df[,c(2:10)], as.numeric)
  
  saveRDS(gamresultsum.df, paste0(interfileFolder, 
                                  '/AllSCDeviation_NormSymptom_Age_Yeo', 
                                  Yeoresolution,'_CV75_TDmatched.rds'))
}


gamresultsum.df.ADHD <- readRDS(paste0(interfileFolder, '/AllSCDeviation_NormSymptom_Age_Yeo', Yeoresolution,'_CV75_ADHD.rds'))
gamresultsum.df.TD <- readRDS(paste0(interfileFolder, '/AllSCDeviation_NormSymptom_Age_Yeo', Yeoresolution,'_CV75_TDmatched_merge.rds'))

gamresultsum.df.ADHD.SC <- gamresultsum.df.ADHD[str_detect(gamresultsum.df.ADHD$parcel, "SC."),]
gamresultsum.df.ADHD.SC$bootstrap_p.fdr <- p.adjust(gamresultsum.df.ADHD.SC$bootstrap_pvalue, method="fdr")
gamresultsum.df.TD.SC <- gamresultsum.df.TD[str_detect(gamresultsum.df.TD$parcel, "SC."),]
sigedge.all <- read.csv(paste0(resultFolder, "/ABCD_sigAgeChangedEdges.csv"))
gamresultsum.df.TD.SC$bootstrap_p.fdr <- p.adjust(gamresultsum.df.TD.SC$bootstrap_pvalue, method="fdr")

interaction.df <- readRDS(paste0(interfileFolder, "/Deviation_TDtest-ADHDMatched_SC", element_num,"_merge.rds"))
interaction.df <- interaction.df[seq(from=2, to=2*element_num, by = 2),]

interaction.df <- interaction.df %>% filter(parcel %in% sigedge.all$x)
interaction.df$bootstrap_pvalue.fdr <- p.adjust(interaction.df$bootstrap_pvalue, method="fdr")
df <- data.frame(parcel=sigedge.all$x, Label=rep(NA, length(sigedge.all$x)), F.TD=rep(NA, length(sigedge.all$x)),
                 F.ADHD=rep(NA, length(sigedge.all$x)),Pfdr.TD=rep(NA, length(sigedge.all$x)), 
                 Pfdr.ADHD=rep(NA, length(sigedge.all$x)), Chisq.interaction=rep(NA, length(sigedge.all$x)),
                 Pfdr.interaction=rep(NA, length(sigedge.all$x)))

for (i in 1:nrow(df)){
  index <- gsub("_deviationZ", "", df$parcel[i])
  index <- as.numeric(gsub("SC.", "", index))
  
  networklabel <- paste0(index2network(index)[2], "-", index2network(index)[1])
  
  df$Label[i] <- networklabel
  df$F.TD[i] <- gamresultsum.df.TD.SC$gamm.smooth.F[gamresultsum.df.TD.SC$parcel==df$parcel[i]]
  df$F.ADHD[i] <- gamresultsum.df.ADHD.SC$gamm.smooth.F[gamresultsum.df.ADHD.SC$parcel==df$parcel[i]]
  df$Pfdr.TD[i] <- gamresultsum.df.TD.SC$bootstrap_p.fdr[gamresultsum.df.TD.SC$parcel==df$parcel[i]]
  df$Pfdr.ADHD[i] <- gamresultsum.df.ADHD.SC$bootstrap_p.fdr[gamresultsum.df.ADHD.SC$parcel==df$parcel[i]]
  
  df$Chisq.interaction[i] <- interaction.df$bootstrap_chisq[interaction.df$parcel==df$parcel[i]]
  df$Pfdr.interaction[i] <- interaction.df$bootstrap_pvalue.fdr[interaction.df$parcel==df$parcel[i]]
  
}


write.csv(df, paste0(resultFolder, "/SCdeviationDevelopment_sig_STats_matched.csv"), row.names = F)

interaction.df2 <- interaction.df %>% select("parcel", "bootstrap_pvalue", "bootstrap_chisq")
for (i in 1:nrow(interaction.df2)){
  parcel <- interaction.df2$parcel[i]
  index <- gsub("SC.", "", parcel)
  index <- gsub("_deviationZ", "", index)
  
  label <- paste0(index2network(as.numeric(index))[1], "-", index2network(as.numeric(index))[2])
  interaction.df2$parcel[i] <- label
}

write.csv(interaction.df2, paste0(resultFolder, "/interaction_Age_Diagnosis_matched_14edges.csv"), row.names = F)

#### Effect size
lmthr <- max(abs(c(gamresultsum.df.ADHD.SC$partialRsq, gamresultsum.df.TD.SC$partialRsq)))
# ADHD
SCrankcorr(gamresultsum.df.ADHD.SC, "partialRsq", Yeoresolution.delLM)
# ADHD
# 15	partialRsq	r=-0.4202567	P=1.764467e-06	

# Matrix
gamresultsum.df.ADHD.SC$sigo <- gamresultsum.df.ADHD.SC$parcel %in% sigedge.all$x
gamresultsum.df.ADHD.SC$sig <- gamresultsum.df.ADHD.SC$parcel %in% interaction.df$parcel[interaction.df$bootstrap_pvalue.fdr < 0.1]
gamresultsum.df.ADHD.SC <- gamresultsum.df.ADHD.SC %>% mutate(
  pvalueo = case_when(sigo==T ~ 0.01,
                      .default = 0.5),
  pvalue = case_when(sig==T ~ 0.01,
                      .default = 0.5)
)

axeslabels=c("VS.P", "SM.A", "VS.C", "SM.B", "DA.A", "DA.B", "FP.C", "DM.C", "VA.A", "TP", "FP.A", "VA.B", "DM.A", "FP.B", "DM.B") 
axeslabelsGap=T
PaletteSet <- list(Name="RdBu", drirection=-1, lmmin = -lmthr, lmmax = lmthr, anglex=45, angley=45,hjustx = 1, hjusty = 1, vjustx = 1, vjusty=0.3)

FigintR2 <- plotmatrix(dataname="gamresultsum.df.ADHD.SC", "partialRsq", ds.resolution=Yeoresolution.delLM, Pvar="pvalue", NAcol="white", 
                       lmthr=lmthr, axeslabels, axeslabelsGap=F, linerange_frame=NA, PaletteSet=PaletteSet, Pvar.noFDR="pvalueo")

print(FigintR2)

ggsave(paste0(FigureFolder, "/SCstrength_Age/Deviationcompare/Matrix", Yeoresolution.delLM,"_DeviationAge_Rsq_ADHD_MatchedPfdr01.tiff"),
       FigintR2, height=14, width=16, units = "cm")


# TD
SCrankcorr(gamresultsum.df.TD.SC, "partialRsq", Yeoresolution.delLM)
# TD
# 15	partialRsq	r=-0.08518825  P=0.3549108	

# scatter plot
df <- SCrankcorr(gamresultsum.df.TD.SC, "partialRsq", Yeoresolution.delLM, T)

scatter <- ggplot(data=df)+
  geom_point(aes(x=SCrank, y=partialRsq, color=partialRsq), size=5)+
  geom_smooth(aes(x=SCrank, y=partialRsq), method ="lm", color="black", linewidth=1.2)+
  scale_color_distiller(type="seq", palette = "RdBu", limits=c(-lmthr, lmthr), direction = -1)+
  labs(x="S-A connectional axis rank", y = expression("Age effect (partial " * italic(R)^2 * ")"))+
  #scale_y_continuous(breaks = c(-0.015, -0.010, -0.005,0, 0.005, 0.010), labels = c(-15,-10, -5, 0, 5, 10))+
  #scale_y_continuous(breaks = c(-0.003,0, 0.003, 0.006), labels = c(-3, 0, 3, 6))+
  scale_x_continuous(breaks = c(0,20,40,60,80,100,120))+
  theme_classic()+
  theme(axis.text=element_text(size=23, color="black"), 
        axis.title =element_text(size=23),aspect.ratio = 1,
        axis.line = element_line(linewidth=0.6),axis.ticks= element_line(linewidth=0.6),
        plot.title = element_text(size=20, hjust = 0.5, vjust=2),
        plot.background=element_rect(fill="transparent"),
        panel.background=element_rect(fill="transparent"),
        legend.position = "none")
ggsave(paste0(FigureFolder, '/SCstrength_Age/Deviationcompare/DeviationRsqTD_SCrankcorr_Matched.tiff'),scatter, width=15, height =15, units = "cm")


# Matrix
lmthr <- max(abs(c(gamresultsum.df.ADHD.SC$partialRsq, gamresultsum.df.TD.SC$partialRsq)))
axeslabels=c("VS.P", "SM.A", "VS.C", "SM.B", "DA.A", "DA.B", "FP.C", "DM.C", "VA.A", "TP", "FP.A", "VA.B", "DM.A", "FP.B", "DM.B") 
axeslabelsGap=T
PaletteSet <- list(Name="RdBu", drirection=-1, lmmin = -lmthr, lmmax = lmthr, anglex=45, angley=45,hjustx = 1, hjusty = 1, vjustx = 1, vjusty=0.3)

FigintR2 <- plotmatrix(dataname="gamresultsum.df.TD.SC", "partialRsq", ds.resolution=Yeoresolution.delLM, Pvar="bootstrap_p.fdr", NAcol="white", lmthr=lmthr, axeslabels, axeslabelsGap=F, linerange_frame=NA, PaletteSet=PaletteSet, Pvar.noFDR=NA)

print(FigintR2)

ggsave(paste0(FigureFolder, "/SCstrength_Age/Deviationcompare/Matrix", Yeoresolution.delLM,"_DeviationAge_Rsq_TD_Matched.tiff"),
       FigintR2, height=14, width=16, units = "cm")


## Interaction effect scatter
gamresultsum.df.TD.SC$if_TD <- "TD"
gamresultsum.df.ADHD.SC$if_TD <- "ADHD"
var.com <- intersect(names(gamresultsum.df.TD.SC), names(gamresultsum.df.ADHD.SC))

gamresult.SC <- rbind(gamresultsum.df.TD.SC[,var.com], gamresultsum.df.ADHD.SC[,var.com])

Yeo_10 <- data.frame(SClabel=paste0("SC.", 1:element_num), rand=rnorm(element_num))
Yeo_10 <- SCrankcorr(Yeo_10, "rand", Yeoresolution.delLM, T)
Yeo_10$SClabel <- paste0("SC.", 1:element_num, "_deviationZ")
gamresult.SC <- gamresult.SC %>% left_join(Yeo_10, join_by(parcel==SClabel))
gamresult.SC$if_TD <- factor(gamresult.SC$if_TD, levels=c("TD", "ADHD"))
intmodel <- lm(partialRsq ~ SCrank*if_TD, data = gamresult.SC)
nullmodel <- lm(partialRsq ~ SCrank+if_TD, data = gamresult.SC)
summary(intmodel)
# Coefficients:
#                  Estimate Std. Error t value Pr(>|t|)    
# (Intercept)      -3.716e-04  5.074e-04  -0.732  0.46471    
# SCrank           -7.029e-06  7.279e-06  -0.966  0.33520    
# if_TDADHD         6.094e-04  7.176e-04   0.849  0.39661    
# SCrank:if_TDADHD -4.027e-05  1.029e-05  -4.205  3.72e-05 ***

anova(nullmodel, intmodel) # F = 17.678, P=3.716e-05 ***


scatter <- ggplot(data=gamresult.SC)+
  geom_point(aes(x=SCrank, y=partialRsq, color=if_TD), size=1.8)+
  geom_smooth(aes(x=SCrank, y=partialRsq, fill=if_TD), method ="lm", color="black", alpha=0.3, linewidth=1.2)+
  scale_color_manual(values = c("gray", "#B2182B"))+
  scale_fill_manual(values = c("gray", "#B2182B"))+
  labs(x="S-A connectional axis rank", y=expression("Age effect (partial " * italic(R)^2 * ")"), color=NULL, fill=NULL)+
  #scale_y_continuous(breaks = c(-0.015, -0.010, -0.005,0, 0.005, 0.010), labels = c(-15,-10, -5, 0, 5, 10))+
  #scale_y_continuous(breaks = c(-0.003,0, 0.003, 0.006), labels = c(-3, 0, 3, 6))+
  scale_x_continuous(breaks = c(0,20,40,60,80,100,120))+
  theme_classic()+
  theme(axis.text=element_text(size=20, color="black"), 
        axis.title =element_text(size=20),aspect.ratio = 0.95,
        axis.line = element_line(linewidth=0.6),axis.ticks= element_line(linewidth=0.6),
        plot.title = element_text(size=20, hjust = 0.5, vjust=2),
        plot.background=element_rect(fill="transparent"),
        panel.background=element_rect(fill="transparent"),
        legend.text = element_text(size=20, color="black"))

ggsave(paste0(FigureFolder, '/SCstrength_Age/Deviationcompare/DeviationRsq_SCrankcorr_Interaction_Matched.tiff'),scatter, width=20.5, height =17, units = "cm")




