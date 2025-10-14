library(ggplot2)
library(tidyverse)
library(mgcv)
library(readxl)
library(openxlsx)
library(parallel)
library(gamlss)
library(scales)
library(MatchIt)
library(tableone)
rm(list=ls())
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

# load data
source(paste0(functionFolder, "/Construct_gamlss_set.R"))
source(paste0(functionFolder, "/plotmatrix.R"))
source(paste0(functionFolder.SCDev, '/SCrankcorr.R'))

deviationdf <- readRDS(paste0(interfileFolder, '/SCdata_Yeo',Yeoresolution, '_CV75_sumSCinvnode.deviations_TDtest_ADHD.rds'))
deviationdf$if_TD <- factor(deviationdf$if_TD, levels = c("TD", "ADHD"))
summary(deviationdf[,c("if_TD", "age", "sex", "meanFD")])
## Match data
SCdata.match.df <- list(); n =0
for (seed in seq(from=500,to=599,by=1)){
  set.seed(seed)
  SCdata.match.m <- matchit(if_TD~ age+sex+meanFD, data = deviationdf, 
                            method = "nearest", ratio=1, m.order="random", replace = F)
  summary(SCdata.match.m)
  tiff(paste0(FigureFolder, "/PSM_scoredistribution_TDtest1_ADHD1_rematch_Normmod.tiff"), width=14, height=10, res=300, unit="cm")
  plot(SCdata.match.m, type = "hist")
  dev.off()
  
  SCdata.match.final <- match.data(SCdata.match.m, drop.unmatched = T)
  n=n+1
  SCdata.match.df[[n]] <- SCdata.match.final
  #saveRDS(SCdata.match.final, paste0(interfileFolder, "/SCdata_Yeo", Yeoresolution, 
  #                                   "_CV75_sumSCinvnode.deviations_TDtest_ADHD_match_seed", seed,".rds"))
}

SCdata.match.TD <- do.call(rbind, lapply(SCdata.match.df, function(x) x[x$if_TD=="TD", ]))
SCdata.match.TD <- as.data.frame(SCdata.match.TD %>% distinct(scanID, .keep_all = T))


## Merge
SCdata.match.final <- rbind(SCdata.match.TD, SCdata.match.ADHD)

saveRDS(SCdata.match.final, paste0(interfileFolder, "/SCdata_Yeo", Yeoresolution, 
                                   "_CV75_sumSCinvnode.deviations_TDtest_ADHD_match_merge.rds"))

SCdata.match.final <- readRDS(paste0(interfileFolder, "/SCdata_Yeo", Yeoresolution, 
                                     "_CV75_sumSCinvnode.deviations_TDtest_ADHD_match_merge.rds"))
print(paste("Within TD group,", sum(table(SCdata.match.final$subID[SCdata.match.final$if_TD=="TD"])>1) / length(unique(SCdata.match.final$subID[SCdata.match.final$if_TD=="TD"])), "subjects have longitudinal measurements in matched sample."))
# Within TD group, 0.209710743801653 subjects have longitudinal measurements in matched sample.
print(paste("Within ADHD group,", sum(table(SCdata.match.final$subID[SCdata.match.final$if_TD=="ADHD"])>1) / length(unique(SCdata.match.final$subID[SCdata.match.final$if_TD=="ADHD"])), "subjects have longitudinal measurements in matched sample."))
# Within ADHD group, 0.175965665236052 subjects have longitudinal measurements in matched sample.

# description
demovar <- c("sex", "age", "meanFD", "distance", "ehi_y_ss_scoreb", "eventname2", 
             "race_ethnicity", "cbcl_scr_syn_attention_t", "cbcl_scr_syn_external_t", "cbcl_scr_dsm5_adhd_t")

tableone.df <- CreateTableOne(demovar, strata="if_TD", data=SCdata.match.final, 
                              factorVars = c("sex", "ehi_y_ss_scoreb", "eventname2", "race_ethnicity"), includeNA = T,
                              test=T)
tableone.df <- print(tableone.df)
write.csv(tableone.df, paste0(resultFolder, "/demoinfo_ADHD_TDtest_matched_Yeo17.csv"), row.names = T)

#### Plot deviation
deviationSCmat.list <- list()
# ADHD
deviationdf.ADHD <- deviationdf %>% filter(if_TD=="ADHD")

deviationSCmat <- deviationdf.ADHD %>% dplyr::select(ends_with("_deviationZ") & starts_with("SC."))
deviationSCmat.mean <- colMeans(deviationSCmat)
deviationSCmat.df <- data.frame(SClabel = names(deviationSCmat), deviationMean = deviationSCmat.mean)
deviationSCmat.list[["ADHD"]] <- deviationSCmat.df

lmthr <- max(abs(deviationSCmat.df$deviationMean))+1e-5

# TD
deviationdf.TD <- SCdata.match.final %>% filter(if_TD=="TD")

deviationSCmat <- deviationdf.TD %>% dplyr::select(ends_with("_deviationZ") & starts_with("SC."))
deviationSCmat.mean <- colMeans(deviationSCmat)
deviationSCmat.df <- data.frame(SClabel = names(deviationSCmat), deviationMean = deviationSCmat.mean)
deviationSCmat.list[["TD"]] <- deviationSCmat.df

axeslabels=c("VS.P", "SM.A", "VS.C", "SM.B", "DA.A", "DA.B", "FP.C", "DM.C", "VA.A", "TP", "FP.A", "VA.B", "DM.A", "FP.B", "DM.B") 
axeslabelsGap=T
PaletteSet <- list(Name="RdBu", drirection=-1, lmmin = -lmthr, lmmax = lmthr, anglex=45, angley=45,hjustx = 1, hjusty = 1, vjustx = 1, vjusty=0.3)

deviationTD <- plotmatrix(dataname="deviationSCmat.df", variable="deviationMean", ds.resolution=Yeoresolution.delLM, Pvar=NA, NAcol="white", lmthr=NA, axeslabels=axeslabels, axeslabelsGap=F, linerange_frame=NA, PaletteSet=PaletteSet, Pvar.noFDR=NA)
deviationTD

ggsave(paste0(FigureFolder, "/Normative_Trajectory/MeanDeviation_TD_MatMatched.tiff"), deviationTD, height = 14, width = 16, units = "cm")

## Correlation to S-A connectional axis
SCrankcorr(deviationSCmat.df, "deviationMean", Yeoresolution.delLM)
#   ds.resolution  Interest.var r.spearman   p.spearman
# 1            15 deviationMean  -0.05807469  0.5286621

# scatter plot
df <- SCrankcorr(deviationSCmat.df, "deviationMean", Yeoresolution.delLM, T)

scatter <- ggplot(data=df)+
  geom_point(aes(x=SCrank, y=deviationMean, color=deviationMean), size=5)+
  geom_smooth(aes(x=SCrank, y=deviationMean), method ="lm", color="black", linewidth=1.2)+
  scale_color_distiller(type="seq", palette = "RdBu", limits=c(-lmthr, lmthr), direction = -1)+
  labs(x="S-A connectional axis rank", y="Average SC deviation")+
  #scale_y_continuous(breaks = c(-0.015, -0.010, -0.005,0, 0.005, 0.010), labels = c(-15,-10, -5, 0, 5, 10))+
  #scale_y_continuous(breaks = c(-0.003,0, 0.003, 0.006), labels = c(-3, 0, 3, 6))+
  scale_x_continuous(breaks = c(0,20,40,60,80,100,120))+
  theme_classic()+
  theme(axis.text=element_text(size=20, color="black"), 
        axis.title =element_text(size=20),aspect.ratio = 0.83,
        axis.line = element_line(linewidth=0.6),axis.ticks= element_line(linewidth=0.6),
        plot.title = element_text(size=20, hjust = 0.5, vjust=2),
        plot.background=element_rect(fill="transparent"),
        panel.background=element_rect(fill="transparent"),
        legend.position = "none")
ggsave(paste0(FigureFolder, '/Normative_Trajectory/AvgDeivationTD_SCrankcorr_matched.tiff'),scatter, width=17, height =15.5, units = "cm")

### Interaction effect
deviationSCmat.list[[1]]$if_TD <- "ADHD"
deviationSCmat.list[[2]]$if_TD <- "TD"

deviationSCmat.all <- do.call(rbind, deviationSCmat.list)
SCrankdf <- SCrankcorr(deviationSCmat.df, "deviationMean", Yeoresolution.delLM, T)
SCrankdf$SClabel <- paste0("SC.", 1:element_num, "_deviationZ")
SCrankdf$deviationMean <- NULL
deviationSCmat.all <- deviationSCmat.all %>% left_join(SCrankdf, by="SClabel")
deviationSCmat.all$if_TD <- factor(deviationSCmat.all$if_TD, levels=c("TD", "ADHD"))
intmodel <- lm(deviationMean ~ SCrank*if_TD, data = deviationSCmat.all)
nullmodel <- lm(deviationMean ~ SCrank+if_TD, data = deviationSCmat.all)
summary(intmodel)
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)      -2.790e-03  7.656e-03  -0.364   0.7159    
# SCrank           -2.656e-05  1.098e-04  -0.242   0.8091    
# if_TDADHD        -2.006e-02  1.083e-02  -1.852   0.0652 .  
# SCrank:if_TDADHD  6.495e-04  1.553e-04   4.182 4.07e-05 ***

anova(nullmodel, intmodel) # F = 17.49, P=4.074e-05


scatter <- ggplot(data=deviationSCmat.all)+
  geom_point(aes(x=SCrank, y=deviationMean, color=if_TD), size=1.8)+
  geom_smooth(aes(x=SCrank, y=deviationMean, fill=if_TD), method ="lm", color="black", alpha=0.3, linewidth=1.2)+
  scale_color_manual(values = c("gray", "#B2182B"))+
  scale_fill_manual(values = c("gray", "#B2182B"))+
  labs(x="S-A connectional axis rank", y="Average SC deviation", color=NULL, fill=NULL)+
  #scale_y_continuous(breaks = c(-0.015, -0.010, -0.005,0, 0.005, 0.010), labels = c(-15,-10, -5, 0, 5, 10))+
  #scale_y_continuous(breaks = c(-0.003,0, 0.003, 0.006), labels = c(-3, 0, 3, 6))+
  scale_x_continuous(breaks = c(0,20,40,60,80,100,120))+
  theme_classic()+
  theme(axis.text=element_text(size=20, color="black"), 
        axis.title =element_text(size=20),aspect.ratio = 0.9,
        axis.line = element_line(linewidth=0.6),axis.ticks= element_line(linewidth=0.6),
        plot.title = element_text(size=20, hjust = 0.5, vjust=2),
        plot.background=element_rect(fill="transparent"),
        panel.background=element_rect(fill="transparent"),
        legend.text = element_text(size=20, color="black"))
ggsave(paste0(FigureFolder, '/Normative_Trajectory/AvgDeivation_SCrankcorr_Interaction_matched.tiff'),scatter, width=20.5, height =18, units = "cm")



