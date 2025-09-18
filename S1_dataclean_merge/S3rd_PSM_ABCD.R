## Propensity Score Matching
rm(list = ls())
library(MatchIt)
library(tidyverse)
library(tableone)

# input Yeo resolution of 7 or 17.
Yeoresolution <- 17
if (Yeoresolution == 7){
  Yeoresolution.delLM = 6
}else if (Yeoresolution == 17){
  Yeoresolution.delLM = 15
}
edgenum <- Yeoresolution.delLM*(Yeoresolution.delLM+1) /2

wd <- getwd()
homepath <- str_split_i(wd, "Normative_model", 1)
demopath <- paste0(homepath, '/Normative_model/demography')
interfileFolder <- paste0(homepath, '/Normative_model/interfileFolder_ABCD')
functionFolder <- paste0(homepath, "/Normative_model/functions")
resultFolder <- paste0(homepath, "/Normative_model/results_ABCD")
FigureFolder <- paste0(homepath, '/Normative_model/Figures_ABCD/Yeo', Yeoresolution,'/CV75')

SCdata <- readRDS(paste0(interfileFolder, "/SCdata_Yeo", Yeoresolution, "_CV75_sumSCinvnode.sum.msmtcsd.combat_TD_ADHDall.rds"))
SCdata$if_TD <- factor(SCdata$if_TD, levels = c(1, 0), labels=c("TD", "ADHD"))

table(SCdata$if_TD, SCdata$eventname)
#      2_year_follow_up_y_arm_1 4_year_follow_up_y_arm_1 baseline_year_1_arm_1
# TD                       2877                      653                  3156
# ADHD                      401                      116                   597

# all
SCdata.match.m <- matchit(if_TD~ age+meanFD, data = SCdata, method = "nearest", ratio=5)
summary(SCdata.match.m)

tiff(paste0(FigureFolder, "/PSM_scoredistribution_TD1_ADHD1_rematch.tiff"), width=14, height=10, res=300, unit="cm")
plot(SCdata.match.m, type = "hist")
dev.off()

SCdata.match.final <- match.data(SCdata.match.m)

# Save out
saveRDS(SCdata.match.final, paste0(interfileFolder, "/SCdata_Yeo", Yeoresolution, "_CV75_sumSCinvnode.sum.msmtcsd.combat_TD_ADHDall_matched.rds"))

print(paste("Within TD group,", sum(table(SCdata.match.final$subID[SCdata.match.final$if_TD=="TD"])>1) / length(unique(SCdata.match.final$subID[SCdata.match.final$if_TD=="TD"])), "subjects have longitudinal measurements in matched sample."))
# Within TD group, 0.259225303690121 subjects have longitudinal measurements in matched sample.

# description
demovar <- c("sex", "age", "meanFD", "distance", "ehi_y_ss_scoreb")

tableone.df <- CreateTableOne(demovar, strata="if_TD", data=SCdata.match.final, 
                              factorVars = c("sex", "ehi_y_ss_scoreb"), includeNA = T,
                              test=T)
tableone.df <- print(tableone.df)
write.csv(tableone.df, paste0(resultFolder, "/demoinfo_ADHD_TD_matched_Yeo17.csv"), row.names = T)





