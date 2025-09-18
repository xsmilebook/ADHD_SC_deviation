library(tidyverse)
library(mgcv)

rm(list = ls())

wd <- getwd()
homepath <- str_split_i(wd, "Normative_model", 1)
demopath <- paste0(homepath, '/Normative_model/demography')
interfileFolder <- paste0(homepath, '/Normative_model/interfileFolder_EFNY')
functionFolder <- paste0(homepath, "/Normative_model/functions")
resultFolder <- paste0(homepath, "/Normative_model/results_EFNY")
Behavior <- read.csv(paste0(demopath, '/basic_demo_PKU6.csv'))
Behavior$sex <- as.factor(Behavior$sex)
Behavior.TD <- Behavior %>% filter(ADHD==0)
Behavior.ADHD <- Behavior %>% filter(ADHD==1)
## Symptom model in TD
summary(Behavior.TD[,c("IA_B", "HI_B", "TO_B", "age", "sex")])
summary(Behavior.ADHD[,c("IA_B", "HI_B", "TO_B", "age", "sex")])

IA_mod <- gam(IA_B ~ s(age, k=3, fx=T) + sex, method = "REML", data=Behavior.TD)
summary(IA_mod)
HI_mod <- gam(HI_B ~ s(age, k=3, fx=T) + sex, method = "REML", data=Behavior.TD)
summary(HI_mod)
TO_mod <- gam(TO_B ~ s(age, k=3, fx=T) + sex, method = "REML", data=Behavior.TD)
summary(TO_mod)

Behavior$IA_B.norm <- Behavior$IA_B - predict(IA_mod, newdata = Behavior)
Behavior$HI_B.norm <- Behavior$HI_B - predict(HI_mod, newdata = Behavior)
Behavior$TO_B.norm <- Behavior$TO_B - predict(TO_mod, newdata = Behavior)

summary(Behavior[,c("IA_B.norm", "HI_B.norm", "TO_B.norm")])

Behavior.TD <- Behavior %>% filter(ADHD==0)
Behavior.ADHD <- Behavior %>% filter(ADHD==1)

corr.test(Behavior.TD$HI_B, Behavior.TD$age)
corr.test(Behavior.TD$HI_B.norm, Behavior.TD$age)

corr.test(Behavior.ADHD$HI_B, Behavior.ADHD$age)
corr.test(Behavior.ADHD$HI_B.norm, Behavior.ADHD$age)

write.csv(Behavior, paste0(demopath, "/basic_demo_PKU6_addSymptom.csv"), row.names = F)

