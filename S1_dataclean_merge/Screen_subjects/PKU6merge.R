rm(list=ls())
library(openxlsx)

PKU6demofolder <- "D:/xuxiaoyu/open_dataset_information/PKU6/demography"

demolist <- read.csv(paste0(PKU6demofolder, "/basic_demo_PKU6.csv"))
nameIDlist <- read.xlsx(paste0(PKU6demofolder, "/20250901_Name_numberList.xlsx"))

names(nameIDlist) <- c("ID", "Name")

demoADHD <- read.xlsx(paste0(PKU6demofolder, "/20250421_ADHD_NC_demo_adhdrs_cantab_forCIBR.xlsx"), sheet = 1)


# list 1
CBCLlist <- read.xlsx(paste0(PKU6demofolder, "/儿童心理行为评定手册--寇彬彬2.xlsx"), sheet=4)
names(CBCLlist) <- gsub("、", "", names(CBCLlist))
CBCLlist <- CBCLlist[,c(6,10:129)]

CBCLlist <- CBCLlist %>% left_join(nameIDlist, by="Name")


# list 2
CBCLlist2 <- read.xlsx(paste0(PKU6demofolder, "/副本儿童心理行为评定手册20172.xlsx"), sheet=4)
names(CBCLlist2) <- gsub("、", "", names(CBCLlist2))
CBCLlist2 <- CBCLlist2[,c(4,7:126)]
CBCLlist2 <- CBCLlist2 %>% left_join(nameIDlist, by="Name")


CBCLlist.all <- rbind(CBCLlist, CBCLlist2) %>% drop_na(ID)

excludeID <- CBCLlist.all$ID[! CBCLlist.all$ID %in% demolist$ID] # n=22

## missing Age or Sex
demoADHD.miss <- demoADHD %>% filter(is.na(`age(m)`) | is.na(sex))
excludeID.miss <- excludeID[excludeID %in% demoADHD.miss$id] # n=9

## Older age
demoADHD.old <- demoADHD %>% filter(`age(m)` > 15.5*12)
excludeID.old <- excludeID[excludeID %in% demoADHD.old$id]

## Missing MRI
excludeID.noMRI <- excludeID[! excludeID %in% MRsublist.df$ID]
excludeID.noMRI <- excludeID.noMRI[! excludeID.noMRI %in% excludeID.miss] # n=11


## Excessive head motion
length(excludeID[excludeID %in% MRsublist.df$ID]) # n = 3



demolist <- demolist %>% left_join(CBCLlist.all, by="ID")

summary(demolist$`10`)


demolist.nameID <- demolist %>% filter(is.na(`10`), ADHD==1) %>% select(ID)
demolist.nameID <- demolist.nameID %>% left_join(nameIDlist, by="ID")

write.csv(demolist, paste0(PKU6demofolder, "/basic_demo_PKU6_addCBCL.csv"), row.names = F)

write.xlsx(demolist.nameID, paste0(PKU6demofolder, "/ADHDlist_lackCBCL.xlsx"))



