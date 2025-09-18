## EFNY data
## This script is to generate a dataframe, in which each column is the strength for an edge.
## For schaefer 400 atlas, 70786 edges left after deleting edges connecting to limbic regions.
library(R.matlab)
library(tidyverse)
library(parallel)
library(openxlsx)
library(rjson)
library(corrplot)
library(RColorBrewer)
library(ggplot2)
rm(list = ls())
# Set atlas
Yeoresolution <- 17
# Set path and load data
wd <- getwd()
homepath <- str_split_i(wd, "Normative_model", 1)
SC_path_EFNY <-'/ibmgpfs/cuizaixu_lab/congjing/brainproject/development/results/defaultatlas'
SC_path_PKU6 <-'/ibmgpfs/cuizaixu_lab/xuxiaoyu/PKU6/SCmat'
SC_path_CCNP <-'/ibmgpfs/cuizaixu_lab/xuxiaoyu/CCNP/processed/SC'
Volume_path_EFNY <-'/ibmgpfs/cuizaixu_lab/congjing/brainproject/development/results/schaefer400_nodevolume'
Volume_path_PKU6 <-'/ibmgpfs/cuizaixu_lab/xuxiaoyu/PKU6/schaefer400_nodevolume'
Volume_path_CCNP <- '/ibmgpfs/cuizaixu_lab/xuxiaoyu/CCNP/processed/schaefer400_nodevolume'

demopath <- paste0(homepath, '/Normative_model/demography')
interfileFolder <- paste0(homepath, '/Normative_model/interfileFolder_EFNYnoCCNP')
functionFolder <- paste0(homepath, "/Normative_model/functions")
resultFolder <- paste0(homepath, "/Normative_model/results_EFNYnoCCNP")
FigureFolder <- paste0(homepath, "/Normative_model/Figures_EFNYnoCCNP")

Behavior <- read.csv(paste0(demopath, '/ALL_SITES_basic_demo.csv'))
Behavior <- Behavior %>% filter(age >= 6.5, age <= 15.5) # N=785
#### import schaefer400 index
schaefer400_index_SA<-read.csv(paste0(interfileFolder, '/schaefer400_index_SA.csv'))
## qsiprep output matrix is in Yeo 7 order, so reorder schaefer400 index to Yeo 7 order
schaefer400_index<-schaefer400_index_SA[order(schaefer400_index_SA$index),]
limbicindex <- which(str_detect(schaefer400_index$label_17network, "Limbic"))
schaefer400_index <- schaefer400_index[-limbicindex, ]
schaefer376_delLM <- schaefer400_index$index
## Rearrange left and right regions.
schaefer400_index$index_7network_LRmixed <- schaefer400_index$index
schaefer400_index$index_7network_LRmixed[str_detect(schaefer400_index$label, "RH_")] <- schaefer400_index$index_7network_LRmixed[str_detect(schaefer400_index$label, "RH_")] - 200
orderYeo_7<-order(schaefer400_index$index_7network_LRmixed)

schaefer400_index$index_17network_LRmixed <- schaefer400_index$index_17network
schaefer400_index$index_17network_LRmixed[str_detect(schaefer400_index$label_17network, "RH_")] <- schaefer400_index$index_17network_LRmixed[str_detect(schaefer400_index$label_17network, "RH_")] - 200
orderYeo_17<-order(schaefer400_index$index_17network_LRmixed)

#### import SC data
#### 376 regions, 377*376/2=70876 SCs
#################################################
colname <- character(length = 70876)
for (i in 1:70876){
  colname[i] <- paste0('SC.', as.character(i))
}

SCdata.sum <- mclapply(1:nrow(Behavior), function(i){
  ID <- Behavior$ID[i]
  site <- Behavior$site[i]
  
  if (site == "EFNY"){
    SCname <- paste0(ID, '_dir-PA_space-T1w_desc-preproc_msmtconnectome.mat')
    SC_file_path <- paste0(SC_path_EFNY, '/', SCname)
    volumefile <- paste0(Volume_path_EFNY, '/', ID, '_Volume7.txt')
  }else if (site == "PKU6"){
    SCname <- paste0(ID, '.mat')
    SC_file_path <- paste0(SC_path_PKU6, '/', SCname)
    volumefile <- paste0(Volume_path_PKU6, '/', ID, '_Volume7.txt')
  }
  
  if (file.exists(SC_file_path)){
    SCmat <- readMat(SC_file_path)
    SCmat <- SCmat$schaefer400.sift.invnodevol.radius2.count.connectivity[schaefer376_delLM, schaefer376_delLM]
    
    #Reorder the nodes
    if (Yeoresolution == 7){
      SCmat <- SCmat[orderYeo_7, orderYeo_7]
    }else if (Yeoresolution == 17){
      SCmat <- SCmat[orderYeo_17, orderYeo_17]
    }else{
      print("Invalid Yeoresolution!")
    }
    
    indexup <- upper.tri(SCmat)
    indexsave <- !indexup ###keep lower triangle and diagonal
    SCdat <- as.data.frame(c(SCmat[indexsave]))
    SCdat <- as.data.frame(t(SCdat), row.names = NULL)
    names(SCdat) <- colname
    row.names(SCdat) <- NULL
    SCdat$ID[1] <- ID
  }
  return(SCdat)
}, mc.cores = 50)
ncoldf <- lapply(SCdata.sum, function(x) ncol(x))
SCdata.df <- do.call(rbind, SCdata.sum)
saveRDS(SCdata.df, paste0(interfileFolder, '/SCdataYeo', Yeoresolution,'.sum.msmtcsd.delLM.rds'))
SCdata.sum.merge <- merge(SCdata.df, Behavior, by="ID")
## calculate CV
meanSC<-colMeans(SCdata.sum.merge[,2:70877])
sd.SC <- mclapply(1:70876, function(x) {
  sd.tmp<-sd(SCdata.sum.merge[,x+1])
  return(sd.tmp)
}, mc.cores = 40)
sd.SC<-as.numeric(sd.SC)
CV.SC<-sd.SC/meanSC
Perct.CV.SC <- quantile(CV.SC, probs=seq(0, 1, 0.25)) # extract 
# 0%         25%        50%        75%       100% 
# 0.3510151  0.9908105  1.2180421  1.5238890 11.0477614
## ID of edges over threshold
deleteindex.delLM <- which(CV.SC>Perct.CV.SC[4])
SCdata.sum.merge[,deleteindex.delLM+1] <- 0
meanSC[deleteindex.delLM] <-0
saveRDS(SCdata.sum.merge, paste0(interfileFolder, '/SCdata.sum.CV75.merge.Yeo', Yeoresolution,'.delLM.rds'))
saveRDS(deleteindex.delLM, paste0(interfileFolder, '/CV75_deleteindex.Yeo', Yeoresolution,'.delLM.rds'))

# Validation Threshold = 25th CV
deleteindex.delLM.25 <- which(CV.SC>Perct.CV.SC[2])
SCdata.sum.merge.CV25 <- SCdata.sum.merge
SCdata.sum.merge.CV25[,deleteindex.delLM.25+1] <-0
meanSC.CV25 <- meanSC; meanSC.CV25[deleteindex.delLM.25] <-0
saveRDS(SCdata.sum.merge.CV25, paste0(interfileFolder, '/SCdata.sum.CV25.merge.Yeo', Yeoresolution,'.delLM.rds'))
saveRDS(deleteindex.delLM.25, paste0(interfileFolder, '/CV25_deleteindex.Yeo', Yeoresolution,'.delLM.rds'))
########################################################################

## plot
SCdata.sum.merge <- readRDS(paste0(interfileFolder, '/SCdata.sum.CV75.merge.Yeo', Yeoresolution,'.delLM.rds'))
SCdata.sum.merge.CV25 <- readRDS(paste0(interfileFolder, '/SCdata.sum.CV25.merge.Yeo', Yeoresolution,'.delLM.rds'))
meanSC<-colMeans(SCdata.sum.merge[,2:70877])
meanSC25<-colMeans(SCdata.sum.merge.CV25[,2:70877])
#376
Matsize<-376
Matrix.376 <- matrix(NA, nrow=Matsize, ncol =Matsize)
indexup <- upper.tri(Matrix.376)
indexsave <- !indexup ###keep lower triangle and diagonal
index <- as.numeric(meanSC)
Matrix.376[indexsave] <- index
Matrix.376[indexup] <- t(Matrix.376)[indexup]
colnames(Matrix.376) <-seq(1, Matsize)
rownames(Matrix.376) <-seq(1, Matsize)
tiff( 
  filename = paste0(FigureFolder, '/SCmatrix/SClog376_CV75_Yeo', Yeoresolution,'.tiff'),
  width = 600, 
  height = 600,
  units = "px",
  bg = "white",
  res = 100)
image(log(Matrix.376), col=rev(COL2(diverging = "RdBu", n=200)), axes = FALSE)
dev.off()

# CV25
index <- as.numeric(meanSC25)
Matrix.376[indexsave] <- index
Matrix.376[indexup] <- t(Matrix.376)[indexup]
colnames(Matrix.376) <-seq(1, Matsize)
rownames(Matrix.376) <-seq(1, Matsize)
tiff( 
  filename = paste0(FigureFolder, '/SCmatrix/SClog376_CV25_Yeo', Yeoresolution,'.tiff'),
  width = 600, 
  height = 600,
  units = "px",
  bg = "white",
  res = 100)
image(log(Matrix.376), col=rev(COL2(diverging = "RdBu", n=200)), axes = FALSE)
dev.off()

## Plot the age distribution for ADHD and TD
################
Behavior$diagnosis <- factor(Behavior$ADHD, levels=c(0, 1), labels = c("TD", "ADHD"))
fillcolor = brewer.pal(3, "Paired")

ggplot(data = Behavior, aes(age, y = ..count.., fill = site)) +
  geom_histogram(binwidth = 1, color = "black", position = "stack",linewidth=0.5) +
  facet_wrap(~diagnosis, ncol = 1, nrow = 2, scales = "free") +
  labs(x = "Age (years)", y = NULL) +
  #scale_x_continuous(limits = c(9, 16), breaks = c(9,10,11,12,13,14,15,16)) +
  #scale_y_continuous(limits = c(0, 50), breaks = c(0,10, 20,30,40, 50)) +
  scale_fill_manual(values = rep(fillcolor,2))+
  geom_hline(aes(yintercept = 10), colour = "red", linetype="dashed")+
  theme_classic()+
  theme(panel.background = element_rect(fill="transparent"),
        plot.background = element_rect(fill = "transparent",colour = NA),aspect.ratio = 0.4,
        plot.title = element_text(color = "black", size = 14, hjust = 0.5),
        strip.text = element_text(size = 14, color = "black"),
        strip.background = element_rect(fill = "transparent",colour = NA),
        axis.title = element_text(color = "black", size = 14),axis.line = element_line(linewidth = 0.4),
        axis.ticks = element_line(linewidth = 0.4),
        axis.text = element_text(color = "black", size = 14), 
        legend.text = element_text(color = "black", size = 12),
        legend.title = element_text(color = "black", size = 12))

ggsave(paste(FigureFolder, '/Fig1_Age_distribution_count.tiff', sep = ''),  width = 20, height = 20, units = "cm")
ggsave(paste(FigureFolder, '/Fig1_Age_distribution_count.svg', sep = ''), width = 20, height = 20, units = "cm")

