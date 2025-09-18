## This script is to generate a dataframe, in which each column is the strength for an edge in large-scale network.
## For schaefer 400 --> Yeo network atlas, elementnum mean edges in large-scale network.
library(R.matlab)
library(ggplot2)
library(tidyverse)
library(parallel)
library(openxlsx)
library(rjson)
library(reshape)
library(psych)
rm(list = ls())
wdpath <- getwd()
homepath <- "/ibmgpfs/cuizaixu_lab/xuhaoshu/ADHD_SC_deviation"
# input Yeo resolution of 7 or 17.
Yeoresolution <- 17
if (Yeoresolution == 7){
  Yeoresolution.delLM = 6
}else if (Yeoresolution == 17){
  Yeoresolution.delLM = 15
}
elementnum <- Yeoresolution.delLM*(Yeoresolution.delLM+1) /2

SC_path <-'/ibmgpfs/cuizaixu_lab/xuxiaoyu/ABCD/processed/qsiPrep/SC_matrix'
Volume_path <-'/ibmgpfs/cuizaixu_lab/xuxiaoyu/ABCD/processed/schaefer400_7_nodevolume'
demopath <- file.path(homepath, "data", 'demography')
interfileFolder <- file.path(homepath, "data", 'interfileFolder', "ABCD")
functionFolder <- file.path(homepath, "functions")
resultFolder <- file.path(homepath, "data", "reports", "results", "ABCD")
FigureFolder <- file.path(homepath, "data", "reports", "figures", "ABCD")

Behavior <- read.csv(paste0(demopath, '/demo_sublist7.csv'))
Behavior <- Behavior %>% distinct(scanID, .keep_all = TRUE)

source(paste0(functionFolder, "/plotmatrix.R"))

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

# filter index of P75th and P25th of CV.
deleteindex75 <- readRDS(paste0(interfileFolder, '/CV75_deleteindex.Yeo', Yeoresolution,'.delLM.rds'))
deleteindex25 <- readRDS(paste0(interfileFolder, '/CV25_deleteindex.Yeo', Yeoresolution,'.delLM.rds'))

# assign each region to Yeo.resolution*Yeo.resolution network.
schaefer400_index.Yeo7 <- schaefer400_index[order(schaefer400_index$index_7network_LRmixed),]
schaefer400_index.Yeo17 <- schaefer400_index[order(schaefer400_index$index_17network_LRmixed),]
schaefer400_index.Yeo7 <- schaefer400_index.Yeo7 %>% mutate(Yeo.resolutionnode = recode_factor(network_label,
                                                                                               "Vis" = 1,
                                                                                               "SomMot" = 2,
                                                                                               "DorsAttn" = 3,
                                                                                               "SalVentAttn" = 4,
                                                                                               "Cont" = 5,
                                                                                               "Default" = 6,
                                                                                               "Limbic" = 0))

summary(schaefer400_index.Yeo7$Yeo.resolutionnode)
schaefer400_index.Yeo17 <- schaefer400_index.Yeo17 %>% 
  mutate(Yeo.resolutionnode = recode_factor(network_label_17network,
                                            "VisCent" = 3,
                                            "VisPeri" = 1,
                                            "SomMotA" =2,
                                            "SomMotB" = 4,
                                            "DorsAttnA" =5,
                                            "DorsAttnB" = 6,
                                            "SalVentAttnA" =9,
                                            "SalVentAttnB" =12,
                                            "ContA" =11,
                                            "ContB" =14,
                                            "ContC" =7,
                                            "DefaultA" = 13,
                                            "DefaultB" = 15,
                                            "DefaultC" = 8,
                                            "TempPar" = 10))
write.csv(schaefer400_index.Yeo17, paste0(resultFolder, "/schaefer376_index_Yeo17.csv"), row.names = F)

summary(schaefer400_index.Yeo17$Yeo.resolutionnode)

if (Yeoresolution == 7){
  Yeo.resolutionnode <- schaefer400_index.Yeo7$Yeo.resolutionnode
  Yeo.resolutionnode <- factor(Yeo.resolutionnode, levels = c(1:Yeoresolution.delLM))
}else if (Yeoresolution == 17){
  Yeo.resolutionnode <- schaefer400_index.Yeo17$Yeo.resolutionnode
  Yeo.resolutionnode <- factor(Yeo.resolutionnode, levels = c(1:Yeoresolution.delLM))
}else{
  print("Invalid Yeoresolution!")
}

# SC 376*376 --> Yeoresolution.delLM*Yeoresolution.delLM
matrixYeo.resolution <- matrix(NA, Yeoresolution.delLM, Yeoresolution.delLM)
matrixYeo.resolution[lower.tri(matrixYeo.resolution, diag = T)] <- c(1:elementnum)
matrixYeo.resolution[upper.tri(matrixYeo.resolution)] <- t(matrixYeo.resolution)[upper.tri(matrixYeo.resolution)]
matrix_SCYeo.resolution <- matrix(NA, 376, 376)
for (x in 1:Yeoresolution.delLM){
  for (y in 1:Yeoresolution.delLM){
    xindex <- which(Yeo.resolutionnode==x)
    yindex <- which(Yeo.resolutionnode==y)
    matrix_SCYeo.resolution[xindex, yindex] <- matrixYeo.resolution[x,y]
  }
}
# an index telling how 376*376 map to 12*12
Yeo.resolution.index <- matrix_SCYeo.resolution[lower.tri(matrix_SCYeo.resolution, diag = T)]
#################################################

#### import SC data
#### Yeoresolution.delLM regions, (Yeoresolution.delLM+1)*Yeoresolution.delLM/2=elementnum SCs
#### extract a dataframe containing elementnum columns, each represents an edge.
#################################################
colname <- character(length = elementnum)
for (i in 1:elementnum){
  colname[i] <- paste0('SC.', as.character(i))
}

if (str_detect(wdpath, "cuizaixu_lab")){
  SCdata.sum <- mclapply(1:nrow(Behavior), function(i){
    scanID <- Behavior$scanID[i]
    siteID <- Behavior$siteID[i]
    eventname <-strsplit(scanID, 'ses-')[[1]][2]
    scanner <- Behavior$scanner[i]
    SCname <- paste0(scanID, '_space-T1w_desc-preproc_msmtconnectome.mat')
    SC_file_path <- paste0(SC_path, '/', eventname, '/', scanner,'/', siteID, '/', SCname)
    volumefile <- paste0(Volume_path, '/', scanID, '.txt')
    if (file.exists(SC_file_path)){
      SCmat <- readMat(SC_file_path)
      # load steamline counts matrix & fiber length matrix
      SCmat_raw <- SCmat$schaefer400.sift.radius2.count.connectivity[schaefer376_delLM, schaefer376_delLM]
      length_raw <- SCmat$schaefer400.radius2.meanlength.connectivity[schaefer376_delLM, schaefer376_delLM]
      if (Yeoresolution == 7){
        SCmat_raw <- SCmat_raw[orderYeo_7, orderYeo_7] # 376*376 nodes sorted by Yeo index
        length_raw <- length_raw[orderYeo_7, orderYeo_7] # 376*376 nodes sorted by Yeo index
      }else if (Yeoresolution == 17){
        SCmat_raw <- SCmat_raw[orderYeo_17, orderYeo_17]
        length_raw <- length_raw[orderYeo_17, orderYeo_17]
      }else{
        print("Invalid Yeoresolution!")
      }
      
      totallength_raw <- length_raw * SCmat_raw
      indexup <- upper.tri(SCmat_raw)
      indexsave <- !indexup
      SCmat_raw <- SCmat_raw[indexsave] # 1*70876 each element represents streamline counts
      SCmat_raw75 <- SCmat_raw25 <- SCmat_raw
      SCmat_raw75[deleteindex75]<-0 # remove top 1/4 inconsistent connetions
      SCmat_raw25[deleteindex25]<-0 # remove top 3/4 inconsistent connetions
      totallength_raw <- totallength_raw[indexsave]
      totallength_raw75 <- totallength_raw25 <- totallength_raw
      totallength_raw75[deleteindex75]<-0
      totallength_raw25[deleteindex25]<-0
      df <- data.frame(
        group = Yeo.resolution.index,
        value75 = SCmat_raw75,
        value25 = SCmat_raw25,
        length75 = totallength_raw75,
        length25 = totallength_raw25
      )
      # compute the sum of streamline counts / length for each fraction, in total of elementnum.
      result <- df %>%
        group_by(group) %>%
        summarise(sum_value75 = sum(value75), sum_value25 = sum(value25), sum_length75=sum(length75), 
                  sum_length25=sum(length25))
      mean_length75 <- (result$sum_length75 / result$sum_value75)[1:elementnum]
      mean_length25 <- (result$sum_length25 / result$sum_value25)[1:elementnum]
      sumSC.raw75 <- result$sum_value75[1:elementnum]
      sumSC.raw25 <- result$sum_value25[1:elementnum]
      ## node volume
      if (file.exists(volumefile)){
        nodevolume <- read_table(volumefile, col_names=F)
        if (nrow(nodevolume)==453){
          nodevolume <- as.numeric(nodevolume$X2[schaefer376_delLM]) # delete limbic regions
          if (Yeoresolution == 7){
            nodevolume <- nodevolume[orderYeo_7] # sorted by Yeo index
          }else if (Yeoresolution == 17){
            nodevolume <- nodevolume[orderYeo_17]
          }else{
            print("Invalid Yeoresolution!")
          }
          
          df2 <- data.frame(
            group = Yeo.resolutionnode,
            value = nodevolume
          )
          result2 <- df2 %>% arrange(group) %>%
            group_by(group) %>% 
            summarise(sum_value = sum(value))
          nodevolume_sum <- result2$sum_value[1:Yeoresolution.delLM] # sum of nodes' volume for each node fraction (Yeo.resolution).
          
          ### Yeo.resolution*Yeo.resolution
          volumeSC <- matrix(NA, Yeoresolution.delLM, Yeoresolution.delLM)
          for (x in 1:Yeoresolution.delLM){
            for (y in 1:Yeoresolution.delLM){
              volumeSC[x,y] <- (nodevolume_sum[x]+nodevolume_sum[y]) / 2
            }
          }
          volumeSC <- volumeSC[lower.tri(volumeSC, diag = T)] # the scale values of node volume for each edge.
          sumSC.invnode75 <- sumSC.raw75 / volumeSC
          sumSC.invnode25 <- sumSC.raw25 / volumeSC
        }else{
          sumSC.invnode75 <- rep(NA, elementnum)
          sumSC.invnode25 <- rep(NA, elementnum)
        }}else{
          sumSC.invnode75 <- rep(NA, elementnum)
          sumSC.invnode25 <- rep(NA, elementnum)
        }
      ###keep lower triangle and diagonal
      SCdat75 <- as.data.frame(sumSC.invnode75)
      SCdat75 <- as.data.frame(t(SCdat75), row.names = NULL)
      names(SCdat75) <- colname
      row.names(SCdat75) <- NULL
      SCdat75$scanID[1] <- scanID
      #SCdata.sum75<-rbind(SCdata.sum75, SCdat75)
      
      SCdat25 <- as.data.frame(sumSC.invnode25)
      SCdat25 <- as.data.frame(t(SCdat25), row.names = NULL)
      names(SCdat25) <- colname
      row.names(SCdat25) <- NULL
      SCdat25$scanID[1] <- scanID
      #SCdata.sum25<-rbind(SCdata.sum25, SCdat25)
      
      # not inverse node volume
      SCdat75_noInvNode <- as.data.frame(sumSC.raw75)
      SCdat75_noInvNode <- as.data.frame(t(SCdat75_noInvNode), row.names = NULL)
      names(SCdat75_noInvNode) <- colname
      row.names(SCdat75_noInvNode) <- NULL
      SCdat75_noInvNode$scanID[1] <- scanID
      
      SCdat25_noInvNode <- as.data.frame(sumSC.raw25)
      SCdat25_noInvNode <- as.data.frame(t(SCdat25_noInvNode), row.names = NULL)
      names(SCdat25_noInvNode) <- colname
      row.names(SCdat25_noInvNode) <- NULL
      SCdat25_noInvNode$scanID[1] <- scanID
    }
    print(i)
    return(data.frame(SCdat75, SCdat25, SCdat75_noInvNode, SCdat25_noInvNode))
  },  mc.cores = 50)
  print("The parallel Mat reading is over!")
  saveRDS(SCdata.sum, paste0(interfileFolder, "/SCdata.sum.list_Yeo", Yeoresolution, ".rds"))
  
  SCdata.sum75 <- do.call(rbind, lapply(SCdata.sum, function(x) as.data.frame(x[,c(1:(elementnum+1))])))
  SCdata.sum25 <- do.call(rbind, lapply(SCdata.sum, function(x) as.data.frame(x[,c((elementnum+2):(elementnum*2+2))])))
  names(SCdata.sum25) <- c(colname, "scanID")
  
  SCdata.sum75_noInvNode <- do.call(rbind, lapply(SCdata.sum, function(x) as.data.frame(x[,c((elementnum*2+3):(elementnum*3+3))])))
  SCdata.sum25_noInvNode <- do.call(rbind, lapply(SCdata.sum, function(x) as.data.frame(x[,c((elementnum*3+4):(elementnum*4+4))])))
  names(SCdata.sum25_noInvNode) <- c(colname, "scanID")
  names(SCdata.sum75_noInvNode) <- c(colname, "scanID")
  
  SCdata.sum75.merge <- merge(SCdata.sum75, Behavior, by="scanID")
  SCdata.sum25.merge <- merge(SCdata.sum25, Behavior, by="scanID")
  SCdata.sum75_noInvNode.merge <- merge(SCdata.sum75_noInvNode, Behavior, by="scanID")
  SCdata.sum25_noInvNode.merge <- merge(SCdata.sum25_noInvNode, Behavior, by="scanID")
  
  # convert variables' types
  SCdata.sum75.merge$subID <- as.factor(SCdata.sum75.merge$subID) ; SCdata.sum75.merge$siteID <- as.factor(SCdata.sum75.merge$siteID)
  SCdata.sum25.merge$subID <- as.factor(SCdata.sum25.merge$subID) ; SCdata.sum25.merge$siteID <- as.factor(SCdata.sum25.merge$siteID)
  SCdata.sum75_noInvNode.merge$subID <- as.factor(SCdata.sum75_noInvNode.merge$subID) ; SCdata.sum75_noInvNode.merge$siteID <- as.factor(SCdata.sum75_noInvNode.merge$siteID)
  SCdata.sum25_noInvNode.merge$subID <- as.factor(SCdata.sum25_noInvNode.merge$subID) ; SCdata.sum25_noInvNode.merge$siteID <- as.factor(SCdata.sum25_noInvNode.merge$siteID)
  
  saveRDS(SCdata.sum75.merge, paste0(interfileFolder, '/SCdata_Yeo',Yeoresolution, '_CV75_sumSCinvnode.sum.msmtcsd.merge.rds'))
  saveRDS(SCdata.sum25.merge, paste0(interfileFolder, '/SCdata_Yeo',Yeoresolution, '_CV25_sumSCinvnode.sum.msmtcsd.merge.rds'))
  saveRDS(SCdata.sum75_noInvNode.merge, paste0(interfileFolder, '/SCdata_Yeo',Yeoresolution, '_CV75_sumSC.sum.msmtcsd.merge.rds'))
  saveRDS(SCdata.sum25_noInvNode.merge, paste0(interfileFolder, '/SCdata_Yeo',Yeoresolution, '_CV25_sumSC.sum.msmtcsd.merge.rds'))
}else{
  SCdata.sum75.merge <- readRDS(paste0(interfileFolder, '/SCdata_Yeo',Yeoresolution, '_CV75_sumSCinvnode.sum.msmtcsd.merge.rds'))
  SCdata.sum25.merge <- readRDS(paste0(interfileFolder, '/SCdata_Yeo',Yeoresolution, '_CV25_sumSCinvnode.sum.msmtcsd.merge.rds'))
}

# plot
## Select a subject randomly
SCdata.tmp <- as.numeric(SCdata.sum75.merge[12, paste0("SC.", 1:elementnum)])
df <- data.frame(SCstrength=SCdata.tmp)

PaletteSet <- list(Name="Blues", direction=1, lmmin = min(SCdata.tmp), lmmax = max(SCdata.tmp), 
                   anglex=45, angley=45,hjustx = 1, hjusty = 1, vjustx = 1, vjusty=0.3)

SCmatPlot <- plotmatrix(dataname="df", variable="SCstrength", 
                          ds.resolution=Yeoresolution.delLM, Pvar=NA, NAcol="white", 
                        lmthr=NA, axeslabels=NULL, axeslabelsGap=F, 
                        linerange_frame=NA, PaletteSet=PaletteSet, Pvar.noFDR=NA)
SCmatPlot

ggsave(paste0(FigureFolder,"/CV75/Matrix_Yeo",Yeoresolution, "_sumSCinvnode_Age8_22/Random_SCmat.tiff"), 
       SCmatPlot, height = 14, width = 16, units = "cm")




