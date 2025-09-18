

index2network <- function(index){
  mattmp <- matrix(NA,15,15)
  mattmp[lower.tri(mattmp, diag=T)] <- 1:120
  axeslabels=c("VS.P", "SM.A", "VS.C", "SM.B", "DA.A", "DA.B", "FP.C", "DM.C", "VA.A", "TP", "FP.A", "VA.B", "DM.A", "FP.B", "DM.B")
  i <- which(mattmp==index, arr.ind = T)
  
  
  return(axeslabels[i])
}

