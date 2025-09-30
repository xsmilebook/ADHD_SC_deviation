
library(tidyverse)
library(mgcv)
library(openxlsx)
library(parallel)
library(gamlss)
library(scales)

execute_boot <- function(n, var){
  # set resolution
  Yeoresolution <- 17
  if (Yeoresolution == 7){
    Yeoresolution.delLM = 6
  }else if (Yeoresolution == 17){
    Yeoresolution.delLM = 15
  }
  element_num <- Yeoresolution.delLM*(Yeoresolution.delLM+1)/2
  
  # input directory
  homepath <- "/ibmgpfs/cuizaixu_lab/xuhaoshu/ADHD_SC_deviation"
  demopath <- file.path(homepath, "data", 'demography')
  interfileFolder <- file.path(homepath, "data", 'interfileFolder', "ABCD")
  functionFolder <- file.path(homepath, "src", "functions")
  resultFolder <- file.path(homepath, "reports", "results", "ABCD")
  functionFolder.SCDev <- file.path(homepath, "src", "gamfunction")
  FigureFolder <- paste0(homepath, '/reports/figures/ABCD/Yeo', Yeoresolution,'/CV75')
  
  # load data
  SCdataTD <- readRDS(paste0(interfileFolder, "/SCdata.TD.trainset_SCYeo", element_num, ".rds"))
  
  source(paste0(functionFolder, "/Compare_distributions_gamlss.R"))
  
  set.seed(925)
  mu.df <- sigma.df <- degree <- 2
  mod.mu.formula <- as.formula(sprintf("%s ~ bs(%s, df = %s, degree = %s) + %s + random(as.factor(%s))", var, "age", 
                                       mu.df, degree, "sex+meanFD", "siteID"))
  con<-gamlss.control(n.cyc=200, trace = TRUE)
  SCdata.sum75.merge.TD <- SCdataTD %>% select(all_of(c(var, "age", "sex", "meanFD", "siteID", "scanID")))
  model.var <- gamlss(mod.mu.formula, sigma.formula =~ bs(age, df = 2, degree =2)+sex+meanFD, nu.formula = ~1, family=GG, 
                            data=SCdata.sum75.merge.TD, control=con)
  
  
  # get input
  n <- as.numeric(n)
  # run bootstrap
  SCdata.sum75.merge.TD$sex <- factor(SCdata.sum75.merge.TD$sex, levels=c(1,2), labels=c("M", "F"))
  SCdata.sum75.merge.TD$siteID <- as.factor(SCdata.sum75.merge.TD$siteID)
  print(n)
  Base.Seed <- 925
  dataname <- "SCdata.sum75.merge.TD"
  smoothvar <- "age"
  model_obj <- model.var
  stratify <- c("sex", "siteID")
  randomvar <- "siteID"
  
  bootstrap.out <- Boot.Function(n, Base.Seed, dataname, smoothvar,randomvar , model_obj, stratify)
  mod.tmp <- bootstrap.out$mod.tmp
  gam.data.subset <- bootstrap.out$gam.data.subset
  gam.data.subset <- as.data.frame(gam.data.subset)
  assign("gam.data.subset", gam.data.subset, envir = .GlobalEnv)
  # estimate quantiles
  quantiles<-c(0.025, 0.5, 0.975)
  n_quantiles <- length(quantiles)
  n_points <- 1000
  n_sites <- length(unique(gam.data.subset$siteID))
  Centiles_male <- array(NA, dim = c(n_sites, n_quantiles, n_points))
  Centiles_female <- array(NA, dim = c(n_sites, n_quantiles, n_points))
  x_female <- seq(min(SCdata.sum75.merge.TD[[smoothvar]]), max(SCdata.sum75.merge.TD[[smoothvar]]), length.out = n_points)
  x_male <- seq(min(SCdata.sum75.merge.TD[[smoothvar]]), max(SCdata.sum75.merge.TD[[smoothvar]]), length.out = n_points)
  FDvalue <- mean(SCdata.sum75.merge.TD$meanFD)
  
  for (j in 1:n_sites){
    for (i in 1:n_quantiles){
      siteID.tmp <- unique(gam.data.subset$siteID)[j]
      command <- "Qua <- getQuantile(mod.tmp, quantile=quantiles[i], term = 'age', fixed.at = list(sex='F', meanFD=FDvalue, randomvar=siteID.tmp), n.points = n_points)"
      eval(parse(text=command))
      Centiles_female[j, i,] <- Qua(x_female)
      
      command <- "Qua <- getQuantile(mod.tmp, quantile=quantiles[i], term = 'age', fixed.at = list(sex='M', meanFD=FDvalue, randomvar=siteID.tmp), n.points = n_points)"
      eval(parse(text=command))
      Centiles_male[j, i,] <- Qua(x_male)
    }
  }
  
  Centiles <- (Centiles_female + Centiles_male) /2
  bootstrap_centiles <- list(Centiles_female=Centiles_female, Centiles_male=Centiles_male, var=var, boostrap_time=n, converged=mod.tmp$converged)
  
  # estimate first derivatives (finite difference)
  eps <- 1e-7
  x0 <- seq(min(SCdata.sum75.merge.TD[[smoothvar]]), max(SCdata.sum75.merge.TD[[smoothvar]]), length.out = n_points)
  x1 <- x0 + eps
  derivative_male <- array(NA, dim = c(n_sites, n_points))
  derivative_female <- array(NA, dim = c(n_sites, n_points))
  
  for (j in 1:n_sites){
    siteID.tmp <- unique(gam.data.subset$siteID)[j]
    command <- "Qua <- getQuantile(mod.tmp, quantile=0.5, term = 'age', fixed.at = list(sex='F', meanFD=FDvalue, randomvar=siteID.tmp), n.points = n_points)"
    eval(parse(text=command))
    derivative_female[j, ] <- (Qua(x1) - Qua(x0)) / eps
    
    command <- "Qua <- getQuantile(mod.tmp, quantile=0.5, term = 'age', fixed.at = list(sex='M', meanFD=FDvalue, randomvar=siteID.tmp), n.points = n_points)"
    eval(parse(text=command))
    derivative_male[j, ] <- (Qua(x1) - Qua(x0)) / eps
  }
  
  derivative <- (derivative_female + derivative_male) /2
  
  bootstrap_centiles[["derivative.P50"]] <- derivative
  bootstrap_centiles[["derivative.P50.female"]] <- derivative_female
  bootstrap_centiles[["derivative.P50.male"]] <- derivative_male
  
  # save out
  if (! dir.exists(paste0(interfileFolder, "/bootstrap/", var))){
    dir.create(paste0(interfileFolder, "/bootstrap/", var))
  }
  
  saveRDS(bootstrap_centiles, paste0(interfileFolder, "/bootstrap/", var,"/centile_bootstrap_", n, ".rds"))
  
  print(paste("Bootsrap", n, "finished!"))
}

