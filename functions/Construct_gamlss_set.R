# GAMLSS functions
library(tidyverse)
library(gamlss)
## Function 1. construct GAMLSS & return model objects and performance.
## This function is used to fit a GAMLSS containing a smooth term.
construct_gamlss <- function(gam.data, dependentvar, smoothterm, covariates,randomvar=NA, mu.df, sigma.df,degree, distribution.fam,IDvar, quantile.vec, stratify=NULL, debug=FALSE){
  # 允许外部提供 gam.data2 作为兜底（兼容旧脚本）
  if (missing(gam.data) || is.null(gam.data)) {
    if (exists("gam.data2", envir = .GlobalEnv)) {
      gam.data <- get("gam.data2", envir = .GlobalEnv)
    } else {
      stop("Input data not provided and 'gam.data2' not found in global env")
    }
  }
  if (is.character(gam.data)) {
    gam.data <- get(gam.data, envir = .GlobalEnv)
  }

  # 数据预处理
  gam.data <- gam.data %>% drop_na(c(dependentvar, smoothterm))
  covariates <- gsub(" ", "", covariates)
  if (!is.na(randomvar)){
    gam.data1 <- gam.data %>% dplyr::select(all_of(c(dependentvar, smoothterm, unlist(strsplit(covariates, "+", fixed=TRUE)), randomvar, IDvar))) %>% drop_na()
    gam.data1$randomvar <- as.factor(gam.data1[[randomvar]])
    gam.data2 <- gam.data1 %>% group_by(randomvar) %>% filter(n() > 30) %>% ungroup()
  }else{
    gam.data2 <- gam.data %>% dplyr::select(all_of(c(dependentvar, smoothterm, unlist(strsplit(covariates, "+", fixed=TRUE)), IDvar))) %>% drop_na()
  }

  con<-gamlss.control(n.cyc=200)
  gam.data2 <- as.data.frame(gam.data2)
  if(nrow(gam.data2) == 0){
    stop("No rows in gam.data2 after filtering")
  }
  required_cols <- unique(c(dependentvar, smoothterm, unlist(strsplit(covariates, "+", fixed=TRUE)), ifelse(is.na(randomvar), NULL, randomvar), IDvar))
  missing_cols <- setdiff(required_cols, colnames(gam.data2))
  if(length(missing_cols) > 0){
    stop(paste("Missing columns:", paste(missing_cols, collapse=",")))
  }

  # 构建模型
  if (! is.na(randomvar)){
    mod.mu.formula <- as.formula(sprintf("%s ~ bs(%s, df = %s, degree=%s) + %s + random(%s)", dependentvar, smoothterm, mu.df, degree, covariates, randomvar))
    mod.mu.null.formula <- as.formula(sprintf("%s ~ %s + random(%s)", dependentvar, covariates, randomvar))
  }else{
    mod.mu.formula <- as.formula(sprintf("%s ~ bs(%s, df = %s, degree=%s) + %s", dependentvar, smoothterm, mu.df, degree, covariates))
    mod.mu.null.formula <- as.formula(sprintf("%s ~ %s", dependentvar, covariates))
  }

  # Build sigma formula
  sigma.formula <- as.formula(sprintf("~ bs(%s, df = %s, degree = %s) + %s", smoothterm, sigma.df, degree, covariates))

  # 拟合主模型
  mod.tmp <- gamlss(mod.mu.formula, sigma.formula = sigma.formula, nu.formula = ~1, family = get(distribution.fam), data = gam.data2, control = con)

  # 拟合空模型（可能失败则设为 NA）
  sigma.null.formula <- as.formula(sprintf("~ %s", covariates))
  mod.null.tmp <- tryCatch(
    gamlss(mod.mu.null.formula, sigma.formula = sigma.null.formula, nu.formula = ~1, family = get(distribution.fam), data = gam.data2, control = con),
    error=function(e) NA
  )

  # performance
  performance.tb <- data.frame(SClabel=rep("SC",1), BIC=rep(0,1), converged = rep(0,1), partialRsq=rep(0,1), Rsq=rep(0,1))
  performance.tb$SClabel[1] <- dependentvar
  performance.tb$BIC[1] <- mod.tmp$sbc
  performance.tb$converged[1] <- mod.tmp$converged
  performance.tb$Rsq[1] <- Rsq(mod.tmp)

  if (length(mod.null.tmp) == 1){
    partialRsq <- NA
  }else{
    sse.model <- sum((mod.tmp$y - fitted(mod.tmp, what="mu"))^2)
    sse.nullmodel <- sum((mod.null.tmp$y - fitted(mod.null.tmp, what="mu"))^2)
    partialRsq <- (sse.nullmodel - sse.model)/sse.nullmodel
  }

  # 一阶导数与方向
  PEF <- getPEF(mod.tmp, term=smoothterm, n.points = 1000, parameter = "mu", type="response", plot = FALSE)
  x_values <- seq(min(gam.data2[[smoothterm]]), max(gam.data2[[smoothterm]]), length.out = 1000)
  PEF_test <- PEF(x_values, deriv=1)
  direction <- sum(PEF_test) / abs(sum(PEF_test))
  performance.tb$partialRsq[1] <- partialRsq * direction

  # 分位数曲线
  n_quantiles <- length(quantile.vec)
  n_points <- 1000
  x.tmp <- seq(min(gam.data2[[smoothterm]]), max(gam.data2[[smoothterm]]), length.out=n_points)
  if (length(stratify)==1){
    centiles_strat <- list()
    for (l in 1:nlevels(gam.data2[[stratify]])){
      centile.tmp <- array(NA, dim=c(n_quantiles, n_points))
      for (q in 1:n_quantiles){
        Qua <- getQuantile(mod.tmp, quantile=quantile.vec[q], term = smoothterm, fixed.at = setNames(list(levels(gam.data2[[stratify]])[l]), stratify), n.points = n_points)
        centile.tmp[q, ] <- Qua(x.tmp)
      }
      centiles_strat[[l]] <- centile.tmp
    }
  }else if (length(stratify)==2){
    stratify.1 <- stratify[1]
    stratify.2 <- stratify[2]
    gam.data2[[stratify.2]] <- droplevels(gam.data2[[stratify.2]])
    centiles_strat <- list()
    for (l in 1:nlevels(gam.data2[[stratify.1]])){
      centile.tmp <- array(NA, dim=c(n_quantiles, n_points, length(levels(gam.data2[[stratify.2]]))))
      for (s in 1:nlevels(gam.data2[[stratify.2]])){
        for (q in 1:n_quantiles){
          fixed_at_list <- setNames(list(levels(gam.data2[[stratify.1]])[l], levels(gam.data2[[stratify.2]])[s]), c(stratify.1, stratify.2))
          Qua <- getQuantile(mod.tmp, quantile=quantile.vec[q], term = smoothterm, fixed.at = fixed_at_list, n.points = n_points)
          centile.tmp[q, ,s] <- Qua(x.tmp)
        }
      }
      # average along the third dimension
      centile.tmp <- apply(centile.tmp, c(1,2), mean)
      centiles_strat[[l]] <- centile.tmp
    }
  }else if (length(stratify)==0){
    centiles_strat <- array(NA, dim=c(n_quantiles, n_points))
    for (q in 1:n_quantiles){
      Qua <- getQuantile(mod.tmp, quantile=quantile.vec[q], term = smoothterm, n.points = n_points)
      centiles_strat[q, ] <- Qua(x.tmp)
    }
  }

  sumlist <- list(performance.tb=performance.tb, mod.tmp=mod.tmp, centiles_strat=centiles_strat)
  return(sumlist)
}