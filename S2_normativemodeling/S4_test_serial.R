rm(list=ls())
library(ggplot2)
library(tidyverse)
library(gamlss)
library(tableone)

# set resolution
Yeoresolution <- 17
Yeoresolution.delLM <- 15
element_num <- Yeoresolution.delLM*(Yeoresolution.delLM+1)/2

# input directory
homepath <- "D:/code/ADHD_SC_deviation"
interfileFolder <- file.path(homepath, "data", 'interfileFolder', "ABCD")
functionFolder <- file.path(homepath, "src", "functions")
resultFolder <- file.path(homepath, "reports", "results", "ABCD")

# source function
source(paste0(functionFolder, "/Construct_gamlss_set.R"))

# Load training data
cat("正在加载训练数据...\n")
SCdata.TD.trainset <- readRDS(paste0(interfileFolder, "/SCdata.TD.trainset_SCYeo", element_num, ".rds"))

# 数据预处理
cat("正在预处理数据...\n")
SCdata.TD.trainset <- as.data.frame(SCdata.TD.trainset)
SCdata.TD.trainset[,c("sex", "siteID")] <- lapply(SCdata.TD.trainset[,c("sex", "siteID")], as.factor)

# 设置模型参数
smoothterm <- "age"
covariates <- "sex+meanFD"
randomvar <- "siteID"
mu.df <- sigma.df <- degree <- 2
distribution.fam <- "GG"
IDvar <- "scanID"
quantile.vec <- c(0.025, 0.5, 0.975)
stratify <- c("sex", "siteID")

# 获取SC变量列表
SC.colnames <- colnames(SCdata.TD.trainset)[grepl("^SC\\.", colnames(SCdata.TD.trainset))]
cat("找到", length(SC.colnames), "个SC变量\n")

# 只测试第一个SC变量
dependentvar <- SC.colnames[1]
cat("测试变量:", dependentvar, "\n")

# 检查数据
cat("数据维度:", dim(SCdata.TD.trainset), "\n")
required_cols <- c(dependentvar, smoothterm, "sex", "meanFD", randomvar, IDvar)
missing_cols <- required_cols[!required_cols %in% colnames(SCdata.TD.trainset)]
if(length(missing_cols) > 0) {
  cat("缺失的列:", paste(missing_cols, collapse=", "), "\n")
  stop("数据中缺少必需的列")
} else {
  cat("所有必需的列都存在\n")
}

# 检查数据的前几行
cat("数据前5行:\n")
print(head(SCdata.TD.trainset[, required_cols], 5))

# 串行调用construct_gamlss函数
cat("开始串行调用construct_gamlss函数...\n")
cat("调用参数:\n")
cat("- 数据框维度:", dim(SCdata.TD.trainset), "\n")
cat("- 依赖变量:", dependentvar, "\n")
cat("- 平滑项:", smoothterm, "\n")
cat("- 协变量:", covariates, "\n")
cat("- 随机变量:", randomvar, "\n")
cat("- 分布族:", distribution.fam, "\n")

tryCatch({
  # 直接调用construct_gamlss函数
  sumlist <- construct_gamlss(SCdata.TD.trainset, dependentvar, smoothterm, covariates, randomvar, 
                              mu.df, sigma.df, degree, distribution.fam, IDvar, quantile.vec, stratify)
  
  cat("construct_gamlss函数调用成功!\n")
  cat("返回的对象类型:", class(sumlist), "\n")
  if(is.list(sumlist)) {
    cat("返回列表的元素名称:", names(sumlist), "\n")
  }
  
}, error = function(e) {
  cat("错误信息:", e$message, "\n")
  cat("错误调用栈:\n")
  traceback()
  
  # 尝试调试：检查全局环境中的对象
  cat("\n全局环境中的对象:\n")
  print(ls(envir = .GlobalEnv))
  
  # 检查是否存在gam.data2
  if(exists("gam.data2", envir = .GlobalEnv)) {
    cat("gam.data2存在于全局环境中\n")
    cat("gam.data2的维度:", dim(get("gam.data2", envir = .GlobalEnv)), "\n")
  } else {
    cat("gam.data2不存在于全局环境中\n")
  }
})

cat("测试完成\n")