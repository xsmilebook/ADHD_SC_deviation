rm(list=ls())
library(ggplot2)
library(tidyverse)
library(gamlss)

# 设置路径
homepath <- "D:/code/ADHD_SC_deviation"
interfileFolder <- file.path(homepath, "data", 'interfileFolder', "ABCD")
functionFolder <- file.path(homepath, "src", "functions")

# 加载函数
source(paste0(functionFolder, "/Construct_gamlss_set.R"))

# 设置参数
Yeoresolution <- 17
Yeoresolution.delLM <- 15
element_num <- Yeoresolution.delLM*(Yeoresolution.delLM+1)/2

# 加载数据
cat("正在加载数据...\n")
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

# 测试单个construct_gamlss调用
cat("开始测试construct_gamlss函数...\n")
dependentvar <- "SC.1"

cat("数据维度:", dim(SCdata.TD.trainset), "\n")
cat("依赖变量:", dependentvar, "\n")
cat("平滑项:", smoothterm, "\n")
cat("协变量:", covariates, "\n")
cat("随机变量:", randomvar, "\n")

# 检查数据中是否包含所需列
required_cols <- c(dependentvar, smoothterm, "sex", "meanFD", randomvar, IDvar)
missing_cols <- required_cols[!required_cols %in% colnames(SCdata.TD.trainset)]
if(length(missing_cols) > 0) {
  cat("缺失的列:", paste(missing_cols, collapse=", "), "\n")
} else {
  cat("所有必需的列都存在\n")
}

# 检查数据的前几行
cat("数据前5行:\n")
print(head(SCdata.TD.trainset[, required_cols], 5))

# 尝试调用construct_gamlss函数
cat("正在调用construct_gamlss函数...\n")
tryCatch({
  sumlist <- construct_gamlss(SCdata.TD.trainset, dependentvar, smoothterm, covariates, randomvar, 
                              mu.df, sigma.df, degree, distribution.fam, IDvar, quantile.vec, stratify)
  cat("construct_gamlss函数调用成功!\n")
  cat("返回的对象结构:\n")
  print(str(sumlist))
}, error = function(e) {
  cat("错误信息:", e$message, "\n")
  cat("错误调用栈:\n")
  print(traceback())
})

cat("测试完成\n")