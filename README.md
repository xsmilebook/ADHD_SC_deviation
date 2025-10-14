# ADHD_SC_divation
# ADHD_SC_deviation

本项目研究ADHD（注意力缺陷多动障碍）患者的结构连接（Structural Connectivity, SC）偏差模式。通过规范建模（Normative Modeling）方法，比较ADHD组与典型发育（TD）组之间的差异，并分析这些偏差与年龄和症状的关系。

## 项目结构

项目分为以下几个主要部分：

### S0_NetworkPlot

网络可视化相关脚本：

- **1st_plotConnectionalAxis_Yeo.R**: 绘制Yeo网络的连接轴
- **DU15_SAaxisrank.m**: DU15网络的S-A轴排序MATLAB脚本
- **Yeo17_SAaxisrank.m**: Yeo17网络的S-A轴排序MATLAB脚本
- **Yeo17_merge.m**: 合并Yeo17网络数据的MATLAB脚本

### S1_dataclean_merge

数据清洗和合并相关脚本：

- **S1st_mergedata_ABCD.R**: 合并ABCD数据集
- **S1st_mergedata_EFNY.R**: 合并EFNY数据集
- **S1st_mergedata_Yeo_sumSC_EFNYnoCCNP_LS.R**: 合并EFNY数据集（不含CCNP）的Yeo网络SC总和
- **S2nd_mergedata_Yeo_sumSC_ABCD.R**: 合并ABCD数据集的Yeo网络SC总和
- **S2nd_mergedata_Yeo_sumSC_EFNYnoCCNP.R**: 合并EFNY数据集（不含CCNP）的Yeo网络SC总和
- **S2nd_mergedata_Yeo_sumSC_PostMedication_EFNYnoCCNP.R**: 合并EFNY数据集（不含CCNP）用药后的Yeo网络SC总和
- **S3rd_PSM_ABCD.R**: ABCD数据集的倾向性评分匹配（PSM）

#### Screen_subjects

被试筛选相关脚本：

- **CCNP_exclusion_demo.Rmd**: CCNP数据集排除标准
- **PKU6_exclusion.Rmd/html**: PKU6数据集排除标准
- **PKU6merge.R**: 合并PKU6数据
- **S1_dataclean_ABCD.Rmd/html**: ABCD数据集清洗
- **S1_dataclean_ABCD_CBCL.Rmd/html**: ABCD数据集CBCL量表数据清洗
- **S1_dataclean_ABCD_nopremature.Rmd**: 排除早产儿的ABCD数据集清洗
- **S1_hbn_exclusion_demo.Rmd/html**: HBN数据集排除标准
- **S1_normalize_symptom_PKU6.R**: PKU6数据集症状标准化
- **cuiBP_exclusion_demo.Rmd/html**: cuiBP数据集排除标准

### S2_normativemodeling

规范建模相关脚本：

- **S1_selectparameters_ABCD.qmd**: ABCD数据集参数选择
- **S1_selectparameters_EFNY_allTDTrain_noCCNP.qmd**: EFNY数据集（不含CCNP）参数选择
- **S2_bootstrap_ABCD.R**: ABCD数据集bootstrap方法
- **S2_bootstrap_ABCD_exe.R**: 执行ABCD数据集bootstrap
- **S2_bootstrap_EFNYnoCCNP.R**: EFNY数据集（不含CCNP）bootstrap方法
- **S2_bootstrap_EFNYnoCCNP_exe.R**: 执行EFNY数据集（不含CCNP）bootstrap
- **S3_summary_nommodel_ABCD.Rmd**: ABCD数据集规范模型结果汇总
- **S4_constructNM_forDeviation_ABCD.R**: 构建ABCD数据集偏差计算的规范模型
- **S4_constructNM_forDeviation_EFNYnoCCNP.R**: 构建EFNY数据集（不含CCNP）偏差计算的规范模型
- **S5_plotNM_eachGroup_ABCD.qmd**: 绘制ABCD数据集各组规范模型结果
- **S5_plotNM_eachGroup_EFNYnoCCNP.qmd**: 绘制EFNY数据集（不含CCNP）各组规范模型结果
- **S6_extremedeviation_symptom.R**: 分析极端偏差与症状的关系
- **V5_plotNM_eachGroup_ABCD.R**: 绘制ABCD数据集各组规范模型结果（R脚本版本）
- **submit_bootstrap.sh**: 提交bootstrap任务的Shell脚本

### S3_directcompare

直接比较分析相关脚本：

- **S1_DeviationChangeByAge_ABCD.Rmd**: ABCD数据集偏差随年龄变化分析
- **S1_DeviationChangeByAge_EFNYnoCCNP_ADHD.R**: EFNY数据集（不含CCNP）ADHD组偏差随年龄变化分析
- **S1_DeviationChangeByAge_EFNYnoCCNP_TD.R**: EFNY数据集（不含CCNP）TD组偏差随年龄变化分析
- **S1_DeviationChangeByAge_EFNYnoCCNP_onlyADHD.Rmd**: EFNY数据集（不含CCNP）仅ADHD组偏差随年龄变化分析
- **S1_compareSCDeviation_adolescens_ABCD.R**: ABCD数据集青少年SC偏差比较
- **S1_compareSCDeviation_adolescens_EFNYnoCCNP.R**: EFNY数据集（不含CCNP）青少年SC偏差比较
- **S1_compareSCDeviation_all_EFNYnoCCNP.R**: EFNY数据集（不含CCNP）全部样本SC偏差比较
- **S2_compareSCDeviation_ADHDall_VS_TD_ABCD.qmd**: ABCD数据集ADHD组与TD组SC偏差比较
- **S2_compareSCDeviation_ADHDall_VS_TD_EFNYnoCCNP.qmd**: EFNY数据集（不含CCNP）ADHD组与TD组SC偏差比较
- **S3_DeviationSimilarity2SAaxis.R**: 偏差与S-A轴相似性分析
- **V1_Deviationchange_TD.R**: TD组偏差变化分析
- **V2_compareSCDeviation_ADHD_TDtest_matched.R**: 匹配后的ADHD组与TD组SC偏差比较
- **test_compareSCDeviation_ADHDall_VS_TD_ABCD.qmd**: ABCD数据集ADHD组与TD组SC偏差比较测试版本

### S4_association_ADHDsymp_SCstrength

ADHD症状与SC强度关联分析相关脚本：

- **S1_ADHDsymp_SCstrengthDeviation_ABCD.Rmd/html**: ABCD数据集ADHD症状与SC强度偏差关联分析
- **S1_ADHDsymp_SCstrengthDeviation_EFNYnoCCNP.Rmd**: EFNY数据集（不含CCNP）ADHD症状与SC强度偏差关联分析
- **S2_ADHDsymp_SCDeviation_Med_ABCD.qmd**: ABCD数据集ADHD症状、SC偏差与用药关系分析
- **S3_delta_deviation_delta_symptom_ABCD.qmd**: ABCD数据集偏差变化与症状变化关系分析
- **S4_ADHDsymp_medication_PKU6.qmd**: PKU6数据集ADHD症状与用药关系分析
- **V2_ADHDsymp_SCDeviation_Med_ABCD_crosssectional.qmd**: ABCD数据集ADHD症状、SC偏差与用药的横断面分析

### functions

通用函数脚本：

- **ComBat_sva.R**: ComBat批次效应校正函数
- **Compare_distributions_gamlss.R**: 使用gamlss比较分布函数
- **Construct_gamlss_set.R**: 构建gamlss模型集函数
- **boot_alignment.R**: bootstrap对齐函数
- **colorbarvalue.R**: 颜色条值函数
- **compute_alignment.R**: 计算对齐函数
- **gam_linearinteraction.R**: GAM线性交互函数
- **gam_smoothinteraction.R**: GAM平滑交互函数
- **gam_varyingfactor.R**: GAM变因子函数
- **gamm_linearinteraction.R**: GAMM线性交互函数
- **gamm_smoothinteraction.R**: GAMM平滑交互函数
- **gamm_tensorinteraction.R**: GAMM张量交互函数
- **gamm_varyingfactor.R**: GAMM变因子函数
- **help.R**: 帮助函数
- **index2network.R**: 索引到网络转换函数
- **linearMixed_idvslope.R**: 线性混合模型ID与斜率函数
- **linearmediation_longitudinal.R**: 纵向线性中介函数
- **lme4symp.R**: lme4症状分析函数
- **lmsymp.R**: 线性模型症状分析函数
- **mediate.R**: 中介分析函数
- **nonlinearmediation_longitudinal.R**: 纵向非线性中介函数
- **parallelmediation.R**: 并行中介分析函数
- **plotmatrix.R**: 矩阵绘图函数

### gamfunction

GAM相关函数脚本：

- **EFA.R**: 探索性因子分析函数
- **PCA.R**: 主成分分析函数
- **SCrankcorr.R**: SC排序相关函数
- **SCrankcorr_mult.R**: SC多重排序相关函数
- **SCrankcorr_nullrho.R**: SC零相关系数函数
- **SCrankcorr_spin.R**: SC自旋相关函数
- **colorbarvalue.R**: 颜色条值函数
- **combat.R**: ComBat批次效应校正函数
- **corrmat.R**: 相关矩阵函数
- **gam_factor_interaction.R**: GAM因子交互函数
- **gam_factor_interaction_disease.R**: GAM疾病因子交互函数
- **gam_varyingcoefficients.R**: GAM变系数函数
- **gamcog.R**: GAM认知函数
- **gamcog_scale.R**: GAM认知尺度函数
- **gamderivatives.R**: GAM导数函数
- **gamderivatives_changed.R**: 修改版GAM导数函数
- **gamestimatesmooth.R**: GAM估计平滑函数
- **gaminteraction.R**: GAM交互函数
- **gamm_factor_interaction.R**: GAMM因子交互函数
- **gamm_varyingcoefficients.R**: GAMM变系数函数
- **gammcog.R**: GAMM认知函数
- **gammestimatesmooth.R**: GAMM估计平滑函数
- **gamminteraction.R**: GAMM交互函数
- **gammsmooth.R**: GAMM平滑函数
- **gamsmooth.R**: GAM平滑函数
- **gamsmooth_nocovariates.R**: 无协变量GAM平滑函数
- **glmcog.R**: GLM认知函数
- **glmlogit.R**: GLM逻辑回归函数
- **linearmediation.R**: 线性中介函数
- **linearmediation_base2year.R**: 基线到年度线性中介函数
- **linearmediation_longitudinal.R**: 纵向线性中介函数
- **nonlinearmediation.R**: 非线性中介函数
- **permspheregam.R**: 球面GAM置换函数
- **permutation_correlation_compare.R**: 置换相关比较函数
- **plotdata_derivatives.R**: 绘制导数数据函数
- **plotdata_generate.R**: 生成绘图数据函数

## 执行顺序

项目分析流程按以下顺序执行：

1. **数据准备阶段** (S1_dataclean_merge)：清洗和合并原始数据
2. **规范建模阶段** (S2_normativemodeling)：构建规范模型并计算偏差
3. **直接比较阶段** (S3_directcompare)：比较不同组别间的偏差差异
4. **关联分析阶段** (S4_association_ADHDsymp_SCstrength)：分析偏差与症状的关系

## 数据集说明

- **ABCD**：青少年脑认知发展研究数据集
- **EFNY**：EFNY数据集（不含CCNP部分）
- **PKU6**：北京大学6号数据集
- **HBN**：健康脑网络数据集
- **CCNP**：中国认知神经影像计划数据集

## 分析方法

本项目主要采用规范建模方法分析ADHD患者的结构连接偏差。通过以下步骤：

1. 使用TD组数据构建规范模型
2. 计算ADHD组相对于规范模型的偏差
3. 比较不同组别间的偏差差异
4. 分析偏差与年龄、症状的关系

## 依赖库

项目依赖的主要R包包括：
- tidyverse
- mgcv
- gamlss
- ggplot2
- openxlsx
- parallel
- psych
- reshape
- lme4
- scales