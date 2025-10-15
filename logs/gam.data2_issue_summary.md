# gam.data2 问题总结

## 1. 问题背景描述
- 错误现象：在串行与并行阶段均出现 `object 'gam.data2' not found`。
- 触发位置：`construct_gamlss` 的导数与分位数计算（`getPEF`/`getQuantile`）阶段使用 `gam.data2` 的列范围与分层信息。
- 并行环境：Windows 下 PSOCK 集群的每个 worker 是独立 R 会话，不继承主进程的全局对象。
- 根因摘要：
  - 计算分位/导数时需要访问名为 `gam.data2` 的数据对象或依赖 `mod.tmp$call$data` 指向的符号；
  - 若未在 worker 环境提供该对象，就会触发 “找不到对象”。

## 2. 解决方案与原理
- 已采用的修复（最小改动，稳定可用）：
  - 在 ABCD 构建脚本中显式赋值并导出：
    ```r
    gam.data2 <- SCdata.TD.trainset
    parallel::clusterExport(cl, c('gam.data2', 'SCdata.TD.trainset', 'smoothterm', 'covariates',
                                  'randomvar', 'mu.df', 'sigma.df', 'degree', 'distribution.fam',
                                  'IDvar', 'quantile.vec', 'stratify', 'construct_gamlss'))
    ```
  - 在 `construct_gamlss` 内保留外部兜底：当未显式传入数据或传入字符串名时，优先从全局环境读取 `gam.data2`，并做必需列/行数校验。
  - 同时移除之前的调试日志与阶段性错误标识，使函数精简且自包含。
- 原理说明：
  - 通过 `clusterExport` 将 `gam.data2` 序列化到各个 worker 的全局环境，使 `getQuantile`/`getPEF` 以及基于 `gam.data2` 的范围、分层计算在子进程中可用；
  - 外部兜底确保老脚本或字符串数据名的调用方式仍然有效。
- 可选的替代方案（按稳健度从高到低）：
  - 完全去除对模型 call 数据的依赖：用“网格 `newdata` + `predict`(mu/sigma/nu) + `qGG/pGG` + 有限差分”替代 `getQuantile/getPEF`。
  - 在函数内部做短暂全局兜底：`assign('gam.data2', gam.data2, envir=.GlobalEnv)` 并用 `on.exit` 自动清理。
  - 用“参数名绑定”：`data_sym <- deparse(substitute(gam.data))` 后从全局用该符号 `get(data_sym)` 作为模型数据。
  - 采用 `future.apply`/`foreach` 并让框架自动传播全局依赖。
- 关键代码片段（示例）：
  ```r
  # 1) 构建脚本中导出
  gam.data2 <- SCdata.TD.trainset
  parallel::clusterExport(cl, c('gam.data2', 'construct_gamlss'))

  # 2) 函数内短暂全局兜底（可选）
  assign('gam.data2', gam.data2, envir = .GlobalEnv)
  on.exit({ if (exists('gam.data2', envir=.GlobalEnv)) remove('gam.data2', envir=.GlobalEnv) }, add = TRUE)

  # 3) 网格预测替代（思路）
  mu <- predict(mod.tmp, newdata=newdata, what='mu', type='response')
  sigma <- predict(mod.tmp, newdata=newdata, what='sigma', type='response')
  nu <- predict(mod.tmp, newdata=newdata, what='nu', type='response')
  curve_p <- qGG(p, mu=mu, sigma=sigma, nu=nu)
  deriv <- (curve_p(x+eps) - curve_p(x)) / eps
  ```
- 适用范围与影响：
  - 当前修复覆盖 ABCD 并行拟合阶段；EFNY 脚本可复用同样模式；
  - 预期结果：不再出现 `object 'gam.data2' not found`，拟合与偏差计算正常；
  - 限制：需保证 `gam.data2` 列完整、行数大于 0；无分层分支中 `Qua` 未定义问题不影响当前脚本。
- 后续建议：
  - 若希望彻底去除全局依赖，建议逐步迁移到“网格预测 + 分布函数”实现；
  - 或将 `clusterExport` 封装为工具函数，统一导出所需对象与函数，减少脚本重复。