# 模块03：目标函数设定与优化原理（对应课程5）

## 目标
通过 e_coli_core 模型理解 FBA 的三要素（稳态约束、通量边界、目标函数），
掌握目标函数查看/更换、FBA 解的简并性、pFBA 与 MOMA 的假设差异及敲除预测，
并用 FVA 与通量取样（flux sampling）从"解空间"角度量化预测的唯一性与柔性。

## 依赖
```bash
pip install cobra pandas matplotlib
```
模型加载见文末"模型加载说明"；本模块附带本地模型文件 `e_coli_core.json`。

## 脚本与运行顺序
| 顺序 | 脚本 | 内容 |
|---|---|---|
| 1 | `01_FBA基础与目标函数设定.py` | FBA 原理；查看/更换目标函数（生物量 / ATPM / EX_ac_e）；解读 objective_value、fluxes、shadow_prices、reduced_costs |
| 2 | `02_FBA与pFBA比较.py` | FVA 展示解的简并性（SUCDi/FRD7 无效循环）；构造"强制循环 FBA 解"；与 pFBA 比较并画条形图 |
| 3 | `03_MOMA敲除模拟.py` | 野生型 pFBA 解作参考；ΔgltA(b0720) 与 Δpgi(b4025) 敲除；FBA vs linear MOMA 生长率与通量分布对比图 |
| 4 | `04_flux_variability_analysis.py` | FVA：野生型 vs Δpgi 关键反应通量区间对比；统计 blocked/fixed/variable 反应数；通量区间对比图 |
| 5 | `05_flux_sampling.py` | OptGPSampler 均匀取样 1000 个可行通量分布；关键反应边际分布直方图 + 两两相关性散点/相关系数；存 png 与 csv |

## 预期输出
- 脚本01：野生型生长率 0.8739 1/h；换目标后 ATPM=175、EX_ac_e=20（此时生物量通量=0）。
- 脚本02：三者生长率同为 0.8739，但"强制循环 FBA" sum(|v|)=1308.3，pFBA=518.4；
  图 `output/02_FBA_vs_pFBA_comparison.png`。
- 脚本03：ΔgltA 两种方法均预测致死（≈0）；Δpgi FBA=0.8632 vs MOMA=0.8426；
  图 `output/03_MOMA_knockout_comparison.png`。
- 脚本04：fraction_of_optimum=0.9 下，野生型与 Δpgi 的 FVA 区间对比；
  Δpgi 解空间收缩（部分反应 variable→fixed / →blocked）；
  图 `output/04_FVA_wildtype_vs_dpgi.png`；表 `output/04_fva_wildtype.csv`、`output/04_fva_dpgi.csv`。
- 脚本05：取样约 1000 点（整体运行约 0.5~2 分钟，取样最耗时）；
  直方图显示"最优解位于取样分布边缘"（最优≠典型）；
  相关系数表给出反应间耦合（FVA 无法提供）；
  图 `output/05_flux_sampling_distributions.png`；表 `output/05_flux_samples.csv`、`output/05_flux_correlation.csv`。

## 注意
- 本环境无 QP 求解器，quadratic MOMA 不可用，脚本03采用 linear MOMA（`moma(..., linear=True)`），教学结论方向一致。
- 脚本05默认 `OptGPSampler(model, processes=1, thinning=100)`；换用 `ACHRSampler` 只需改一行（两者接口一致）。

## 代码规范（第二轮统一）
- 变量名/函数名使用英文 snake_case；print 与 docstring 为英文；解释性注释保持中文逐行/逐块；每个脚本开头保留中文块注释（目的、用法）。
- 图表标注一律英文。

## 模型加载说明（统一策略）
- 各脚本内置 `load_ecoli_core()`：优先读取本模块目录下的本地模型文件
  （`e_coli_core.json`，BiGG 当前版，生物量反应 ID 为 `BIOMASS_Ecoli_core_w_GAM`），
  本地文件缺失时才尝试 `cobra.io.load_model("e_coli_core")` 在线下载
  （cobra 0.32 在线下载可能被 BiGG 301 重定向打断），两者都失败则给出明确报错。
- 生物量反应 ID 由 `find_biomass_rxn()` 自动探测（遍历 id 含 "biomass" 的反应，
  兜底取当前目标反应），兼容 `Biomass_Ecoli_core`（textbook 版）与
  `BIOMASS_Ecoli_core_w_GAM`（BiGG 当前版）两套命名。

## 对应课次
课程5：目标函数与优化原理（FBA / pFBA / MOMA / FVA / Flux Sampling）
