# 酵母单细胞代谢建模（GIMME）分析流程的 Python 复现

本工程用 **纯 Python** 完整复现了
[hongzhonglu/yeast_GEM_multi_omics_analysis](https://github.com/hongzhonglu/yeast_GEM_multi_omics_analysis/tree/main/single_cell)
中 `single_cell/` 目录的单细胞建模分析流程。原流程的核心建模步骤依赖
MATLAB（RAVEN + COBRA Toolbox），本复现将其全部移植到 Python，并已端到端
运行验证。

## 生物学背景

- 模型：yeast-GEM v9.0.0（酿酒酵母基因组尺度代谢模型，4130 反应 /
  2805 代谢物 / 1162 基因，目标函数为生长反应 `r_2111`）
- 数据：GSE102475（Gasch 2017），BY4741 酵母单细胞 RNA-seq 标准化表达矩阵，
  6697 基因 × 163 细胞（80 个 0.7 M NaCl 胁迫细胞 + 83 个对照细胞）
- 方法：**GIMME**（Becker & Palsson 2008）——以表达阈值 20 为界，在保持
  最大生长速率 90% 的约束下，最小化低表达反应的加权使用，从而提取每个
  细胞特异的代谢子网络

## 流程步骤与原仓库文件对应关系

| 步骤 | 本工程脚本 | 原仓库文件 | 说明 |
|------|-----------|-----------|------|
| 0 | `code/gimme.py` | `code/get_GPR.m`、`code/singlecellmodel.m` 中调用的 COBRA Toolbox 函数（`findUsedGenesLevels`、`selectGeneFromGPR`、`GIMME`） | GIMME 算法核心移植 |
| 1 | `code/step1_build_single_cell_models.py` | `code/singlecellmodel.m` | 构建 163 个单细胞模型（SBML） |
| 2 | `code/step2_get_pfba_flux.py` | `code/get_pfba_flux.py` | 逐模型 pFBA，生成通量矩阵 `pfba_flux.xlsx` |
| 3 | `code/step3_growth_violin.py` | `code/Grow bar.py` | 两组生长速率小提琴图 |
| 4 | `code/step4_t_test.py` | `code/t_test.py` | 逐反应独立双样本 t 检验 + 饼图 |
| 5 | `code/step5_umap_analysis.py` | `code/umap_analysis.py` | UMAP 可视化 + Leiden/层次聚类 |
| 6 | `code/step6_classification.py` | `code/classification.py` | 随机森林分类（6 种特征表示） |

## GIMME 算法（Python 移植要点）

`code/gimme.py` 严格对应 COBRA Toolbox 的实现，分四个函数：

1. **`get_parsed_gpr`**（对应 `get_GPR.m`）：把 grRule 字符串按 `' or '` 拆分、
   去括号后按 `' and '` 拆分，得到「OR 子句 × AND 基因」结构；
2. **`find_used_genes_levels`**（对应 `findUsedGenesLevels.m`）：表达数据映射到
   模型基因，缺失记 -1，重复基因取均值；
3. **`select_gene_from_gpr`**（对应 `selectGeneFromGPR.m`，minSum=false）：
   AND→min，OR→max，得到反应表达值 `expressionRxns`；
4. **`gimme`**（对应 `GIMME.m`）：
   - 将模型转为不可逆形式（可逆反应拆成正/反两个方向）；
   - 最大化目标函数得 `maxObjective`；
   - 对 `-1 < expr < threshold` 的反应设惩罚 `threshold - expr`，
     约束 `目标 ≥ 0.9 × maxObjective`，最小化总惩罚（稀疏 LP，**GLPK** 求解）；
   - `expr > threshold` 或无数据（-1）的反应保留；低表达但携带通量的保留；
     其余删除，同时清理孤儿代谢物。

> 求解器说明：原 MATLAB 脚本用 `changeCobraSolver('glpk')` 指定 GLPK。
> 这一步至关重要——GIMME 的惩罚 LP 存在大量退化解，内点法（如 HiGHS）
> 给出的解带有数值尘埃（如某反应通量为 -1e-8），会把「负荷反应」误判为
> 无活性从而杀死个别细胞的模型；GLPK 单纯形顶点解则给出干净的 0 通量，
> 与参考实现行为一致。本复现通过 optlang 的 GLPK 接口求解该 LP（每个进程
> 预构建一次问题，每个细胞仅替换目标函数表达式）。

> 效率说明：原 MATLAB 实现对每个细胞都重建不可逆模型并两次调用 GLPK。
> 本复现中不可逆模型的稀疏 LP 结构只构建一次（所有细胞共享），每个细胞
> 仅需一次 LP 求解，数学上与原版完全等价。

### 与原仓库的一处重要差异说明

原 `singlecellmodel.m` 中有一行被注释掉的代码：

```matlab
% model.lb(find(strcmp('r_0446', model.rxns()))) = -1000;
```

但仓库发布的单细胞模型（`output/scmodel/*.xml`）中 `r_0446` 的下界实际
**为 -1000**，且参考通量表中的最大生长速率 0.087859 恰好等于放开该边界后
全模型的最大生长（原边界下为 0.0859）。这说明作者实际运行时该修改是生效
的。为了复现已发表的结果，`step1` 中启用了该边界修改（代码中有注释标明）。
如需严格按上传的脚本运行，删除 `step1_build_single_cell_models.py` 中的
对应一行即可。

## 运行方法

```bash
pip install -r requirements.txt

cd code
python3 step1_build_single_cell_models.py --n-jobs 8   # 构建 163 个单细胞模型
python3 step2_get_pfba_flux.py --n-jobs 8              # pFBA 通量矩阵
python3 step3_growth_violin.py                          # 生长小提琴图
python3 step4_t_test.py                                 # t 检验
python3 step5_umap_analysis.py                          # UMAP + 聚类
python3 step6_classification.py                         # 随机森林分类
```

Leiden 聚类需要额外安装 `igraph` 和 `leidenalg`；不需要时可用
`python3 step5_umap_analysis.py --skip-clustering` 跳过。

## 复现验证

本工程已端到端完整运行（163 个细胞全部成功建模并通过 pFBA），关键结果与
原仓库一致：

### 模型构建

- 163 个单细胞模型全部可以生长（FBA 生长速率 0.0849-0.0879 h⁻¹，与参考
  模型范围一致）；细胞 `BY4741_Stressed_..._C50_BY1` 与原仓库发布的对应
  SBML 对比：反应集合重合度 99.4%（差异全部为低表达的边缘反应，GIMME LP
  存在多重最优解，不同求解器的选择略有不同）；BY1 的 pFBA 生长速率
  0.084897 h⁻¹ 与参考值 0.084897 h⁻¹ 精确吻合。

### 生长速率（step3）

| 组别 | 细胞数 | 平均生长速率 (h⁻¹) |
|------|-------|-------------------|
| Stressed（胁迫） | 80 | 0.08615 |
| Unstressed（对照） | 83 | 0.08697 |

与原结果一致：对照组平均生长速率略高于胁迫组。

### t 检验（step4）

非零通量反应 1803 个，其中 p ≤ 0.05 的 452 个。

### 随机森林分类（step6）

| 特征表示 | 测试准确率 | 参考结果 |
|----------|-----------|---------|
| 转录组（rf trans） | **1.0000** | 1.00 |
| 通量（rf flux） | **0.7959** | 0.73-0.76 |
| 转录组 + KernelPCA | 1.0000 | — |
| 转录组 + UMAP | 0.9592 | — |
| 通量 + KernelPCA | 0.7959 | — |
| 通量 + UMAP | 0.5918 | — |

转录组特征 100% 区分两组细胞、通量特征约 80% 的核心结论与原文一致。

### 两个已知的技术说明

1. **pFBA 数值细节**：2 个细胞（`Stressed_BY40`、`Unstressed_BY17`）的模型
   FBA 生长正常，但 GLPK 单纯形在 pFBA 第二阶段（最小化总通量）的 LP 上
   数值失败；这两个细胞改用 `glpk_exact`（精确有理单纯形，同属 GLPK）求解，
   生长速率与其余 161 个细胞的结果完全在同一范围（0.08498 / 0.08786）。
2. **Leiden 聚类**：原仓库 `umap_analysis.py` 中该段代码为注释状态，本复现
   按其原样移植（全连接余弦相似度图 + resolution=1.01）；在该参数下本数据
   得到的是每细胞一个社区，说明该注释段本身是探索性代码，不影响主流程结论。
