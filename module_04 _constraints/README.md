# 模块04：GEMs运行所需的约束（对应课程6）

## 目标
理解约束（培养基、摄取速率、通量上下界、方向性、维持能耗）如何定义可行域并决定 FBA 预测；
通过单因素扫描与二维表型相面，体会"约束越多，预测越具体"。

## 依赖
```bash
pip install cobra pandas matplotlib
```
模型加载见文末"模型加载说明"；本模块附带本地模型文件 `e_coli_core.json`。

## 脚本与运行顺序
| 顺序 | 脚本 | 内容 |
|---|---|---|
| 1 | `01_培养基与摄取约束.py` | 查看 model.medium；葡萄糖/氧气/氨摄取下界扫描；厌氧模拟与发酵产物谱对比 |
| 2 | `02_通量上下界与方向约束.py` | 反应边界与不可逆性；ATPM 维持消耗下界扫描（含致死临界）；收紧 CYTBD 上界演示溢流代谢 |
| 3 | `03_组合约束与表型相面.py` | 葡萄糖×氧气 21×26 网格手写 FBA 扫描（546 次）；生长率与乙酸分泌双热力图 |

## 预期输出
- 脚本01：葡萄糖 10→0 生长率 0.8739→0；O2=0 厌氧生长 0.2117，分泌乙酸 8.50、乙醇 8.28、甲酸 17.8；
  图 `output/01_medium_and_secretion_comparison.png`。
- 脚本02：ATPM 下界 0→55 生长率 0.9166→0.6366，下界 180 时 LP infeasible（致死）；
  CYTBD 上界收紧时乙酸溢出；图 `output/02_bound_constraint_effects.png`。
- 脚本03：相面统计——respiratory growth 289 点、O2-limited fermentation 221 点、no growth 36 点；
  图 `output/03_glucose_oxygen_phenotype_phase_plane.png`。

## 注意
- 葡萄糖=0 时因 ATPM 维持约束 LP 不可行，脚本按"不生长"处理并静音相应警告。

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
课程6：GEMs运行所需的约束
