# -*- coding: utf-8 -*-
"""
step4_t_test.py
===============
对应原仓库 single_cell/code/t_test.py。

对胁迫 (前 80 行) 与对照 (后 83 行) 两组间同一反应的 pFBA 通量做
独立双样本 t 检验, 输出:
  - ../output/stress_flux_mean.xlsx    胁迫组各反应平均通量
  - ../output/unstress_flux_mean.xlsx  对照组各反应平均通量
  - ../output/t-test.tif               p 值 <= 0.05 / > 0.05 反应比例饼图
"""

import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

# 反应子系统注释 (原脚本读入, 用于后续通路分析)
pathway = pd.read_table(os.path.join(BASE, "data", "uniqueSubsystems.tsv"))
pathway.index = pathway.loc[:, "ID"]

t_data = pd.read_excel(os.path.join(BASE, "output", "pfba_flux.xlsx"))
t_data = t_data.iloc[:, 3:]  # 去掉索引列与 metnum/20、rxnnum/20

# 删除全零反应列
co = t_data.columns.tolist()
for z in co:
    if sum(t_data.loc[:, z]) == 0:
        t_data = t_data.drop(z, axis=1)

stress = t_data.loc[0:79, :]
unstress = t_data.loc[80:, :]
none_zero_l = stress.columns.tolist()

rxn_p, p_xiao, p_da = {}, {}, {}
for te in none_zero_l:
    sample1 = np.asarray(stress.loc[:, te], dtype=float)
    sample2 = np.asarray(unstress.loc[:, te], dtype=float)
    r = stats.ttest_ind(sample1, sample2)
    rxn_p[te] = r.pvalue
    if r.pvalue <= 0.05:
        p_xiao[te] = r.pvalue
    else:
        p_da[te] = r.pvalue

stress_mean = stress.mean()
unstress_mean = unstress.mean()
stress_mean.to_excel(os.path.join(BASE, "output", "stress_flux_mean.xlsx"))
unstress_mean.to_excel(os.path.join(BASE, "output", "unstress_flux_mean.xlsx"))

# 饼图
name_l = ["p value > 0.05", "p value ≤ 0.05"]
values = [len(none_zero_l) - len(p_xiao), len(p_xiao)]
font1 = {"family": "Times New Roman", "weight": "normal", "size": 23}
font2 = {"family": "Times New Roman", "weight": "normal", "size": 20}
font3 = {"family": "Times New Roman", "weight": "normal", "size": 12}

plt.subplots(figsize=(8, 6))
plt.pie(values,
        labeldistance=1,
        wedgeprops={"linewidth": 3, "edgecolor": "white"},
        colors=["#64F0E1", "#F09165"],
        textprops=font2)
plt.title("t-test of nonzero flux", fontdict=font1)
plt.legend(name_l, bbox_to_anchor=(0.9, 0.9), prop=font3)
plt.savefig(os.path.join(BASE, "output", "t-test.tif"))
print(f"t-test.tif 已保存: 非零通量反应 {len(none_zero_l)} 个, "
      f"其中 p<=0.05 的 {len(p_xiao)} 个")
