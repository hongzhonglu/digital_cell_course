# -*- coding: utf-8 -*-
"""
step3_growth_violin.py
======================
对应原仓库 single_cell/code/Grow bar.py。

绘制胁迫/对照两组单细胞模型生长速率 (r_2111 的 pFBA 通量) 的小提琴图,
输出 ../output/Grow violin.jpg。

注: 原脚本使用了 pandas 2.x 中已移除的 DataFrame.append,
本复现改为 pd.concat (行为等价)。
"""

import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

data = pd.read_excel(os.path.join(BASE, "output", "pfba_flux.xlsx"))

# 前 80 行为胁迫细胞, 后 83 行为对照细胞
a = pd.DataFrame({"group": np.repeat("Stress", 80), "value": data.loc[0:79, "r_2111"].values})
b = pd.DataFrame({"group": np.repeat("Unstressed", 83), "value": data.loc[80:, "r_2111"].values})
df = pd.concat([b, a], ignore_index=True)

font3 = {"family": "Arial", "weight": "normal", "size": 26}

plt.subplots(figsize=(8, 6))
sns.violinplot(x="group", y="value", data=df)
plt.xticks(fontproperties="Arial", size=26)
plt.yticks(fontproperties="Arial", size=24)
plt.xlabel("")
plt.ylabel("Growth rate(h-1)", fontdict=font3)
plt.tight_layout()
plt.savefig(os.path.join(BASE, "output", "Grow violin.jpg"), dpi=800)
print("Grow violin.jpg 已保存")
print("Stress 组均值 %.5f, Unstressed 组均值 %.5f (h-1)"
      % (a.value.mean(), b.value.mean()))
