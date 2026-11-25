# -*- coding: utf-8 -*-
"""
step5_umap_analysis.py
======================
对应原仓库 single_cell/code/umap_analysis.py。

对 pFBA 通量矩阵做标准化 + UMAP 降维 (n_neighbors=100, min_dist=0.5,
n_components=2), 绘制胁迫/对照两组散点图, 输出 ../output/umap12.jpg。
"""

import argparse
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import umap
from sklearn.preprocessing import StandardScaler

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def umapfigure(embedding, a, b, outpath):
    font3 = {"family": "Arial", "weight": "normal", "size": 26}
    plt.subplots(figsize=(8, 6))
    plt.scatter(embedding[80:163, a - 1], embedding[80:163, b - 1],
                s=20, color="#3174A1", alpha=0.6, label="Unstressed")
    plt.scatter(embedding[0:80, a - 1], embedding[0:80, b - 1],
                s=20, color="#E1812C", alpha=0.6, label="Stress")
    ax = plt.gca()
    for side in ("top", "right", "bottom", "left"):
        ax.spines[side].set_linewidth(1)
    plt.xlabel("Dimension {}".format(a), fontdict=font3)
    plt.ylabel("Dimension {}".format(b), fontdict=font3)
    plt.legend(loc="best", shadow=False, scatterpoints=1, prop=font3)
    plt.xticks(fontproperties="Arial", size=24)
    plt.yticks(fontproperties="Arial", size=26)
    plt.tight_layout()
    plt.savefig(outpath, dpi=600)
    plt.close()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--skip-clustering", action="store_true")
    args = ap.parse_args()

    umap_data = pd.read_excel(os.path.join(BASE, "output", "pfba_flux.xlsx"))
    umap_data = umap_data.iloc[:, 1:]
    umap_data.loc[:, "metnum/20"] = umap_data.loc[:, "metnum/20"] * 20
    umap_data.loc[:, "rxnnum/20"] = umap_data.loc[:, "rxnnum/20"] * 20

    reducer = umap.UMAP(n_neighbors=100, min_dist=0.5, n_components=2,
                        metric="euclidean")
    scaled = StandardScaler().fit_transform(umap_data)
    embedding = reducer.fit_transform(scaled)
    umapfigure(embedding, 1, 2, os.path.join(BASE, "output", "umap12.jpg"))
    print("umap12.jpg 已保存")
if __name__ == "__main__":
    main()
