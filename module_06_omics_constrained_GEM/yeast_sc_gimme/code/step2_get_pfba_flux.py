# -*- coding: utf-8 -*-
"""
step2_get_pfba_flux.py
======================
对应原仓库 single_cell/code/get_pfba_flux.py。

对 output/scmodel/ 中的 163 个单细胞模型逐一运行 pFBA
(parsimonious FBA: 在保持最优目标值的前提下最小化总通量),
将每个细胞模型在全模型反应集合上的通量对齐成矩阵, 并附加两个模型大小
特征 (rxnnum/20, metnum/20), 输出 ../output/pfba_flux.xlsx。

与原脚本一致:
  - 反应行顺序按字母序排列文件名 (Stressed 80 个在前, Unstressed 83 个在后);
  - 单细胞模型中不存在的反应通量记为 0;
  - pFBA 不可行的模型跳过并打印提示。

用法:  python3 step2_get_pfba_flux.py [--n-jobs 8]
"""

import argparse
import os
import time
from multiprocessing import Pool

import numpy as np
import pandas as pd
import cobra

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SCDIR = os.path.join(BASE, "output", "scmodel")

_RXN_IDS = None


def _pfba_one(fname):
    path = os.path.join(SCDIR, fname)
    name = fname[:-4] if fname.endswith(".xml") else fname
    try:
        m = cobra.io.read_sbml_model(path)
        m.solver = "glpk"
        try:
            sol = cobra.flux_analysis.pfba(m)
        except Exception:
            # GLPK 单纯形在个别模型 pFBA 第二阶段 LP 上数值失败,
            # 回退到精确有理单纯形 (同属 GLPK 家族, 结果一致)
            m.solver = "glpk_exact"
            sol = cobra.flux_analysis.pfba(m)
        if sol.status != "optimal":
            return name, None, None, None
        flux = np.zeros(len(_RXN_IDS))
        for rid, v in sol.fluxes.items():
            j = _IDX.get(rid)
            if j is not None:
                flux[j] = v
        return name, flux, len(m.reactions), len(m.metabolites)
    except Exception:
        return name, None, None, None


# 反应 id -> 列索引 (worker 内构建)
_IDX = None


def _init():
    global _RXN_IDS, _IDX
    full = cobra.io.read_sbml_model(os.path.join(BASE, "data", "yeast-GEM.xml"))
    _RXN_IDS = [r.id for r in full.reactions]
    _IDX = {r: i for i, r in enumerate(_RXN_IDS)}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-jobs", type=int, default=8)
    ap.add_argument("--scdir", type=str, default=None,
                    help="单细胞模型 SBML 目录 (默认 ../output/scmodel)")
    args = ap.parse_args()

    global SCDIR
    if args.scdir:
        SCDIR = os.path.abspath(args.scdir)

    files = sorted(f for f in os.listdir(SCDIR) if f.endswith(".xml"))
    print(f"共 {len(files)} 个单细胞模型")
    assert files[0].startswith("BY4741_Stressed") and files[-1].startswith("BY4741_Unstressed"), \
        "文件排序异常: 应为 80 个 Stressed 在前, 83 个 Unstressed 在后"

    t0 = time.time()
    results = []
    with Pool(args.n_jobs, initializer=_init) as pool:
        for k, res in enumerate(pool.imap(_pfba_one, files)):
            results.append(res)
            if (k + 1) % 20 == 0 or k + 1 == len(files):
                print(f"[{k + 1}/{len(files)}]")

    _init()
    pca_data = pd.DataFrame(index=_RXN_IDS, columns=[r[0] for r in results], dtype=float)
    rxnnum, metnum, keep_cols = [], [], []
    for name, flux, nr, nm in results:
        if flux is None:
            print(f"{name}: pFBA 不可行, 跳过")
            continue
        pca_data[name] = flux
        rxnnum.append(nr)
        metnum.append(nm)
        keep_cols.append(name)
    pca_data = pca_data[keep_cols].T

    # 与原脚本一致: 附加模型大小特征 (除以 20 缩放)
    pca_data.insert(0, "rxnnum/20", [r / 20 for r in rxnnum])
    pca_data.insert(0, "metnum/20", [m / 20 for m in metnum])
    pca_data.to_excel(os.path.join(BASE, "output", "pfba_flux.xlsx"))
    print(f"pfba_flux.xlsx 已保存, 形状 {pca_data.shape}, 耗时 {time.time() - t0:.0f}s")


if __name__ == "__main__":
    main()
