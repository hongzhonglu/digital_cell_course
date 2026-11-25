# -*- coding: utf-8 -*-
"""
step1_build_single_cell_models.py
=================================
对应原仓库 single_cell/code/singlecellmodel.m (MATLAB/RAVEN+COBRA)。

流程:
  1. 读取 yeast-GEM.xml (yeastGEM v9.0.0, 4130 反应);
  2. 读取 GSE102475 (Gasch 2017) 单细胞 RNA-seq 标准化表达矩阵
     (6697 基因 x 163 细胞: 前 80 个为 NaCl 胁迫细胞, 后 83 个为对照细胞);
  3. 对每个细胞, 用 GIMME 算法 (threshold=20, obj_frac=0.9) 提取
     上下文特异模型;
  4. 每个单细胞模型写出为 SBML 文件到 ../output/scmodel/<细胞名>.xml。

用法:  python3 step1_build_single_cell_models.py [--n-jobs 8] [--cells 10]
"""

import argparse
import os
import shutil
import sys
import tempfile
import time
from multiprocessing import Pool

import numpy as np
import pandas as pd
import cobra

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from gimme import (get_parsed_gpr_all, prepare_gimme, build_single_cell_model,
                   _build_glpk_problem)

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(BASE, "data")
OUT = os.path.join(BASE, "output", "scmodel")

# 全局变量 (每个 worker 进程初始化一次)
_MODEL = None
_PREPARED = None
_PARSED_GPR = None
_EXPR_GENES = None
_EXPR_MATRIX = None
_CELL_NAMES = None


def _init_worker():
    global _MODEL, _PREPARED, _PARSED_GPR, _EXPR_GENES, _EXPR_MATRIX, _CELL_NAMES
    _MODEL = cobra.io.read_sbml_model(os.path.join(DATA, "yeast-GEM.xml"))
    _MODEL.solver = "glpk"
    # 注意: 原仓库 singlecellmodel.m 中下面这行被注释掉了:
    #     % model.lb(find(strcmp('r_0446', model.rxns()))) = -1000;
    # 但仓库发布的单细胞模型 (output/scmodel/*.xml) 中 r_0446 的下界实际为
    # -1000, 且参考通量表中的最大生长速率 0.087859 恰好等于放开该边界后的
    # 全模型最大生长 —— 说明作者实际运行时该修改是生效的。为复现已发表的
    # 结果, 此处启用该边界修改。
    _MODEL.reactions.get_by_id("r_0446").lower_bound = -1000
    _PREPARED = prepare_gimme(_MODEL)
    # 预构建 GLPK 惩罚 LP 问题 (每进程一次); 原 MATLAB 脚本用
    # changeCobraSolver('glpk'), GLPK 单纯形顶点解可避免内点法在退化解上
    # 的数值尘埃 (该尘埃会把负荷反应误判为无活性从而杀死模型)。
    _build_glpk_problem(_PREPARED, 0.9)
    _PARSED_GPR = get_parsed_gpr_all(_MODEL)
    expr = pd.read_excel(os.path.join(DATA, "GSE102475_GASCH_NaCl-scRNAseq_NormData.xlsx"))
    _EXPR_GENES = expr.iloc[:, 0].astype(str).values
    _CELL_NAMES = expr.columns[1:].tolist()
    _EXPR_MATRIX = expr.iloc[:, 1:].values.astype(float)


def _write_sbml_verified(cs_model, outfile):
    """先写本地临时文件, 再复制到目标位置并逐字节校验。

    持久盘 (FUSE 挂载) 对 libsbml 的写入模式可能出现空洞/截断,
    顺序写的小文件 + 复制 + 校验可规避。
    """
    tmp = tempfile.NamedTemporaryFile(suffix=".xml", delete=False, dir="/tmp")
    tmp.close()
    try:
        cobra.io.write_sbml_model(cs_model, tmp.name)
        with open(tmp.name, "rb") as f:
            ref = f.read()
        for _ in range(8):
            try:
                shutil.copyfile(tmp.name, outfile)
                with open(outfile, "rb") as f:
                    got = f.read()
                if got == ref:
                    return
                os.unlink(outfile)
            except OSError:
                pass  # FUSE 抖动 / 挂载瞬时消失, 重试
            time.sleep(1.0)
        raise IOError(f"SBML 写出校验失败: {outfile}")
    finally:
        try:
            os.unlink(tmp.name)
        except OSError:
            pass


def _build_one(cell_idx):
    """构建单个细胞模型并写出 SBML; 返回 (细胞名, 反应数, 代谢物数, 耗时)。"""
    name = _CELL_NAMES[cell_idx]
    outfile = os.path.join(OUT, name + ".xml")
    if os.path.exists(outfile):
        return name, -1, -1, 0.0  # 断点续跑: 已存在则跳过
    t0 = time.time()
    cs_model = build_single_cell_model(
        _PREPARED, _MODEL, _PARSED_GPR, _EXPR_GENES,
        _EXPR_MATRIX[:, cell_idx],
        threshold=20, obj_frac=0.9)   # threshold=20 与 singlecellmodel.m 一致
    cs_model.id = name
    _write_sbml_verified(cs_model, outfile)
    return name, len(cs_model.reactions), len(cs_model.metabolites), time.time() - t0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-jobs", type=int, default=8)
    ap.add_argument("--cells", type=int, default=None,
                    help="只构建前 N 个细胞 (调试用)")
    ap.add_argument("--scdir", type=str, default=None,
                    help="单细胞模型 SBML 输出目录 (默认 ../output/scmodel)")
    args = ap.parse_args()

    global OUT
    if args.scdir:
        OUT = os.path.abspath(args.scdir)
    os.makedirs(OUT, exist_ok=True)
    _init_worker()
    n_cells = len(_CELL_NAMES) if args.cells is None else min(args.cells, len(_CELL_NAMES))
    print(f"共 {n_cells} 个细胞, GIMME threshold=20, obj_frac=0.9; "
          f"全模型最大生长 = {_PREPARED.max_objective[0]:.5f} h-1")

    t0 = time.time()
    stats = []
    with Pool(args.n_jobs, initializer=_init_worker) as pool:
        for k, res in enumerate(pool.imap_unordered(_build_one, range(n_cells))):
            stats.append(res)
            if (k + 1) % 10 == 0 or k + 1 == n_cells:
                print(f"[{k + 1}/{n_cells}] 最近完成: {res[0]} "
                      f"({res[1]} rxns, {res[2]} mets, {res[3]:.1f}s)")

    df = pd.DataFrame(stats, columns=["cell", "n_rxns", "n_mets", "seconds"])
    df.to_csv(os.path.join(BASE, "output", "step1_model_stats.csv"), index=False)
    print(f"全部完成, 总耗时 {time.time() - t0:.0f}s; "
          f"模型大小 {df.n_rxns.min()}-{df.n_rxns.max()} rxns")


if __name__ == "__main__":
    main()
