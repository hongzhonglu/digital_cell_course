# -*- coding: utf-8 -*-
"""
gimme.py — GIMME 算法的纯 Python 移植
=====================================

对应原仓库 single_cell/code 中的 MATLAB 函数:

1. ``get_GPR.m``                -> :func:`get_parsed_gpr` / :func:`get_parsed_gpr_all`
   将 grRules 字符串按 ' or ' 拆分、去括号后再按 ' and ' 拆分,
   得到 "OR 子句列表, 每个子句是 AND 基因列表" 的结构。

2. ``findUsedGenesLevels.m``    -> :func:`find_used_genes_levels`
   (COBRA Toolbox) 将表达数据映射到模型基因上; 无数据的基因记为 -1,
   表达数据中重复的基因取均值。

3. ``selectGeneFromGPR.m``      -> :func:`select_gene_from_gpr`
   (COBRA Toolbox, minSum=false 模式) 由 GPR 计算反应表达值:
   AND -> min, OR -> max; 无 GPR/无数据的反应记为 -1。

4. ``GIMME.m``                  -> :func:`prepare_gimme` + :func:`gimme`
   (COBRA Toolbox, Becker & Palsson 2008) 上下文特异模型提取:
   - 将模型转为不可逆形式 (可逆反应拆成正/反两个方向);
   - 先最大化目标函数得到 maxObjective;
   - 对表达值满足 -1 < expr < threshold 的不可逆反应施加惩罚
     (threshold - expr) * v, 约束目标 >= obj_frac * maxObjective
     (obj_frac 默认 0.9), 最小化总惩罚;
   - 反应活性: expr > threshold 或 expr == -1 -> 保留 (activity 1);
     低表达但 GIMME 解中携带通量 -> 保留 (activity 2);
     否则从模型中删除。

为兼顾 163 个细胞的批量计算, LP 部分基于稀疏矩阵 + HiGHS 求解器实现,
不可逆模型结构只构建一次 (见 :func:`prepare_gimme`), 数学上与
COBRA Toolbox 的 GIMME.m 完全等价。
"""

from __future__ import annotations

from dataclasses import dataclass

import os
import tempfile

import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy.optimize import linprog

import cobra
from cobra.util.array import create_stoichiometric_matrix


# ---------------------------------------------------------------------------
# 1. get_GPR.m 的移植
# ---------------------------------------------------------------------------
def get_parsed_gpr(gr_rule: str) -> list[list[str]]:
    """解析单个反应的 grRule 字符串。

    与 MATLAB ``get_GPR.m`` 完全一致的简化解析:
    先按 ``' or '`` 拆分, 去掉括号, 再按 ``' and '`` 拆分。

    Returns
    -------
    list[list[str]]
        OR 子句列表, 每个子句是 AND 连接的基因列表。
        空规则返回 ``[['']]`` (与 MATLAB 版本行为一致)。
    """
    if not gr_rule or not gr_rule.strip():
        return [[""]]
    clauses = []
    for or_part in gr_rule.split(" or "):
        cleaned = or_part.replace("(", "").replace(")", "")
        genes = [g.strip() for g in cleaned.split(" and ")]
        clauses.append(genes)
    return clauses


def get_parsed_gpr_all(model: cobra.Model) -> list[list[list[str]]]:
    """对整个模型调用 :func:`get_parsed_gpr` (对应 MATLAB 中的 parsedGPR 细胞数组)。"""
    return [get_parsed_gpr(r.gene_reaction_rule) for r in model.reactions]


# ---------------------------------------------------------------------------
# 2. findUsedGenesLevels.m 的移植
# ---------------------------------------------------------------------------
def find_used_genes_levels(model_genes,
                           expr_genes,
                           expr_values) -> tuple[list[str], np.ndarray]:
    """将表达数据映射到模型基因。

    - 模型基因在表达数据中找不到 -> 表达值记为 -1 (无数据标记);
    - 表达数据中出现多次的基因 -> 取均值。
    """
    expr_genes = np.asarray(expr_genes).astype(str)
    expr_values = np.asarray(expr_values, dtype=float)

    uniq, inv = np.unique(expr_genes, return_inverse=True)
    sums = np.zeros(len(uniq))
    cnts = np.zeros(len(uniq))
    np.add.at(sums, inv, expr_values)
    np.add.at(cnts, inv, 1.0)
    lookup = dict(zip(uniq, sums / cnts))

    gene_expr = np.array([lookup.get(g, -1.0) for g in model_genes], dtype=float)
    return list(model_genes), gene_expr


# ---------------------------------------------------------------------------
# 3. selectGeneFromGPR.m 的移植 (minSum = false)
# ---------------------------------------------------------------------------
def select_gene_from_gpr(model: cobra.Model,
                         gene_names,
                         gene_exp: np.ndarray,
                         parsed_gpr: list[list[list[str]]],
                         min_sum: bool = False) -> np.ndarray:
    """由 GPR 规则把基因表达映射为反应表达值。

    minSum=false (原仓库用法): AND -> min, OR -> max。
    无 GPR 或无任何基因数据的反应 -> -1。
    与原 MATLAB 实现一致: min/max 在包含 -1 (无数据) 的值上直接计算。
    """
    exp_map = dict(zip(gene_names, gene_exp))
    expression_col = np.full(len(model.reactions), -1.0)

    for i, clauses in enumerate(parsed_gpr):
        cur_expr = []
        for genes in clauses:
            if len(genes) >= 1:
                vals = [exp_map[g] for g in genes if g in exp_map]
                if vals:
                    if min_sum:
                        cur_expr.append(float(np.sum(vals)))
                    else:
                        cur_expr.append(float(np.min(vals)))  # AND -> min
        if cur_expr:
            if min_sum:
                expression_col[i] = float(np.min(cur_expr))
            else:
                expression_col[i] = float(np.max(cur_expr))  # OR -> max
    return expression_col


# ---------------------------------------------------------------------------
# 4. GIMME.m 的移植 (稀疏 LP 实现)
# ---------------------------------------------------------------------------
@dataclass
class GimmeLP:
    """预处理好的不可逆模型 LP 结构 (对所有细胞共享)。"""
    S: sp.csc_matrix          # 不可逆模型化学计量矩阵 (m x n_irrev)
    lb: np.ndarray            # 不可逆反应下界
    ub: np.ndarray            # 不可逆反应上界
    irrev2rev: np.ndarray     # 不可逆反应 j -> 原反应索引
    obj_fwd: list[int]        # 目标反应在不可逆模型中的正向副本索引
    max_objective: np.ndarray  # 各目标的最大值 (对全模型只求一次)
    rxn_ids: list[str]        # 原模型反应 id
    glpk: object = None       # 可选: 预构建的 optlang GLPK 问题
    glpk_vars: list = None    # GLPK 变量列表
    glpk_mul_terms: list = None  # 目标函数各项 (coef 占位)


def _max_objective_lp(c, S, lb, ub):
    """求解初始 FBA (linprog 约定: 最小化 c^T x), 多求解器回退。

    某些平台/scipy 版本上 HiGHS 默认变体对该 LP 报数值失败 (status=4),
    依次尝试 HiGHS 各变体, 最后回退到 GLPK (与参考实现一致的求解器)。
    """
    for method in ("highs", "highs-ds", "highs-ipm"):
        res = linprog(c, A_eq=S, b_eq=np.zeros(S.shape[0]),
                      bounds=list(zip(lb, ub)), method=method)
        if res.status == 0:
            return -res.fun
    # GLPK 回退
    from optlang.glpk_interface import (Model as GLPKModel, Variable,
                                        Constraint, Objective)
    import optlang
    n = S.shape[1]
    glp = GLPKModel()
    vars_ = [Variable("x%d" % j, lb=float(lb[j]), ub=float(ub[j]))
             for j in range(n)]
    glp.add(vars_)
    Scoo = S.tocoo()
    Srow = [dict() for _ in range(S.shape[0])]
    for r_, c_, v_ in zip(Scoo.row, Scoo.col, Scoo.data):
        Srow[r_][c_] = float(v_)
    for d in Srow:
        e = (optlang.symbolics.add([v_ * vars_[j_] for j_, v_ in d.items()])
             if d else 0)
        glp.add([Constraint(e, lb=0, ub=0)])
    nz = np.nonzero(c)[0]
    glp.objective = Objective(
        optlang.symbolics.add([float(c[j]) * vars_[j] for j in nz]),
        direction="min")
    if glp.optimize() != "optimal":
        raise RuntimeError("GIMME 初始 FBA 求解失败 (所有求解器均失败)")
    # objective.value 等价于 linprog 的 res.fun (最小化值), 取负得最大值
    return -glp.objective.value


def prepare_gimme(model: cobra.Model) -> GimmeLP:
    """构建不可逆模型 LP 结构并求解目标函数最大值 (对应 GIMME.m 前半部分)。

    convertToIrreversible 规则:
      - 可逆反应 (lb < 0 且 ub > 0): 正向副本界 [max(lb,0), ub],
        反向副本 S 列取负, 界 [max(-ub,0), -lb];
      - 纯负向反应 (ub <= 0): 翻转为正向;
      - 其余保持不变。
    """
    S0 = create_stoichiometric_matrix(model, array_type="lil")
    lb0 = np.array([r.lower_bound for r in model.reactions], dtype=float)
    ub0 = np.array([r.upper_bound for r in model.reactions], dtype=float)
    rxn_ids = [r.id for r in model.reactions]
    n = len(rxn_ids)

    cols = []          # 不可逆模型各列 (coo 三元组)
    irrev_lb, irrev_ub, irrev2rev = [], [], []

    S0c = S0.tocoo()
    # 按列组织原 S
    S_cols = [None] * n
    S_csr = S0.tocsr()
    S_csc = S0.tocsc()

    for i in range(n):
        col = S_csc.getcol(i)
        rows, vals = col.indices, col.data
        lb, ub = lb0[i], ub0[i]
        if lb < 0 and ub > 0:
            # 正向副本
            cols.append((rows, vals, len(irrev_lb)))
            irrev_lb.append(max(lb, 0.0)); irrev_ub.append(ub); irrev2rev.append(i)
            # 反向副本 (S 取负)
            cols.append((rows, -vals, len(irrev_lb)))
            irrev_lb.append(max(-ub, 0.0)); irrev_ub.append(-lb); irrev2rev.append(i)
        elif lb < 0 and ub <= 0:
            # 纯负向: 翻转
            cols.append((rows, -vals, len(irrev_lb)))
            irrev_lb.append(-ub); irrev_ub.append(-lb); irrev2rev.append(i)
        else:
            cols.append((rows, vals, len(irrev_lb)))
            irrev_lb.append(lb); irrev_ub.append(ub); irrev2rev.append(i)

    n_irrev = len(irrev_lb)
    r_all = np.concatenate([c[0] for c in cols]) if cols else np.array([], dtype=int)
    v_all = np.concatenate([c[1] for c in cols]) if cols else np.array([], dtype=float)
    c_all = np.concatenate([np.full(len(c[0]), c[2]) for c in cols]) if cols else np.array([], dtype=int)
    S = sp.coo_matrix((v_all, (r_all, c_all)),
                      shape=(S0.shape[0], n_irrev)).tocsc()

    lb = np.array(irrev_lb)
    ub = np.array(irrev_ub)
    irrev2rev = np.asarray(irrev2rev)

    # 目标反应及其在不可逆模型中的正向副本
    obj_idx = [i for i, r in enumerate(model.reactions)
               if abs(r.objective_coefficient) > 0]
    if not obj_idx:
        raise ValueError("模型没有设定目标函数")
    obj_fwd = []
    for oi in obj_idx:
        copies = np.where((irrev2rev == oi))[0]
        # 正向副本: 对应原反应方向 (S 列未取负的那个, 即第一个)
        obj_fwd.append(int(copies[0]))

    # 目标函数最大值 (对不可逆全模型最大化, 每个细胞相同, 只求一次)
    max_obj = np.zeros(len(obj_fwd))
    for k, ji in enumerate(obj_fwd):
        c = np.zeros(n_irrev)
        c[ji] = 1.0
        max_obj[k] = _max_objective_lp(-c, S, lb, ub)

    return GimmeLP(S=S, lb=lb, ub=ub, irrev2rev=irrev2rev, obj_fwd=obj_fwd,
                   max_objective=max_obj, rxn_ids=rxn_ids)


def _build_glpk_problem(prepared: GimmeLP, obj_frac: float):
    """基于不可逆 LP 结构预构建 optlang/GLPK 问题 (每个进程一次)。

    原 MATLAB 脚本用 ``changeCobraSolver('glpk')`` 指定 GLPK 求解;
    GIMME 的惩罚 LP 存在大量退化解, 单纯形顶点解能提供干净的 0 通量,
    与参考实现的行为一致。
    """
    from optlang.glpk_interface import (Model as GLPKModel, Variable,
                                        Constraint, Objective)
    import optlang

    n = prepared.S.shape[1]
    lb = prepared.lb.copy()
    for ji, fmax in zip(prepared.obj_fwd, prepared.max_objective):
        lb[ji] = obj_frac * fmax
    glp = GLPKModel()
    vars_ = [Variable("x%d" % j, lb=float(lb[j]), ub=float(prepared.ub[j]))
             for j in range(n)]
    glp.add(vars_)
    Srow = [dict() for _ in range(prepared.S.shape[0])]
    Scoo = prepared.S.tocoo()
    for r_, c_, v_ in zip(Scoo.row, Scoo.col, Scoo.data):
        Srow[r_][c_] = float(v_)
    cons = []
    for d in Srow:
        if d:
            e = optlang.symbolics.add(
                [v_ * vars_[j_] for j_, v_ in d.items()])
        else:
            e = 0
        cons.append(Constraint(e, lb=0, ub=0))
    glp.add(cons)
    glp.objective = Objective(0, direction="min")
    prepared.glpk = glp
    prepared.glpk_vars = vars_


def _solve_penalty_lp(prepared: GimmeLP, penalty: np.ndarray,
                      solver: str, obj_frac: float = 0.9
                      ) -> tuple[int, np.ndarray]:
    """求解 GIMME 惩罚 LP, 返回 (status, x)。status=0 表示最优。"""
    n = len(penalty)
    if solver == "glpk":
        import optlang
        from optlang.glpk_interface import Objective
        if prepared.glpk is None:
            _build_glpk_problem(prepared, obj_frac)
        glp = prepared.glpk
        vars_ = prepared.glpk_vars
        terms = [float(penalty[j]) * vars_[j]
                 for j in range(n) if penalty[j] != 0]
        glp.objective = Objective(
            optlang.symbolics.add(terms) if terms else 0, direction="min")
        status = glp.optimize()
        if status != "optimal":
            return 1, np.zeros(n)
        x = np.array([glp.primal_values["x%d" % j] for j in range(n)])
        return 0, x
    else:
        lb = prepared.lb.copy()
        for ji, fmax in zip(prepared.obj_fwd, prepared.max_objective):
            lb[ji] = obj_frac * fmax
        res = linprog(penalty, A_eq=prepared.S,
                      b_eq=np.zeros(prepared.S.shape[0]),
                      bounds=list(zip(lb, prepared.ub)),
                      method="highs")
        if res.status != 0:
            return 1, np.zeros(n)
        return 0, res.x


def gimme(prepared: GimmeLP,
          model: cobra.Model,
          expression_rxns: np.ndarray,
          threshold: float = 20,
          obj_frac: float = 0.9,
          solver: str = "glpk") -> cobra.Model:
    """GIMME 主算法 (Becker & Palsson 2008)。

    Parameters
    ----------
    prepared : GimmeLP
        :func:`prepare_gimme` 预处理结果。
    model : cobra.Model
        原始模型 (用于生成删减后的组织特异模型)。
    expression_rxns : np.ndarray
        与 model.reactions 一一对应的反应表达值, 无数据为 -1。
    threshold : float
        表达阈值 (原脚本 options.threshold = 20)。
    obj_frac : float
        目标函数需保持的最小分数 (COBRA 默认 0.9, 原脚本未修改)。
    solver : str
        惩罚 LP 的求解器: "glpk" (默认, 与参考实现的
        changeCobraSolver('glpk') 一致, 单纯形顶点解) 或 "highs"。

    Returns
    -------
    cobra.Model
        上下文特异单细胞模型 (已删除无活性反应及孤儿代谢物)。
    """
    expression_rxns = np.asarray(expression_rxns, dtype=float)
    n_rxns = len(model.reactions)
    expr_irrev = expression_rxns[prepared.irrev2rev]
    n_irrev = len(expr_irrev)

    # 惩罚向量: -1 < expr < threshold -> threshold - expr
    penalty = np.zeros(n_irrev)
    mask = (expr_irrev > -1) & (expr_irrev < threshold)
    penalty[mask] = threshold - expr_irrev[mask]

    status, fluxes = _solve_penalty_lp(prepared, penalty, solver, obj_frac)
    if status != 0:
        # 与 MATLAB 版一致: 求解失败时通量取 0
        fluxes = np.zeros(n_irrev)

    # 反应活性判定
    activity_irrev = np.zeros(n_irrev, dtype=int)
    activity_irrev[(expr_irrev > threshold) | (expr_irrev == -1)] = 1
    low = (expr_irrev <= threshold) & (expr_irrev != -1)
    activity_irrev[low & (fluxes > 0)] = 2

    # 映射回可逆模型
    activity = np.zeros(n_rxns, dtype=int)
    np.maximum.at(activity, prepared.irrev2rev, activity_irrev)

    # 删除无活性反应
    remove_ids = [prepared.rxn_ids[i] for i in range(n_rxns) if activity[i] == 0]
    tissue_model = model.copy()
    tissue_model.remove_reactions(remove_ids)

    # 删除不再参与任何反应的代谢物 (与原 MATLAB 流程写出的 SBML 一致)。
    # 注意: 直接在深拷贝模型上连续调用 remove_metabolites 会污染
    # cobrapy/optlang 的内部状态 (导致后续模型求解异常), 因此先经 SBML
    # 序列化得到干净的新对象再修剪。
    tmp = tempfile.NamedTemporaryFile(suffix=".xml", delete=False)
    tmp.close()
    try:
        cobra.io.write_sbml_model(tissue_model, tmp.name)
        fresh = cobra.io.read_sbml_model(tmp.name)
    finally:
        os.unlink(tmp.name)
    orphan_mets = [m for m in fresh.metabolites if len(m.reactions) == 0]
    fresh.remove_metabolites(orphan_mets)
    return fresh


def build_single_cell_model(prepared: GimmeLP,
                            model: cobra.Model,
                            parsed_gpr: list[list[list[str]]],
                            expr_genes,
                            cell_expr_values: np.ndarray,
                            threshold: float = 20,
                            obj_frac: float = 0.9,
                            solver: str = "glpk") -> cobra.Model:
    """单细胞模型构建全流程 (对应 singlecellmodel.m 的循环体)。

    findUsedGenesLevels -> selectGeneFromGPR -> GIMME
    """
    gene_id, gene_expr = find_used_genes_levels(
        [g.id for g in model.genes], expr_genes, cell_expr_values)
    expression_rxns = select_gene_from_gpr(
        model, gene_id, gene_expr, parsed_gpr, min_sum=False)
    return gimme(prepared, model, expression_rxns, threshold, obj_frac, solver)
