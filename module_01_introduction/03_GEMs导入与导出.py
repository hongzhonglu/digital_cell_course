# -*- coding: utf-8 -*-
# =====================================================================
# 03_GEMs导入与导出.py
# ---------------------------------------------------------------------
# 模块一《入门知识介绍》（对应课程3）第 3 个上机脚本。
#
# 目的：掌握基因组尺度代谢模型（GEMs）的常见文件格式与导入导出：
#       1) SBML 格式：read_sbml_model / write_sbml_model
#       2) JSON 格式：load_json_model / save_json_model
#       3) 模型内容转 pandas.DataFrame 并导出 CSV（反应表、代谢物表）
#
# 用法：python 03_GEMs导入与导出.py
# 依赖：cobra、pandas（pip install cobra pandas）
# =====================================================================
"""Read/write GEMs in SBML and JSON formats, and export model content to CSV tables."""

import os

import pandas as pd

from cobra.io import (load_json_model, read_sbml_model, save_json_model,
                      write_sbml_model)
from model_loader import load_ecoli_core

# 输出目录：与本脚本同级的 output/ 文件夹，先确保其存在
OUTPUT_DIR = os.path.join(os.getcwd(), "module_01_introduction/output")
os.makedirs(OUTPUT_DIR, exist_ok=True)
print("Output directory:", OUTPUT_DIR)





# ----------------------------------------------------------------------
# 1. 加载模型（在线或本地缓存，见 model_loader.py）
# ----------------------------------------------------------------------
model = load_ecoli_core()
print("\nModel loaded:", model.id,
      f"({len(model.reactions)} reactions, {len(model.metabolites)} metabolites, "
      f"{len(model.genes)} genes)")

# ----------------------------------------------------------------------
# 2. SBML 格式导出与再导入
#    SBML (Systems Biology Markup Language) 是系统生物学模型的标准 XML 格式
#    （GEM 社区常用 sbml 级别 3 + fbc 插件存储上下界、GPR 等信息）
# ----------------------------------------------------------------------
sbml_path = os.path.join(OUTPUT_DIR, "e_coli_core.xml")   # SBML 输出路径
write_sbml_model(model, sbml_path)                         # 把模型写成 SBML 文件
print(f"\nSBML file written: {sbml_path} ({os.path.getsize(sbml_path)} bytes)")

model_from_sbml = read_sbml_model(sbml_path)               # 再从 SBML 读回来
print("Re-imported from SBML:", model_from_sbml.id,
      f"({len(model_from_sbml.reactions)} reactions)")     # 反应数应与原模型一致

# ----------------------------------------------------------------------
# 3. JSON 格式导出与再导入
#    JSON 文件便于程序读取，也是 Escher 可视化工具使用的格式；
#    此外还有 MATLAB 的 .mat 格式，可用 load_matlab_model / save_matlab_model 读写
# ----------------------------------------------------------------------
json_path = os.path.join(OUTPUT_DIR, "e_coli_core.json")
save_json_model(model, json_path)                          # 把模型写成 JSON 文件
print(f"\nJSON file written: {json_path} ({os.path.getsize(json_path)} bytes)")

model_from_json = load_json_model(json_path)               # 再从 JSON 读回来
print("Re-imported from JSON:", model_from_json.id,
      f"({len(model_from_json.reactions)} reactions)")

# 验证：重新导入的模型求解结果应与原模型一致
growth = model_from_json.optimize().objective_value
print(f"Max growth rate of the re-imported model: {growth:.4f} 1/h")

# ----------------------------------------------------------------------
# 4. 模型内容转 pandas 表格
# ----------------------------------------------------------------------
def reactions_to_dataframe(model):
    """Collect all reactions of ``model`` into a pandas.DataFrame.

    Parameters
    ----------
    model : cobra.Model
        The metabolic model.

    Returns
    -------
    pandas.DataFrame
        One row per reaction with columns: reaction id, name, equation,
        lower/upper bounds and GPR rule.
    """
    rows = []
    for rxn in model.reactions:                    # 遍历所有反应
        rows.append({
            "reaction_id": rxn.id,               # 反应编号
            "reaction_name": rxn.name,           # 反应名称
            "equation": rxn.reaction,            # 化学计量反应式
            "lower_bound": rxn.lower_bound,      # 通量下界
            "upper_bound": rxn.upper_bound,      # 通量上界
            "gpr": rxn.gene_reaction_rule,       # 基因-蛋白-反应关联规则
        })
    return pd.DataFrame(rows)                      # 列表 of dict -> DataFrame


def metabolites_to_dataframe(model):
    """Collect all metabolites of ``model`` into a pandas.DataFrame.

    Parameters
    ----------
    model : cobra.Model
        The metabolic model.

    Returns
    -------
    pandas.DataFrame
        One row per metabolite with columns: id, name, formula, compartment.
    """
    rows = []
    for met in model.metabolites:                  # 遍历所有代谢物
        rows.append({
            "metabolite_id": met.id,             # 代谢物编号（含区室后缀）
            "metabolite_name": met.name,         # 代谢物名称
            "formula": met.formula,              # 化学式
            "compartment": met.compartment,      # 区室（c/e/p）
        })
    return pd.DataFrame(rows)


df_reactions = reactions_to_dataframe(model)
df_metabolites = metabolites_to_dataframe(model)

print("\nFirst 5 rows of the reaction table:")
print(df_reactions.head().to_string(index=False))     # to_string 去掉行索引更美观
print("\nFirst 5 rows of the metabolite table:")
print(df_metabolites.head().to_string(index=False))

# ----------------------------------------------------------------------
# 5. 表格导出为 CSV，可用 Excel / R / pandas 进一步分析
# ----------------------------------------------------------------------
reactions_csv = os.path.join(OUTPUT_DIR, "e_coli_core_反应表.csv")
metabolites_csv = os.path.join(OUTPUT_DIR, "e_coli_core_代谢物表.csv")
# encoding="utf-8-sig" 让 Excel 打开 UTF-8 文件不乱码
df_reactions.to_csv(reactions_csv, index=False, encoding="utf-8-sig")
df_metabolites.to_csv(metabolites_csv, index=False, encoding="utf-8-sig")
print(f"\nCSV exported: {reactions_csv} ({len(df_reactions)} rows)")
print(f"CSV exported: {metabolites_csv} ({len(df_metabolites)} rows)")

# 顺便把本次 FBA 通量也存成 CSV
flux_csv = os.path.join(OUTPUT_DIR, "e_coli_core_FBA通量.csv")
solution = model.optimize()
solution.fluxes.to_frame("flux").to_csv(flux_csv, encoding="utf-8-sig")
print(f"CSV exported: {flux_csv} ({len(solution.fluxes)} rows)")

print("\n=== Script 03_GEMs导入与导出.py finished ===")
