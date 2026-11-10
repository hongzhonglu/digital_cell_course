# -*- coding: utf-8 -*-
"""
03_模型分析_底物利用.py
=======================
课程7&8：对"自动构建产物模型"做底物（碳源）利用能力分析。

Substrate (carbon source) utilization analysis of an automatically
reconstructed GEM — the typical downstream analysis after gapseq/AlphaGEM.

Model source (two options, handled automatically):
  A) download iAF692 (Methylobacterium extorquens AM1, a representative
     "non-model" small model) from the BiGG database;
  B) if the download fails (no network), fall back to COBRApy's built-in
     e_coli_core as the demo model.

Analysis idea: set each carbon-containing exchange reaction as the "sole
carbon source" in turn (close all other carbon exchanges, open this one),
run FBA, and call the substrate usable if growth rate exceeds a threshold.
Results are saved to substrate_utilization.csv.
"""

import os                      # 路径与文件存在性判断
import re                      # 分子式中"碳元素"的正则识别
import pandas as pd            # 结果表格
from cobra.io import load_json_model  # 读 json 模型

# ---------------- 1. Load the model (try A then B) ----------------
__file__='module_05_GEM_reconstruction_for_non_model_organism'
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))      # 脚本所在目录
MODEL_JSON = os.path.join(SCRIPT_DIR, "iAF692.json")         # 本地下载的 iAF692
BIGG_URL = "http://bigg.ucsd.edu/static/models/iAF692.json"  # BiGG 下载地址

model = None
if not os.path.exists(MODEL_JSON):
    # 尝试方案 A：从 BiGG 下载 iAF692 并缓存到本地
    try:
        import requests  # 发起 HTTP 请求
        print("[Plan A] Downloading iAF692 from BiGG ...")
        response = requests.get(BIGG_URL, timeout=60)
        response.raise_for_status()
        with open(MODEL_JSON, "wb") as f:
            f.write(response.content)
        print("         Download succeeded, cached as iAF692.json")
    except Exception as exc:  # 网络失败等
        print(f"[Plan A] Download failed ({exc}); falling back to Plan B.")

if os.path.exists(MODEL_JSON):
    model = load_json_model(MODEL_JSON)      # 方案 A：读入 iAF692
    print(f"[Model] iAF692: {len(model.reactions)} reactions / {len(model.genes)} genes")
else:
    # 方案 B：用 e_coli_core 充当"自动构建模型"演示算法
    from cobra.io import load_model
    model = load_model("textbook")
    print(f"[Model] Fallback to e_coli_core: {len(model.reactions)} reactions")

# 兼容求解器：优先 glpk（免费），避免依赖 cplex
try:
    model.solver = "glpk"
except Exception:
    pass

# ---------------- 2. Find candidate carbon-source exchange reactions ----------------
# 交换反应：model.exchanges；候选碳源：交换代谢物分子式含碳元素。
# 注意：分子式中 Co（钴）、Cu（铜）、Cl（氯）、Ca（钙）等也含字母 "C"，
# 因此必须用正则识别"碳元素 C"（大写 C 后面不跟小写字母），否则会误判。
def contains_carbon(formula):
    """Check whether a molecular formula contains carbon (excluding Co/Cu/Cl/Ca etc.)."""
    if not formula:
        return False
    return re.search(r"C(?![a-z])", formula) is not None

# 明确排除的无机/非营养性含碳交换：CO2、碳酸氢盐等
EXCLUDE_EXCHANGES = {"EX_co2_e", "EX_hco3_e", "EX_co_e"}

carbon_exchanges = []   # 候选碳源交换反应列表
for rxn in model.exchanges:
    if rxn.id in EXCLUDE_EXCHANGES:
        continue
    # 取该交换反应涉及的代谢物，检查是否真正含碳元素
    mets = list(rxn.metabolites.keys())
    if not mets:
        continue
    met = mets[0]
    if contains_carbon(met.formula or ""):
        carbon_exchanges.append(rxn)

print(f"[Candidate carbon sources] {len(carbon_exchanges)} carbon-containing exchange reactions")

# ---------------- 2b. Identify "essential supplements" (vitamins/amino acids for auxotrophy) ----------------
# iAF692 等模型存在营养缺陷：如烟酸(nac)、对氨基苯甲酸(4abz)，
# 它们虽含碳，却是"培养基补充成分"——缺少即无法生长，但并不充当主要碳源。
# 判据（两轮）：
#   第一轮：缺省培养基下单独关闭该交换即不生长 → "疑似必需成分"；
#   第二轮：疑似成分若能"单独作为唯一碳源支持生长"，说明它其实是主要碳源
#           （如 iAF692 的甲醇——该菌为甲基营养菌），否则判定为补充成分。
SUPPLEMENT_UPTAKE = 1.0   # 给补充成分设置的小摄取上限（仅满足生物量微量需求）

# 第一轮：找出关闭后致死（疑似必需）的成分
suspect_supplements = []
for rxn in carbon_exchanges:
    with model:
        rxn.lower_bound = 0.0                      # 单独关闭该成分
        growth = model.optimize().objective_value
        if growth is None or growth < 0.001:
            suspect_supplements.append(rxn)

# 第二轮：疑似成分逐一测试"能否独立供碳"，区分真碳源 vs 补充成分
essential_supplements = []
for rxn in suspect_supplements:
    with model:
        for other in carbon_exchanges:             # 先关闭全部候选碳源
            other.lower_bound = 0.0
        for other in suspect_supplements:          # 其他疑似成分小量开放
            if other is not rxn:
                other.lower_bound = -SUPPLEMENT_UPTAKE
        rxn.lower_bound = -10.0                    # 当前成分作为唯一碳源
        growth = model.optimize().objective_value
        if growth is None or growth < 0.001:
            essential_supplements.append(rxn)      # 不能独立供碳 → 补充成分

print(f"[Essential supplements] {[r.id for r in essential_supplements]}")
print(f"[Sole-carbon-capable suspects] "
      f"{[r.id for r in suspect_supplements if r not in essential_supplements]}")
# 从候选碳源中剔除补充成分
carbon_exchanges = [r for r in carbon_exchanges if r not in essential_supplements]

# ---------------- 3. Set each candidate as sole carbon source and run FBA ----------------
GROWTH_THRESHOLD = 0.001   # 生长率阈值（1/h），高于此判定为"可利用"
UPTAKE_BOUND = 10.0        # 设定的摄取速率上限 (mmol/gDW/h)

results = []               # 收集 (碳源, 生长率, 是否可利用)

with model:                # with 上下文：退出时自动恢复模型所有修改
    # 关闭所有候选碳源的摄取（lower bound 置 0），构建"无碳"基底
    for rxn in carbon_exchanges:
        rxn.lower_bound = 0.0
    # 必需补充成分保持小量开放（满足营养缺陷，但不作为主要碳源）
    for rxn in essential_supplements:
        rxn.lower_bound = -SUPPLEMENT_UPTAKE

    for rxn in carbon_exchanges:
        rxn.lower_bound = -UPTAKE_BOUND    # 打开当前碳源摄取
        solution = model.optimize()        # FBA
        growth = solution.objective_value if solution.status == "optimal" else 0.0
        growth = growth if growth is not None else 0.0
        usable = growth > GROWTH_THRESHOLD
        results.append({
            "exchange_id": rxn.id,                          # 交换反应编号
            "metabolite": list(rxn.metabolites)[0].name,    # 代谢物名称
            "growth_rate": round(growth, 4),                # 生长率
            "usable": usable,                               # 是否可利用
        })
        rxn.lower_bound = 0.0              # 关闭，测试下一个碳源

# ---------------- 4. Save and print results ----------------
result_df = pd.DataFrame(results).sort_values("growth_rate", ascending=False)
out_csv = os.path.join(SCRIPT_DIR, "substrate_utilization.csv")
result_df.to_csv(out_csv, index=False, encoding="utf-8-sig")   # utf-8-sig 便于 Excel 打开中文

usable_df = result_df[result_df["usable"]]
print(f"\n[Result] usable substrates: {len(usable_df)} / {len(result_df)}:")
print(usable_df[["exchange_id", "metabolite", "growth_rate"]].to_string(index=False))
print(f"\nSaved: {out_csv}")
