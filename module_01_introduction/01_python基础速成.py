# -*- coding: utf-8 -*-
# =====================================================================
# 01_python基础速成.py
# ---------------------------------------------------------------------
# 模块一《入门知识介绍》（对应课程3）第 1 个上机脚本。
#
# 目的：用最贴近生物信息学的例子，快速复习/学习本课程需要的 Python 基础：
#       变量与数值运算、列表(list)、字典(dict)、循环、函数、numpy 基础。
#
# 用法：python 01_python基础速成.py
# 依赖：仅 numpy（pip install numpy）
# =====================================================================
"""Quick refresher of the Python basics needed for this course (bioinformatics examples)."""

# ----------------------------------------------------------------------
# 1. 变量与基本数据类型
# ----------------------------------------------------------------------
# 变量无需声明类型，直接赋值；# 后为行尾注释，解释本行含义
gene_name = "b0351"          # 字符串：一个大肠杆菌基因的 BiGG 编号（pfkA）
expression_level = 245.6     # 浮点数：假设的转录组表达量（FPKM）
replicate_count = 3          # 整数：生物学重复数
is_essential = True          # 布尔值：该基因是否为必需基因（预测）

# type() 查看变量类型；print() 打印到屏幕，方便课堂逐行观察
print("Gene ID:", gene_name, "| type:", type(gene_name).__name__)
print("Expression level:", expression_level, "| type:", type(expression_level).__name__)
print("Replicate count:", replicate_count, "| type:", type(replicate_count).__name__)
print("Is essential:", is_essential, "| type:", type(is_essential).__name__)

# 字符串的常用操作：切片与格式化输出（f-string）
print("First 4 characters of the gene ID:", gene_name[:4])            # 切片取前 4 个字符
print(f"Gene {gene_name} has expression level {expression_level:.1f} FPKM")  # f-string 格式化，保留 1 位小数

# ----------------------------------------------------------------------
# 2. 列表 list：有序的集合，常用于存放基因列表、反应列表等
# ----------------------------------------------------------------------
# 糖酵解途径中的几个关键基因（大肠杆菌 BiGG 编号）
glycolysis_genes = ["b0351", "b3916", "b1676", "b1852"]  # pfkA, pfkB, fbaA, eno 等
print("\nGlycolysis gene list:", glycolysis_genes)
print("List length (number of genes):", len(glycolysis_genes))  # len() 求元素个数
print("First gene:", glycolysis_genes[0])                       # 索引从 0 开始
print("Last gene:", glycolysis_genes[-1])                       # 负索引 -1 表示倒数第一个

glycolysis_genes.append("b2925")     # append() 在列表末尾添加一个基因（pgk）
print("List after append:", glycolysis_genes)

# 列表推导式：把每个基因名前面加上 "Eco_" 前缀，一步生成新列表
prefixed_genes = ["Eco_" + g for g in glycolysis_genes]
print("List with prefix added:", prefixed_genes)

# 切片取子列表：前 3 个基因
print("First 3 genes (slice):", glycolysis_genes[:3])

# ----------------------------------------------------------------------
# 3. 字典 dict：键-值对，常用于“基因 -> 表达量”这类映射
# ----------------------------------------------------------------------
# 构造一个“基因表达量字典”：键是基因编号，值是表达量
expression_dict = {
    "b0351": 245.6,   # pfkA
    "b3916": 12.3,    # pfkB
    "b1676": 560.1,   # fbaA
    "b1852": 430.8,   # eno
}
print("\nExpression dict:", expression_dict)
print("Expression of b0351:", expression_dict["b0351"])         # 通过键取值
print("All genes in dict:", list(expression_dict.keys()))       # keys() 取所有键
print("All expression values:", list(expression_dict.values())) # values() 取所有值

# get() 方法在键不存在时返回默认值，避免 KeyError 报错
print("Expression of a missing gene (default 0):", expression_dict.get("b0000", 0))

# ----------------------------------------------------------------------
# 4. 循环 for：批量处理生物数据
# ----------------------------------------------------------------------
print("\nPrint expression per gene (for loop over dict):")
threshold = 100.0                                # 设定“高表达”阈值
high_expression_genes = []                       # 空列表，收集高表达基因
for gene, expr in expression_dict.items():       # items() 同时取键和值
    print(f"  gene {gene}: {expr} FPKM")
    if expr > threshold:                         # if 条件判断
        high_expression_genes.append(gene)       # 满足条件则加入列表
print("High-expression genes (> 100 FPKM):", high_expression_genes)

# ----------------------------------------------------------------------
# 5. 函数 def：把一段可复用的分析逻辑封装起来
# ----------------------------------------------------------------------
def compute_mean_expression(expr_dict):
    """Return the mean expression value over all genes in ``expr_dict``.

    Parameters
    ----------
    expr_dict : dict
        Mapping of gene ID (str) to expression level (float).

    Returns
    -------
    float
        Mean expression level; 0.0 if the dict is empty.
    """
    if len(expr_dict) == 0:            # 防御空字典，避免除以 0
        return 0.0
    return sum(expr_dict.values()) / len(expr_dict)   # 总和 / 个数


def find_top_gene(expr_dict):
    """Return the gene with the highest expression level.

    Parameters
    ----------
    expr_dict : dict
        Mapping of gene ID to expression level.

    Returns
    -------
    tuple (str, float)
        (top gene ID, its expression level).
    """
    top_gene = max(expr_dict, key=expr_dict.get)      # max 按字典的值比较
    return top_gene, expr_dict[top_gene]


mean_expr = compute_mean_expression(expression_dict)   # 调用自定义函数
top_gene, top_value = find_top_gene(expression_dict)
print(f"\nMean expression: {mean_expr:.2f} FPKM")
print(f"Top expressed gene: {top_gene} ({top_value} FPKM)")

# ----------------------------------------------------------------------
# 6. numpy 基础：数值计算的利器，FBA 通量结果就是 numpy 数组
# ----------------------------------------------------------------------
import numpy as np                                   # 惯例把 numpy 缩写为 np

# 用 numpy 数组存放 5 个样本的葡萄糖摄取速率（mmol/gDW/h，模拟 FBA 结果）
uptake_rates = np.array([10.0, 8.5, 6.0, 4.2, 2.0])
print("\nGlucose uptake rates array:", uptake_rates)
print("Array shape:", uptake_rates.shape)            # shape 返回数组维度
print("Mean:", np.mean(uptake_rates))                # 向量化统计函数
print("Standard deviation:", np.std(uptake_rates))
print("Max / Min:", np.max(uptake_rates), np.min(uptake_rates))

# 向量化运算：不用写循环，整组数据同时乘一个系数（例如单位换算）
growth_rates = 0.087 * uptake_rates                  # 假设生长速率与摄取速率成正比
print("Converted growth rates:", np.round(growth_rates, 3))  # np.round 保留 3 位小数

# 布尔索引：筛选出摄取速率大于 5 的样本
fast_samples = uptake_rates[uptake_rates > 5]
print("Samples with uptake > 5:", fast_samples)

# 生成等差序列：后续模块做“表型相平面图”扫描摄取速率时会用到
scan_points = np.linspace(0, 10, 6)                  # 0 到 10 均匀取 6 个点
print("Evenly spaced scan points:", scan_points)

print("\n=== Script 01_python基础速成.py finished ===")