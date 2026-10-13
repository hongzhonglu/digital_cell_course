# 模块01：入门知识与 COBRApy 基础

**对应课次：课程 3《入门知识介绍》**

## 模块目标
1. 以生物数据为例复习 Python 基础（变量、列表、字典、循环、函数、numpy）。
2. 认识 COBRApy 的核心对象 `Model / Reaction / Metabolite / Gene`，学会查看反应式、通量上下界与 GPR 规则。
3. 运行一次 FBA 并解读目标值与关键通量。
4. 掌握 GEMs 的导入导出（SBML / JSON / CSV 表格）。

## 代码规范说明（第二轮完善）
- 变量名/函数名为英文 snake_case；print 输出与函数 docstring 为英文；解释性注释保持中文；每个脚本开头有中文块注释说明目的与用法。
- 每个脚本均可独立 `python xxx.py` 运行，关键中间结果均有 print。

## 依赖安装
```bash
pip install cobra numpy pandas
```

## 文件清单与运行顺序
| 顺序 | 文件 | 内容 |
|---|---|---|
| 1 | `01_python基础速成.py` | Python 基础速成（结合基因/表达量数据） |
| 2 | `02_cobrapy框架入门.py` | 遍历 e_coli_core：反应/代谢物/基因/GPR、`reactions.query` 按名称搜索、单基因敲除演示（`gene.knock_out()`）、FBA 结果筛选（`solution.to_frame()`）、`model.summary()` |
| 3 | `03_GEMs导入与导出.py` | SBML/JSON 读写（并说明 .mat 格式的读写函数），模型转 pandas 表格并导出 CSV |
| — | `model_loader.py` | 模型加载辅助：优先在线 `load_model("e_coli_core")`，离线自动回退本地缓存 |
| — | `e_coli_core.xml.gz` | e_coli_core 模型本地缓存（BiGG 官方文件） |

运行方式（按顺序执行，均支持逐行调试）：
```bash
python 01_python基础速成.py
python 02_cobrapy框架入门.py
python 03_GEMs导入与导出.py
```

## 预期输出（关键数字）
- 模型规模：**95 反应、72 代谢物、137 基因**。
- FBA 最优解：**最大生长速率 ≈ 0.8739 1/h**，葡萄糖摄取 -10.0、耗氧 -21.8、CO₂ 释放 22.8 mmol/gDW/h。
- 敲除演示：敲除 gapA（b1779）后生长速率降为 0（GAPD 为糖酵解必需反应）。
- `03_GEMs导入与导出.py` 在 `output/` 目录生成：`e_coli_core.xml`（SBML）、`e_coli_core.json`、`e_coli_core_反应表.csv`（95 行，英文列名）、`e_coli_core_代谢物表.csv`（72 行）、`e_coli_core_FBA通量.csv`。

## 说明
- 课堂若无法联网，`model_loader.py` 会自动使用同目录的 `e_coli_core.xml.gz`，结果与在线下载完全一致。
- SBML 为 GEM 交换的标准 XML 格式（level 3 + fbc 插件）；JSON 便于程序读取且为 Escher 可视化所用格式；MATLAB 用户可用 `load_matlab_model` / `save_matlab_model` 读写 .mat 文件。
