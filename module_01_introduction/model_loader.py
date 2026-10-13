# -*- coding: utf-8 -*-
# =====================================================================
# model_loader.py —— 模型加载辅助模块
# ---------------------------------------------------------------------
# 目的：提供 load_ecoli_core()：优先用 cobra.io.load_model("e_coli_core")
#       从 BiGG 数据库在线下载；若当前环境无网络，则自动回退到本文件夹内
#       预先下载好的本地缓存文件 e_coli_core.xml.gz，保证课堂演示离线可运行。
#
# 用法：from model_loader import load_ecoli_core
# =====================================================================
"""Helper that loads the e_coli_core model online, with an offline local-cache fallback."""

import gzip
import os

from cobra.io import read_sbml_model


def load_ecoli_core():
    """Load the E. coli core metabolic model e_coli_core.

    Returns
    -------
    cobra.Model
        The e_coli_core model (95 reactions, 72 metabolites, 137 genes).
    """
    try:
        # 标准用法：从 BiGG Models 在线数据库下载并解析模型
        from cobra.io import load_model
        model = load_model("e_coli_core")
        print("Loaded e_coli_core from the BiGG online database")
        return model
    except Exception as err:                       # 捕获网络异常等一切错误
        # 离线回退：读取与本脚本同目录的本地缓存文件
        local_gz = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "e_coli_core.xml.gz")
        print(f"Online download failed ({type(err).__name__}), "
              f"using local cache instead: {local_gz}")
        with gzip.open(local_gz, "rt", encoding="utf-8") as fh:  # 解压并按文本读取
            return read_sbml_model(fh)             # 解析 SBML 得到 Model 对象
