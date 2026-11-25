# -*- coding: utf-8 -*-
"""
step6_classification.py
=======================
对应原仓库 single_cell/code/classification.py。

以 80 个胁迫细胞 (label 0) 与 83 个对照细胞 (label 1) 为样本, 用随机森林
分类器比较六种特征表示对细胞状态的区分能力:
  - rf_umap_flux / rf_umap_trans : UMAP (5 维) 降维后的通量 / 转录组
  - rf_flux / rf_trans           : 原始通量 / 转录组
  - rf_kpca_flux / rf_kpca_trans : KernelPCA 降维后的通量 / 转录组

每种表示用 5 折交叉验证做网格搜索, 汇总结果输出 ../output/result.xlsx,
各分项结果分别输出 rf_*.xlsx。
"""

import os
import time

import pandas as pd
import umap
from sklearn.decomposition import KernelPCA
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import accuracy_score
from sklearn.model_selection import GridSearchCV, train_test_split
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def _count_confusion(train_label, pre_train, test_label, pre_test):
    """统计 预测/真实 组合计数 (与原脚本一致):
    ss=胁迫->胁迫, s_ns=胁迫->对照, ns_s=对照->胁迫, ns_ns=对照->对照"""
    ss = s_ns = ns_s = ns_ns = 0
    for true, pred in zip(list(train_label) + list(test_label),
                          list(pre_train) + list(pre_test)):
        if true == 0 and pred == 0:
            ss += 1
        elif true == 1 and pred == 1:
            ns_ns += 1
        elif true == 0 and pred == 1:
            s_ns += 1
        else:
            ns_s += 1
    return ss, s_ns, ns_ns, ns_s


def _make_result(method, datatype, train_acc, test_acc, counts):
    ss, s_ns, ns_ns, ns_s = counts
    result = pd.DataFrame(index=["value_{}_{}".format(method, datatype)],
                          columns=["test_acc", "train_acc", "ss", "s_ns", "ns_ns", "ns_s"])
    result.iloc[0] = [test_acc, train_acc, ss, s_ns, ns_ns, ns_s]
    return result


def rf_kpca(data, label, method, datatype):
    train_data, test_data, train_label, test_label = train_test_split(
        data, label, random_state=1, train_size=0.7, test_size=0.3)

    pipeline = Pipeline([("kpca", KernelPCA()),
                         ("rf", RandomForestClassifier(random_state=42))])
    param_grid = {
        "kpca__n_components": [5, 10, 100, 500],
        "kpca__kernel": ["linear", "poly", "rbf"],
        "rf__n_estimators": [10, 100, 500],
        "rf__max_depth": [None, 10, 20],
        "rf__min_samples_split": [2, 5, 10],
        "rf__min_samples_leaf": [1, 2, 4],
    }
    grid_search = GridSearchCV(pipeline, param_grid, cv=5, n_jobs=-1)
    grid_search.fit(train_data, train_label)
    best_model = grid_search.best_estimator_
    pre_train = best_model.predict(train_data)
    pre_test = best_model.predict(test_data)
    train_acc = accuracy_score(train_label, pre_train)
    test_acc = accuracy_score(test_label, pre_test)
    print("  test_acc: %.4f, train_acc: %.4f" % (test_acc, train_acc))
    print("  Best parameters:", grid_search.best_params_)
    print("  Best CV score: {:.2f}".format(grid_search.best_score_))
    counts = _count_confusion(train_label, pre_train, test_label, pre_test)
    return _make_result(method, datatype, train_acc, test_acc, counts)


def rf(data, label, method, datatype):
    train_data, test_data, train_label, test_label = train_test_split(
        data, label, random_state=1, train_size=0.7, test_size=0.3)
    param_grid = {
        "n_estimators": [10, 100, 500],
        "max_depth": [None, 10, 20],
        "min_samples_split": [2, 5, 10],
        "min_samples_leaf": [1, 2, 4],
    }
    grid_search = GridSearchCV(RandomForestClassifier(random_state=42),
                               param_grid, cv=5, n_jobs=-1)
    grid_search.fit(train_data, train_label)
    print("  Best parameters:", grid_search.best_params_)
    best_rf = RandomForestClassifier(**grid_search.best_params_, random_state=42)
    best_rf.fit(train_data, train_label)
    pre_train = best_rf.predict(train_data)
    pre_test = best_rf.predict(test_data)
    train_acc = accuracy_score(train_label, pre_train)
    test_acc = accuracy_score(test_label, pre_test)
    print("  test_acc: %.4f, train_acc: %.4f" % (test_acc, train_acc))
    counts = _count_confusion(train_label, pre_train, test_label, pre_test)
    return _make_result(method, datatype, train_acc, test_acc, counts)


def main():
    outdir = os.path.join(BASE, "output")

    data = pd.read_excel(os.path.join(outdir, "pfba_flux.xlsx"))
    data = data.iloc[:, 1:]
    data.loc[:, "metnum/20"] = data.loc[:, "metnum/20"] * 20
    data.loc[:, "rxnnum/20"] = data.loc[:, "rxnnum/20"] * 20

    label = [0] * 80 + [1] * 83  # 0 = 胁迫, 1 = 对照

    reducer = umap.UMAP(n_components=5)
    scaled_flux = StandardScaler().fit_transform(data)
    embedding_flux = reducer.fit_transform(scaled_flux)

    transdata = pd.read_excel(
        os.path.join(BASE, "data", "GSE102475_GASCH_NaCl-scRNAseq_NormData.xlsx"),
        index_col=0)
    transdata = transdata.T

    t0 = time.time()
    print("rf_umap flux:")
    rf_umap_flux = rf(embedding_flux, label, "umap", "flux")
    rf_umap_flux.to_excel(os.path.join(outdir, "rf_umap_flux.xlsx"))

    print("rf_umap trans:")
    scaled_trans = StandardScaler().fit_transform(transdata)
    embedding_trans = reducer.fit_transform(scaled_trans)
    rf_umap_trans = rf(embedding_trans, label, "umap", "trans")
    rf_umap_trans.to_excel(os.path.join(outdir, "rf_umap_trans.xlsx"))

    print("rf flux:")
    rf_flux = rf(data, label, "off", "flux")
    rf_flux.to_excel(os.path.join(outdir, "rf_flux.xlsx"))

    print("rf trans:")
    rf_trans = rf(transdata, label, "off", "trans")
    rf_trans.to_excel(os.path.join(outdir, "rf_trans.xlsx"))

    print("rf_kpca flux:")
    rf_kpca_flux = rf_kpca(data, label, "kpca", "flux")
    rf_kpca_flux.to_excel(os.path.join(outdir, "rf_kpca_flux.xlsx"))

    print("rf_kpca trans:")
    rf_kpca_trans = rf_kpca(transdata, label, "kpca", "trans")
    rf_kpca_trans.to_excel(os.path.join(outdir, "rf_kpca_trans.xlsx"))

    result = pd.concat([rf_kpca_flux, rf_kpca_trans, rf_umap_flux,
                        rf_umap_trans, rf_flux, rf_trans], axis=0)
    result.to_excel(os.path.join(outdir, "result.xlsx"))
    print("result.xlsx 已保存, 总耗时 %.0fs" % (time.time() - t0))
    print(result)


if __name__ == "__main__":
    main()
