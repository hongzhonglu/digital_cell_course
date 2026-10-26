---
title: "在 PyCharm 中配置 Python 与安装软件包指南"
author: "AI Assistant"
date: "`r Sys.Date()`"
output: html_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

# 在 PyCharm 中安装 Python、Numpy、Cobrapy 等软件包

在 PyCharm 中安装 Python 解释器以及第三方软件包（如 numpy、cobrapy 等）非常直观。你可以通过 PyCharm 内置的图形化界面（Python Packages 或 Settings）来完成，也可以直接在终端使用 pip 安装。

以下是具体的步骤指南：

---

## 第一步：检查或安装 Python 解释器

在安装第三方包之前，确保 PyCharm 已经配置好了 Python 解释器（即 Python 环境）。

1. 打开 PyCharm，进入设置界面：
   * **Windows/Linux**: 点击菜单栏的 `File` -> `Settings`
   * **macOS**: 点击菜单栏的 `PyCharm` -> `Settings` (或 `Preferences`)
2. 在左侧菜单中导航到：**`Project: 你的项目名`** -> **`Python Interpreter`**。
3. 如果右侧显示了 Python 版本（例如 Python 3.10 或 3.11），说明已经配置好，可以直接进行第二步。
4. 如果显示为空或需要新建环境：
   * 点击右上角的齿轮图标 (`⚙️`) -> **`Add Interpreter`** -> **`Add Local Interpreter...`**
   * 选择 **`System Interpreter`**（使用你电脑上已经安装的 Python）或者新建一个 Virtualenv 虚拟环境，点击 **`OK`** 保存。

> **如何验证**：在 PyCharm 下方的 **Terminal（终端）** 标签页中输入 `python --version` 或 `python3 --version`，如果能正确输出版本号，说明解释器配置成功。

---

## 第二步：安装 Numpy、Cobrapy 等软件包

你可以选择以下任意一种方法进行安装：

### 方法一：通过图形化界面安装（最推荐、简单直观）

1. 在刚才的 **`Settings`** 窗口中，继续停留在 **`Python Interpreter`** 页面。
2. 你会看到当前环境中已安装的包列表。点击列表上方的 **`+` (加号)** 图标。
3. 在弹出的搜索框中，输入你想安装的包名（例如：`numpy`）。
4. 选中搜索结果中的 `numpy`，然后点击右下角的 **`Install Package`** 按钮。
5. 等待下方日志显示安装成功后，用同样的步骤搜索并安装 `cobra`（注意：cobrapy 在 PyPI 上的包名通常是 `cobra`）。

> **如何验证**：安装完成后，返回 Python Interpreter 列表，检查列表中是否出现了 `numpy` 和 `cobra`。

### 方法二：通过 PyCharm 底部终端（Terminal）使用 pip 安装

如果你习惯用命令行，或者图形化安装遇到网络超时，可以直接用命令行：

1. 点击 PyCharm 窗口最下方的 **`Terminal`** 标签页，打开内置终端。
2. 输入以下命令并回车安装所需的包：
   ```bash
   pip install numpy cobra
   ```
   *(注：生物信息学中的 cobrapy 对应的安装包名称是 `cobra`，而不是 `cobrapy`)*
3. 如果下载速度慢，可以使用国内的清华镜像源加速：
   ```bash
   pip install numpy cobra -i https://pypi.tuna.tsinghua.edu.cn/simple
   ```

> **如何验证**：在终端输入 `pip list`，如果在打印出的列表中看到了 `numpy` 和 `cobra`，说明安装成功。或者在代码中输入 `import numpy as np` 和 `import cobra`，如果没有红色波浪线报错，则代表环境配置完全正常。