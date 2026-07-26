# Delta Dependency：跨实体瘤类型的旁系同源基因合成致死候选优先排序

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.21502030.svg)](https://doi.org/10.5281/zenodo.21502030)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

**[English](README.md)** | **中文**

---

## 概述

Delta Dependency（DD）是一个简单、可解释的发现阶段指标，用于从癌症依赖图谱（DepMap）CRISPR 筛选数据中优先排序基于旁系同源基因的合成致死（SL）候选靶点。DD 测量驱动基因突变型与野生型细胞系之间 Chronos 基因效应评分的偏移量，按癌症谱系分别计算。

**核心发现（详见论文全文）：**
- **金标准评估（证据分级、引用逐条核验）：** 完整 12 对策展集 AUROC = 0.682；去除 2 对 DepMap 来源标签后 0.728；仅用 DepMap 之前证据 0.774；方向严格重标后不变（0.682）。两条可评估的 Tier A 记录（方向性外部实验证据）均排在全部未标记对照之上（Tier A AUROC = 1.000）
- **DD vs. 已发表方法（背景参照，CV3 类框架）：** 无需训练 AUROC = 0.682；已发表 CV3 最佳 0.790（SLMGAE）——评估框架类似但不完全相同
- **头对头比较（同一 75 对测试集、6 个阳性）：** 可解释组合分（0.841）在留一交叉验证下超过全部四种多特征分类器（SVM-RBF 0.744、RF 0.617、SVM-Linear 0.217、LR 0.138；仅用 DD 0.551）——小样本结果附明确检验效能说明
- DD + ≥30% 序列一致性过滤 → AUROC = 1.000

## 方法原理

```
DD(D, P, c) = G̅(P, D-WT, c) − G̅(P, D-MUT, c)
```

其中 `G̅` 为 Chronos 基因效应评分均值，`D` 为驱动基因，`P` 为候选旁系同源基因，`c` 为癌种。正值 DD 表示旁系同源基因在突变型细胞中更为必需——符合旁系同源补偿机制。

## 目录结构

```
paralog_sl_predictor/
├── README.md / README_CN.md   # 英文/中文说明
├── LICENSE                    # MIT 许可证
├── requirements.txt           # Python 依赖
├── Dockerfile                 # 可复现环境
│
├── manuscript.tex/.pdf        # 主文稿（LaTeX + PDF）
├── supplementary.tex/.pdf     # 补充材料
├── cover_letter.md            # 投稿信
│
├── config.py                  # 驱动基因、已知SL对、参数配置
├── data_loader.py             # DepMap 数据加载工具
├── main.py                    # 主分析管线
├── pcs.py                     # DD/PCS 计算引擎
├── pancancer.py               # 泛癌分析（23种实体瘤）
├── prism_analysis.py          # PRISM 药物敏感性分析
├── msi_analysis.py            # MSI 分层分析
├── mutation_type_analysis.py  # 截短突变 vs 错义突变分析
├── protein_features.py        # UniProt 蛋白特征提取
├── alphafold_analysis.py      # 结构相似性分析
│
├── R_fig1.R ~ R_fig4.R        # 主图生成（R）
├── R_figS1.R ~ R_figS9.R      # 补充图生成（R）
│
├── data/                      # 输入数据（大文件通过 .gitignore 排除）
│   ├── README.md              # 数据下载说明
│   └── cptac_cache/           # CPTAC 蛋白丰度缓存（7个队列）
│
├── output/                    # 分析输出
│   ├── figures/               # PDF 图表（主图 + 补充图）
│   └── tables/                # TSV 数据表（S1-S6）
│
└── R_package/                 # paralogSL R 包（v1.0.0）
    ├── DESCRIPTION
    ├── R/, man/, data/
    └── README.md
```

## 快速开始

### Python 分析管线

```bash
# 1. 安装依赖
pip install -r requirements.txt

# 2. 下载 DepMap 数据（见 data/README.md）

# 3. 运行主分析
python main.py

# 4. 运行泛癌分析
python pancancer.py

# 5. 生成图表（R）
Rscript R_fig1.R
Rscript R_fig2.R
Rscript R_fig3.R
Rscript R_fig4.R
```

### R 包安装与使用

```r
# 安装
devtools::install_github("tjogzt/paralogSL")

# 快速分析
library(paralogSL)
result <- compute_dd(dep_matrix, driver_gene = "ARID1A",
                     paralog_gene = "ARID1B",
                     mut_lines = mut_ids, wt_lines = wt_ids)

# result$DD = 0.182, result$p_value = 1.4e-26
```

### Docker 环境

```bash
docker build -t paralog-sl .
docker run -v $(pwd)/data:/app/data paralog-sl python main.py
```

## 数据来源

| 数据 | 来源 | 链接 |
|------|------|------|
| DepMap 26Q1 | CRISPR 依赖性评分、突变、表达、拷贝数 | https://depmap.org/portal/download/ |
| CPTAC | 7 个队列蛋白质组学数据 | 通过 cBioPortal API（已缓存） |
| PRISM | 1,482 种化合物药物敏感性 | https://depmap.org/portal/download/ |
| 处理后表格 | 补充表 S1-S6 | `output/tables/` |

## 分析管线流程

```
DepMap 基因效应矩阵
    ↓ build_mutation_matrix()     构建突变矩阵
    ↓ compute_dd()               计算 DD（按驱动基因×旁系同源×癌种）
    ↓ compute_auroc()            与金标准对照评估
    ↓ 正交验证层
    │   ├── CPTAC 蛋白质共变
    │   ├── PRISM 药物敏感性
    │   ├── TCGA 生存分析
    │   └── MSI/突变类型分层
    ↓ compute_therapeutic_window() 治疗窗排序
    ↓ 结构可药性分析
    ↓ 输出：候选靶点排名表 + 图表
```

## 可复现性

- **随机种子**：所有分析使用 `set.seed(42)`（R）或 `random_state=42`（Python）
- **Docker**：提供 Dockerfile 保证环境一致性
- **数据版本**：DepMap 26Q1（2026年第一季度发布）

## R 包

独立 R 包仓库：**[tjogzt/paralogSL](https://github.com/tjogzt/paralogSL)** [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.21502114.svg)](https://doi.org/10.5281/zenodo.21502114)

| 功能 | 函数 |
|------|------|
| DD 计算 | `compute_dd()` |
| PCS 计算 | `compute_pcs()` |
| AUROC 评估 | `compute_auroc()` |
| 治疗窗分类 | `compute_therapeutic_window()` |
| 泛癌可视化 | `plot_pancancer_summary()` |
| 治疗窗气泡图 | `plot_therapeutic_window()` |

内置数据集：`known_sl_pairs`（12对金标准）、`solid_tumor_summary`（23种癌型 AUROC）、`benchmark_methods`（8种已发表方法 CV3 值）

## 引用

```
Mo Q, Zhu T. Delta Dependency Prioritizes Paralog-Based Synthetic Lethality
Candidates Across Solid Tumor Types. Genome Biology (2026).
DOI: 10.5281/zenodo.21502031
```

```bibtex
@article{Zhu2026,
  title   = {Delta Dependency Prioritizes Paralog-Based Synthetic Lethality
             Candidates Across Solid Tumor Types},
  author  = {Zhu, Tao},
  journal = {Genome Biology},
  year    = {2026},
  doi     = {10.5281/zenodo.21502031},
}
```

## 许可证

MIT 许可证。详见 [LICENSE](LICENSE)。
