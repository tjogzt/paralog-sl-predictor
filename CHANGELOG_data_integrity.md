# 数据诚信整改变更对照（2026-07-26，commit 7b5de99）

> 用途：(1) 投稿前终审图 3 与 TCGA 段落的对照依据；(2) 审稿阶段如被问及数字
> 来源时的可追溯记录。本轮所有变动均可由 `verify_all.sh` 一键复现，每项数字
> 的 artifact 出处见下表；闭环审计 155/155 通过
> （`output/manuscript_number_audit.tsv`）。

## 一、背景

全面代码清查发现三类违反"禁止模拟/随机/硬编码数据"规则的遗留问题：

1. **静默 fallback 分支**：artifact 拉取失败时用内置数值顶替的 else 分支
   （R_fig1/R_fig2/R_fig4/R_figS2，及 tcga_survival.py）；
2. **被取代的旧版绘图文件**：`validation_viz.py` 的 fig1–fig6 函数（手稿图
   早已改由 R 脚本生成，旧函数仍产出与现行数字矛盾的图）；
3. **隐性硬编码 panel**：数值与 artifact 一致但写死在绘图脚本中、未从
   artifact 读取、审计无法覆盖（R_fig3 四个 panel）。

## 二、文稿数字变动对照（需要作者确认接受的部分）

### 2.1 TCGA 生存分析（Results 段落、图 3c、Methods）——本轮唯一实质科学变动

| 项目 | 旧值（伪造，已删除） | 新值（真实 Cox PH 重跑） | artifact 出处 |
|---|---|---|---|
| 数据源（Methods） | "UCSC Xena browser" | cBioPortal API，study `brca_tcga_pan_can_atlas_2018` | `tcga_survival.py` |
| 样本量 | n=1,082 | **n=1,069**（151 死亡事件） | `output/tcga_survival_associations.csv` |
| 分析方法 | 未说明（实为中位数比值 + Mann-Whitney，忽略删失；CI 无计算来源） | Cox 比例风险模型（statsmodels PHReg），中位数分高低表达，Wald 95% CI | 同上 |
| 显著基因 | BRCA2 (HR=1.116, p=0.032)；ATR (HR=1.112, p=0.039) | **仅 ARID1B**（HR=1.613，95% CI 1.163–2.238，p=0.004） | 同上 |
| BRCA2 | 显著 | 不显著（HR=1.245，p=0.181） | 同上 |
| ATR | 显著 | 不显著（HR=0.980，p=0.902） | 同上 |
| 多重校正表述 | "两基因经 Bonferroni 均不再显著" | ARID1B 通过 8 个关键基因校正（α=0.00625），但不通过全部 32 基因校正（文稿如实说明） | 同上 |
| 分析基因数 | 8（实为 10 个 fallback 值中的 8 个展示） | 32（16 个优先旁系对的所有旁系基因），图 3c 展示 8 个关键基因 | 同上 |

**叙事影响**：旧版是"两个弱效应基因且过不了校正"（自我削弱的负面结果）；
新版是"全流程头号候选 ARID1A→ARID1B 的旁系基因是唯一显著的生存关联基因"
（HR=1.613，61% 风险升高），与主候选形成跨层一致的故事线。

### 2.2 共突变分析（图 3d、caption、新增 Methods 段落）

| 旁系对 | 旧 OR（无来源，已删除） | 新 OR（Fisher 精确检验，1,208 细胞系） | p 值（新） |
|---|---|---|---|
| PIK3CA/PIK3CB | 5.157 | **6.836** [3.285–14.226] | 3.3×10⁻⁶ |
| EP300/CREBBP | 4.548 | **6.111** [3.744–9.973] | 2.4×10⁻¹¹ |
| ARID1A/ARID1B | 6.147 | **5.722** [3.646–8.978] | 2.0×10⁻¹² |
| BRCA1/BRCA2 | 3.422 | **4.528** [1.290–15.895] | 0.040 |
| SMARCA4/SMARCA2 | 4.753 | **1.623** [0.627–4.203] | 0.372（ns） |

- caption 更正："range 3.4–6.1，all p<10⁻⁴" → "range 1.6–6.8，4/5 名义显著"；
- 核心结论"所有关键对 OR>1（共突变而非互斥）"**仍然成立**；
- artifact：`output/cooccurrence_analysis.csv`（新脚本 `cooccurrence_analysis.py`，
  复用主管线 driver-rule 突变定义：TSG=LikelyLoF，癌基因=Hotspot）；
- Methods 新增"Mutational co-occurrence"段落。

### 2.3 无数字变动、仅改为从 artifact 读取的 panel

| 位置 | 内容 | 现状 |
|---|---|---|
| 图 3a（MSI） | 0.767/0.712/0.838/0.556，n=14/45/17/11 | 数值不变，改读 `msi_key_numbers_min3.json` |
| 图 3b（突变类型） | 0.388/0.020/0.464/0.150/0.080/0.000/−0.136/0.000 | 数值不变，改读 `muttype_{ovarian,colorectal,breast}_results.csv` |
| 图 S9 | r=0.88/0.91/0.41，n=50/13/37 | **经重算全部为真**（0.882/0.906/0.413，p=2.5×10⁻¹⁷），新落盘 `kmer_nw_correlation.csv` 并纳入审计 |

## 三、代码层处置（无文稿数字影响）

| 处置 | 位置 |
|---|---|
| 静默 fallback 改为 `stop()` 硬报错（缺 artifact 即失败，绝不顶替） | R_fig1.R、R_fig2.R（含 runif 伪造 CPTAC 热图分支）、R_fig4.R（4 处）、R_figS2.R |
| 删除六个遗留伪造图函数（旧假 HR、"estimated" CNV R²、无源转移矩阵与候选表），保留 main.py 使用的验证套件 | validation_viz.py（624→177 行） |
| cBioPortal 请求加自动重试（偶发 SSL 故障），修正已废弃的 `clinicalDataType=SURVIVAL` 枚举为 `PATIENT` | tcga_survival.py |

## 四、验证状态

- `audit_manuscript_numbers.py` 扩充 +30 断言（TCGA×11、共突变×14、S9×7），
  **155/155 全部匹配**；
- `verify_all.sh` 五步全绿，31 个 pytest 全过；
- manuscript.pdf（27 页）/ supplementary.pdf（16 页）重编零错误，PDF 内旧值
  （1.116、1.112、n=1,082）经全文检索确认已清除；
- 文档层清扫（cover_letter.md、highlights、README、R 包）：零残留；R 包源码
  无随机/模拟模式，无需更新。
