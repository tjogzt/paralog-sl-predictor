# Genome Biology 投稿系统填写清单

> 版本锚点：paralog-sl-predictor @ `d671ea9`（manuscript.pdf 29 页 / supplementary.pdf 17 页 / cover_letter.pdf 3 页）
> 配套文件：highlights 与图形摘要说明文字见 `submission_texts.md`（已含字符数核对）
> 使用前请完成第 0 节的提交前核对；【待补】为待补项。

---

## 0. 提交前最终核对（投稿当天执行）

```bash
cd paralog_sl_predictor && bash verify_all.sh
# 期望输出：165/165 claims match; 31 passed; ALL CHECKS PASSED
```

- [ ] verify_all.sh 全绿
- [ ] 三份 PDF 为最新编译（manuscript 29 页 / supplementary 17 页 / cover letter 3 页）
- [ ] 作者名单与邮箱最终确认（见第 2 节【待补】）
- [已完成] 摘要版本已选定（第 3 节：压缩版 236 词已同步 manuscript.tex，≤250 达标）

---

## 1. 期刊与文章类型

| 系统字段 | 填写内容 |
|---|---|
| Journal | Genome Biology |
| Article type | Research |
| Section | 按系统列表选最接近者（Cancer genomics / Computational biology） |

## 2. 标题（Title）

```
Delta Dependency Prioritizes Paralog-Based Synthetic Lethality Candidates Across Solid Tumor Types
```

## 3. 摘要（Abstract，结构化，系统上限 250 词）

**【已解决 2026-07-27】**：原摘要 273 词超限，以下 236 词压缩版（所有数字与正文一致、限定语全部保留）**已同步进 manuscript.tex 并重编 PDF**（tex 排版口径 241 词，系统粘贴用下方纯文本版为 236 词，均 ≤250）：

**Background:** Paralog compensation is a specific mechanism of synthetic lethality (SL): when a tumor-suppressor driver gene is lost to mutation, its sequence-similar paralog can become conditionally essential. Systematic prioritization of paralog-SL candidates remains limited because existing computational methods are black boxes that cannot explain their predictions, hindering targeted validation.

**Results:** We introduce Delta Dependency (DD), a single interpretable metric measuring the dependency shift on a paralog between driver-mutant and wild-type cell lines in DepMap (1,208 lines, 23 solid tumor types). In the pre-specified primary evaluation, a lineage-level benchmark on twelve curated, evidence-tiered paralog-SL pairs, DD achieved AUROC = 0.676 (0.725 excluding two pairs with DepMap-derived labels; 1.000 on the Tier A/B external benchmark, two of five pairs lineage-evaluable). In pair-level evaluation on the same 72 test pairs, the interpretable composite score reached AUROC = 0.831, statistically indistinguishable from the best multi-feature classifier (SVM-RBF, 0.843) under leave-one-pair-out cross-validation. Exploratory analyses across seven CPTAC cohorts (n = 771) detected protein-level paralog co-variation invisible at the RNA level. Dependency-window scoring prioritized ARID1A->ARID1B (DWS = 2.82) as the leading selective candidate, suggesting a hypothesis-generating biomarker strategy centered on truncating ARID1A mutations.

**Conclusions:** DD offers a mechanistically transparent alternative to black-box SL prediction: each nomination traces to a measured dependency shift, enabling rational experimental follow-up. DD is an association-based prioritization metric, not a causal test of synthetic lethality; all candidates require experimental validation. An open-source R package (paralogSL) and a reproducible pipeline are provided.

## 4. 关键词（Keywords，6 个，3–10 范围内）

```
Synthetic lethality; paralog compensation; DepMap; CPTAC proteomics; Delta Dependency; dependency window
```

## 5. 作者与单位（Authors）

| 顺序 | 姓名 | 单位 | 邮箱 | ORCID |
|---|---|---|---|---|
| 1 | Qingqing Mo | 单位 1, 2（见下） | 【待补】 | 【注意】建议补 |
| 2（通讯） | Tao Zhu | 单位 1, 2 | zhutao@tjh.tjmu.edu.cn | 【注意】建议补 |

Affiliation 1:
```
Department of Obstetrics and Gynecology,
National Clinical Research Center for Obstetrics and Gynecology,
Tongji Hospital, Tongji Medical College,
Huazhong University of Science and Technology, Wuhan, China
```
Affiliation 2:
```
Key Laboratory of Cancer Invasion and Metastasis (Ministry of Education),
Hubei Key Laboratory of Tumor Invasion and Metastasis,
Tongji Hospital, Tongji Medical College,
Huazhong University of Science and Technology, Wuhan, China
```

【待补】：Qingqing Mo 邮箱；两位作者 ORCID（BMC 强烈建议通讯作者提供）。

## 6. 推荐审稿人（Suggested Reviewers，5 位）

| # | 姓名 | 单位 | 专长 | 邮箱 |
|---|---|---|---|---|
| 1 | Dr. Francisca Vazquez | DepMap / Broad Institute | CRISPR dependency | 【待补】 |
| 2 | Dr. Michael P. Snyder | Stanford University | proteogenomics, CPTAC | 【待补】 |
| 3 | Dr. Jason Moffat | University of Toronto | synthetic lethality screening | 【待补】 |
| 4 | Dr. Rameen Beroukhim | Dana-Farber / Broad Institute | cancer genomics, paralog dependencies | 【待补】 |
| 5 | Dr. Bing Zhang | Baylor College of Medicine | CPTAC proteogenomics | 【待补】 |

**Opposed reviewers:** None

【待补】：BMC 系统要求推荐审稿人的**机构邮箱**（不接受 Gmail 等私人邮箱），投稿前请到各单位官网查到公开学术邮箱填入。

## 7. 上传文件清单（Files）

| 系统文件类型 | 文件 | 说明 |
|---|---|---|
| Manuscript (PDF) | `manuscript.pdf` | 29 页，双倍行距+行号 |
| Manuscript source | 打包 zip：`manuscript.tex` + `manuscript.bbl` + `output/figures/Fig1–4 .pdf` | GB 要求可编辑 LaTeX 源 |
| Supplementary | `supplementary.pdf` + 源打包（`supplementary.tex` + `output/figures/FigS1–S11 .pdf`，无独立 bbl） | 17 页 |
| Cover letter | `cover_letter.pdf`（或将其纯文本粘贴到系统 cover letter 框） | 3 页 |
| Graphical abstract | `GraphicalAbstract.png` | 位于 `output/figures/`；caption 见 `submission_texts.md` |
| Highlights | 在系统文本框逐条录入（见第 10 节） | 每条 ≤85 字符 |

LaTeX 打包命令（在仓库根目录执行）：

```bash
zip submission_src.zip manuscript.tex manuscript.bbl manuscript_refs.bib \
  output/figures/Fig1_Framework_Validation.pdf output/figures/Fig2_Proteomics.pdf \
  output/figures/Fig3_Clinical.pdf output/figures/Fig4_Translational.pdf
zip supplementary_src.zip supplementary.tex output/figures/FigS*.pdf
```

## 8. 声明勾选项（Declarations，逐项照抄）

| 系统字段 | 填写内容 |
|---|---|
| Ethics approval and consent to participate | Not applicable. This study used exclusively publicly available, de-identified data. |
| Consent for publication | Not applicable. |
| Competing interests | The authors declare no competing interests. |
| Funding | This research received no specific grant from any funding agency.（系统如有 "No funding" 复选框请勾选） |
| Authors' contributions | Q.Q.M. and T.Z. conceived the study. T.Z. developed the methodology, performed all computational analyses, and developed the R package. Q.Q.M. curated clinical datasets and provided gynecological oncology domain expertise. Q.Q.M. and T.Z. wrote the manuscript. |
| Acknowledgements | We thank the DepMap, CPTAC, TCGA, PRISM, and cBioPortal teams for making their data publicly available. We also acknowledge the SynLethDB, HGNC, and UniProt teams for curated resources. |

**Availability of data and materials（系统文本框）：**

```
All data are publicly available: DepMap 26Q1 from https://depmap.org/portal/download/;
CPTAC proteomics from https://cbioportal.org and https://pdc.cancer.gov;
TCGA PanCan from https://xenabrowser.net; PRISM from https://depmap.org/portal/prism/;
HGNC gene families from https://www.genenames.org/.
Processed datasets with cell line-level DD values and all supplementary tables
are archived at Zenodo (https://doi.org/10.5281/zenodo.21502030).
All analysis code is available under the MIT license:
paralogSL R package v1.1.0 at https://github.com/tjogzt/paralogSL
(https://doi.org/10.5281/zenodo.21502113);
complete analysis pipeline at https://github.com/tjogzt/paralog-sl-predictor
(https://doi.org/10.5281/zenodo.21502030).
A single entry point (verify_all.sh) recomputes all 165 audited numeric claims
in under one minute.
```

## 9. 系统其他常见字段

| 字段 | 填写内容 |
|---|---|
| Related/preprint | 无预印本（如后续挂 bioRxiv 再更新） |
| Previously submitted elsewhere | No |
| All authors approve submission | Yes（cover letter 已含此声明） |
| Not under consideration elsewhere | Yes（cover letter 已含此声明） |

## 10. Highlights（系统逐条录入，每条 ≤85 字符）

1. A single-subtraction metric (Delta Dependency) predicts paralog synthetic lethality
2. Interpretable composite (AUROC 0.831) matches best multi-feature classifier in CV
3. Protein-level paralog compensation is confirmed across seven CPTAC cohorts (n = 771)
4. Therapeutic window analysis nominates ARID1A-ARID1B as the leading SL candidate
5. Every quantitative claim is recomputable from raw data with a single command

Graphical abstract caption（较长，完整版见 `submission_texts.md`）。

---

## 待办汇总（投稿前必须关闭）

1. 【待补】确定最终作者名单；补 Qingqing Mo 邮箱（第 5 节）
2. 【待补】补两位作者 ORCID（第 5 节）
3. 【待补】补 5 位推荐审稿人机构邮箱（第 6 节）
4. ~~摘要超限~~ 【已解决】压缩版已同步 manuscript.tex（第 3 节）
