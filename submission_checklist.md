# Genome Biology 投稿系统填写清单

> 版本锚点：【待补——最终 commit 后回填 commit hash】（manuscript.pdf 46 页 / supplementary.pdf 26 页 / cover_letter.pdf 2 页，2026-08-02 编译）
> 配套文件：highlights 与图形摘要说明文字见 `submission_texts.md`（已含字符数核对）
> 使用前请完成第 0 节的提交前核对；【待补】为待补项。

---

## 0. 提交前最终核对（投稿当天执行）

```bash
cd paralog_sl_predictor && bash verify_all.sh
# 期望输出：420/420 claims match; 31 passed; ALL CHECKS PASSED（实测 ~30 s）
```

- [ ] verify_all.sh 全绿
- [ ] 四份 PDF 为最新编译（manuscript 46 页 / supplementary 26 页 / cover letter 2 页 / submission_texts 1 页）
- [ ] 作者名单与邮箱最终确认（见第 2 节【待补】）
- [已完成] 摘要版本已选定（第 3 节：当前版纯文本 237 词已同步 manuscript.tex，≤250 达标）

---

## 1. 期刊与文章类型

| 系统字段 | 填写内容 |
|---|---|
| Journal | Genome Biology |
| Article type | Research |
| Section | 按系统列表选最接近者（Cancer genomics / Computational biology） |

## 2. 标题（Title）

```
Delta Dependency Prioritizes Candidate Paralog Dependencies Across Solid Tumor Types
```

## 3. 摘要（Abstract，结构化，系统上限 250 词）

**【已同步 2026-08-01】**：以下摘要已与 manuscript.tex 当前版本**逐字一致**（纯文本化：去 LaTeX 符号，∪ 写作 Tier A + B，→ 写作 ->），纯文本 237 词，低于 250 上限。本轮已写入 pair-clustered bootstrap CI（0.253-0.933）与 CPTAC 肿瘤纯度稳健性表述，并删除"invisible at RNA level"旧结论。

**Background:** Paralog compensation is a specific synthetic-lethal (SL) mechanism: when a tumor-suppressor driver is lost to mutation, its sequence-similar paralog can become conditionally essential. Systematic prioritization of paralog-SL candidates remains limited; existing models integrate multiple features, whereas a measurable dependency-shift statistic offers a complementary, transparent strategy.

**Results:** We introduce Delta Dependency (DD), an interpretable metric measuring the dependency shift on a paralog between driver-mutant and wild-type cell lines in DepMap (1,208 lines, 23 solid tumor types). On twelve curated, evidence-tiered pairs in the pre-specified primary benchmark (three gynecological lineages), signed DD achieved AUROC = 0.629 (pair-clustered bootstrap 95% CI 0.253-0.933; 0.613 excluding two DepMap-labelled pairs; 0.525 for the eight pairs with pre-DepMap experimental evidence; 1.000 on the Tier A + B benchmark, both evaluable entries being the same pair). At pair level on the same 72 pairs, the interpretable composite score reached AUROC = 0.831, indistinguishable from the best classifier (SVM-RBF, 0.841; six positives, so descriptive). Exploratory analysis of seven CPTAC cohorts (n = 771) detected protein-level paralog co-variation that is robust to tumor-purity adjustment. Dependency-window scoring recovered the established ARID1A->ARID1B pair (DWS = 2.82) as the leading selective candidate and nominated KMT2D->KMT2C for experimental follow-up.

**Conclusions:** DD offers a mechanistically transparent complement to multi-feature SL prediction: each nomination traces to a measured dependency shift. DD is an association-based prioritization metric, not a causal test of synthetic lethality; all candidates require experimental validation. An open-source R package (paralogSL) and reproducible pipeline are provided.

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

| # | 姓名 | 单位 | 邮箱 |
|---|---|---|---|
| 1 | Dr. Francisca Vazquez | DepMap / Broad Institute | vazquez@broadinstitute.org |
| 2 | Dr. Michael P. Snyder | Stanford University | mpsnyder@stanford.edu |
| 3 | Dr. Jason Moffat | University of Toronto / SickKids | jason.moffat@sickkids.ca |
| 4 | Dr. Rameen Beroukhim | Dana-Farber / Broad Institute | rameen_beroukhim@dfci.harvard.edu |
| 5 | Dr. Bing Zhang | Baylor College of Medicine | bing.zhang@bcm.edu |

**Opposed reviewers:** None

【已核实 2026-07-27】邮箱均取自各单位官网：Vazquez（[target-discovery.depmap.org](https://target-discovery.depmap.org/)）、Snyder（[med.stanford.edu/snyderlab](https://med.stanford.edu/snyderlab/about.html)）、Moffat（[moleculargenetics.utoronto.ca](https://moleculargenetics.utoronto.ca/faculty/jason-moffat)）、Beroukhim（[ogephd.hms.harvard.edu](https://ogephd.hms.harvard.edu/people/rameen-beroukhim)）、Zhang（[bcm.edu/people-search](https://www.bcm.edu/people-search/bing-zhang-33575)）。全部为机构邮箱，已同步至 cover_letter.md。

## 7. 上传文件清单（Files）

| 系统文件类型 | 文件 | 说明 |
|---|---|---|
| Manuscript (PDF) | `manuscript.pdf` | 46 页，双倍行距+行号 |
| Manuscript source | 打包 zip：`manuscript.tex` + `output/figures/Fig1–4 .pdf`（参考文献已内嵌 `thebibliography`，单文件即可编译，无需 bbl/bib） | GB 要求可编辑 LaTeX 源 |
| Supplementary | `supplementary.pdf` + 源打包（`supplementary.tex` + `output/figures/FigS1–S10` 组合图 PDF，注意排除 `_panel_` 单图） | 26 页 |
| Cover letter | `cover_letter.pdf`（或将其纯文本粘贴到系统 cover letter 框） | 2 页 |
| Graphical abstract | `GraphicalAbstract.png` | 位于 `output/figures/`；caption 见 `submission_texts.md` |
| Highlights | 在系统文本框逐条录入（见第 10 节） | 每条 ≤85 字符 |

LaTeX 打包命令（在仓库根目录执行；参考文献内嵌于 `manuscript.tex`，无需 bbl）：

```bash
zip submission_src.zip manuscript.tex \
  output/figures/Fig1_Framework_Validation.pdf output/figures/Fig2_Proteomics.pdf \
  output/figures/Fig3_Clinical.pdf output/figures/Fig4_Translational.pdf
zip supplementary_src.zip supplementary.tex \
  output/figures/FigS1_CellLine_Landscape.pdf output/figures/FigS2_Evaluation_Robustness.pdf \
  output/figures/FigS3_CNV_Independence.pdf output/figures/FigS4_CPTAC_PerCohort.pdf \
  output/figures/FigS5_MutationType.pdf output/figures/FigS6_Cooccurrence_TMB.pdf \
  output/figures/FigS7_PRISM_Selectivity.pdf output/figures/FigS8_TherapeuticWindow.pdf \
  output/figures/FigS9_Sequence_Structure_Descriptors.pdf
```

> 提示：投稿包 `submission_2026-07-31/` 中已按此结构备妥 `manuscript_source/` 与 `supplementary_source/`，可直接压缩上传。

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
paralogSL R package v1.1.1 at https://github.com/tjogzt/paralogSL
(https://doi.org/10.5281/zenodo.21502113);
complete analysis pipeline at https://github.com/tjogzt/paralog-sl-predictor
(https://doi.org/10.5281/zenodo.21502030).
A single entry point (verify_all.sh) recomputes all 420 audited numeric claims
in ~30 s.
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
3. Protein-level paralog co-variation is detectable across seven CPTAC cohorts (n = 771)
4. Dependency-window analysis nominates ARID1A-ARID1B as leading selective candidate
5. Every quantitative claim is recomputable from raw data with a single command

（与 `submission_texts.md` 逐字一致；2026-07-31 逐条复核字符数 1:83 2:81 3:85 4:81 5:76，均 ≤85）

Graphical abstract caption（较长，完整版见 `submission_texts.md`）。

---

## 待办汇总（投稿前必须关闭）

1. 【待补】确定最终作者名单；补 Qingqing Mo 邮箱（第 5 节）
2. 【待补】补两位作者 ORCID（第 5 节）
3. ~~补 5 位推荐审稿人机构邮箱~~ 【已解决 2026-07-27】官网核实，见第 6 节
4. ~~摘要超限~~ 【已解决】压缩版已同步 manuscript.tex（第 3 节）
