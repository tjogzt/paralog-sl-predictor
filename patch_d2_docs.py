#!/usr/bin/env python3
"""D2 sync: README.md, README_CN.md, cover_letter.md, submission_texts.md."""
from pathlib import Path

ROOT = Path(__file__).parent

PATCHES = {
 "README.md": [
  (r"""AUROC = 0.682 on the full 12-pair curated set; 0.728 excluding two DepMap-derived pairs; 0.774 on pre-DepMap evidence only; unchanged (0.682) under direction-strict relabelling. Both lineage-evaluable Tier A pairs (directional external evidence) rank above all unlabeled controls (Tier A AUROC = 1.000)""",
   r"""AUROC = 0.676 on the full 12-pair curated set; 0.725 excluding two DepMap-derived pairs; 0.774 on pre-DepMap evidence only; unchanged (0.676) under direction-strict relabelling. Both lineage-evaluable positives on the Tier A∪B external benchmark rank above all unlabeled controls (AUROC = 1.000; 2 of 5 benchmark pairs lineage-evaluable)"""),
  (r"""AUROC = 0.682 without training; best published CV3 result 0.790 (SLMGAE)""",
   r"""AUROC = 0.676 without training; best published CV3 result 0.790 (SLMGAE)"""),
  (r"""**Head-to-head (identical 75-pair test set, 6 positives):** the interpretable composite score (0.841) outperforms all four multi-feature classifiers under leave-one-pair-out CV (SVM-RBF 0.744, RF 0.617, SVM-Linear 0.217, LR 0.138; DD alone 0.551)""",
   r"""**Head-to-head (identical 72-pair test set, 6 positives):** the interpretable composite score (0.831) matches the best of four multi-feature classifiers under leave-one-pair-out CV (SVM-RBF 0.843, RF 0.722, SVM-Linear 0.114, LR 0.136; DD alone 0.566)"""),
 ],
 "README_CN.md": [
  (r"""完整 12 对策展集 AUROC = 0.682；去除 2 对 DepMap 来源标签后 0.728；仅用 DepMap 之前证据 0.774；方向严格重标后不变（0.682）。两条可评估的 Tier A 记录（方向性外部实验证据）均排在全部未标记对照之上（Tier A AUROC = 1.000）""",
   r"""完整 12 对策展集 AUROC = 0.676；去除 2 对 DepMap 来源标签后 0.725；仅用 DepMap 之前证据 0.774；方向严格重标后不变（0.676）。Tier A∪B 外部基准上两条可评估阳性记录均排在全部未标记对照之上（AUROC = 1.000；5 对基准对中 2 对可做谱系级评估）"""),
  (r"""无需训练 AUROC = 0.682；已发表 CV3 最佳 0.790（SLMGAE）""",
   r"""无需训练 AUROC = 0.676；已发表 CV3 最佳 0.790（SLMGAE）"""),
  (r"""**头对头比较（同一 75 对测试集、6 个阳性）：** 可解释组合分（0.841）在留一交叉验证下超过全部四种多特征分类器（SVM-RBF 0.744、RF 0.617、SVM-Linear 0.217、LR 0.138；仅用 DD 0.551）""",
   r"""**头对头比较（同一 72 对测试集、6 个阳性）：** 可解释组合分（0.831）在留一交叉验证下与最优多特征分类器持平（SVM-RBF 0.843、RF 0.722、SVM-Linear 0.114、LR 0.136；仅用 DD 0.566）"""),
 ],
 "cover_letter.md": [
  (r"""DD achieves AUROC = 0.682 — robust to label quality (0.728 excluding two DepMap-derived pairs; unchanged under a direction-strict relabelling) — and both lineage-evaluable pairs with direct directional experimental evidence rank above all unlabeled controls (Tier A AUROC = 1.000). In a head-to-head benchmark on the same 75 pairs (6 known positives), the interpretable composite score (AUROC = 0.841) outperforms four multi-feature classifiers under leave-one-pair-out cross-validation (SVM-RBF 0.744, random forest 0.617; linear classifiers unstable at this sample size).""",
   r"""DD achieves AUROC = 0.676 — robust to label quality (0.725 excluding two DepMap-derived pairs; unchanged under a direction-strict relabelling) — and both lineage-evaluable positives on the Tier A∪B external benchmark (direct or genotype-conditional experimental evidence) rank above all unlabeled controls (AUROC = 1.000; 2 of 5 benchmark pairs lineage-evaluable). In a head-to-head benchmark on the same 72 pairs (6 known positives), the interpretable composite score (AUROC = 0.831) matches the best of four multi-feature classifiers under leave-one-pair-out cross-validation (SVM-RBF 0.843, random forest 0.722; linear classifiers unstable at this sample size)."""),
 ],
 "submission_texts.md": [
  (r"""DD nominates priority paralog-SL candidates without any training (AUROC 0.682 on the tiered gold-standard set; 1.000 on the two pairs with directional external experimental evidence). Bottom: exploratory support across genomic (23 solid tumor lineages; AUROC > 0.7 in 10 of 14 evaluable lineages)""",
   r"""DD nominates priority paralog-SL candidates without any training (AUROC 0.676 on the tiered gold-standard set; 1.000 on the Tier A∪B external benchmark, 2 of 5 pairs lineage-evaluable). Bottom: exploratory support across genomic (23 solid tumor lineages; AUROC > 0.7 in 7 of 8 evaluable lineages)"""),
  (r"""2. Interpretable composite score (AUROC 0.841) outperforms four multi-feature classifiers in cross-validation""",
   r"""2. Interpretable composite score (AUROC 0.831) matches best multi-feature classifier in CV"""),
  (r"""1: 83  2: 85  3: 84  4: 79  5: 76 -->""",
   r"""1: 83  2: 84  3: 84  4: 79  5: 76 -->"""),
 ],
}

errors = []
for fname, pairs in PATCHES.items():
    p = ROOT / fname
    text = p.read_text()
    for old, new in pairs:
        n = text.count(old)
        if n != 1:
            errors.append(f"{fname}: {n} occurrences: {old[:70]!r}")
        else:
            text = text.replace(old, new, 1)
    if not errors:
        p.write_text(text)

if errors:
    print("FAILED:")
    for e in errors:
        print(" ", e)
    raise SystemExit(1)
print("OK: README.md, README_CN.md, cover_letter.md, submission_texts.md synced")
