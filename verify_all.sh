#!/usr/bin/env bash
# verify_all.sh — one-command reproduction of every quantitative claim in the
# manuscript (headline metrics, ML benchmark, regression controls, and the
# 237-claim closed-loop number audit covering manuscript + supplementary) plus
# the test suite. Each script ends with an automatic claims check that
# compares its recomputed numbers against the values written in the
# manuscript and exits non-zero on mismatch.
#
# Usage:
#   ./verify_all.sh            # fast: reuses cached data slices when present
#   VERIFY_FULL=1 ./verify_all.sh   # full: rebuilds every cache from raw data
#
# Requirements: DepMap CSVs under data/ (see README "Reproducibility
# checklist"), Python 3 with pandas/numpy/scipy/statsmodels/scikit-learn.
set -euo pipefail
cd "$(dirname "$0")"

PY=${PYTHON:-python}
CACHE=output/cache
NEED_CACHE=("$CACHE/dep_slice.pkl" "$CACHE/mut_matrix.pkl" \
            "$CACHE/expr_slice.pkl" "$CACHE/cnv_slice.pkl")

echo "==> [1/5] Headline metrics (compute_headline_metrics.py)"
"$PY" compute_headline_metrics.py

echo "==> [2/5] ML benchmark (ml_benchmark.py)"
"$PY" ml_benchmark.py

echo "==> [3/5] Regression controls (regression_controls.py)"
FULL=${VERIFY_FULL:-0}
CACHE_OK=1
for f in "${NEED_CACHE[@]}"; do [ -f "$f" ] || CACHE_OK=0; done
if [ "$FULL" = "1" ] || [ "$CACHE_OK" = "0" ]; then
  echo "    full recompute (VERIFY_FULL=$FULL, cache complete=$CACHE_OK)"
  "$PY" regression_controls.py --stage all
else
  echo "    using cached slices in $CACHE (set VERIFY_FULL=1 to rebuild)"
  "$PY" regression_controls.py --stage analyze
fi

echo "==> [4/5] Manuscript number audit (audit_manuscript_numbers.py)"
"$PY" audit_manuscript_numbers.py

echo "==> [5/5] Test suite (pytest)"
"$PY" -m pytest tests/ -q

echo
echo "==> Session information"
"$PY" -c "import sys, platform, numpy, scipy, pandas, sklearn, statsmodels; \
print('python    :', sys.version.split()[0]); \
print('platform  :', platform.platform()); \
print('numpy     :', numpy.__version__); \
print('scipy     :', scipy.__version__); \
print('pandas    :', pandas.__version__); \
print('sklearn   :', sklearn.__version__); \
print('statsmodels:', statsmodels.__version__)"

echo
echo "ALL CHECKS PASSED — see output/headline_metrics.json,"
echo "output/ml_benchmark.json, output/regression_controls.json and"
echo "output/manuscript_number_audit.tsv for the per-claim comparison tables."
