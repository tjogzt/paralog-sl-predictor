#!/usr/bin/env bash
# verify_all.sh — one-command reproduction of every quantitative claim in the
# manuscript (headline metrics, ML benchmark, regression controls) plus the
# test suite. Each script ends with an automatic claims check that compares
# its recomputed numbers against the values written in manuscript.tex and
# exits non-zero on mismatch.
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

echo "==> [1/4] Headline metrics (compute_headline_metrics.py)"
"$PY" compute_headline_metrics.py

echo "==> [2/4] ML benchmark (ml_benchmark.py)"
"$PY" ml_benchmark.py

echo "==> [3/4] Regression controls (regression_controls.py)"
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

echo "==> [4/4] Test suite (pytest)"
"$PY" -m pytest tests/ -q

echo
echo "ALL CHECKS PASSED — see output/headline_metrics.json,"
echo "output/ml_benchmark.json and output/regression_controls.json"
echo "for the per-claim comparison tables."
