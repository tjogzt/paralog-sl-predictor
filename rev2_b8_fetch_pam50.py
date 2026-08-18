#!/usr/bin/env python3
"""
rev2_b8_fetch_pam50.py  (Stage-4 revision, item B8, data fetch)
===============================================================
Fetches the sample-level PAM50 subtype call (PAM50_CALL_RNASEQ) for
brca_tcga_pan_can_atlas_2018 from the cBioPortal API so the BRCA ARID1B
Cox model can be refit with a PAM50 covariate (rev2_b8_tcga_ph_retest.R).
Saved to output/revision_stage4/rev2_b8_pam50_calls.csv.
If the fetch fails, this script exits non-zero and B8 records
NOT COMPUTABLE for the PAM50-adjusted model.
"""

import sys
import time
from pathlib import Path

import pandas as pd
import requests

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "output" / "revision_stage4"
OUT.mkdir(parents=True, exist_ok=True)

BASE = "https://www.cbioportal.org/api"
STUDY = "brca_tcga_pan_can_atlas_2018"


def req(method, url, retries=4, backoff=3.0, **kw):
    last = None
    for att in range(retries):
        try:
            r = requests.request(method, url, **kw)
            if r.status_code < 500:
                return r
            last = RuntimeError(f"HTTP {r.status_code}")
        except Exception as e:  # noqa: BLE001
            last = e
        time.sleep(backoff * (2 ** att))
    raise RuntimeError(f"request failed: {last}")


def main():
    url = f"{BASE}/studies/{STUDY}/clinical-data"
    params = {"clinicalDataType": "SAMPLE", "projection": "DETAILED",
              "pageSize": 100000}
    r = req("GET", url, params=params, timeout=180)
    r.raise_for_status()
    rows = [i for i in r.json() if i.get("clinicalAttributeId") == "PAM50_CALL_RNASEQ"]
    if not rows:
        raise RuntimeError("PAM50_CALL_RNASEQ attribute not present in SAMPLE clinical data")
    df = pd.DataFrame(rows)[["sampleId", "patientId", "value"]].rename(
        columns={"sampleId": "sample_id", "patientId": "patient_id", "value": "pam50"})
    out = OUT / "rev2_b8_pam50_calls.csv"
    df.to_csv(out, index=False)
    print(f"fetched {len(df)} PAM50 calls -> {out}")
    print(df["pam50"].value_counts().to_string())


if __name__ == "__main__":
    try:
        main()
    except Exception as e:  # noqa: BLE001
        print(f"PAM50 FETCH FAILED: {e}", file=sys.stderr)
        sys.exit(1)
