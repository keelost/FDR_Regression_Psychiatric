#!/usr/bin/env python3
"""
Build a merged long-format summary table across all targets.

For each target and analysis type (SNP/fdrreg, MAGMA, MetaXcan per brain region,
SMultiXcan), count entries with qval < 0.01 and fdr_theoretical < 0.01, plus totals.
Risk loci counts (threshold 0.01) are attached to the SNP row only.
"""

import os
import re
import pandas as pd

BASE = os.environ.get("FDRREG_RESULTS_DIR", "data/pipeline")
OUT_DIR = os.path.join(BASE, "01.extra.analysis", "09.summary_tables")
OUT_FILE = os.path.join(OUT_DIR, "all_targets_fdr_summary.csv")

THR = 0.01

# Column names differ across analysis types.
# (analysis_type, qval_col, fdrthe_col)
QVAL_COLS = {
    "SNP": ("qval_bh", "fdr_theoretical"),
    "MAGMA": ("qval", "FDR.the"),
    "MetaXcan": ("qval", "FDR_theoretical"),
    "SMultiXcan": ("qval", "FDR.the"),
}


def is_target_dir(name, path):
    """A target dir does not start with 'NN.' (those are analysis dirs)."""
    if not os.path.isdir(path):
        return False
    if re.match(r"^\d+\.", name):
        return False
    return True


def count_file(path, sep, qcol, fcol):
    r"""Return (n_qval_sig, n_qval_total, n_fdrthe_sig, n_fdrthe_total).

    sep: ',' for csv, r'\s+' for whitespace-delimited txt.
    Missing files or columns yield None so we can flag them downstream.
    """
    if not os.path.isfile(path):
        return None
    try:
        df = pd.read_csv(path, sep=sep, engine="python")
    except Exception as e:
        print(f"[WARN] failed to read {path}: {e}")
        return None

    res = {}
    for label, col in (("q", qcol), ("f", fcol)):
        if col not in df.columns:
            print(f"[WARN] column '{col}' missing in {path}")
            res[label + "_sig"] = None
            res[label + "_tot"] = None
            continue
        vals = pd.to_numeric(df[col], errors="coerce")
        res[label + "_tot"] = int(vals.notna().sum())
        res[label + "_sig"] = int((vals < THR).sum())
    return (res["q_sig"], res["q_tot"], res["f_sig"], res["f_tot"])


def get_risk_loci(target):
    """Extract n_risk_loci at threshold 0.01 for qval_bh and fdr_theoretical."""
    path = os.path.join(
        BASE, "01.extra.analysis", "05.risk_loci", target,
        f"{target}.risk_loci_summary.csv",
    )
    out = {"qval_bh": None, "fdr_theoretical": None}
    if not os.path.isfile(path):
        return out
    try:
        rl = pd.read_csv(path)
    except Exception as e:
        print(f"[WARN] failed to read {path}: {e}")
        return out
    rl["threshold"] = pd.to_numeric(rl["threshold"], errors="coerce")
    for metric in ("qval_bh", "fdr_theoretical"):
        sub = rl[(rl["metric"] == metric) & (rl["threshold"] == THR)]
        if len(sub):
            out[metric] = int(sub["n_risk_loci"].iloc[0])
    return out


def build_row(target, analysis, region, counts, risk=None):
    q_sig, q_tot, f_sig, f_tot = counts if counts else (None, None, None, None)
    row = {
        "target": target,
        "analysis_type": analysis,
        "region": region if region else "",
        "n_qval_sig_lt0.01": q_sig,
        "n_qval_total": q_tot,
        "n_fdrthe_sig_lt0.01": f_sig,
        "n_fdrthe_total": f_tot,
        "n_risk_loci_qval_bh_0.01": risk["qval_bh"] if risk else "",
        "n_risk_loci_fdrthe_0.01": risk["fdr_theoretical"] if risk else "",
    }
    return row


def process_target(target):
    rows = []
    tdir = os.path.join(BASE, target)

    # 1. SNP (fdrreg) -- comma separated; risk loci attached here.
    snp_path = os.path.join(
        tdir, "02.fdrreg_results", "fdr_values_per_snp_nolasso.csv")
    qc, fc = QVAL_COLS["SNP"]
    counts = count_file(snp_path, ",", qc, fc)
    risk = get_risk_loci(target)
    rows.append(build_row(target, "SNP", "", counts, risk))

    # 2. MAGMA -- whitespace separated.
    magma_path = os.path.join(
        tdir, "05.magma_fdrreg", f"{target}.gene.bio.fdrreg.txt")
    qc, fc = QVAL_COLS["MAGMA"]
    counts = count_file(magma_path, r"\s+", qc, fc)
    rows.append(build_row(target, "MAGMA", "", counts))

    # 3. MetaXcan -- one row per brain region.
    mx_root = os.path.join(tdir, "07.metaxcan_fdrreg", "01.fdrreg_results")
    qc, fc = QVAL_COLS["MetaXcan"]
    if os.path.isdir(mx_root):
        for region in sorted(os.listdir(mx_root)):
            rdir = os.path.join(mx_root, region)
            if not os.path.isdir(rdir):
                continue
            rpath = os.path.join(rdir, f"{region}.gene.bio.fdrreg.txt")
            counts = count_file(rpath, r"\s+", qc, fc)
            rows.append(build_row(target, "MetaXcan", region, counts))

    # 4. SMultiXcan -- whitespace separated.
    smx_path = os.path.join(
        tdir, "09.smultixcan_fdrreg", f"{target}.gene.bio.fdrreg.txt")
    qc, fc = QVAL_COLS["SMultiXcan"]
    counts = count_file(smx_path, r"\s+", qc, fc)
    rows.append(build_row(target, "SMultiXcan", "", counts))

    return rows


def main():
    os.makedirs(OUT_DIR, exist_ok=True)
    targets = sorted(
        d for d in os.listdir(BASE)
        if is_target_dir(d, os.path.join(BASE, d))
    )
    print(f"Found {len(targets)} targets: {targets}")

    all_rows = []
    for t in targets:
        all_rows.extend(process_target(t))

    cols = [
        "target", "analysis_type", "region",
        "n_qval_sig_lt0.01", "n_qval_total",
        "n_fdrthe_sig_lt0.01", "n_fdrthe_total",
        "n_risk_loci_qval_bh_0.01", "n_risk_loci_fdrthe_0.01",
    ]
    summary = pd.DataFrame(all_rows, columns=cols)
    summary.to_csv(OUT_FILE, index=False)
    print(f"Saved summary with {len(summary)} rows to:\n{OUT_FILE}")


if __name__ == "__main__":
    main()
