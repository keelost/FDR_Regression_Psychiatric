#!/usr/bin/env python3
"""
Add a panel_variant_id column (GTEx varID) to a GWAS file by mapping its
rsID column through the master rsid -> varID table.
Rows whose rsID is not in the map are dropped (they cannot match the models).
Allele orientation is left to S-PrediXcan, which aligns a1/a2 against the
model ref/eff alleles, so no manual z-score flipping is done here.
"""

import argparse
import sys

import pandas as pd


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gwas", required=True, help="Input GWAS file")
    ap.add_argument("--map", required=True, help="rsid_to_varid.tsv from script one")
    ap.add_argument("--snp_column", default="snpid", help="rsID column in GWAS")
    ap.add_argument("--sep", default=r"\s+", help="GWAS field separator regex")
    ap.add_argument("--out", required=True, help="Output GWAS file")
    args = ap.parse_args()

    # Load GWAS
    gwas = pd.read_csv(args.gwas, sep=args.sep, engine="python")
    if args.snp_column not in gwas.columns:
        sys.exit("Column '{}' not found in GWAS".format(args.snp_column))
    n_in = len(gwas)

    # Load mapping (only need rsid -> varID)
    m = pd.read_csv(args.map, sep="\t", usecols=["rsid", "varID"])
    m = m.rename(columns={"rsid": args.snp_column, "varID": "panel_variant_id"})

    # Merge
    merged = gwas.merge(m, on=args.snp_column, how="inner")
    n_out = len(merged)

    # Strip any stray whitespace from object (string) columns to avoid
    # broken field alignment in the output file.
    for col in merged.select_dtypes(include=["object"]).columns:
        merged[col] = merged[col].astype(str).str.strip()

    # Write strictly tab-separated output with no extra whitespace.
    merged.to_csv(args.out, sep="\t", index=False, na_rep="NA")

    print("Input GWAS rows: {}".format(n_in))
    print("Mapped (kept) rows: {}".format(n_out))
    print("Unmapped (dropped) rows: {}".format(n_in - n_out))
    print("Output written to: {}".format(args.out))


if __name__ == "__main__":
    main()
