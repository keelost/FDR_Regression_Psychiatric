#!/usr/bin/env python3
"""
Build a master rsid -> varID mapping from all GTEx v8 mashr .db models.
The mapping is consistent across tissues (same GTEx genotype reference),
so taking the union over all models gives the most complete table.
rsIDs that map to more than one varID are ambiguous and are dropped.
"""

import glob
import os
import sqlite3
import sys

# Directory containing all mashr_*.db files
DB_DIR = os.environ.get("FDRREG_V8_MODEL_DIR", "")
OUT_FILE = os.environ.get("FDRREG_MAP_FILE", "rsid_to_varid.tsv")


def main():
    db_files = sorted(glob.glob(os.path.join(DB_DIR, "mashr_*.db")))
    if not db_files:
        sys.exit("No mashr_*.db files found in {}".format(DB_DIR))

    # rsid -> set of (varID, ref_allele, eff_allele)
    mapping = {}

    for db in db_files:
        con = sqlite3.connect(db)
        cur = con.cursor()
        cur.execute(
            "SELECT DISTINCT rsid, varID, ref_allele, eff_allele FROM weights"
        )
        for rsid, varid, ref, eff in cur.fetchall():
            if rsid is None or varid is None:
                continue
            mapping.setdefault(rsid, set()).add((varid, ref, eff))
        con.close()
        print("Processed {}".format(os.path.basename(db)))

    n_total = len(mapping)
    n_ambiguous = sum(1 for v in mapping.values() if len(v) > 1)

    with open(OUT_FILE, "w") as out:
        out.write("rsid\tvarID\tref_allele\teff_allele\n")
        for rsid, entries in mapping.items():
            if len(entries) != 1:
                continue  # drop ambiguous rsIDs (multi-allelic / conflicting)
            varid, ref, eff = next(iter(entries))
            out.write("{}\t{}\t{}\t{}\n".format(rsid, varid, ref, eff))

    print("Total unique rsIDs: {}".format(n_total))
    print("Ambiguous rsIDs dropped: {}".format(n_ambiguous))
    print("Mapping written to: {}".format(OUT_FILE))


if __name__ == "__main__":
    main()
