"""Emit metadata rows for SAbDab summary entries that are absent from antibodies.json.

These are structures SAbDab lists but never produced a Chothia file for (e.g. 9PY5),
so build_db.py could not include them. We reuse build_db's summary parsing and metadata
cleaning here, then hand off to number_missing.js which computes the CDRs client-side
from RCSB sequences. Output: data/missing_meta.json (metadata only, cdr_* left null).
"""
import csv
import json
import os
import sys

HERE = os.path.dirname(__file__)
sys.path.insert(0, HERE)
import build_db  # noqa: E402

JSON_PATH = os.path.join(HERE, "..", "docs", "data", "antibodies.json")
OUT_PATH = os.path.join(HERE, "missing_meta.json")

CDR_COLS = ["cdr_h1", "cdr_h2", "cdr_h3", "cdr_l1", "cdr_l2", "cdr_l3"]


def norm_antigen(value):
    """Summary antigen_chain uses '|' as separator; the JSON schema uses ';'."""
    value = build_db.normalize_na(value)
    return value.replace("|", ";") if value else None


def main():
    existing = {
        (r["pdb"], (r["hchain"] or "").upper(), (r["lchain"] or "").upper())
        for r in json.load(open(JSON_PATH))
    }

    rows = []
    seen = set()
    with open(build_db.SUMMARY_PATH) as f:
        for raw in csv.DictReader(f):
            pdb = build_db.normalise_pdb(raw["PDB"])
            hchain = raw["Hchain"].strip().upper()
            lchain = raw["Lchain"].strip().upper()
            key = (pdb, hchain, lchain)
            if key in existing or key in seen:
                continue
            seen.add(key)

            meta = build_db.build_metadata_row(raw)
            row = {
                "pdb": pdb,
                "hchain": None if hchain in build_db.NA_VALUES else hchain,
                "lchain": None if lchain in build_db.NA_VALUES else lchain,
                "antigen_chain": norm_antigen(raw.get("antigen_chain")),
            }
            for col in build_db.METADATA_COLUMNS:
                row[col] = meta[col]
            row["resolution"] = meta["resolution"]
            for c in CDR_COLS:
                row[c] = None
            rows.append(row)

    with open(OUT_PATH, "w") as f:
        json.dump(rows, f)
    print(f"Wrote {len(rows)} missing metadata rows to {OUT_PATH}")


if __name__ == "__main__":
    main()
