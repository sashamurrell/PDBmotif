"""Parse SAbDab summary + Chothia-renumbered structures into a SQLite database.

For each structure file, reads the REMARK 5 PAIRED_HL lines to find H/L/antigen chain
triples, extracts the 6 CDR sequences per H/L pairing from the Chothia-numbered ATOM
records, and joins with the matching summary row(s) for metadata.
"""

import csv
import os
import sqlite3

from cdr_ranges import CDR_RANGES

HERE = os.path.dirname(__file__)
SUMMARY_PATH = os.path.join(HERE, "summary_all.tsv")
STRUCTURES_DIR = os.path.join(HERE, "structures")
DB_PATH = os.path.join(HERE, "antibodies.db")

AA3to1 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q",
    "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
    "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
    "TYR": "Y", "VAL": "V", "MSE": "M", "SEC": "U", "PYL": "O",
}

NA_VALUES = {"", "NA", "None"}


def parse_paired_hl(lines):
    """Return list of dicts: {hchain, lchain, agchains} from REMARK 5 PAIRED_HL/SINGLE lines.

    PAIRED_HL covers Fab/Fv/scFv with both H and L chains; SINGLE covers VHH/nanobody
    (HCHAIN only) or VL-only constructs (LCHAIN only) -- the missing chain defaults to "NA".
    """
    pairs = []
    for line in lines:
        tokens = line.split()
        if len(tokens) < 3 or tokens[2] not in ("PAIRED_HL", "SINGLE"):
            continue
        fields = {}
        for token in tokens[3:]:
            if "=" in token:
                key, val = token.split("=", 1)
                fields[key] = val
        agchain = fields.get("AGCHAIN", "NONE")
        agchains = [] if agchain == "NONE" else agchain.split(";")
        pairs.append({
            "hchain": fields.get("HCHAIN", "NA"),
            "lchain": fields.get("LCHAIN", "NA"),
            "agchains": agchains,
        })
    return pairs


def chain_residues(atom_lines, chain_id):
    """Return ordered list of (resnum_int, inscode, resname1) for a chain's ATOM records,
    one entry per residue (deduplicated across atoms of the same residue).

    Matches chain id case-insensitively: REMARK 5 lines always report HCHAIN/LCHAIN in
    uppercase, but the actual ATOM chain id can be lowercase (e.g. PDB 11OV uses chains
    'h'/'l'), so a case-sensitive match would silently find zero residues.
    """
    seen = set()
    residues = []
    for line in atom_lines:
        if line[21].upper() != chain_id.upper():
            continue
        resnum_str = line[22:27].strip()  # e.g. "100A" or "52"
        resname3 = line[17:20].strip()
        key = (line[21], resnum_str)
        if key in seen:
            continue
        seen.add(key)
        # split numeric part from insertion code
        i = len(resnum_str)
        while i > 0 and not resnum_str[i - 1].isdigit():
            i -= 1
        resnum = int(resnum_str[:i])
        inscode = resnum_str[i:]
        aa = AA3to1.get(resname3, "X")
        residues.append((resnum, inscode, aa))
    return residues


def extract_cdr(residues, cdr_name):
    start, end = CDR_RANGES[cdr_name]
    seq = "".join(aa for resnum, inscode, aa in residues if start <= resnum <= end)
    return seq or None


def load_summary():
    """Return dict keyed by (pdb, hchain.upper(), lchain.upper()) -> row dict.

    Chain letters are normalized to uppercase because REMARK 5 lines in the structure
    files always report HCHAIN/LCHAIN in uppercase regardless of the summary's casing
    (the summary uses case to distinguish same-letter scFv H/L domains, e.g.
    Hchain='A'/Lchain='a' for a single physical chain 'A' holding both domains).
    """
    by_key = {}
    with open(SUMMARY_PATH) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            pdb = row["pdb"].strip().lower()
            hchain = row["Hchain"].strip().upper()
            lchain = row["Lchain"].strip().upper()
            by_key[(pdb, hchain, lchain)] = row
    return by_key


def coerce_resolution(value):
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def coerce_bool(value):
    return value.strip().lower() == "true"


def normalize_na(value):
    value = (value or "").strip()
    return None if value in NA_VALUES else value


# SAbDab's species fields sometimes record the recombinant expression host (e.g.
# "escherichia coli k-12" for a humanized Fab actually expressed there for
# crystallography) rather than the antibody's real immunological origin, alongside
# outright data artifacts (plant/algae/insect names, "synthetic construct", etc).
# Only keep comma-separated parts whose first word is a known antibody-bearing
# vertebrate genus/family/order/class, built from what's actually observed in this
# dataset -- not an exhaustive taxonomy, so a genuinely new vertebrate species absent
# from this list would also get dropped.
PLAUSIBLE_SPECIES_TOKENS = {
    "homo", "mus", "rattus", "oryctolagus", "camelus", "lama", "vicugna", "camelidae",
    "gallus", "bos", "ovis", "sus", "canis", "equus", "mesocricetus", "cricetulus",
    "macaca", "pan", "platyrrhini", "mammalia", "chiloscyllium", "ginglymostoma",
    "okamejei", "orectolobus", "squalus", "myodes", "ondatra", "murinae",
}


def clean_species(value):
    """Drop non-vertebrate / expression-host parts, keeping genuine co-annotated origin
    species when present (e.g. "homo sapiens, synthetic construct" -> "homo sapiens").
    Returns None if no part is a recognizable vertebrate (no real origin recorded)."""
    if not value:
        return None
    kept = [p.strip() for p in value.split(",") if p.strip().split()[0].lower() in PLAUSIBLE_SPECIES_TOKENS]
    return ", ".join(kept) if kept else None


METADATA_COLUMNS = [
    "antigen_type", "antigen_name", "organism", "heavy_species", "light_species",
    "antigen_species", "method", "engineered", "scfv", "heavy_subclass",
    "light_subclass", "light_ctype", "date", "compound", "pmid",
]


def build_metadata_row(summary_row):
    if summary_row is None:
        return {col: None for col in METADATA_COLUMNS} | {"resolution": None}
    meta = {col: normalize_na(summary_row.get(col)) for col in METADATA_COLUMNS}
    meta["heavy_species"] = clean_species(meta["heavy_species"])
    meta["light_species"] = clean_species(meta["light_species"])
    meta["engineered"] = coerce_bool(summary_row.get("engineered", ""))
    meta["scfv"] = coerce_bool(summary_row.get("scfv", ""))
    meta["resolution"] = coerce_resolution(summary_row.get("resolution"))
    return meta


def create_schema(conn):
    columns = ["pdb TEXT", "hchain TEXT", "lchain TEXT", "antigen_chain TEXT"]
    columns += [f"{c} TEXT" if c not in ("engineered", "scfv", "resolution") else
                (f"{c} INTEGER" if c in ("engineered", "scfv") else f"{c} REAL")
                for c in METADATA_COLUMNS + ["resolution"]]
    columns += [f"cdr_{name} TEXT" for name in CDR_RANGES]
    conn.execute(f"CREATE TABLE antibodies ({', '.join(columns)})")
    for name in CDR_RANGES:
        conn.execute(f"CREATE INDEX idx_cdr_{name} ON antibodies(cdr_{name})")
    for col in ("heavy_species", "light_species", "organism", "antigen_type", "method"):
        conn.execute(f"CREATE INDEX idx_{col} ON antibodies({col})")


def main():
    summary_by_key = load_summary()
    print(f"Loaded {len(summary_by_key)} summary rows")

    pdb_files = sorted(f for f in os.listdir(STRUCTURES_DIR) if f.endswith(".pdb"))
    print(f"Found {len(pdb_files)} structure files")

    if os.path.exists(DB_PATH):
        os.remove(DB_PATH)
    conn = sqlite3.connect(DB_PATH)
    create_schema(conn)

    n_rows = 0
    n_parse_failures = 0
    n_no_summary_match = 0
    n_scfv_same_chain = 0

    insert_cols = (
        ["pdb", "hchain", "lchain", "antigen_chain"]
        + METADATA_COLUMNS + ["resolution"]
        + [f"cdr_{name}" for name in CDR_RANGES]
    )
    placeholders = ", ".join("?" for _ in insert_cols)
    insert_sql = f"INSERT INTO antibodies ({', '.join(insert_cols)}) VALUES ({placeholders})"

    for filename in pdb_files:
        pdb_id = filename[:-4]
        path = os.path.join(STRUCTURES_DIR, filename)
        try:
            with open(path, errors="replace") as f:
                lines = f.readlines()
        except OSError:
            n_parse_failures += 1
            continue

        remark_lines = [l for l in lines if l.startswith("REMARK   5")]
        atom_lines = [l for l in lines if l.startswith("ATOM")]
        pairs = parse_paired_hl(remark_lines)
        if not pairs:
            n_parse_failures += 1
            continue

        for pair in pairs:
            hchain, lchain = pair["hchain"].upper(), pair["lchain"].upper()
            # scFv constructs sometimes fuse H and L domains into one physical chain
            # (HCHAIN==LCHAIN). Chothia renumbering doesn't reliably restart at the
            # domain boundary in that case, so we can't safely tell which residues
            # belong to which domain -- skip CDR extraction rather than risk silently
            # wrong sequences. Metadata is still recorded and searchable.
            same_chain_scfv = hchain != "NA" and lchain != "NA" and hchain == lchain
            cdrs = {}
            try:
                if same_chain_scfv:
                    n_scfv_same_chain += 1
                    cdrs = {name: None for name in CDR_RANGES}
                else:
                    if hchain != "NA":
                        h_residues = chain_residues(atom_lines, hchain)
                        for name in ("h1", "h2", "h3"):
                            cdrs[name] = extract_cdr(h_residues, name)
                    else:
                        cdrs["h1"] = cdrs["h2"] = cdrs["h3"] = None
                    if lchain != "NA":
                        l_residues = chain_residues(atom_lines, lchain)
                        for name in ("l1", "l2", "l3"):
                            cdrs[name] = extract_cdr(l_residues, name)
                    else:
                        cdrs["l1"] = cdrs["l2"] = cdrs["l3"] = None
            except Exception:
                n_parse_failures += 1
                continue

            summary_row = summary_by_key.get((pdb_id, hchain, lchain))
            if summary_row is None:
                n_no_summary_match += 1
            meta = build_metadata_row(summary_row)

            values = (
                [pdb_id, hchain, lchain, ";".join(pair["agchains"]) or None]
                + [meta[c] for c in METADATA_COLUMNS] + [meta["resolution"]]
                + [cdrs[name] for name in CDR_RANGES]
            )
            conn.execute(insert_sql, values)
            n_rows += 1

    conn.commit()
    conn.close()

    print(f"\nWrote {n_rows} rows to {DB_PATH}")
    print(f"Structure parse failures: {n_parse_failures}")
    print(f"Rows with no matching summary row: {n_no_summary_match}")
    print(f"scFv same-physical-chain rows (CDRs skipped, metadata kept): {n_scfv_same_chain}")


if __name__ == "__main__":
    main()
