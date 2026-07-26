"""Export the antibodies SQLite database to the static JSON file used by the viewer."""

import json
import os
import sqlite3

HERE = os.path.dirname(__file__)
DB_PATH = os.path.join(HERE, "antibodies.db")
JSON_PATH = os.path.join(HERE, "..", "docs", "data", "antibodies.json")

COLUMNS = [
    "pdb", "hchain", "lchain", "antigen_chain", "antigen_type", "antigen_name",
    "organism", "heavy_species", "light_species", "antigen_species", "method",
    "engineered", "scfv", "heavy_subclass", "light_subclass", "light_ctype",
    "resolution", "date", "compound", "pmid",
    "cdr_h1", "cdr_h2", "cdr_h3", "cdr_l1", "cdr_l2", "cdr_l3",
]


def main():
    conn = sqlite3.connect(DB_PATH)
    conn.row_factory = sqlite3.Row
    rows = conn.execute(f"SELECT {', '.join(COLUMNS)} FROM antibodies").fetchall()
    conn.close()

    data = [dict(row) for row in rows]
    with open(JSON_PATH, "w") as f:
        json.dump(data, f, separators=(",", ":"))

    print(f"Wrote {len(data)} rows to {JSON_PATH}")


if __name__ == "__main__":
    main()
