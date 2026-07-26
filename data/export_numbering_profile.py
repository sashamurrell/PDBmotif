"""Export 128-column IMGT position-specific profiles (log-odds + occupancy) built
from ANARCI's bundled germline sequences, for the client-side numbering engine.

The profiles are the only data the JS aligner needs -- no HMMER, no germline files.
For each chain type (H, K, L) we aggregate that type's V and J germlines into a
per-column amino-acid frequency, then store log2 odds vs a uniform background plus
the fraction of germlines occupying each column (used to weight deletion penalties).
"""
import json
import math
import os

from anarci.germlines import all_germlines

AAS = "ACDEFGHIKLMNPQRSTVWY"
NCOL = 128
BG = 0.05
PSEUDO = 0.5

HERE = os.path.dirname(__file__)
OUT = os.path.join(HERE, "..", "docs", "data", "numbering_profile.json")


def build_profile(chain_type):
    counts = [{a: 0.0 for a in AAS} for _ in range(NCOL)]
    occ = [0] * NCOL
    n = 0
    for region in ("V", "J"):
        for genes in all_germlines[region][chain_type].values():
            for seq in genes.values():
                n += 1
                for c in range(NCOL):
                    aa = seq[c]
                    if aa in counts[c]:
                        counts[c][aa] += 1
                        occ[c] += 1
    logodds = []
    for c in range(NCOL):
        tot = sum(counts[c].values())
        row = []
        for a in AAS:
            f = (counts[c][a] + PSEUDO) / (tot + PSEUDO * 20)
            row.append(round(math.log2(f / BG), 4))
        logodds.append(row)
    occupancy = [round(o / n, 4) for o in occ]
    return {"logodds": logodds, "occupancy": occupancy}


def main():
    data = {
        "alphabet": AAS,
        "ncol": NCOL,
        "unknown_score": -4.0,
        "profiles": {ct: build_profile(ct) for ct in ("H", "K", "L")},
    }
    with open(OUT, "w") as f:
        json.dump(data, f, separators=(",", ":"))
    size = os.path.getsize(OUT)
    print(f"Wrote {OUT} ({size/1024:.1f} KB)")


if __name__ == "__main__":
    main()
