"""Download the SAbDab summary CSV covering every antibody structure in the PDB."""

import os

import requests

SUMMARY_URL = "https://sabdab.opig.stats.ox.ac.uk/api/download/search-summary"
OUT_PATH = os.path.join(os.path.dirname(__file__), "summary_all.csv")


def main():
    print(f"Downloading {SUMMARY_URL}")
    resp = requests.post(SUMMARY_URL, json={}, timeout=120)
    resp.raise_for_status()
    with open(OUT_PATH, "w") as f:
        f.write(resp.text)
    n_rows = resp.text.count("\n") - 1
    print(f"Wrote {OUT_PATH} ({n_rows} rows)")


if __name__ == "__main__":
    main()
