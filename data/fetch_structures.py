"""Download Chothia-renumbered structure files for every unique PDB id in the summary.

Resumable: skips files that already exist and look valid (contain an ATOM record).
Uses a small thread pool to stay polite to the OPIG server.
"""

import csv
import os
import time
from concurrent.futures import ThreadPoolExecutor, as_completed

import requests

HERE = os.path.dirname(__file__)
SUMMARY_PATH = os.path.join(HERE, "summary_all.csv")
OUT_DIR = os.path.join(HERE, "structures")
URL_TEMPLATE = "https://opig.stats.ox.ac.uk/webapps/abdb/entries/{pdb}/structure/chothia/{pdb}.pdb"
MAX_WORKERS = 8
TIMEOUT = 30


def normalise_pdb(raw):
    """Strip the pdb_0000 prefix used in the new SAbDab ID format."""
    raw = raw.strip().lower()
    if raw.startswith("pdb_"):
        raw = raw[4:].lstrip("0") or raw[4:]
    return raw


def unique_pdb_ids():
    with open(SUMMARY_PATH) as f:
        reader = csv.DictReader(f)
        return sorted({normalise_pdb(row["PDB"]) for row in reader if row.get("PDB", "").strip()})


def is_valid(path):
    if not os.path.isfile(path) or os.path.getsize(path) == 0:
        return False
    with open(path) as f:
        return "ATOM" in f.read()


def download_one(pdb_id, session):
    out_path = os.path.join(OUT_DIR, f"{pdb_id}.pdb")
    if is_valid(out_path):
        return pdb_id, "skipped"
    url = URL_TEMPLATE.format(pdb=pdb_id)
    try:
        resp = session.get(url, timeout=TIMEOUT)
        if resp.status_code != 200 or "ATOM" not in resp.text:
            return pdb_id, f"failed (status={resp.status_code})"
        with open(out_path, "w") as f:
            f.write(resp.text)
        return pdb_id, "downloaded"
    except requests.RequestException as exc:
        return pdb_id, f"failed ({exc})"


def main():
    os.makedirs(OUT_DIR, exist_ok=True)
    pdb_ids = unique_pdb_ids()
    print(f"{len(pdb_ids)} unique PDB ids to fetch")

    counts = {"downloaded": 0, "skipped": 0}
    failures = []
    session = requests.Session()

    with ThreadPoolExecutor(max_workers=MAX_WORKERS) as pool:
        futures = {pool.submit(download_one, pdb_id, session): pdb_id for pdb_id in pdb_ids}
        done = 0
        start = time.time()
        for future in as_completed(futures):
            pdb_id, status = future.result()
            done += 1
            if status.startswith("failed"):
                failures.append((pdb_id, status))
            else:
                counts[status] += 1
            if done % 250 == 0 or done == len(pdb_ids):
                elapsed = time.time() - start
                print(f"  {done}/{len(pdb_ids)} processed "
                      f"({counts['downloaded']} downloaded, {counts['skipped']} skipped, "
                      f"{len(failures)} failed) [{elapsed:.0f}s]")

    print(f"\nDone. downloaded={counts['downloaded']} skipped={counts['skipped']} "
          f"failed={len(failures)}")
    if failures:
        log_path = os.path.join(HERE, "fetch_failures.log")
        with open(log_path, "w") as f:
            for pdb_id, status in failures:
                f.write(f"{pdb_id}\t{status}\n")
        print(f"Failure list written to {log_path}")


if __name__ == "__main__":
    main()
