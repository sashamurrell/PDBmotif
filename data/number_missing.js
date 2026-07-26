/*
 * Compute CDRs for the SAbDab entries missing from antibodies.json (data/missing_meta.json,
 * produced by prepare_missing.py), using the same client-side engine the app ships
 * (docs/numbering.js). Fetches each structure's sequences from RCSB, numbers the H and L
 * chains, and merges the completed rows into docs/data/antibodies.json.
 *
 *   node data/number_missing.js
 */
const fs = require("fs");
const path = require("path");
const https = require("https");

const REPO = path.join(__dirname, "..");
const Numbering = require(path.join(REPO, "docs/numbering.js"));
Numbering.setProfile(JSON.parse(fs.readFileSync(path.join(REPO, "docs/data/numbering_profile.json"))));

const JSON_PATH = path.join(REPO, "docs/data/antibodies.json");
const META_PATH = path.join(__dirname, "missing_meta.json");
const CONCURRENCY = 10;

function fetchText(url) {
  return new Promise((resolve, reject) => {
    https.get(url, (res) => {
      if (res.statusCode === 404) { res.resume(); return resolve(null); }
      if (res.statusCode !== 200) { res.resume(); return reject(new Error("HTTP " + res.statusCode)); }
      let d = ""; res.on("data", (c) => (d += c)); res.on("end", () => resolve(d));
    }).on("error", reject);
  });
}

// Parse an RCSB FASTA into { authChainId: sequence }. Headers look like
// ">7T76_1|Chains A, G[auth I]|...": the author id (used by SAbDab) is the [auth X] value.
function parseFasta(text) {
  const chains = {};
  if (!text) return chains;
  let header = null, buf = [];
  const flush = () => {
    if (!header) return;
    const seq = buf.join("");
    const parts = header.split("|");
    if (parts.length >= 2 && /Chain/.test(parts[1])) {
      const field = parts[1].replace(/^[^A-Za-z]*Chains?\s*/, "");
      for (let token of field.split(",")) {
        token = token.trim();
        const m = token.match(/(\w+)\s*\[auth\s+(\w+)\]/);
        const cid = m ? m[2] : token;
        if (cid) chains[cid] = seq;
      }
    }
  };
  for (const line of text.split("\n")) {
    if (line.startsWith(">")) { flush(); header = line; buf = []; }
    else buf.push(line.trim());
  }
  flush();
  return chains;
}

// Resolve a SAbDab chain id to its RCSB sequence. SAbDab suffixes scFv/single-chain
// constructs (e.g. VH="A1", VL="A2" both on physical chain "A"); RCSB only lists "A",
// so fall back to the digit-stripped id and let the aligner pick out the right domain.
function chainSeq(fasta, chainId) {
  if (!chainId) return null;
  return fasta[chainId] || fasta[chainId.replace(/\d+$/, "")] || null;
}

function matchCount(seq, ct) {
  try { return Numbering.alignProfile(seq, ct).filter((e) => e[0][1] === "m").length; }
  catch (e) { return 0; }
}

// Number a (possibly scFv) sequence for the requested domain. Accept only a confident,
// complete domain: >=70 profile-match columns AND all three CDRs resolved. This gates out
// non-antibody chains and partial/garbage fits (e.g. the second domain of a diabody).
function cdrsFor(seq, wantHeavy) {
  if (!seq) return null;
  const score = wantHeavy ? matchCount(seq, "H")
                          : Math.max(matchCount(seq, "K"), matchCount(seq, "L"));
  if (score < 70) return null;
  const r = Numbering.numberChain(seq, wantHeavy);
  if (!r) return null;
  const names = wantHeavy ? ["h1", "h2", "h3"] : ["l1", "l2", "l3"];
  const cdrs = Numbering.extractCdrs(r.numbering, names);
  if (names.some((n) => !cdrs[n])) return null;
  return cdrs;
}

async function fetchFasta(pdb) {
  const url = `https://www.rcsb.org/fasta/entry/${pdb.toUpperCase()}`;
  for (let attempt = 0; attempt < 3; attempt++) {
    try { return parseFasta(await fetchText(url)); }
    catch (e) { await new Promise((r) => setTimeout(r, 300 * (attempt + 1))); }
  }
  return null;
}

async function processPdb(pdb, rows) {
  const fasta = await fetchFasta(pdb);
  if (!fasta) return { numbered: 0 }; // gave up -> leave CDRs null, still keep metadata rows
  let numbered = 0;
  for (const row of rows) {
    const h = cdrsFor(chainSeq(fasta, row.hchain), true);
    if (h) { row.cdr_h1 = h.h1; row.cdr_h2 = h.h2; row.cdr_h3 = h.h3; numbered++; }
    const l = cdrsFor(chainSeq(fasta, row.lchain), false);
    if (l) { row.cdr_l1 = l.l1; row.cdr_l2 = l.l2; row.cdr_l3 = l.l3; numbered++; }
  }
  return { numbered };
}

async function main() {
  const existing = JSON.parse(fs.readFileSync(JSON_PATH));
  const existingKeys = new Set(existing.map((r) => `${r.pdb}|${r.hchain || ""}|${r.lchain || ""}`));
  const meta = JSON.parse(fs.readFileSync(META_PATH));

  // Group missing rows by PDB (one fetch per structure).
  const byPdb = new Map();
  for (const row of meta) {
    const key = `${row.pdb}|${row.hchain || ""}|${row.lchain || ""}`;
    if (existingKeys.has(key)) continue; // idempotent: don't re-add
    if (!byPdb.has(row.pdb)) byPdb.set(row.pdb, []);
    byPdb.get(row.pdb).push(row);
  }
  const pdbs = [...byPdb.keys()];
  console.log(`Numbering ${pdbs.length} structures (${meta.length} rows)...`);

  let done = 0, chainsNumbered = 0;
  for (let i = 0; i < pdbs.length; i += CONCURRENCY) {
    const batch = pdbs.slice(i, i + CONCURRENCY);
    const results = await Promise.all(batch.map((p) => processPdb(p, byPdb.get(p))));
    for (const r of results) chainsNumbered += r.numbered;
    done += batch.length;
    if (done % 100 < CONCURRENCY) process.stdout.write(`\r  ${done}/${pdbs.length} structures  (${chainsNumbered} chains numbered)`);
  }
  console.log(`\r  ${done}/${pdbs.length} structures  (${chainsNumbered} chains numbered)`);

  // Merge and write.
  const added = [];
  for (const rows of byPdb.values()) for (const row of rows) added.push(row);
  const merged = existing.concat(added);
  fs.writeFileSync(JSON_PATH, JSON.stringify(merged));

  const withCdr = added.filter((r) => r.cdr_h3 || r.cdr_l3).length;
  console.log(`Added ${added.length} rows (${withCdr} with at least one CDR set).`);
  console.log(`antibodies.json now has ${merged.length} rows.`);
}

main();
