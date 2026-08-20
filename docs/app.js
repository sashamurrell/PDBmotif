const AA3TO1 = {
  ALA:"A",ARG:"R",ASN:"N",ASP:"D",CYS:"C",GLN:"Q",GLU:"E",GLY:"G",HIS:"H",
  ILE:"I",LEU:"L",LYS:"K",MET:"M",PHE:"F",PRO:"P",SER:"S",THR:"T",TRP:"W",
  TYR:"Y",VAL:"V",MSE:"M",SEC:"U",PYL:"O",
};

const EPITOPE_DISTANCE = 4.5;

const FILTER_SELECTS = {
  heavy_species: "heavy-species-select",
  light_species: "light-species-select",
  organism: "organism-select",
  antigen_type: "antigen-type-select",
  method: "method-select",
};

const RESULT_LIMIT = 500;

let viewer = null;
let allData = null;
let lastResults = [];
let lastAllResults = [];
let lastTotalCount = 0;
let lastCdr = "h3";
let lastMotif = "";
let profileReady = false;

function $(id) { return document.getElementById(id); }

async function loadData() {
  const resp = await fetch("data/antibodies.json");
  allData = await resp.json();
  fetch("data/numbering_profile.json")
    .then((r) => r.json())
    .then((p) => { Numbering.setProfile(p); profileReady = true; });
  buildFilterOptions();
  $("heavy-species-select").value = "homo sapiens";
  $("light-species-select").value = "homo sapiens";
  $("complex-checkbox").checked = true;
  $("collapse-duplicates-checkbox").checked = true;
}

function buildFilterOptions() {
  for (const [col, selectId] of Object.entries(FILTER_SELECTS)) {
    const values = [...new Set(allData.map(r => r[col]).filter(Boolean))].sort();
    const select = $(selectId);
    for (const value of values) {
      const opt = document.createElement("option");
      opt.value = value;
      opt.textContent = value;
      select.appendChild(opt);
    }
  }
}

// Convert SQL LIKE pattern (with _ and % wildcards) to a JS RegExp.
function motifToRegex(motif) {
  const escaped = motif.replace(/[.+^${}()|[\]\\]/g, "\\$&");
  return new RegExp(escaped.replace(/%/g, ".*").replace(/_/g, "."));
}

function isUnpublished(row) {
  if (row.pmid === "TBD") return true;
  if (!row.pmid && row.date) {
    const yr = row.date.split("/").pop();
    return yr === "25" || yr === "26";
  }
  return false;
}

// High-confidence antigen-name markers -> source organism + category. Many antigens are
// deposited as synthetic-peptide constructs with a null antigen_species, so the pathogen can
// only be recovered from the protein name itself (e.g. "gp41 MPER Peptide" is HIV with no
// species field). Used both to categorise and, when antigen_species is missing, to prefix the
// target with the organism. Ordered — first match wins. Kept high-precision: nothing non-HIV is
// named gp120/gp160/MPER, etc.
// `speciesRe` recognises a deposited antigen_species token that already names this pathogen, so
// the more specific real species is kept (e.g. "simian immunodeficiency virus" stays SIV, not
// "HIV-1"); otherwise the species field is the expression host/scaffold and the target uses the
// marker organism.
const ANTIGEN_MARKERS = [
  { re: /\bgp41\b|\bgp120\b|\bgp160\b|\bmper\b|envelope glycoprotein gp/i, organism: "HIV-1", category: "virus", speciesRe: /immunodeficiency|\bhiv\b|\bsiv\b|lentivir/i },
  { re: /neuraminidase|hemagglutinin/i, organism: "Influenza", category: "virus", speciesRe: /influenza [ab]|influenza virus|orthomyxo/i },
  { re: /spike glycoprotein|spike protein s/i, organism: "Coronavirus", category: "virus", speciesRe: /coronavirus|sars|\bmers\b|betacoronavirus/i },
  { re: /circumsporozoite/i, organism: "Plasmodium", category: "parasite", speciesRe: /plasmodium|malaria/i },
];
function antigenMarker(row) {
  const name = cleanAntigenName(row.antigen_name);
  return ANTIGEN_MARKERS.find((m) => m.re.test(name)) || null;
}

// Categorise the antigen as self (host), a pathogen class, or a small molecule.
// A high-confidence pathogen protein marker (gp120, spike, circumsporozoite, …) wins over
// everything, because antigen_species is frequently the *expression host* — e.g. an HIV gp120
// immunogen produced in human cells is tagged "homo sapiens" and would otherwise mis-fire the
// self test (and scaffold-fusion antigens mis-fire bacteria/fungal). Otherwise: self, then the
// keyword scans over `antigen_species + antigen_name` (so null-species pathogens like "Endogenous
// retrovirus …" still resolve). "influenza\b" avoids matching the bacterium H. influenzae, and
// "adenovirus" (not bare "adeno") avoids adenosine/adenocarcinoma. Bacteria uses a "bacter"
// catch-all plus a curated genus list from the data. Falls back to "small molecule" for
// hapten-only antigens, "other" for any remaining real antigen, and "" (apo) when there's none.
function categorizeAntigen(row) {
  const marker = antigenMarker(row);
  if (marker) return marker.category;
  const sp = (row.antigen_species || "").split("|")[0].trim().toLowerCase();
  const hsp = (row.heavy_species || "").toLowerCase().trim();
  const lsp = (row.light_species || "").toLowerCase().trim();
  if (sp && ((hsp && sp === hsp) || (lsp && sp === lsp))) return "self";
  const text = `${sp} ${cleanAntigenName(row.antigen_name)}`.toLowerCase();
  if (/virus|phage|viral|coronavirus|hiv|siv|influenza\b|hepatitis|dengue|zika|rabies|ebola|herpes|adenovirus|retrovirus|sars|mers|norovirus|rotavirus/.test(text)) return "virus";
  if (/plasmodium|trypanosoma|leishmania|toxoplasma|malaria|babesia|theileria|cryptosporidium|giardia|schistosoma/.test(text)) return "parasite";
  if (/bacter|bacill|coccus|streptococ|staphyloco|streptomyces|escherichia|salmonella|shigella|klebsiella|pseudomonas|mycoplasma|clostridi|listeria|yersinia|neisseria|bordetella|vibrio|thermotoga|thermus|aquifex|borreli|porphyromonas|shewanella|ruegeria|amycolatopsis|saccharopolyspora|dickeya|anaplasma|haemophilus|enterococc|burkholderia|francisella|brucella|legionella|chlamydia|treponema|leptospira|rickettsia|coxiella|serratia|prevotella|lactobacill|lactococc|bacteroides|deinococcus|geobacill|paenibacill|xanthomonas|agrobacter|rhizobium|moraxella/.test(text)) return "bacteria";
  if (/aspergillus|candida|cryptococcus|saccharomyces|fungal|fungus|pichia|histoplasma|blastomyces|coccidioides|pneumocystis|malassezia|trichophyton/.test(text)) return "fungal";
  if (/synthetic|artificial|de novo|designed/.test(text)) return "synthetic";
  const at = (row.antigen_type || "").toUpperCase();
  if (!/PROTEIN|PEPTIDE|SUGAR|RNA|DNA/.test(at) && /HAPTEN/.test(at)) return "small molecule";
  if (/PROTEIN|PEPTIDE|SUGAR|RNA|DNA/.test(at) || cleanAntigenName(row.antigen_name)) return "other";
  return "";
}

// Components of antigen_name (joined by "|") that are crystallisation/buffer additives rather
// than the antigen itself: ions ("X ION"), buffers, cryo-solvents/precipitants, reductants, and
// the N-linked glycan sugars on glycoproteins (IUPAC names always contain "pyranose"/"furanose").
// These are dropped from the displayed target. Ions are matched by a trailing " ION" so genuine
// proteins such as "…ligand-gated ion channel" are preserved; everything else is matched by
// EXACT name so genuine haptens are kept — e.g. "Glycerol-3-phosphate dehydrogenase" survives
// GLYCEROL, and "…BENZIMIDAZOLE…" survives IMIDAZOLE.
const ADDITIVE_NAMES = new Set([
  // buffers
  "IMIDAZOLE", "TRIS", "BIS-TRIS", "BIS-TRIS PROPANE", "HEPES", "MES", "MOPS", "PIPES",
  "2-AMINO-2-HYDROXYMETHYL-PROPANE-1,3-DIOL", "2-(N-MORPHOLINO)-ETHANESULFONIC ACID",
  // cryo-solvents / precipitants
  "GLYCEROL", "1,2-ETHANEDIOL", "ETHYLENE GLYCOL", "DI(HYDROXYETHYL)ETHER", "TRIETHYLENE GLYCOL",
  "TETRAETHYLENE GLYCOL", "PENTAETHYLENE GLYCOL", "HEXAETHYLENE GLYCOL", "POLYETHYLENE GLYCOL",
  "2-METHYL-2,4-PENTANEDIOL", "(4R)-2-METHYLPENTANE-2,4-DIOL", "(4S)-2-METHYLPENTANE-2,4-DIOL",
  "DIMETHYL SULFOXIDE", "SUCROSE", "TREHALOSE", "BETAINE", "GLYCINE BETAINE", "UREA", "GUANIDINE",
  "FORMAMIDE", "ETHANOL", "METHANOL", "ISOPROPYL ALCOHOL", "2-PROPANOL", "1,4-DIOXANE", "ACETONE", "ACETONITRILE",
  // acids / reductants / misc additives
  "ACETIC ACID", "FORMIC ACID", "CITRIC ACID", "SUCCINIC ACID", "TARTARIC ACID", "MALONIC ACID", "MALIC ACID",
  "EDTA", "ETHYLENEDIAMINETETRAACETIC ACID", "DITHIOTHREITOL", "1,4-DITHIOTHREITOL",
  "2,3-DIHYDROXY-1,4-DITHIOBUTANE", "2-MERCAPTOETHANOL", "BETA-MERCAPTOETHANOL",
  "SPERMIDINE", "SPERMINE", "PUTRESCINE", "BENZAMIDINE", "PHENOL", "WATER",
]);
function isSolventComponent(part) {
  return / ion$/i.test(part) || /pyranose|furanose/i.test(part) || ADDITIVE_NAMES.has(part.toUpperCase());
}
function cleanAntigenName(name) {
  if (!name) return "";
  const parts = name.split("|").map(s => s.trim()).filter(Boolean);
  return parts.filter(p => !isSolventComponent(p)).join(", ");
}

function formatTarget(row) {
  const name = cleanAntigenName(row.antigen_name);
  if (categorizeAntigen(row) === "self") return name || "n/a";
  const marker = antigenMarker(row);
  let species;
  if (marker) {
    // A pathogen marker fired: keep the deposited species only if it already names that pathogen
    // (more specific), otherwise it's the expression host/scaffold — use the marker organism.
    const tokens = (row.antigen_species || "").split("|").map((s) => s.trim()).filter(Boolean);
    species = tokens.find((t) => marker.speciesRe.test(t)) || marker.organism;
  } else {
    species = row.antigen_species ? row.antigen_species.split("|")[0].trim() : "";
  }
  // Avoid "HIV-1 HIV-1 GP120 …": if the name already mentions the organism, don't prefix it.
  if (species && name) {
    return name.toLowerCase().includes(species.toLowerCase()) ? name : `${species} ${name}`;
  }
  return name || species || (row.antigen_chain == null ? "n/a" : "");
}

function highlightMotif(sequence, motif) {
  if (!sequence) return "";
  if (!motif) return sequence;
  const idx = sequence.indexOf(motif);
  if (idx === -1) return sequence;
  return (
    sequence.slice(0, idx) +
    "<mark>" + sequence.slice(idx, idx + motif.length) + "</mark>" +
    sequence.slice(idx + motif.length)
  );
}

function dedupeByPdb(results) {
  const byPdb = new Map();
  for (const row of results) {
    if (byPdb.has(row.pdb)) byPdb.get(row.pdb).count += 1;
    else byPdb.set(row.pdb, { row, count: 1 });
  }
  return Array.from(byPdb.values());
}

function renderResults(data, cdr, motif) {
  lastAllResults = data.results;
  lastTotalCount = data.total_count;
  lastCdr = cdr;
  lastMotif = motif;
  renderTable();
}

function renderTable() {
  const cdr = lastCdr;
  const motif = lastMotif;
  const collapse = $("collapse-duplicates-checkbox").checked;

  const entries = collapse
    ? dedupeByPdb(lastAllResults)
    : lastAllResults.map((row) => ({ row, count: 1 }));
  lastResults = entries.map((e) => e.row);

  const summary = $("results-summary");
  summary.textContent = collapse
    ? `Showing ${entries.length} unique structures (from ${lastAllResults.length} of ${lastTotalCount} chain matches)`
    : `Showing ${lastAllResults.length} of ${lastTotalCount} results`;

  const body = $("results-body");
  body.innerHTML = "";
  entries.forEach(({ row, count }, i) => {
    const tr = document.createElement("tr");
    tr.dataset.index = i;
    const cdrSeq = row[`cdr_${cdr}`] || "";
    const dupBadge = count > 1 ? `<span class="dup-badge">&times;${count}</span>` : "";
    tr.innerHTML = `
      <td><a href="https://www.rcsb.org/structure/${row.pdb}" target="_blank" rel="noopener">${row.pdb}</a>${dupBadge}${isUnpublished(row) ? '<span class="tbp-badge">TBP</span>' : ""}</td>
      <td>${row.hchain}/${row.lchain}</td>
      <td>${highlightMotif(cdrSeq, motif)}</td>
      <td>${row.heavy_species || ""}</td>
      <td>${categorizeAntigen(row)}</td>
      <td>${formatTarget(row)}</td>
      <td>${row.resolution ?? ""}</td>
      <td>${row.method || ""}</td>
    `;
    tr.addEventListener("click", () => selectRow(i, cdr));
    body.appendChild(tr);
  });
}

function runSearch() {
  if (!allData) return;

  const cdr = $("cdr-select").value.toLowerCase();
  const motif = $("motif-input").value.trim().toUpperCase();

  let results = allData;

  if (motif) {
    const re = motifToRegex(`%${motif}%`);
    results = results.filter(r => re.test(r[`cdr_${cdr}`] || ""));
  }

  for (const [col, selectId] of Object.entries(FILTER_SELECTS)) {
    const value = $(selectId).value;
    if (value) results = results.filter(r => r[col] === value);
  }

  const resMax = parseFloat($("resolution-max-input").value);
  if (!isNaN(resMax)) results = results.filter(r => r.resolution != null && r.resolution <= resMax);

  if ($("engineered-checkbox").checked) results = results.filter(r => r.engineered === 1);
  if ($("complex-checkbox").checked) results = results.filter(r => r.antigen_chain);

  const total_count = results.length;
  const limited = results.slice(0, RESULT_LIMIT);

  lastMotif = motif;
  renderResults({ results: limited, total_count }, cdr, motif);
}


// Parse the _atom_site loop of an mmCIF string and return per-atom metadata
// ({ comp, seq, icode }) in file order. Returns null if the loop or the required
// columns can't be found. Used to recover author numbering (see fixCifResidueNumbering).
function parseCifAtomMeta(text) {
  const cols = [];
  let icodeCol = -1, compCol = -1, seqCol = -1, sawHeader = false, inData = false;
  const meta = [];
  for (const raw of text.split("\n")) {
    const t = raw.trim();
    if (t.startsWith("_atom_site.")) {
      if (t === "_atom_site.pdbx_PDB_ins_code") icodeCol = cols.length;
      if (t === "_atom_site.label_comp_id") compCol = cols.length;
      if (t === "_atom_site.auth_seq_id") seqCol = cols.length;
      cols.push(t);
      sawHeader = true;
      continue;
    }
    if (sawHeader && !inData) {
      inData = true;
      if (icodeCol === -1 || compCol === -1 || seqCol === -1) return null;
    }
    if (!inData) continue;
    if (raw.startsWith("ATOM") || raw.startsWith("HETATM")) {
      const f = t.split(/\s+/);
      let ic = f[icodeCol];
      if (ic === "?" || ic === "." || ic == null) ic = "";
      meta.push({ comp: f[compCol], seq: f[seqCol], icode: ic });
    } else if (meta.length && (t === "#" || t === "" || t === "loop_" || t.startsWith("_"))) {
      break; // end of the atom_site data block
    }
  }
  return meta.length ? meta : null;
}

// 3Dmol's mmCIF parser ignores pdbx_PDB_ins_code, so every residue in an insertion-coded
// loop (e.g. a CDR-H3 numbered 100A..100T, as in 7T77) is loaded with the same resi (100)
// and no icode. Two consequences: hover labels can't distinguish them, and — because 3Dmol's
// cartoon builder groups residues by resi — the whole block collapses to one residue and the
// ribbon is never drawn through the loop. Re-parse the raw CIF, give every residue a unique
// resi (per chain, in file order) so the cartoon renders, restore the icode, and stash the
// true author label (e.g. "100A") on each atom for the hover readout. Atoms are matched
// positionally (both are in atom_site order); we bail without mutating if counts or residue
// names don't line up, and no-op for PDB input whose parser already reads column-27.
function fixCifResidueNumbering(viewer, text, fmt) {
  if (fmt !== "cif") return false;
  const meta = parseCifAtomMeta(text);
  if (!meta) return false;
  const atoms = viewer.selectedAtoms({});
  if (atoms.length !== meta.length) return false;
  for (let i = 0; i < atoms.length; i++) {
    if (atoms[i].resn !== meta[i].comp) return false; // order mismatch — leave atoms untouched
  }
  const counters = {};
  let prevKey = null, resi = 0;
  for (let i = 0; i < atoms.length; i++) {
    const a = atoms[i], m = meta[i];
    const key = `${a.chain}|${m.seq}|${m.icode}`;
    if (key !== prevKey) {
      counters[a.chain] = (counters[a.chain] || 0) + 1;
      resi = counters[a.chain];
      prevKey = key;
    }
    a.resLabel = `${m.seq}${m.icode}`; // true author number+insertion code for display
    a.icode = m.icode;
    a.resi = resi;                     // unique per residue so the cartoon renders
  }
  return true;
}

// Build residue groups for a chain from the already-loaded 3Dmol model.
// Works for PDB and mmCIF since 3Dmol has already parsed the atoms.
// Uses sequential scan rather than resi+icode key so that structures where
// multiple residues share the same auth_seq_id (CIF insertion-coded loops
// stored with pdbx_PDB_ins_code="?") are correctly split into separate groups.
function getChainResidueGroups(viewer, chainId) {
  const atoms = viewer.selectedAtoms({ chain: chainId });
  const groups = [];
  let prevResi = null, prevIcode = null, prevResn = null;
  let seenAtomNames = new Set();
  for (const atom of atoms) {
    const icode = atom.icode || "";
    const contextChanged = atom.resi !== prevResi || icode !== prevIcode || atom.resn !== prevResn;
    // Split when an atom name repeats within the same residue context: handles CIF structures
    // where consecutive identical residues share the same auth_seq_id (pdbx_PDB_ins_code="?").
    const repeatedAtom = !contextChanged && atom.atom != null && seenAtomNames.has(atom.atom);
    if (contextChanged || repeatedAtom) {
      // aa is null for non-standard residues (e.g. TYS sulfotyrosine) — they are kept in the
      // group array so their atoms can be highlighted, but excluded from sequence matching.
      groups.push({ aa: AA3TO1[atom.resn] || null, idx: groups.length, serials: [] });
      prevResi = atom.resi; prevIcode = icode; prevResn = atom.resn;
      seenAtomNames = new Set();
    }
    if (atom.atom != null) seenAtomNames.add(atom.atom);
    groups[groups.length - 1].serials.push(atom.serial);
  }
  return groups;
}

// Find CDR residue groups by matching the known CDR sequence in the chain's sequence.
// Non-standard residues (aa: null) are excluded from matching but remain in chainGroups
// so expandedSerials can include them in the highlighted region.
function findCdrGroupsBySeq(chainGroups, cdrSeq) {
  if (!cdrSeq || !chainGroups.length) return null;
  const aaGroups = chainGroups.filter(g => g.aa);
  const chainSeq = aaGroups.map(g => g.aa).join("");
  const idx = chainSeq.indexOf(cdrSeq);
  if (idx === -1) return null;
  return aaGroups.slice(idx, idx + cdrSeq.length);
}

// Fuzzy match: like findCdrGroupsBySeq but allows up to maxGaps residues from cdrSeq
// to be absent from the chain (disordered/unmodeled residues in the RCSB structure).
// Returns an array of length cdrSeq.length with a chain group at each matched position
// and null at gap positions, preserving the 1:1 index correspondence with cdrSeq so
// motif slicing still works correctly.
function findCdrGroupsFuzzy(chainGroups, cdrSeq, maxGaps = 4) {
  if (!cdrSeq || !chainGroups.length) return null;
  const aaGroups = chainGroups.filter(g => g.aa);
  const chainSeq = aaGroups.map(g => g.aa).join("");
  const n = cdrSeq.length;
  const m = chainSeq.length;
  for (let start = 0; start <= m - (n - maxGaps); start++) {
    const result = new Array(n).fill(null);
    let ci = start, gaps = 0, ok = true;
    for (let qi = 0; qi < n; qi++) {
      if (ci >= m) { gaps += (n - qi); break; }
      if (chainSeq[ci] === cdrSeq[qi]) {
        result[qi] = aaGroups[ci++];
      } else {
        if (++gaps > maxGaps) { ok = false; break; }
      }
    }
    if (ok && gaps <= maxGaps) return result;
  }
  return null;
}

// Scan every chain in the loaded model and return the first one containing cdrSeq.
// Handles PDB/CIF chain ID differences — SAbDab renames chains to single letters
// but RCSB CIF files may use multi-character author chain IDs.
// Falls back to a fuzzy match when residues are disordered/unmodeled in RCSB.
function findCdrInModel(viewer, cdrSeq) {
  if (!cdrSeq) return null;
  const allChainIds = [...new Set(viewer.selectedAtoms({}).map(a => a.chain))];
  for (const chainId of allChainIds) {
    const groups = getChainResidueGroups(viewer, chainId);
    const cdrGroups = findCdrGroupsBySeq(groups, cdrSeq);
    if (cdrGroups) return { chainId, cdrGroups, allGroups: groups };
  }
  for (const chainId of allChainIds) {
    const groups = getChainResidueGroups(viewer, chainId);
    const cdrGroups = findCdrGroupsFuzzy(groups, cdrSeq);
    if (cdrGroups) return { chainId, cdrGroups, allGroups: groups };
  }
  return null;
}

// Collect atom serials for a set of CDR groups, expanding to include any non-standard
// residues (aa: null, e.g. TYS sulfotyrosine) that sit between the first and last
// matched group in the full chain group array.
function expandedSerials(cdrGroups, allGroups) {
  const valid = cdrGroups.filter(g => g);
  if (!valid.length) return [];
  const minIdx = Math.min(...valid.map(g => g.idx));
  const maxIdx = Math.max(...valid.map(g => g.idx));
  return allGroups.slice(minIdx, maxIdx + 1).flatMap(g => g.serials);
}

// Create the 3Dmol viewer on first use.
function ensureViewer() {
  if (!viewer) {
    viewer = $3Dmol.createViewer($("viewer-container"), { backgroundColor: "#0e1117" });
    $("viewer-container").addEventListener("mouseleave", () => {
      viewer.setHover(null);
      viewer.render();
    });
  }
  return viewer;
}

// Attach the residue-label-on-hover behaviour to the current model.
function setupViewerHover() {
  viewer.setHoverable(
    {},
    true,
    (atom) => {
      if (atom.label) return;
      const label = atom.resLabel != null ? atom.resLabel : `${atom.resi}${(atom.icode || "").trim()}`;
      atom.label = viewer.addLabel(`${atom.resn} ${label} (chain ${atom.chain})`, {
        position: atom,
        backgroundColor: "black",
        backgroundOpacity: 0.7,
        fontColor: "white",
        fontSize: 12,
        borderThickness: 0,
      });
      viewer.render();
    },
    (atom) => {
      if (!atom.label) return;
      viewer.removeLabel(atom.label);
      delete atom.label;
      viewer.render();
    }
  );
}

// Fetch a structure from RCSB, preferring PDB format and falling back to mmCIF.
async function fetchStructure(pdb) {
  let resp = await fetch(`https://files.rcsb.org/download/${pdb}.pdb`);
  if (resp.ok) return { text: await resp.text(), fmt: "pdb" };
  resp = await fetch(`https://files.rcsb.org/download/${pdb}.cif`);
  if (resp.ok) return { text: await resp.text(), fmt: "cif" };
  return null;
}

// Annotate an arbitrary structure (RCSB fetch or user upload) that isn't in the DB:
// number every chain client-side, colour heavy/light chains, and highlight all 6 CDRs.
function annotateStructure(structData, structFmt, label) {
  if (!profileReady) { setViewerNote("Numbering data still loading — try again in a moment."); return; }
  $("viewer-placeholder").style.display = "none";
  setViewerNote(null);
  renderCdrReadout(null);
  ensureViewer();
  viewer.clear();
  viewer.addModel(structData, structFmt);
  fixCifResidueNumbering(viewer, structData, structFmt);
  setupViewerHover();
  viewer.setStyle({}, { cartoon: { color: "#888888" } });

  const chainIds = [...new Set(viewer.selectedAtoms({}).map((a) => a.chain))];
  const detected = [];
  for (const chainId of chainIds) {
    const groups = getChainResidueGroups(viewer, chainId);
    const seq = groups.filter((g) => g.aa).map((g) => g.aa).join("");
    const res = Numbering.numberChainAuto(seq);
    if (!res) continue;
    detected.push({ chainId, chainType: res.chainType });
    viewer.setStyle({ chain: chainId }, { cartoon: { color: res.isHeavy ? "#7fb3ff" : "#8fe3a0" } });
    for (const [cdrName, cdrSeq] of Object.entries(res.cdrs)) {
      if (!cdrSeq) continue;
      const cdrGroups = findCdrGroupsBySeq(groups, cdrSeq) || findCdrGroupsFuzzy(groups, cdrSeq);
      if (cdrGroups) {
        viewer.setStyle({ serial: expandedSerials(cdrGroups, groups) }, { cartoon: { color: "#ff9f43" } });
        // Show the CDR-H3 sequence readout (the annotate path clears it on load but
        // previously never populated it — only the search-results path did).
        if (res.isHeavy && cdrName === "h3") renderCdrReadout(cdrName, cdrSeq, cdrGroups);
      }
    }
  }

  if (detected.length) {
    const summary = detected.map((d) => `${d.chainId} (${d.chainType})`).join(", ");
    setViewerNote(`${label} — numbered ${detected.length} chain(s): ${summary}. CDRs highlighted.`);
  } else {
    setViewerNote(`${label} — no antibody variable domains detected.`);
  }
  viewer.zoomTo();
  viewer.render();
}

function setViewerNote(text) {
  const note = $("viewer-note");
  if (!text) {
    note.classList.add("hidden");
    note.textContent = "";
  } else {
    note.textContent = text;
    note.classList.remove("hidden");
  }
}

// Show which residues of the selected CDR are resolved in the loaded structure.
// cdrGroups is the per-position array from findCdrInModel: index i corresponds to
// cdrSeq[i], and a null entry means that residue is present in the sequence but not
// modeled in the coordinates (disordered). Because our CDRs are numbered from the
// sequence, a disordered residue is still listed here — this readout surfaces the gap.
function renderCdrReadout(cdrName, cdrSeq, cdrGroups) {
  const el = $("cdr-readout");
  if (!cdrSeq || !cdrGroups) { el.classList.add("hidden"); el.innerHTML = ""; return; }
  const positions = cdrSeq.split("").map((aa, i) => ({ aa, resolved: !!cdrGroups[i] }));
  const nRes = positions.filter((p) => p.resolved).length;
  const nTot = positions.length;
  const seqHtml = positions
    .map((p) => `<span class="${p.resolved ? "res" : "dis"}">${p.aa}</span>`)
    .join("");
  const msg = nRes === nTot
    ? `all ${nTot} residues resolved`
    : `${nRes} of ${nTot} residues resolved · ${nTot - nRes} disordered`;
  el.innerHTML =
    `<span class="cdr-label">CDR-${cdrName.toUpperCase()}</span>` +
    `<span class="cdr-seq">${seqHtml}</span>` +
    `<span class="cdr-msg">${msg}</span>`;
  el.classList.remove("hidden");
}

async function selectRow(index, cdr) {
  document.querySelectorAll("#results-body tr").forEach((tr) => tr.classList.remove("selected"));
  const tr = document.querySelector(`#results-body tr[data-index="${index}"]`);
  if (tr) tr.classList.add("selected");

  const row = lastResults[index];
  $("viewer-placeholder").style.display = "none";
  setViewerNote(null);
  renderCdrReadout(null);

  ensureViewer();

  const struct = await fetchStructure(row.pdb);
  if (!struct) return;
  const structData = struct.text;
  const structFmt = struct.fmt;

  viewer.clear();
  viewer.addModel(structData, structFmt);
  fixCifResidueNumbering(viewer, structData, structFmt);
  setupViewerHover();
  viewer.setStyle({}, { cartoon: { color: "#888888" } });

  const cdrSeq = row[`cdr_${cdr}`] || "";
  const isHeavyCdr = cdr[0] === "h";

  // Find antibody chains by sequence rather than by chain ID — SAbDab renames chains
  // to single letters (H, L) but RCSB CIF files may use multi-character author IDs.
  const cdrFind = findCdrInModel(viewer, cdrSeq);
  const partnerSeq = isHeavyCdr ? (row.cdr_l3 || row.cdr_l1 || "") : (row.cdr_h3 || row.cdr_h1 || "");
  const partnerFind = findCdrInModel(viewer, partnerSeq);

  // If the selected CDR wasn't found (e.g. disordered residues), fall back to another
  // CDR of the same chain just to identify the chain for colouring.
  let chainFind = cdrFind;
  if (!chainFind) {
    const fallbackSeqs = isHeavyCdr
      ? [row.cdr_h1, row.cdr_h2, row.cdr_h3]
      : [row.cdr_l1, row.cdr_l2, row.cdr_l3];
    for (const seq of fallbackSeqs) {
      if (seq && seq !== cdrSeq) {
        chainFind = findCdrInModel(viewer, seq);
        if (chainFind) break;
      }
    }
  }

  const heavyChainId = isHeavyCdr ? chainFind?.chainId : partnerFind?.chainId;
  const lightChainId = isHeavyCdr ? partnerFind?.chainId : chainFind?.chainId;
  if (heavyChainId) viewer.setStyle({ chain: heavyChainId }, { cartoon: { color: "#7fb3ff" } });
  if (lightChainId) viewer.setStyle({ chain: lightChainId }, { cartoon: { color: "#8fe3a0" } });

  if (cdrFind) {
    const cdrSerials = expandedSerials(cdrFind.cdrGroups, cdrFind.allGroups);
    viewer.setStyle({ serial: cdrSerials }, { cartoon: { color: "#ff9f43" } });

    // Show the resolved/disordered breakdown and tint residues flanking each gap so
    // the disorder location is visible on the 3D loop (the missing atom can't be drawn).
    renderCdrReadout(cdr, cdrSeq, cdrFind.cdrGroups);
    const flankSerials = [];
    for (let i = 0; i < cdrFind.cdrGroups.length; i++) {
      if (cdrFind.cdrGroups[i] === null) {
        for (const j of [i - 1, i + 1]) {
          const g = cdrFind.cdrGroups[j];
          if (g) flankSerials.push(...g.serials);
        }
      }
    }
    if (flankSerials.length) viewer.setStyle({ serial: flankSerials }, { cartoon: { color: "#a34d00" } });

    let matchedSerials = [];
    if (lastMotif) {
      const matchIdx = cdrSeq.indexOf(lastMotif);
      if (matchIdx !== -1) {
        matchedSerials = expandedSerials(cdrFind.cdrGroups.slice(matchIdx, matchIdx + lastMotif.length), cdrFind.allGroups);
        if (matchedSerials.length) {
          viewer.setStyle({ serial: matchedSerials }, { cartoon: { color: "#e63946" }, stick: { colorscheme: "redCarbon", radius: 0.3 } });
        }
      }
    }

    const epitopeBasisSel = matchedSerials.length ? { serial: matchedSerials } : { serial: cdrSerials };
    const dbAntigenChains = (row.antigen_chain || "").split(";").filter(Boolean);
    let antigenChainIds = dbAntigenChains;

    // If db antigen chain IDs match nothing in this model, fall back to all non-antibody chains
    if (antigenChainIds.length && !viewer.selectedAtoms({ chain: antigenChainIds }).length) {
      const abChainIds = new Set([heavyChainId, lightChainId].filter(Boolean));
      const allChainIds = [...new Set(viewer.selectedAtoms({}).map(a => a.chain))];
      antigenChainIds = allChainIds.filter(c => !abChainIds.has(c));
    }

    if (antigenChainIds.length) {
      const contactAtoms = viewer.selectedAtoms({ chain: antigenChainIds, within: { distance: EPITOPE_DISTANCE, sel: epitopeBasisSel } });
      const contactGroupKeys = new Set(contactAtoms.map((a) => `${a.chain}:${a.resi}:${a.icode || ""}`));
      const allAntigenAtoms = viewer.selectedAtoms({ chain: antigenChainIds });
      const epitopeSerials = allAntigenAtoms
        .filter((a) => contactGroupKeys.has(`${a.chain}:${a.resi}:${a.icode || ""}`))
        .map((a) => a.serial);
      if (epitopeSerials.length) {
        viewer.setStyle(
          { serial: epitopeSerials },
          { cartoon: { color: "#ff6b9d" }, stick: { colorscheme: "magentaCarbon", radius: 0.22 } }
        );
      } else {
        const basis = matchedSerials.length ? "motif" : "CDR loop";
        setViewerNote(`No binding to ${basis} within ${EPITOPE_DISTANCE}Å.`);
      }
    } else {
      setViewerNote("No antigen present in this structure.");
    }

    viewer.zoomTo({ serial: cdrSerials });
  } else {
    viewer.zoomTo();
  }
  viewer.render();
}

$("search-button").addEventListener("click", runSearch);
$("motif-input").addEventListener("keydown", (e) => { if (e.key === "Enter") runSearch(); });
$("collapse-duplicates-checkbox").addEventListener("change", renderTable);

// Annotate an arbitrary PDB ID (not necessarily in the DB) via client-side numbering.
async function annotateFromInput() {
  const id = $("annotate-input").value.trim();
  if (!id) return;
  setViewerNote(`Fetching ${id.toUpperCase()}…`);
  $("viewer-placeholder").style.display = "none";
  try {
    const struct = await fetchStructure(id);
    if (!struct) { setViewerNote(`Could not fetch ${id.toUpperCase()} from RCSB.`); return; }
    annotateStructure(struct.text, struct.fmt, id.toUpperCase());
  } catch (e) {
    setViewerNote(`Error loading ${id.toUpperCase()}: ${e.message}`);
  }
}

$("annotate-button").addEventListener("click", annotateFromInput);
$("annotate-input").addEventListener("keydown", (e) => { if (e.key === "Enter") annotateFromInput(); });
$("upload-button").addEventListener("click", () => $("upload-file-input").click());
$("upload-file-input").addEventListener("change", (e) => {
  const file = e.target.files[0];
  if (!file) return;
  const reader = new FileReader();
  reader.onload = () => {
    const fmt = /\.cif$/i.test(file.name) ? "cif" : "pdb";
    annotateStructure(reader.result, fmt, file.name);
  };
  reader.readAsText(file);
  e.target.value = ""; // allow re-uploading the same file
});

loadData();
