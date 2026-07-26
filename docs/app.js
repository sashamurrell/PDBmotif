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

function categorizeAntigen(row) {
  const ag = (row.antigen_species || "").split("|")[0].trim().toLowerCase();
  if (!ag) return "";
  const hsp = (row.heavy_species || "").toLowerCase().trim();
  const lsp = (row.light_species || "").toLowerCase().trim();
  if ((hsp && ag === hsp) || (lsp && ag === lsp)) return "self";
  if (/virus|phage|viral|coronavirus|hiv|siv|influenza|hepatitis|dengue|zika|rabies|ebola|herpes|adeno|retrovirus/.test(ag)) return "virus";
  if (/plasmodium|trypanosoma|leishmania|toxoplasma/.test(ag)) return "parasite";
  if (/bacteri|bacillus|streptococ|staphyloco|mycobacter|escherichia|salmonella|clostridium|pseudomonas|klebsiella|bordetella|helicobacter/.test(ag)) return "bacteria";
  if (/aspergillus|candida|cryptococcus|saccharomyces/.test(ag)) return "fungal";
  if (/synthetic|artificial/.test(ag)) return "synthetic";
  return "other";
}

// SAbDab lists every resolved antigen component in antigen_name, joined by "|",
// including the N-linked glycans on a glycoprotein (HIV Env, SARS spike, etc).
// Drop those glycan entries — their IUPAC names always contain "pyranose"/"furanose" —
// so the target shows the actual protein antigen. If a structure's only antigen is a
// sugar (hapten-like), keep it rather than returning nothing.
function cleanAntigenName(name) {
  if (!name) return "";
  const parts = name.split("|").map(s => s.trim()).filter(Boolean);
  const nonSugar = parts.filter(p => !/pyranose|furanose/i.test(p));
  return (nonSugar.length ? nonSugar : parts).join(", ");
}

function formatTarget(row) {
  const name = cleanAntigenName(row.antigen_name);
  if (categorizeAntigen(row) === "self") return name || "n/a";
  const species = row.antigen_species
    ? row.antigen_species.split("|")[0].trim()
    : "";
  if (species && name) return `${species} ${name}`;
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
      const icode = (atom.icode || "").trim();
      atom.label = viewer.addLabel(`${atom.resn} ${atom.resi}${icode} (chain ${atom.chain})`, {
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
    for (const cdrSeq of Object.values(res.cdrs)) {
      if (!cdrSeq) continue;
      const cdrGroups = findCdrGroupsBySeq(groups, cdrSeq) || findCdrGroupsFuzzy(groups, cdrSeq);
      if (cdrGroups) {
        viewer.setStyle({ serial: expandedSerials(cdrGroups, groups) }, { cartoon: { color: "#ff9f43" } });
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
