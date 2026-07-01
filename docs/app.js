// Mirrors data/cdr_ranges.py -- keep in sync if the Chothia boundaries there change.
const CDR_RANGES = {
  h1: [26, 32], h2: [52, 56], h3: [95, 102],
  l1: [24, 34], l2: [50, 56], l3: [89, 97],
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

function $(id) { return document.getElementById(id); }

async function loadData() {
  const resp = await fetch("data/antibodies.json");
  allData = await resp.json();
  buildFilterOptions();
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
      <td>${row.target_category || ""}</td>
      <td>${row.target_antigen || row.antigen_name || (row.target_category == null ? "n/a" : "")}</td>
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

  const total_count = results.length;
  const limited = results.slice(0, RESULT_LIMIT);

  lastMotif = motif;
  renderResults({ results: limited, total_count }, cdr, motif);
}

function resiRange(start, end) {
  const out = [];
  for (let i = start; i <= end; i++) out.push(i);
  return out;
}

function parseChainResidueGroups(pdbText, chainId) {
  const indexByKey = new Map();
  const groups = [];
  for (const line of pdbText.split("\n")) {
    if (!line.startsWith("ATOM")) continue;
    if (line[21].toUpperCase() !== chainId.toUpperCase()) continue;
    const resnumStr = line.slice(22, 27).trim();
    const serial = parseInt(line.slice(6, 11).trim(), 10);
    let groupIdx = indexByKey.get(resnumStr);
    if (groupIdx === undefined) {
      let i = resnumStr.length;
      while (i > 0 && !/[0-9]/.test(resnumStr[i - 1])) i--;
      groupIdx = groups.length;
      indexByKey.set(resnumStr, groupIdx);
      groups.push({ resnum: parseInt(resnumStr.slice(0, i), 10), serials: [] });
    }
    groups[groupIdx].serials.push(serial);
  }
  return groups;
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

async function selectRow(index, cdr) {
  document.querySelectorAll("#results-body tr").forEach((tr) => tr.classList.remove("selected"));
  const tr = document.querySelector(`#results-body tr[data-index="${index}"]`);
  if (tr) tr.classList.add("selected");

  const row = lastResults[index];
  $("viewer-placeholder").style.display = "none";
  setViewerNote(null);

  if (!viewer) {
    viewer = $3Dmol.createViewer($("viewer-container"), { backgroundColor: "#0e1117" });
    $("viewer-container").addEventListener("mouseleave", () => {
      viewer.setHover(null);
      viewer.render();
    });
  }

  const resp = await fetch(`https://files.rcsb.org/download/${row.pdb}.pdb`);
  if (!resp.ok) return;
  const pdbData = await resp.text();

  viewer.clear();
  viewer.addModel(pdbData, "pdb");
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
  viewer.setStyle({}, { cartoon: { color: "#888888" } });
  if (row.hchain && row.hchain !== "NA") {
    viewer.setStyle({ chain: row.hchain }, { cartoon: { color: "#7fb3ff" } });
  }
  if (row.lchain && row.lchain !== "NA") {
    viewer.setStyle({ chain: row.lchain }, { cartoon: { color: "#8fe3a0" } });
  }

  const cdrChain = cdr[0] === "h" ? row.hchain : row.lchain;
  if (cdrChain && cdrChain !== "NA") {
    const [start, end] = CDR_RANGES[cdr];
    const sel = { chain: cdrChain, resi: resiRange(start, end) };
    viewer.setStyle(sel, { cartoon: { color: "#ff9f43" } });

    let matchedSerials = [];
    if (lastMotif) {
      const cdrSeq = row[`cdr_${cdr}`] || "";
      const matchIdx = cdrSeq.indexOf(lastMotif);
      if (matchIdx !== -1) {
        const cdrGroups = parseChainResidueGroups(pdbData, cdrChain).filter((g) => g.resnum >= start && g.resnum <= end);
        matchedSerials = cdrGroups.slice(matchIdx, matchIdx + lastMotif.length).flatMap((g) => g.serials);
        if (matchedSerials.length) {
          viewer.setStyle({ serial: matchedSerials }, { cartoon: { color: "#e63946" }, stick: { colorscheme: "redCarbon", radius: 0.3 } });
        }
      }
    }

    const epitopeBasisSel = matchedSerials.length ? { serial: matchedSerials } : sel;

    const antigenChains = (row.antigen_chain || "").split(";").filter(Boolean);
    const contactBasis = matchedSerials.length ? "the matched motif" : "the CDR loop";
    if (antigenChains.length) {
      const contactAtoms = viewer.selectedAtoms({ chain: antigenChains, within: { distance: EPITOPE_DISTANCE, sel: epitopeBasisSel } });
      const contactGroupKeys = new Set(contactAtoms.map((a) => `${a.chain}:${a.resi}:${a.icode || ""}`));
      const allAntigenAtoms = viewer.selectedAtoms({ chain: antigenChains });
      const epitopeSerials = allAntigenAtoms
        .filter((a) => contactGroupKeys.has(`${a.chain}:${a.resi}:${a.icode || ""}`))
        .map((a) => a.serial);
      if (epitopeSerials.length) {
        viewer.setStyle(
          { serial: epitopeSerials },
          { cartoon: { color: "#ff6b9d" }, stick: { colorscheme: "magentaCarbon", radius: 0.22 } }
        );
      } else {
        setViewerNote(`No antigen contacts within ${EPITOPE_DISTANCE}Å of ${contactBasis}.`);
      }
    } else {
      setViewerNote("No antigen present in this structure.");
    }

    viewer.zoomTo(sel);
  } else {
    viewer.zoomTo();
  }
  viewer.render();
}

$("search-button").addEventListener("click", runSearch);
$("motif-input").addEventListener("keydown", (e) => { if (e.key === "Enter") runSearch(); });
$("collapse-duplicates-checkbox").addEventListener("change", renderTable);

loadData();
