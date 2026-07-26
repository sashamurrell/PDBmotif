/*
 * Client-side antibody Chothia numbering.
 *
 * A faithful JS port of ANARCI's Chothia renumbering (schemes.py), with HMMER's
 * profile alignment replaced by a position-specific profile built from ANARCI's
 * germline sequences (data/export_numbering_profile.py -> data/numbering_profile.json).
 *
 * numberChain(seq, "H"|"L") returns { numbering, chainType } where numbering is a
 * list of [[num, ins], aa]. extractCdrs() reads CDR strings off it using the same
 * Chothia ranges as the Python pipeline, so results stay consistent with the DB.
 */
(function (global) {
  "use strict";

  // Chothia CDR residue-number ranges (inclusive), matching data/cdr_ranges.py.
  const CDR_RANGES = {
    h1: [26, 32], h2: [52, 56], h3: [93, 102],
    l1: [24, 34], l2: [50, 56], l3: [89, 97],
  };

  // ANARCI insertion alphabet: A..Z, AA..ZZ, then " " at the final index.
  // Python indexes this with -1 to mean the blank (no insertion) code.
  const ALPHABET = (function () {
    const a = [];
    for (let i = 0; i < 26; i++) a.push(String.fromCharCode(65 + i));
    for (let i = 0; i < 26; i++) a.push(String.fromCharCode(65 + i).repeat(2));
    a.push(" ");
    return a;
  })();
  function alpha(i) { return i < 0 ? " " : ALPHABET[i]; }

  const NEG = -Infinity;

  // Profile data (loaded from JSON via setProfile).
  let PROFILE = null;      // { alphabet, ncol, unknown_score, profiles: {H,K,L} }
  let AA_INDEX = null;     // residue -> column index into logodds rows
  function setProfile(p) {
    PROFILE = p;
    AA_INDEX = {};
    for (let i = 0; i < p.alphabet.length; i++) AA_INDEX[p.alphabet[i]] = i;
  }

  // ---- profile alignment (semi-global; free query overhang at both ends) ----
  // Returns an ANARCI state vector: [ [ [col, type], seqIndex ] , ... ]
  // col is 1-based IMGT position; type is "m"|"i"|"d"; seqIndex is 0-based or null.
  function alignProfile(query, chainType, gapOpen, gapExt, delCost) {
    gapOpen = gapOpen == null ? -10.0 : gapOpen;
    gapExt = gapExt == null ? -1.0 : gapExt;
    delCost = delCost == null ? 16.0 : delCost;

    const prof = PROFILE.profiles[chainType];
    const logodds = prof.logodds, occ = prof.occupancy;
    const unknown = PROFILE.unknown_score;
    const n = query.length, m = PROFILE.ncol;

    // M/Ix/Iy score matrices ((n+1) x (m+1)).
    const M = [], Ix = [], Iy = [];
    for (let i = 0; i <= n; i++) {
      M.push(new Float64Array(m + 1).fill(NEG));
      Ix.push(new Float64Array(m + 1).fill(NEG));
      Iy.push(new Float64Array(m + 1).fill(NEG));
    }
    for (let i = 0; i <= n; i++) Ix[i][0] = 0.0;            // free leading query gap
    for (let j = 0; j <= m; j++) {
      Iy[0][j] = -delCost * (j >= 1 ? occ[j - 1] : 0);
      if (j > 0) Iy[0][j] += Iy[0][j - 1] !== NEG ? Iy[0][j - 1] : 0;
    }

    for (let i = 1; i <= n; i++) {
      const qi = query[i - 1];
      const aidx = AA_INDEX[qi];
      for (let j = 1; j <= m; j++) {
        const s = aidx === undefined ? unknown : logodds[j - 1][aidx];
        let best = M[i - 1][j - 1];
        if (Ix[i - 1][j - 1] > best) best = Ix[i - 1][j - 1];
        if (Iy[i - 1][j - 1] > best) best = Iy[i - 1][j - 1];
        if (0.0 > best) best = 0.0;                          // local domain start
        M[i][j] = s + best;
        Ix[i][j] = Math.max(M[i - 1][j] + gapOpen, Ix[i - 1][j] + gapExt);
        const dcost = -delCost * occ[j - 1];
        Iy[i][j] = Math.max(M[i][j - 1] + dcost, Iy[i][j - 1] + dcost);
      }
    }

    // End anchor: best match cell reaching into the FW4 columns (>=100).
    let bestVal = NEG, bi = n, bj = m;
    for (let i = 1; i <= n; i++) {
      for (let j = 100; j <= m; j++) {
        if (M[i][j] > bestVal) { bestVal = M[i][j]; bi = i; bj = j; }
      }
    }

    const sv = [];
    let i = bi, j = bj, state = "M";
    while (i > 0 && j > 0) {
      if (state === "M") {
        sv.push([[j, "m"], i - 1]);
        // tie-break order matches Python max(): M, Ix, Iy, START.
        let lbl = "M", val = M[i - 1][j - 1];
        if (Ix[i - 1][j - 1] > val) { lbl = "Ix"; val = Ix[i - 1][j - 1]; }
        if (Iy[i - 1][j - 1] > val) { lbl = "Iy"; val = Iy[i - 1][j - 1]; }
        if (0.0 > val) { lbl = "START"; val = 0.0; }
        i -= 1; j -= 1;
        if (lbl === "START") break;
        state = lbl;
      } else if (state === "Ix") {
        sv.push([[j < m ? j + 1 : j, "i"], i - 1]);
        state = (M[i - 1][j] + gapOpen >= Ix[i - 1][j] + gapExt) ? "M" : "Ix";
        i -= 1;
      } else { // Iy (deletion)
        const dcost = -delCost * occ[j - 1];
        sv.push([[j, "d"], null]);
        state = (M[i][j - 1] + dcost >= Iy[i][j - 1] + dcost) ? "M" : "Iy";
        j -= 1;
      }
    }
    sv.reverse();
    return sv;
  }

  // ---- smooth_insertions (port of schemes.py) ----
  const ENFORCED_PATTERNS = [
    [[25, "m"], [26, "m"], [27, "m"], [28, "i"]],
    [[38, "i"], [38, "m"], [39, "m"], [40, "m"]],
    [[54, "m"], [55, "m"], [56, "m"], [57, "i"]],
    [[65, "i"], [65, "m"], [66, "m"], [67, "m"]],
    [[103, "m"], [104, "m"], [105, "m"], [106, "i"]],
    [[117, "i"], [117, "m"], [118, "m"], [119, "m"]],
  ];

  function smoothInsertions(stateVector) {
    const stateBuffer = [];
    const sv = [];
    let reg = -2;
    for (const entry of stateVector) {
      const stateId = entry[0][0], stateType = entry[0][1], si = entry[1];
      if (stateId < 23) { stateBuffer.push(entry); reg = -1; continue; }
      else if (stateId >= 25 && stateId < 28) { stateBuffer.push(entry); reg = 0; continue; }
      else if (stateId > 37 && stateId <= 40) { stateBuffer.push(entry); reg = 1; continue; }
      else if (stateId >= 54 && stateId < 57) { stateBuffer.push(entry); reg = 2; continue; }
      else if (stateId > 64 && stateId <= 67) { stateBuffer.push(entry); reg = 3; continue; }
      else if (stateId >= 103 && stateId < 106) { stateBuffer.push(entry); reg = 4; continue; }
      else if (stateId > 116 && stateId <= 119) { stateBuffer.push(entry); reg = 5; continue; }
      else if (stateBuffer.length !== 0) {
        const nins = stateBuffer.filter(s => s[0][1] === "i").length;
        if (nins > 0) {
          if (reg === -1) { // FW1: only adjust if >= as many N-terminal deletions as insertions
            let ntDels = stateBuffer[0][0][0] - 1;
            for (const s of stateBuffer) {
              if (s[0][1] === "d" || s[1] === null) ntDels += 1;
              else break;
            }
            if (ntDels >= nins) {
              let newStates = stateBuffer.filter(s => s[0][1] === "m").map(s => s[0]);
              const first = newStates[0][0];
              const buf = stateBuffer.filter(s => s[0][1] !== "d");
              const add = buf.length - newStates.length;
              const prefix = [];
              for (let k = first - add; k < first; k++) prefix.push([k, "m"]);
              newStates = prefix.concat(newStates);
              for (let k = 0; k < buf.length; k++) sv.push([newStates[k], buf[k][1]]);
            } else {
              for (const s of stateBuffer) sv.push(s);
            }
          } else {
            const buf = stateBuffer.filter(s => s[0][1] !== "d");
            const pat = ENFORCED_PATTERNS[reg];
            let newStates;
            if (reg % 2) { // N-terminal framework
              const head = [];
              for (let k = 0; k < Math.max(0, buf.length - 3); k++) head.push(pat[0]);
              newStates = head.concat(pat.slice(Math.max(4 - buf.length, 1)));
            } else { // C-terminal framework
              const tail = [];
              for (let k = 0; k < Math.max(0, buf.length - 3); k++) tail.push(pat[2]);
              newStates = pat.slice(0, 3).concat(tail);
            }
            for (let k = 0; k < buf.length; k++) sv.push([newStates[k], buf[k][1]]);
          }
        } else {
          for (const s of stateBuffer) sv.push(s);
        }
        sv.push(entry);
        stateBuffer.length = 0;
      } else {
        sv.push(entry);
      }
    }
    return sv;
  }

  // ---- _number_regions (port of schemes.py) ----
  function numberRegions(sequence, stateVector, stateString, regionString,
                         regionIndexDict, relsIn, nRegions, excludeDeletions) {
    const sv = smoothInsertions(stateVector);
    const rels = Object.assign({}, relsIn); // mutated below; keep caller's copy clean
    const regions = [];
    for (let r = 0; r < nRegions; r++) regions.push([]);

    let insertion = -1;
    let previousStateId = 1;
    let previousStateType = "d";
    let startIndex = null, endIndex = null;
    let region = null;

    for (const entry of sv) {
      const stateId = entry[0][0], stateType = entry[0][1], si = entry[1];

      if (stateType !== "i" || region === null) {
        region = regionIndexDict[regionString[stateId - 1]];
      }

      if (stateType === "m") {
        if (stateString[stateId - 1] === "I") {
          if (previousStateType !== "d") insertion += 1;
          rels[region] -= 1;
        } else {
          insertion = -1;
        }
        regions[region].push([[stateId + rels[region], alpha(insertion)], sequence[si]]);
        previousStateId = stateId;
        if (startIndex === null) startIndex = si;
        endIndex = si;
        previousStateType = stateType;
      } else if (stateType === "i") {
        insertion += 1;
        regions[region].push([[previousStateId + rels[region], alpha(insertion)], sequence[si]]);
        if (startIndex === null) startIndex = si;
        endIndex = si;
        previousStateType = stateType;
      } else { // deletion
        previousStateType = stateType;
        if (stateString[stateId - 1] === "I") {
          rels[region] -= 1;
          continue;
        }
        insertion = -1;
        previousStateId = stateId;
      }

      if (insertion >= 25 && excludeDeletions.indexOf(region) !== -1) insertion = 0;
      if (insertion >= 25) throw new Error("Too many insertions for numbering scheme to handle");
    }
    return { regions, startIndex, endIndex };
  }

  // ---- get_cdr3_annotations (chothia only) ----
  function getCdr3Annotations(length, chainType) {
    const az = "ABCDEFGHIJKLMNOPQRSTUVWXYZ";
    let ordered, insertions, anchor;
    if (chainType === "heavy") {
      insertions = Math.max(length - 10, 0);
      ordered = [[100, " "], [99, " "], [98, " "], [97, " "], [96, " "], [95, " "],
                 [101, " "], [102, " "], [94, " "], [93, " "]];
      anchor = 100;
    } else {
      insertions = Math.max(length - 9, 0);
      ordered = [[95, " "], [94, " "], [93, " "], [92, " "], [91, " "],
                 [96, " "], [97, " "], [90, " "], [89, " "]];
      anchor = 95;
    }
    if (insertions >= 27) throw new Error("Too many insertions for numbering scheme to handle");
    const base = ordered.slice(Math.max(0, ordered.length - length));
    const ins = [];
    for (let a = 0; a < insertions; a++) ins.push([anchor, az[a]]);
    const all = base.concat(ins);
    all.sort((x, y) => x[0] - y[0] || (x[1] < y[1] ? -1 : x[1] > y[1] ? 1 : 0));
    return all;
  }

  // ---- gap_missing ----
  function gapMissing(numbering) {
    const flat = [];
    for (const region of numbering) for (const e of region) flat.push(e);
    const num = [[[0, " "], "-"]];
    for (const [p, a] of flat) {
      if (p[0] > num[num.length - 1][0][0] + 1) {
        for (let k = num[num.length - 1][0][0] + 1; k < p[0]; k++) num.push([[k, " "], "-"]);
      }
      num.push([p, a]);
    }
    return num.slice(1);
  }

  // Build [[num, ins], aa] region from an annotations list + the region's residues.
  function applyAnnotations(annotations, regionResidues) {
    const out = [];
    for (let i = 0; i < regionResidues.length; i++) {
      out.push([annotations[i], regionResidues[i][1]]);
    }
    return out;
  }

  // ---- number_chothia_heavy ----
  const H_STATE = "XXXXXXXXXIXXXXXXXXXXXXXXXXXXXXIIIIXXXXXXXXXXXXXXXXXXXXXXXIXIIXXXXXXXXXXXIXXXXXXXXXXXXXXXXXXIIIXXXXXXXXXXXXXXXXXXIIIXXXXXXXXXXXXX";
  const H_REGION = "11111111112222222222222333333333333333444444444444444455555555555666666666666666666666666666666666666666777777777777788888888888";
  const H_RIDX = { "1": 0, "2": 1, "3": 2, "4": 3, "5": 4, "6": 5, "7": 6, "8": 7 };
  const H_RELS = { 0: 0, 1: -1, 2: -1, 3: -5, 4: -5, 5: -8, 6: -12, 7: -15 };

  function numberChothiaHeavy(stateVector, sequence) {
    const { regions, startIndex, endIndex } = numberRegions(
      sequence, stateVector, H_STATE, H_REGION, H_RIDX, H_RELS, 8, [0, 2, 4, 6]);
    const numbering = [[], regions[1], [], regions[3], [], regions[5], [], regions[7]];

    // Region 0: insertions placed at Chothia 6.
    let insertions = regions[0].filter(x => x[0][1] !== " ").length;
    if (insertions) {
      const start = regions[0][0][0][0];
      const ann = [];
      for (let k = start; k < 7; k++) ann.push([k, " "]);
      for (let k = 0; k < insertions; k++) ann.push([6, alpha(k)]);
      ann.push([7, " "], [8, " "], [9, " "]);
      numbering[0] = applyAnnotations(ann, regions[0]);
    } else {
      numbering[0] = regions[0];
    }

    // CDR1 (region 2): insertions on 31.
    let length = regions[2].length;
    insertions = Math.max(length - 11, 0);
    let ann;
    if (insertions) {
      ann = [];
      for (let k = 23; k < 32; k++) ann.push([k, " "]);
      for (let k = 0; k < insertions; k++) ann.push([31, alpha(k)]);
      ann.push([32, " "], [33, " "]);
    } else {
      ann = [];
      for (let k = 23; k < 32; k++) ann.push([k, " "]);
      ann = ann.slice(0, length - 2).concat([[32, " "], [33, " "]].slice(0, length));
    }
    numbering[2] = applyAnnotations(ann, regions[2]);

    // CDR2 (region 4): insertions on 52.
    length = regions[4].length;
    insertions = Math.max(length - 8, 0);
    ann = [[50, " "], [51, " "], [52, " "]].slice(0, Math.max(0, length - 5));
    for (let k = 0; k < insertions; k++) ann.push([52, alpha(k)]);
    ann = ann.concat([[53, " "], [54, " "], [55, " "], [56, " "], [57, " "]]
      .slice(Math.abs(Math.min(0, length - 5))));
    numbering[4] = applyAnnotations(ann, regions[4]);

    // CDR3 (region 6): insertions on 100.
    length = regions[6].length;
    if (length > 36) return { numbering: [], startIndex, endIndex };
    numbering[6] = applyAnnotations(getCdr3Annotations(length, "heavy"), regions[6]);

    return { numbering: gapMissing(numbering), startIndex, endIndex };
  }

  // ---- number_chothia_light ----
  const L_STATE = "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXIIIIIIXXXXXXXXXXXXXXXXXXXXXXIIIIIIIXXXXXXXXIXXXXXXXIIXXXXXXXXXXXXXXXXXXXXXXXXXXXIIIIXXXXXXXXXXXXXXX";
  const L_REGION = "11111111111111111111111222222222222222223333333333333333444444444445555555555555555555555555555555555555666666666666677777777777";
  const L_RIDX = { "1": 0, "2": 1, "3": 2, "4": 3, "5": 4, "6": 5, "7": 6 };
  const L_RELS = { 0: 0, 1: 0, 2: -6, 3: -6, 4: -13, 5: -16, 6: -20 };

  function numberChothiaLight(stateVector, sequence) {
    const { regions, startIndex, endIndex } = numberRegions(
      sequence, stateVector, L_STATE, L_REGION, L_RIDX, L_RELS, 7, [1, 3, 4, 5]);
    const numbering = [regions[0], [], regions[2], [], regions[4], [], regions[6]];

    // CDR1 (region 1): insertions on 30.
    let length = regions[1].length;
    let insertions = Math.max(length - 11, 0);
    let ann = [[24, " "], [25, " "], [26, " "], [27, " "], [28, " "], [29, " "], [30, " "]]
      .slice(0, Math.max(0, length));
    for (let k = 0; k < insertions; k++) ann.push([30, alpha(k)]);
    ann = ann.concat([[31, " "], [32, " "], [33, " "], [34, " "]]
      .slice(Math.abs(Math.min(0, length - 11))));
    numbering[1] = applyAnnotations(ann, regions[1]);

    // CDR2 (region 3): insertions on 52.
    length = regions[3].length;
    insertions = Math.max(length - 4, 0);
    if (insertions > 0) {
      ann = [[51, " "], [52, " "]];
      for (let k = 0; k < insertions; k++) ann.push([52, alpha(k)]);
      ann.push([53, " "], [54, " "]);
      numbering[3] = applyAnnotations(ann, regions[3]);
    } else {
      numbering[3] = regions[3];
    }

    // FW3 (region 4): insertions on 68.
    length = regions[4].length;
    insertions = Math.max(length - 34, 0);
    if (insertions > 0) {
      ann = [];
      for (let k = 55; k < 69; k++) ann.push([k, " "]);
      for (let k = 0; k < insertions; k++) ann.push([68, alpha(k)]);
      for (let k = 69; k < 89; k++) ann.push([k, " "]);
      numbering[4] = applyAnnotations(ann, regions[4]);
    } else if (length === 33) {
      ann = [];
      for (let k = 55; k < 68; k++) ann.push([k, " "]);
      for (let k = 69; k < 89; k++) ann.push([k, " "]);
      numbering[4] = applyAnnotations(ann, regions[4]);
    } else {
      numbering[4] = regions[4];
    }

    // CDR3 (region 5): insertions on 95.
    length = regions[5].length;
    if (length > 35) return { numbering: [], startIndex, endIndex };
    numbering[5] = applyAnnotations(getCdr3Annotations(length, "light"), regions[5]);

    return { numbering: gapMissing(numbering), startIndex, endIndex };
  }

  // ---- top-level ----
  function countMatches(sv) {
    let c = 0;
    for (const e of sv) if (e[0][1] === "m") c++;
    return c;
  }

  // Number a chain sequence. `heavy` true -> H; false -> pick better of K/L.
  // Returns { numbering, chainType } or null if numbering failed.
  function numberChain(seq, heavy) {
    if (!PROFILE) throw new Error("Numbering profile not loaded (call setProfile)");
    seq = seq.toUpperCase().replace(/[^A-Z]/g, "");
    try {
      if (heavy) {
        const sv = alignProfile(seq, "H");
        const { numbering } = numberChothiaHeavy(sv, seq);
        return numbering.length ? { numbering, chainType: "H" } : null;
      }
      const svk = alignProfile(seq, "K");
      const svl = alignProfile(seq, "L");
      const useK = countMatches(svk) >= countMatches(svl);
      const sv = useK ? svk : svl;
      const { numbering } = numberChothiaLight(sv, seq);
      return numbering.length ? { numbering, chainType: useK ? "K" : "L" } : null;
    } catch (e) {
      return null;
    }
  }

  // Minimum profile-match columns for a sequence to be accepted as an antibody
  // variable domain. Real V-domains match ~83-115 of the 128 columns; non-antibody
  // chains (antigens) force-fit at most ~16, and multidomain fusion constructs where
  // the domain is a tiny fraction score ~40-60. 70 cleanly separates the two.
  const MIN_DOMAIN_MATCHES = 70;

  // Auto-detect chain type (H vs K vs L) for an unknown chain and number it.
  // Returns { numbering, chainType, isHeavy, cdrs } or null if it isn't a confident
  // antibody variable domain.
  function numberChainAuto(seq) {
    if (!PROFILE) throw new Error("Numbering profile not loaded (call setProfile)");
    seq = seq.toUpperCase().replace(/[^A-Z]/g, "");
    if (seq.length < 60) return null; // too short to be a variable domain
    let best = null, bestScore = -1;
    for (const ct of ["H", "K", "L"]) {
      let sv;
      try { sv = alignProfile(seq, ct); } catch (e) { continue; }
      const score = countMatches(sv);
      if (score > bestScore) { bestScore = score; best = { ct, sv }; }
    }
    if (!best || bestScore < MIN_DOMAIN_MATCHES) return null;
    const isHeavy = best.ct === "H";
    try {
      const { numbering } = isHeavy
        ? numberChothiaHeavy(best.sv, seq)
        : numberChothiaLight(best.sv, seq);
      if (!numbering.length) return null;
      const names = isHeavy ? ["h1", "h2", "h3"] : ["l1", "l2", "l3"];
      return { numbering, chainType: best.ct, isHeavy, cdrs: extractCdrs(numbering, names) };
    } catch (e) {
      return null;
    }
  }

  // Extract CDR strings for the given names from a numbering.
  function extractCdrs(numbering, names) {
    const out = {};
    for (const name of names) {
      const [lo, hi] = CDR_RANGES[name];
      let s = "";
      for (const [pos, aa] of numbering) {
        if (pos[0] >= lo && pos[0] <= hi && aa !== "-" && aa != null) s += aa;
      }
      out[name] = s || null;
    }
    return out;
  }

  const api = { setProfile, numberChain, numberChainAuto, extractCdrs, alignProfile, CDR_RANGES };
  if (typeof module !== "undefined" && module.exports) module.exports = api;
  else global.Numbering = api;
})(typeof globalThis !== "undefined" ? globalThis : this);
