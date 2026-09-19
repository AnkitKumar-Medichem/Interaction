import { PredictionResult, PredictionMethod, FunctionalGroupReactivity } from "./gemini";

/**
 * ============================================================================
 *  GENERAL-PURPOSE FUNCTIONAL-GROUP REACTION ENGINE
 * ============================================================================
 *  Domain-agnostic: works for any organic / organometallic / inorganic species
 *  given as a SMILES string or a common chemical name.
 *
 *  How it works (pipeline):
 *   1. SMILES  -> molecular graph (own parser: brackets, charges, aromatic and
 *      Kekule input, ring closures, salts / multi-component input).
 *   2. Graph-based functional-group perception (topology, not text regexes), so
 *      the result does not depend on how a SMILES happens to be written.
 *   3. Per-group reactivity profile for 7 conditions (acid, base, neutral
 *      hydrolysis, photolysis, thermal, oxidative, cross-interaction), adjusted
 *      by local context (aryl vs alkyl ester, ring strain, amine class, ...).
 *   4. Reaction templates that EDIT the graph, so every predicted product is a
 *      valid, atom-balanced structure (formula and mass shift are reported).
 *   5. Co-reactant analysis: complementary reactive classes are matched
 *      (nucleophile/electrophile, acid/base, oxidant, carbonyl/amine ...).
 *   6. Ranking: heuristic likelihood (from the same vulnerability grades shown
 *      in the functional-group table) and Boltzmann weights at 298.15 K.
 *
 *  Public API (unchanged): PHARMA_COMPOUNDS (deprecated alias of COMPOUND_LIBRARY), lookupCompoundSmiles,
 *  detectFunctionalGroupsDetailed, detectFunctionalGroups, generateComputationalPrediction.
 *  New (optional): registerCompound(names, smiles), parseSmiles / writeSmiles / canonicalSmiles / formulaOf.
 *  Unresolvable input never falls back to a placeholder compound: an empty, clearly explained result is returned.
 *
 *  Limitations (by design, no external dependencies): stereochemistry and
 *  isotopes are parsed but not propagated to products; free energies are
 *  reaction-class estimates on a compressed kcal/mol scale, not QM values;
 *  names are descriptive, not systematic IUPAC.
 * ============================================================================
 */

type Vuln = "Critical" | "High" | "Moderate" | "Low" | "Resistant";
type Cond = "Oxidation" | "Acidic Hydrolysis" | "Basic Hydrolysis" | "Hydrolysis" | "Photodegradation" | "Thermal Degradation";

/* ------------------------------------------------------------------ */
/*  Element data                                                       */
/* ------------------------------------------------------------------ */
const MASS: Record<string, number> = {
  H: 1.007825, Li: 7.016003, B: 11.009305, C: 12, N: 14.003074, O: 15.994915, F: 18.998403, Na: 22.98977,
  Mg: 23.985042, Al: 26.981538, Si: 27.976927, P: 30.973762, S: 31.972071, Cl: 34.968853, K: 38.963707,
  Ca: 39.962591, Ti: 47.947946, V: 50.94396, Cr: 51.940508, Mn: 54.938045, Fe: 55.934936, Co: 58.933195,
  Ni: 57.935343, Cu: 62.929598, Zn: 63.929142, Ga: 68.925581, Ge: 73.921178, As: 74.921596, Se: 79.916522,
  Br: 78.918338, Rb: 84.911789, Sr: 87.905612, Zr: 89.904704, Mo: 97.905408, Pd: 105.903486, Ag: 106.905093,
  Cd: 113.903358, Sn: 119.902197, Sb: 120.903812, I: 126.904473, Cs: 132.905447, Ba: 137.905247,
  Ce: 139.905439, W: 183.950931, Pt: 194.964791, Au: 196.966569, Hg: 201.970643, Pb: 207.976652, Bi: 208.980399,
};
const VALENCE: Record<string, number[]> = {
  H: [1], B: [3], C: [4], N: [3, 5], O: [2], F: [1], Si: [4], P: [3, 5], S: [2, 4, 6],
  Cl: [1], Se: [2, 4, 6], As: [3, 5], Br: [1], I: [1, 3, 5],
};
const ORGANIC = new Set(["B", "C", "N", "O", "P", "S", "F", "Cl", "Br", "I"]);
const AROMATIC_ALLOWED = new Set(["B", "C", "N", "O", "P", "S", "Se", "As"]);
const METALS = new Set([
  "Li", "Na", "K", "Rb", "Cs", "Mg", "Ca", "Sr", "Ba", "Al", "Ga", "Zn", "Cd", "Hg", "Sn", "Pb", "Bi", "Sb",
  "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Ag", "Au", "Pd", "Pt", "Zr", "Mo", "W", "Ce",
]);
const HALOGENS = new Set(["F", "Cl", "Br", "I"]);

/* ------------------------------------------------------------------ */
/*  Molecular graph                                                    */
/* ------------------------------------------------------------------ */
interface Atom { el: string; arom: boolean; q: number; h: number; alive: boolean; bracket: boolean; }

export class Mol {
  atoms: Atom[] = [];
  adj: Map<number, number>[] = [];
  add(el: string, o: Partial<Atom> = {}): number {
    this.atoms.push({ el, arom: false, q: 0, h: 0, alive: true, bracket: false, ...o });
    this.adj.push(new Map());
    return this.atoms.length - 1;
  }
  bond(a: number, b: number, o = 1): void { this.adj[a].set(b, o); this.adj[b].set(a, o); }
  unbond(a: number, b: number): void { this.adj[a].delete(b); this.adj[b].delete(a); }
  order(a: number, b: number): number { return this.adj[a].get(b) ?? 0; }
  remove(i: number): void {
    for (const j of [...this.adj[i].keys()]) this.unbond(i, j);
    this.atoms[i].alive = false;
  }
  nb(i: number): number[] { return [...this.adj[i].keys()].filter((j) => this.atoms[j].alive); }
  deg(i: number): number { return this.nb(i).length; }
  el(i: number): string { return this.atoms[i].el; }
  alive(): number[] { const r: number[] = []; this.atoms.forEach((a, i) => { if (a.alive) r.push(i); }); return r; }
  clone(): Mol {
    const c = new Mol();
    c.atoms = this.atoms.map((a) => ({ ...a }));
    c.adj = this.adj.map((m) => new Map(m));
    return c;
  }
  /** neighbours of i excluding one atom */
  nbx(i: number, ex: number): number[] { return this.nb(i).filter((j) => j !== ex); }
}

/** Append a copy of `b` into a clone of `a`; returns merged mol and index offset of b's atoms. */
function mergeMol(a: Mol, b: Mol): { mol: Mol; off: number } {
  const mol = a.clone();
  const off = mol.atoms.length;
  b.atoms.forEach((at) => { mol.atoms.push({ ...at }); mol.adj.push(new Map()); });
  b.adj.forEach((mp, i) => mp.forEach((o, j) => mol.adj[i + off].set(j + off, o)));
  return { mol, off };
}

function components(m: Mol): number[][] {
  const seen = new Set<number>(); const out: number[][] = [];
  for (const s of m.alive()) {
    if (seen.has(s)) continue;
    const comp: number[] = []; const st = [s]; seen.add(s);
    while (st.length) {
      const u = st.pop()!; comp.push(u);
      for (const v of m.nb(u)) if (!seen.has(v)) { seen.add(v); st.push(v); }
    }
    out.push(comp.sort((x, y) => x - y));
  }
  return out;
}

/** default implicit hydrogens of an un-bracketed (organic subset) atom */
function implicitH(m: Mol, i: number): number {
  const a = m.atoms[i];
  let sum = 0; let aromNb = 0;
  for (const j of m.nb(i)) { const o = m.order(i, j); if (o === 1.5) { sum += 1; aromNb++; } else sum += o; }
  if (a.arom) {
    if (a.el === "C") return Math.max(0, 4 - (sum + (aromNb > 0 ? 1 : 0)));
    if (a.el === "N" || a.el === "O" || a.el === "S" || a.el === "Se") return 0;
    sum += 1;
  }
  const vals = VALENCE[a.el];
  if (!vals) return 0;
  for (const v of vals) if (v >= sum) return v - sum;
  return 0;
}

/** basic valence sanity check used to reject invalid predicted products */
function valenceOK(m: Mol, i: number): boolean {
  const a = m.atoms[i];
  if (!a.alive || a.arom) return true;
  const base = VALENCE[a.el];
  if (!base) return true;
  let sum = a.h;
  for (const j of m.nb(i)) { const o = m.order(i, j); sum += o === 1.5 ? 1 : o; }
  const allowed = base.map((v) => {
    if (a.el === "C" || a.el === "Si") return v - Math.abs(a.q);
    if (a.el === "B") return v - a.q;
    return v + a.q;
  }).filter((v) => v >= 0);
  return allowed.includes(sum) || (a.el === "H" && sum <= 1);
}
function molValid(m: Mol): boolean { return m.alive().every((i) => valenceOK(m, i)); }

/* ------------------------------------------------------------------ */
/*  Rings + aromaticity                                                */
/* ------------------------------------------------------------------ */
function bondInRing(m: Mol, a: number, b: number): boolean {
  const seen = new Set<number>([a]); const st = [a];
  while (st.length) {
    const u = st.pop()!;
    for (const v of m.nb(u)) {
      if ((u === a && v === b) || (u === b && v === a)) continue;
      if (v === b) return true;
      if (!seen.has(v)) { seen.add(v); st.push(v); }
    }
  }
  return false;
}

function findRings(m: Mol, maxSize = 8): number[][] {
  const out: number[][] = []; const seen = new Set<string>();
  for (const s of m.alive()) {
    const path = [s]; const onPath = new Set<number>([s]);
    const dfs = (u: number): void => {
      if (out.length > 400) return;
      for (const v of m.nb(u)) {
        if (v === s && path.length >= 3) {
          const key = [...path].sort((x, y) => x - y).join(",");
          if (!seen.has(key)) { seen.add(key); out.push([...path]); }
          continue;
        }
        if (v <= s || onPath.has(v) || path.length >= maxSize) continue;
        path.push(v); onPath.add(v); dfs(v); path.pop(); onPath.delete(v);
      }
    };
    dfs(s);
  }
  return out;
}

function ringSizeOfBond(rings: number[][], a: number, b: number): number {
  let best = 0;
  for (const r of rings) {
    for (let p = 0; p < r.length; p++) {
      const x = r[p], y = r[(p + 1) % r.length];
      if ((x === a && y === b) || (x === b && y === a)) { if (!best || r.length < best) best = r.length; }
    }
  }
  return best;
}
function ringSizeOfAtom(rings: number[][], a: number): number {
  let best = 0;
  for (const r of rings) if (r.includes(a) && (!best || r.length < best)) best = r.length;
  return best;
}

function isHuckelRing(m: Mol, ring: number[], rings: number[][]): boolean {
  if (ring.length !== 5 && ring.length !== 6) return false;
  const inRing = new Set(ring);
  let e = 0;
  for (const i of ring) {
    const a = m.atoms[i];
    if (a.arom) {
      if (a.el === "C") e += 1;
      else if (a.el === "O" || a.el === "S" || a.el === "Se") e += 2;
      else if (a.el === "N") e += (a.h > 0 || m.deg(i) + a.h >= 3) && a.q === 0 ? 2 : 1;
      else e += 1;
      continue;
    }
    let dbl = -1; let trip = false;
    for (const j of m.nb(i)) { const o = m.order(i, j); if (o === 2) dbl = j; if (o === 3) trip = true; }
    if (trip) return false;
    if (dbl >= 0) {
      if (inRing.has(dbl)) { e += 1; continue; }
      if (m.el(dbl) === "C" && ringSizeOfAtom(rings, dbl) > 0 && m.el(i) === "C") { e += 1; continue; }
      return false;
    }
    if (a.q === 0 && (a.el === "O" || a.el === "S" || (a.el === "N" && m.deg(i) + a.h === 3)) && ring.length === 5) { e += 2; continue; }
    if (a.el === "C" && a.q === -1 && ring.length === 5) { e += 2; continue; }
    return false;
  }
  return e === 6;
}

function perceiveAromaticity(m: Mol): void {
  const rings = findRings(m, 6).filter((r) => r.length === 5 || r.length === 6);
  const toSet: number[][] = [];
  for (const r of rings) {
    if (r.every((i) => m.atoms[i].arom)) continue;
    if (isHuckelRing(m, r, rings)) toSet.push(r);
  }
  for (const r of toSet) {
    for (const i of r) m.atoms[i].arom = true;
    for (let p = 0; p < r.length; p++) {
      const a = r[p], b = r[(p + 1) % r.length];
      if (m.order(a, b) > 0) m.bond(a, b, 1.5);
    }
  }
}

/* ------------------------------------------------------------------ */
/*  SMILES parser                                                      */
/* ------------------------------------------------------------------ */
export function parseSmiles(input: string): Mol | null {
  const s = (input || "").trim().split(/\s+/)[0];
  if (!s || s.length > 4000) return null;
  const m = new Mol();
  const stack: number[] = [];
  const ringOpen = new Map<number, { atom: number; order: number | null }>();
  let prev = -1; let pending: number | null = null; let i = 0;
  const link = (a: number, b: number, ex: number | null): void => {
    const o = ex !== null ? ex : (m.atoms[a].arom && m.atoms[b].arom ? 1.5 : 1);
    m.bond(a, b, o);
  };
  while (i < s.length) {
    const ch = s[i];
    if (ch === "(") { if (prev < 0) return null; stack.push(prev); i++; continue; }
    if (ch === ")") { if (!stack.length) return null; prev = stack.pop()!; i++; continue; }
    if (ch === ".") { prev = -1; pending = null; i++; continue; }
    if (ch === "$") return null;
    if ("-=#:/\\".includes(ch)) { pending = ch === "=" ? 2 : ch === "#" ? 3 : ch === ":" ? 1.5 : 1; i++; continue; }
    if ((ch >= "0" && ch <= "9") || ch === "%") {
      let num: number;
      if (ch === "%") { num = parseInt(s.substr(i + 1, 2), 10); i += 3; } else { num = parseInt(ch, 10); i++; }
      if (prev < 0 || Number.isNaN(num)) return null;
      const op = ringOpen.get(num);
      if (op) {
        if (op.atom === prev) return null;
        link(op.atom, prev, pending !== null ? pending : op.order);
        ringOpen.delete(num);
      } else ringOpen.set(num, { atom: prev, order: pending });
      pending = null; continue;
    }
    let idx = -1;
    if (ch === "[") {
      const end = s.indexOf("]", i);
      if (end < 0) return null;
      const body = s.slice(i + 1, end);
      const mm = /^(\d+)?([A-Za-z][a-z]?)(@{1,2}(?:TH\d|AL\d|SP\d|TB\d+|OH\d+)?)?(?:H(\d*))?(\+\d+|-\d+|\++|-+)?(?::\d+)?$/.exec(body);
      if (!mm) return null;
      const raw = mm[2];
      const arom = raw[0] === raw[0].toLowerCase();
      const el = raw[0].toUpperCase() + raw.slice(1);
      if (!(el in MASS)) return null;
      if (arom && !AROMATIC_ALLOWED.has(el)) return null;
      const h = body.includes("H") && mm[4] !== undefined ? (mm[4] === "" ? 1 : parseInt(mm[4], 10)) : 0;
      let q = 0; const qs = mm[5];
      if (qs) {
        if (/^[+-]\d+$/.test(qs)) q = (qs[0] === "+" ? 1 : -1) * parseInt(qs.slice(1), 10);
        else q = (qs[0] === "+" ? 1 : -1) * qs.length;
      }
      idx = m.add(el, { arom, q, h, bracket: true });
      i = end + 1;
    } else {
      let sym = ch; let arom = false;
      if (ch === "C" && s[i + 1] === "l") { sym = "Cl"; i += 2; }
      else if (ch === "B" && s[i + 1] === "r") { sym = "Br"; i += 2; }
      else if (ORGANIC.has(ch)) { i++; }
      else if ("bcnops".includes(ch)) { sym = ch.toUpperCase(); arom = true; i++; }
      else return null;
      idx = m.add(sym, { arom });
    }
    if (prev >= 0) link(prev, idx, pending);
    pending = null; prev = idx;
  }
  if (ringOpen.size || stack.length || m.atoms.length === 0) return null;

  // aromatic bonds outside rings are plain single bonds (e.g. biphenyl link)
  for (const a of m.alive()) for (const b of m.nb(a)) if (a < b && m.order(a, b) === 1.5 && !bondInRing(m, a, b)) m.bond(a, b, 1);
  for (const a of m.alive()) {
    if (m.atoms[a].arom && !m.nb(a).some((b) => m.order(a, b) === 1.5)) return null; // stray aromatic atom => not a real SMILES
  }
  // implicit hydrogens
  for (const a of m.alive()) if (!m.atoms[a].bracket) m.atoms[a].h = implicitH(m, a);
  // fold explicit [H] atoms
  for (const a of m.alive()) {
    const at = m.atoms[a];
    if (at.el === "H" && at.q === 0 && m.deg(a) === 1) {
      const n = m.nb(a)[0];
      if (m.el(n) !== "H") { m.atoms[n].h += 1; m.remove(a); }
    }
  }
  perceiveAromaticity(m);
  return m;
}

/* ------------------------------------------------------------------ */
/*  Canonical ranking + SMILES writer                                  */
/* ------------------------------------------------------------------ */
function rankKeys(keys: Map<number, string>): Map<number, number> {
  const uniq = [...new Set(keys.values())].sort();
  const idx = new Map<string, number>(); uniq.forEach((k, n) => idx.set(k, n));
  const out = new Map<number, number>();
  keys.forEach((k, i) => out.set(i, idx.get(k)!));
  return out;
}

function canonicalRanks(m: Mol): Map<number, number> {
  const atoms = m.alive();
  const keys = new Map<number, string>();
  for (const i of atoms) {
    const a = m.atoms[i];
    keys.set(i, `${a.el}|${a.arom ? 1 : 0}|${a.q}|${a.h}|${m.deg(i)}`);
  }
  let ranks = rankKeys(keys);
  let count = new Set(ranks.values()).size;
  for (let it = 0; it < atoms.length + 2; it++) {
    const nk = new Map<number, string>();
    for (const i of atoms) {
      const nbs = m.nb(i).map((j) => `${ranks.get(j)}/${m.order(i, j)}`).sort().join(",");
      nk.set(i, `${ranks.get(i)}:${nbs}`);
    }
    const nr = rankKeys(nk);
    const nc = new Set(nr.values()).size;
    ranks = nr;
    if (nc === count) break;
    count = nc;
  }
  return ranks;
}

function atomToken(m: Mol, i: number): string {
  const a = m.atoms[i];
  const sym = a.arom ? a.el.toLowerCase() : a.el;
  if (ORGANIC.has(a.el) && a.q === 0 && implicitH(m, i) === a.h) return sym;
  const hs = a.h > 0 ? "H" + (a.h > 1 ? a.h : "") : "";
  const qs = a.q === 0 ? "" : (a.q > 0 ? "+" : "-") + (Math.abs(a.q) > 1 ? Math.abs(a.q) : "");
  return `[${sym}${hs}${qs}]`;
}

function bondSymbol(m: Mol, u: number, v: number): string {
  const o = m.order(u, v);
  const both = m.atoms[u].arom && m.atoms[v].arom;
  if (o === 2) return "="; if (o === 3) return "#";
  if (o === 1.5) return both ? "" : ":";
  return both ? "-" : "";
}

function writeComponent(m: Mol, comp: number[], ranks: Map<number, number>): string {
  const rk = (i: number): number => ranks.get(i)!;
  const start = comp.reduce((b, i) => (rk(i) < rk(b) || (rk(i) === rk(b) && i < b) ? i : b), comp[0]);
  const visited = new Set<number>();
  const children = new Map<number, number[]>();
  const openAt = new Map<number, { d: number; sym: string }[]>();
  const closeAt = new Map<number, number[]>();
  const closureKeys = new Set<string>();
  let counter = 0;
  const dfs = (u: number, parent: number): void => {
    visited.add(u); children.set(u, []);
    const nbs = m.nb(u).sort((x, y) => rk(x) - rk(y) || x - y);
    for (const v of nbs) {
      if (v === parent) continue;
      if (!visited.has(v)) { children.get(u)!.push(v); dfs(v, u); }
      else {
        const key = u < v ? `${u}-${v}` : `${v}-${u}`;
        if (closureKeys.has(key)) continue;
        closureKeys.add(key);
        const d = ++counter;
        if (!openAt.has(v)) openAt.set(v, []);
        openAt.get(v)!.push({ d, sym: bondSymbol(m, v, u) });
        if (!closeAt.has(u)) closeAt.set(u, []);
        closeAt.get(u)!.push(d);
      }
    }
  };
  dfs(start, -1);
  const dig = (d: number): string => (d < 10 ? String(d) : `%${d}`);
  const emit = (u: number): string => {
    let out = atomToken(m, u);
    for (const r of openAt.get(u) ?? []) out += r.sym + dig(r.d);
    for (const d of closeAt.get(u) ?? []) out += dig(d);
    const ch = children.get(u)!;
    ch.forEach((v, n) => {
      const sub = bondSymbol(m, u, v) + emit(v);
      out += n < ch.length - 1 ? `(${sub})` : sub;
    });
    return out;
  };
  return emit(start);
}

/** canonical SMILES of one connected component, computed in isolation (independent of any other component) */
function writeIsolated(m: Mol, comp: number[]): string {
  const sub = subMol(m, comp);
  return writeComponent(sub, sub.alive(), canonicalRanks(sub));
}
export function writeSmiles(m: Mol): string {
  return components(m).map((c) => writeIsolated(m, c)).sort().join(".");
}
export function canonicalSmiles(smiles: string): string | null {
  const m = parseSmiles(smiles);
  return m ? writeSmiles(m) : null;
}

/* ------------------------------------------------------------------ */
/*  Formula + mass                                                     */
/* ------------------------------------------------------------------ */
export function formulaOf(m: Mol, only?: number[]): { formula: string; mono: number } {
  const cnt: Record<string, number> = {};
  let mono = 0;
  for (const i of only ?? m.alive()) {
    const a = m.atoms[i];
    cnt[a.el] = (cnt[a.el] || 0) + 1;
    if (a.h) cnt.H = (cnt.H || 0) + a.h;
    mono += (MASS[a.el] || 0) + a.h * MASS.H;
  }
  const keys = Object.keys(cnt).sort((x, y) => {
    if (cnt.C) { if (x === "C") return -1; if (y === "C") return 1; if (x === "H") return -1; if (y === "H") return 1; }
    return x < y ? -1 : 1;
  });
  return { formula: keys.map((k) => k + (cnt[k] > 1 ? cnt[k] : "")).join(""), mono };
}

/* ------------------------------------------------------------------ */
/*  General chemistry compound library (names -> canonical SMILES)     */
/* ------------------------------------------------------------------ */
// [aliases separated by "|" (first = primary key), SMILES, display name, category]
type LibRow = [string, string, string, string];
const LIB_ROWS: LibRow[] = [
  // ---- Solvents & simple neutral molecules
  ["water|h2o", "O", "Water", "Solvent"],
  ["methanol|meoh|methyl alcohol", "CO", "Methanol", "Alcohol"],
  ["ethanol|etoh|ethyl alcohol", "CCO", "Ethanol", "Alcohol"],
  ["propan-2-ol|isopropanol|isopropyl alcohol|2-propanol|ipa", "CC(C)O", "Propan-2-ol", "Alcohol"],
  ["1-propanol|n-propanol|propanol", "CCCO", "Propan-1-ol", "Alcohol"],
  ["1-butanol|n-butanol|butanol", "CCCCO", "Butan-1-ol", "Alcohol"],
  ["tert-butanol|t-butanol|2-methyl-2-propanol", "CC(C)(C)O", "tert-Butanol", "Alcohol"],
  ["cyclohexanol", "OC1CCCCC1", "Cyclohexanol", "Alcohol"],
  ["benzyl alcohol", "OCc1ccccc1", "Benzyl alcohol", "Alcohol"],
  ["ethylene glycol|1,2-ethanediol", "OCCO", "Ethylene glycol", "Polyol"],
  ["propylene glycol|1,2-propanediol", "CC(O)CO", "Propylene glycol", "Polyol"],
  ["glycerol|glycerin|glycerine", "OCC(O)CO", "Glycerol", "Polyol"],
  ["diethylene glycol", "OCCOCCO", "Diethylene glycol", "Polyether"],
  ["polyethylene glycol|peg", "OCCOCCOCCO", "Polyethylene glycol (triethylene glycol model)", "Polyether"],
  ["acetone|propanone", "CC(C)=O", "Acetone", "Ketone"],
  ["2-butanone|mek|methyl ethyl ketone", "CCC(C)=O", "Butan-2-one", "Ketone"],
  ["cyclohexanone", "O=C1CCCCC1", "Cyclohexanone", "Ketone"],
  ["acetophenone", "CC(=O)c1ccccc1", "Acetophenone", "Ketone"],
  ["benzophenone", "O=C(c1ccccc1)c1ccccc1", "Benzophenone", "Ketone"],
  ["acetonitrile|methyl cyanide", "CC#N", "Acetonitrile", "Nitrile"],
  ["benzonitrile", "N#Cc1ccccc1", "Benzonitrile", "Nitrile"],
  ["dimethylformamide|dmf|n,n-dimethylformamide", "CN(C)C=O", "N,N-Dimethylformamide", "Amide"],
  ["dimethylacetamide|dmac|n,n-dimethylacetamide", "CN(C)C(C)=O", "N,N-Dimethylacetamide", "Amide"],
  ["n-methyl-2-pyrrolidone|nmp|n-methylpyrrolidone", "CN1CCCC1=O", "N-Methyl-2-pyrrolidone", "Lactam"],
  ["dimethyl sulfoxide|dmso", "CS(C)=O", "Dimethyl sulfoxide", "Sulfoxide"],
  ["dimethyl sulfone", "CS(C)(=O)=O", "Dimethyl sulfone", "Sulfone"],
  ["tetrahydrofuran|thf", "C1CCOC1", "Tetrahydrofuran", "Ether"],
  ["1,4-dioxane|dioxane", "C1COCCO1", "1,4-Dioxane", "Ether"],
  ["diethyl ether|ether|ethoxyethane", "CCOCC", "Diethyl ether", "Ether"],
  ["methyl tert-butyl ether|mtbe", "COC(C)(C)C", "Methyl tert-butyl ether", "Ether"],
  ["anisole|methoxybenzene", "COc1ccccc1", "Anisole", "Ether"],
  ["ethyl acetate|etoac", "CCOC(C)=O", "Ethyl acetate", "Ester"],
  ["methyl acetate", "COC(C)=O", "Methyl acetate", "Ester"],
  ["butyl acetate|n-butyl acetate", "CCCCOC(C)=O", "Butyl acetate", "Ester"],
  ["tert-butyl acetate", "CC(=O)OC(C)(C)C", "tert-Butyl acetate", "Ester"],
  ["phenyl acetate", "CC(=O)Oc1ccccc1", "Phenyl acetate", "Ester"],
  ["vinyl acetate", "CC(=O)OC=C", "Vinyl acetate", "Ester"],
  ["methyl benzoate", "COC(=O)c1ccccc1", "Methyl benzoate", "Ester"],
  ["ethyl formate", "CCOC=O", "Ethyl formate", "Ester"],
  ["ethyl acrylate", "C=CC(=O)OCC", "Ethyl acrylate", "Michael acceptor"],
  ["methyl methacrylate|mma", "COC(=O)C(C)=C", "Methyl methacrylate", "Michael acceptor"],
  ["dimethyl carbonate", "COC(=O)OC", "Dimethyl carbonate", "Carbonate"],
  ["diethyl phthalate", "CCOC(=O)c1ccccc1C(=O)OCC", "Diethyl phthalate", "Ester"],
  ["triacetin|glyceryl triacetate", "CC(=O)OCC(COC(C)=O)OC(C)=O", "Triacetin", "Ester"],
  ["acetylsalicylic acid|aspirin|2-acetoxybenzoic acid", "CC(=O)Oc1ccccc1C(=O)O", "2-Acetoxybenzoic acid (acetylsalicylic acid)", "Aryl ester / acid"],
  ["dichloromethane|dcm|methylene chloride", "ClCCl", "Dichloromethane", "Alkyl halide"],
  ["chloroform|trichloromethane", "ClC(Cl)Cl", "Chloroform", "Alkyl halide"],
  ["carbon tetrachloride|tetrachloromethane", "ClC(Cl)(Cl)Cl", "Carbon tetrachloride", "Alkyl halide"],
  ["1,2-dichloroethane|dce", "ClCCCl", "1,2-Dichloroethane", "Alkyl halide"],
  ["hexane|n-hexane", "CCCCCC", "Hexane", "Hydrocarbon"],
  ["heptane|n-heptane", "CCCCCCC", "Heptane", "Hydrocarbon"],
  ["cyclohexane", "C1CCCCC1", "Cyclohexane", "Hydrocarbon"],
  ["methane", "C", "Methane", "Hydrocarbon"],
  ["ethylene|ethene", "C=C", "Ethene", "Alkene"],
  ["propene|propylene", "CC=C", "Propene", "Alkene"],
  ["isobutene|isobutylene", "CC(C)=C", "Isobutene", "Alkene"],
  ["1-hexene", "CCCCC=C", "1-Hexene", "Alkene"],
  ["cyclohexene", "C1=CCCCC1", "Cyclohexene", "Alkene"],
  ["styrene", "C=Cc1ccccc1", "Styrene", "Alkene"],
  ["1,3-butadiene|butadiene", "C=CC=C", "1,3-Butadiene", "Diene"],
  ["cyclopentadiene", "C1=CCC=C1", "Cyclopentadiene", "Diene"],
  ["acetylene|ethyne", "C#C", "Ethyne", "Alkyne"],
  ["phenylacetylene", "C#Cc1ccccc1", "Phenylacetylene", "Alkyne"],
  ["benzene", "c1ccccc1", "Benzene", "Aromatic"],
  ["toluene|methylbenzene", "Cc1ccccc1", "Toluene", "Aromatic"],
  ["o-xylene|xylene", "Cc1ccccc1C", "o-Xylene", "Aromatic"],
  ["ethylbenzene", "CCc1ccccc1", "Ethylbenzene", "Aromatic"],
  ["cumene|isopropylbenzene", "CC(C)c1ccccc1", "Cumene", "Aromatic"],
  ["naphthalene", "c1ccc2ccccc2c1", "Naphthalene", "Aromatic"],
  ["biphenyl", "c1ccc(cc1)-c1ccccc1", "Biphenyl", "Aromatic"],
  ["chlorobenzene", "Clc1ccccc1", "Chlorobenzene", "Aryl halide"],
  ["bromobenzene", "Brc1ccccc1", "Bromobenzene", "Aryl halide"],
  ["iodobenzene", "Ic1ccccc1", "Iodobenzene", "Aryl halide"],
  ["fluorobenzene", "Fc1ccccc1", "Fluorobenzene", "Aryl halide"],
  ["benzotrifluoride|trifluorotoluene", "FC(F)(F)c1ccccc1", "Benzotrifluoride", "Aryl halide"],
  ["1-chloro-2,4-dinitrobenzene|cdnb", "Clc1ccc(cc1[N+]([O-])=O)[N+]([O-])=O", "1-Chloro-2,4-dinitrobenzene", "Aryl halide"],
  ["nitrobenzene", "[O-][N+](=O)c1ccccc1", "Nitrobenzene", "Nitro"],
  ["4-nitrophenol|p-nitrophenol", "Oc1ccc(cc1)[N+]([O-])=O", "4-Nitrophenol", "Phenol"],
  ["nitromethane", "C[N+]([O-])=O", "Nitromethane", "Nitro"],
  ["nitroethane", "CC[N+]([O-])=O", "Nitroethane", "Nitro"],
  // ---- Carboxylic acids & derivatives
  ["acetic acid|ethanoic acid|aa", "CC(=O)O", "Acetic acid", "Carboxylic acid"],
  ["formic acid|methanoic acid", "OC=O", "Formic acid", "Carboxylic acid"],
  ["propionic acid|propanoic acid", "CCC(=O)O", "Propanoic acid", "Carboxylic acid"],
  ["butyric acid|butanoic acid", "CCCC(=O)O", "Butanoic acid", "Carboxylic acid"],
  ["benzoic acid", "OC(=O)c1ccccc1", "Benzoic acid", "Carboxylic acid"],
  ["salicylic acid|2-hydroxybenzoic acid", "OC(=O)c1ccccc1O", "Salicylic acid", "Carboxylic acid / phenol"],
  ["gentisic acid|2,5-dihydroxybenzoic acid", "OC(=O)c1cc(O)ccc1O", "Gentisic acid", "Carboxylic acid / phenol"],
  ["gallic acid", "OC(=O)c1cc(O)c(O)c(O)c1", "Gallic acid", "Carboxylic acid / phenol"],
  ["phthalic acid", "OC(=O)c1ccccc1C(O)=O", "Phthalic acid", "Carboxylic acid"],
  ["cinnamic acid", "OC(=O)C=Cc1ccccc1", "Cinnamic acid", "Michael acceptor"],
  ["oxalic acid", "OC(=O)C(O)=O", "Oxalic acid", "Carboxylic acid"],
  ["malonic acid", "OC(=O)CC(O)=O", "Malonic acid", "Carboxylic acid"],
  ["succinic acid", "OC(=O)CCC(O)=O", "Succinic acid", "Carboxylic acid"],
  ["adipic acid", "OC(=O)CCCCC(O)=O", "Adipic acid", "Carboxylic acid"],
  ["maleic acid", "OC(=O)C=CC(O)=O", "Maleic acid", "Michael acceptor"],
  ["fumaric acid", "OC(=O)C=CC(O)=O", "Fumaric acid", "Michael acceptor"],
  ["acetoacetic acid|3-oxobutanoic acid", "CC(=O)CC(O)=O", "Acetoacetic acid", "Keto acid"],
  ["pyruvic acid", "CC(=O)C(O)=O", "Pyruvic acid", "Keto acid"],
  ["lactic acid", "CC(O)C(O)=O", "Lactic acid", "Hydroxy acid"],
  ["glycolic acid", "OCC(O)=O", "Glycolic acid", "Hydroxy acid"],
  ["tartaric acid", "OC(C(O)C(O)=O)C(O)=O", "Tartaric acid", "Hydroxy acid"],
  ["citric acid", "OC(=O)CC(O)(CC(O)=O)C(O)=O", "Citric acid", "Hydroxy acid"],
  ["ascorbic acid|vitamin c", "OCC(O)C1OC(=O)C(O)=C1O", "Ascorbic acid", "Enediol lactone"],
  ["trifluoroacetic acid|tfa", "OC(=O)C(F)(F)F", "Trifluoroacetic acid", "Carboxylic acid"],
  ["chloroacetic acid", "OC(=O)CCl", "Chloroacetic acid", "Carboxylic acid"],
  ["trichloroacetic acid|tca", "OC(=O)C(Cl)(Cl)Cl", "Trichloroacetic acid", "Carboxylic acid"],
  ["lauric acid|dodecanoic acid", "CCCCCCCCCCCC(O)=O", "Lauric acid", "Fatty acid"],
  ["palmitic acid|hexadecanoic acid", "CCCCCCCCCCCCCCCC(O)=O", "Palmitic acid", "Fatty acid"],
  ["stearic acid|octadecanoic acid", "CCCCCCCCCCCCCCCCCC(O)=O", "Stearic acid", "Fatty acid"],
  ["oleic acid", "CCCCCCCCC=CCCCCCCCC(O)=O", "Oleic acid", "Fatty acid"],
  ["linoleic acid", "CCCCCC=CCC=CCCCCCCCC(O)=O", "Linoleic acid", "Fatty acid"],
  ["methyl oleate", "CCCCCCCCC=CCCCCCCCC(=O)OC", "Methyl oleate", "Ester"],
  ["sodium acetate", "CC([O-])=O.[Na+]", "Sodium acetate", "Carboxylate salt"],
  ["sodium benzoate", "[O-]C(=O)c1ccccc1.[Na+]", "Sodium benzoate", "Carboxylate salt"],
  ["sodium citrate|trisodium citrate", "OC(CC([O-])=O)(CC([O-])=O)C([O-])=O.[Na+].[Na+].[Na+]", "Trisodium citrate", "Carboxylate salt"],
  ["magnesium stearate", "[Mg+2].[O-]C(=O)CCCCCCCCCCCCCCCCC.[O-]C(=O)CCCCCCCCCCCCCCCCC", "Magnesium stearate", "Metal carboxylate"],
  ["calcium stearate", "[Ca+2].[O-]C(=O)CCCCCCCCCCCCCCCCC.[O-]C(=O)CCCCCCCCCCCCCCCCC", "Calcium stearate", "Metal carboxylate"],
  ["sodium stearate", "[Na+].[O-]C(=O)CCCCCCCCCCCCCCCCC", "Sodium stearate", "Metal carboxylate"],
  ["sodium stearyl fumarate", "CCCCCCCCCCCCCCCCCCOC(=O)C=CC([O-])=O.[Na+]", "Sodium stearyl fumarate", "Ester / carboxylate"],
  ["acetic anhydride", "CC(=O)OC(C)=O", "Acetic anhydride", "Anhydride"],
  ["succinic anhydride", "O=C1CCC(=O)O1", "Succinic anhydride", "Anhydride"],
  ["maleic anhydride", "O=C1C=CC(=O)O1", "Maleic anhydride", "Anhydride"],
  ["phthalic anhydride", "O=C1OC(=O)c2ccccc12", "Phthalic anhydride", "Anhydride"],
  ["acetyl chloride", "CC(Cl)=O", "Acetyl chloride", "Acyl halide"],
  ["benzoyl chloride", "O=C(Cl)c1ccccc1", "Benzoyl chloride", "Acyl halide"],
  ["benzoyl peroxide", "O=C(OOC(=O)c1ccccc1)c1ccccc1", "Benzoyl peroxide", "Peroxide"],
  ["gamma-butyrolactone|butyrolactone|gbl", "O=C1CCCO1", "gamma-Butyrolactone", "Lactone"],
  ["delta-valerolactone|valerolactone", "O=C1CCCCO1", "delta-Valerolactone", "Lactone"],
  ["beta-propiolactone|propiolactone", "O=C1CCO1", "beta-Propiolactone", "Lactone"],
  ["epsilon-caprolactone|caprolactone", "O=C1CCCCCO1", "epsilon-Caprolactone", "Lactone"],
  ["coumarin", "O=c1ccc2ccccc2o1", "Coumarin", "Lactone"],
  ["caprolactam|epsilon-caprolactam", "O=C1CCCCCN1", "epsilon-Caprolactam", "Lactam"],
  ["2-pyrrolidone|pyrrolidone|gamma-butyrolactam", "O=C1CCCN1", "2-Pyrrolidone", "Lactam"],
  ["2-azetidinone|beta-lactam", "O=C1CCN1", "2-Azetidinone (beta-lactam)", "Lactam"],
  ["n-vinylpyrrolidone|nvp|polyvinylpyrrolidone|pvp|povidone", "C=CN1CCCC1=O", "N-Vinylpyrrolidone (PVP repeat-unit model)", "Lactam"],
  ["succinimide", "O=C1CCC(=O)N1", "Succinimide", "Imide"],
  ["phthalimide", "O=C1NC(=O)c2ccccc12", "Phthalimide", "Imide"],
  ["maleimide", "O=C1C=CC(=O)N1", "Maleimide", "Imide"],
  ["acetamide", "CC(N)=O", "Acetamide", "Amide"],
  ["formamide", "NC=O", "Formamide", "Amide"],
  ["benzamide", "NC(=O)c1ccccc1", "Benzamide", "Amide"],
  ["acetanilide", "CC(=O)Nc1ccccc1", "Acetanilide", "Anilide"],
  ["n-(4-hydroxyphenyl)acetamide|4-acetamidophenol|paracetamol|acetaminophen", "CC(=O)Nc1ccc(O)cc1", "N-(4-Hydroxyphenyl)acetamide", "Anilide / phenol"],
  ["urea|carbamide", "NC(N)=O", "Urea", "Urea"],
  ["thiourea", "NC(N)=S", "Thiourea", "Thiocarbonyl"],
  ["ethyl carbamate|urethane", "CCOC(N)=O", "Ethyl carbamate", "Carbamate"],
  ["di-tert-butyl dicarbonate|boc anhydride|boc2o", "CC(C)(C)OC(=O)OC(=O)OC(C)(C)C", "Di-tert-butyl dicarbonate", "Anhydride"],
  ["tert-butyl carbamate|boc-nh2", "CC(C)(C)OC(N)=O", "tert-Butyl carbamate", "Carbamate"],
  // ---- Aldehydes / carbonyl / quinones
  ["formaldehyde|methanal", "C=O", "Formaldehyde", "Aldehyde"],
  ["acetaldehyde|ethanal", "CC=O", "Acetaldehyde", "Aldehyde"],
  ["propionaldehyde|propanal", "CCC=O", "Propanal", "Aldehyde"],
  ["benzaldehyde", "O=Cc1ccccc1", "Benzaldehyde", "Aldehyde"],
  ["glyoxal", "O=CC=O", "Glyoxal", "Aldehyde"],
  ["acrolein|propenal", "C=CC=O", "Acrolein", "Michael acceptor"],
  ["methyl vinyl ketone", "CC(=O)C=C", "Methyl vinyl ketone", "Michael acceptor"],
  ["vanillin", "COc1cc(C=O)ccc1O", "Vanillin", "Aldehyde / phenol"],
  ["benzoquinone|1,4-benzoquinone|p-benzoquinone", "O=C1C=CC(=O)C=C1", "1,4-Benzoquinone", "Quinone"],
  // ---- Carbohydrates
  ["glucose|dextrose|d-glucose", "OCC1OC(O)C(O)C(O)C1O", "Glucose (pyranose)", "Reducing sugar"],
  ["fructose|levulose", "OCC1(O)OCC(O)C(O)C1O", "Fructose (pyranose)", "Reducing sugar"],
  ["galactose", "OCC1OC(O)C(O)C(O)C1O", "Galactose (pyranose)", "Reducing sugar"],
  ["lactose|lactose monohydrate", "OCC1OC(OC2C(CO)OC(O)C(O)C2O)C(O)C(O)C1O", "Lactose", "Reducing disaccharide"],
  ["maltose", "OCC1OC(OC2C(CO)OC(O)C(O)C2O)C(O)C(O)C1O", "Maltose", "Reducing disaccharide"],
  ["sucrose|table sugar", "OCC1OC(OC2(CO)OC(CO)C(O)C2O)C(O)C(O)C1O", "Sucrose", "Non-reducing disaccharide"],
  ["mannitol", "OCC(O)C(O)C(O)C(O)CO", "Mannitol", "Polyol"],
  ["sorbitol", "OCC(O)C(O)C(O)C(O)CO", "Sorbitol", "Polyol"],
  ["xylitol", "OCC(O)C(O)C(O)CO", "Xylitol", "Polyol"],
  ["cellulose|microcrystalline cellulose|mcc|starch|corn starch", "OCC1OC(O)C(O)C(O)C1O", "Glucose repeat unit (cellulose / starch model)", "Carbohydrate"],
  // ---- Amines / amino acids
  ["ammonia", "N", "Ammonia", "Amine"],
  ["methylamine", "CN", "Methylamine", "Amine"],
  ["dimethylamine", "CNC", "Dimethylamine", "Amine"],
  ["trimethylamine", "CN(C)C", "Trimethylamine", "Amine"],
  ["ethylamine", "CCN", "Ethylamine", "Amine"],
  ["diethylamine", "CCNCC", "Diethylamine", "Amine"],
  ["triethylamine|tea|et3n", "CCN(CC)CC", "Triethylamine", "Amine"],
  ["diisopropylethylamine|dipea|hunig's base", "CCN(C(C)C)C(C)C", "N,N-Diisopropylethylamine", "Amine"],
  ["ethylenediamine", "NCCN", "Ethylenediamine", "Amine"],
  ["ethanolamine|2-aminoethanol", "NCCO", "Ethanolamine", "Amino alcohol"],
  ["diethanolamine", "OCCNCCO", "Diethanolamine", "Amino alcohol"],
  ["triethanolamine", "OCCN(CCO)CCO", "Triethanolamine", "Amino alcohol"],
  ["tris|tromethamine|tris(hydroxymethyl)aminomethane", "NC(CO)(CO)CO", "Tris(hydroxymethyl)aminomethane", "Amino alcohol"],
  ["meglumine", "CNCC(O)C(O)C(O)C(O)CO", "N-Methylglucamine (meglumine)", "Amino polyol"],
  ["piperidine", "C1CCNCC1", "Piperidine", "Amine"],
  ["morpholine", "C1COCCN1", "Morpholine", "Amine"],
  ["piperazine", "C1CNCCN1", "Piperazine", "Amine"],
  ["pyrrolidine", "C1CCNC1", "Pyrrolidine", "Amine"],
  ["dbu|1,8-diazabicyclo[5.4.0]undec-7-ene", "C1CCC2=NCCCN2CC1", "DBU", "Amidine base"],
  ["benzylamine", "NCc1ccccc1", "Benzylamine", "Amine"],
  ["aniline|aminobenzene", "Nc1ccccc1", "Aniline", "Aromatic amine"],
  ["n-methylaniline", "CNc1ccccc1", "N-Methylaniline", "Aromatic amine"],
  ["n,n-dimethylaniline", "CN(C)c1ccccc1", "N,N-Dimethylaniline", "Aromatic amine"],
  ["4-aminophenol|p-aminophenol", "Nc1ccc(O)cc1", "4-Aminophenol", "Aromatic amine / phenol"],
  ["4-nitroaniline|p-nitroaniline", "Nc1ccc(cc1)[N+]([O-])=O", "4-Nitroaniline", "Aromatic amine"],
  ["guanidine", "NC(N)=N", "Guanidine", "Guanidine"],
  ["hydrazine", "NN", "Hydrazine", "Hydrazine"],
  ["phenylhydrazine", "NNc1ccccc1", "Phenylhydrazine", "Hydrazine"],
  ["hydroxylamine", "NO", "Hydroxylamine", "Hydroxylamine"],
  ["acetone oxime", "CC(C)=NO", "Acetone oxime", "Oxime"],
  ["benzaldehyde oxime", "ON=Cc1ccccc1", "Benzaldehyde oxime", "Oxime"],
  ["n-nitrosodimethylamine|ndma", "CN(C)N=O", "N-Nitrosodimethylamine", "N-Nitrosamine"],
  ["n-nitrosodiethylamine|ndea", "CCN(CC)N=O", "N-Nitrosodiethylamine", "N-Nitrosamine"],
  ["azobenzene", "c1ccc(cc1)N=Nc1ccccc1", "Azobenzene", "Azo"],
  ["phenyl azide", "[N-]=[N+]=Nc1ccccc1", "Phenyl azide", "Azide"],
  ["methyl isocyanate", "CN=C=O", "Methyl isocyanate", "Isocyanate"],
  ["phenyl isothiocyanate", "S=C=Nc1ccccc1", "Phenyl isothiocyanate", "Isothiocyanate"],
  ["glycine", "NCC(O)=O", "Glycine", "Amino acid"],
  ["alanine", "CC(N)C(O)=O", "Alanine", "Amino acid"],
  ["lysine", "NCCCCC(N)C(O)=O", "Lysine", "Amino acid"],
  ["glutamic acid", "NC(CCC(O)=O)C(O)=O", "Glutamic acid", "Amino acid"],
  ["cysteine", "NC(CS)C(O)=O", "Cysteine", "Amino acid"],
  ["methionine", "CSCCC(N)C(O)=O", "Methionine", "Amino acid"],
  ["phenylalanine", "NC(Cc1ccccc1)C(O)=O", "Phenylalanine", "Amino acid"],
  ["tyrosine", "NC(Cc1ccc(O)cc1)C(O)=O", "Tyrosine", "Amino acid"],
  ["tryptophan", "NC(Cc1c[nH]c2ccccc12)C(O)=O", "Tryptophan", "Amino acid"],
  ["proline", "OC(=O)C1CCCN1", "Proline", "Amino acid"],
  ["cetyltrimethylammonium bromide|ctab", "CCCCCCCCCCCCCCCC[N+](C)(C)C.[Br-]", "Cetyltrimethylammonium bromide", "Quaternary ammonium"],
  ["benzalkonium chloride", "CCCCCCCCCCCC[N+](C)(C)Cc1ccccc1.[Cl-]", "Benzalkonium chloride (C12)", "Quaternary ammonium"],
  ["tetramethylammonium hydroxide|tmah", "C[N+](C)(C)C.[OH-]", "Tetramethylammonium hydroxide", "Strong base"],
  // ---- Phenols / oxygenates / aromatics
  ["phenol|hydroxybenzene", "Oc1ccccc1", "Phenol", "Phenol"],
  ["catechol|1,2-dihydroxybenzene", "Oc1ccccc1O", "Catechol", "Phenol"],
  ["resorcinol", "Oc1cccc(O)c1", "Resorcinol", "Phenol"],
  ["hydroquinone|1,4-dihydroxybenzene", "Oc1ccc(O)cc1", "Hydroquinone", "Phenol"],
  ["p-cresol|4-methylphenol", "Cc1ccc(O)cc1", "p-Cresol", "Phenol"],
  ["butylated hydroxytoluene|bht", "Cc1cc(c(O)c(c1)C(C)(C)C)C(C)(C)C", "Butylated hydroxytoluene", "Hindered phenol"],
  ["indole", "c1ccc2[nH]ccc2c1", "Indole", "Heteroaromatic"],
  ["pyridine", "c1ccncc1", "Pyridine", "Heteroaromatic"],
  ["pyrrole", "c1cc[nH]c1", "Pyrrole", "Heteroaromatic"],
  ["furan", "c1ccoc1", "Furan", "Heteroaromatic"],
  ["thiophene", "c1ccsc1", "Thiophene", "Heteroaromatic"],
  ["imidazole", "c1c[nH]cn1", "Imidazole", "Heteroaromatic"],
  ["caffeine|1,3,7-trimethylxanthine", "Cn1cnc2c1c(=O)n(C)c(=O)n2C", "1,3,7-Trimethylxanthine", "Heteroaromatic"],
  ["benzyl chloride", "ClCc1ccccc1", "Benzyl chloride", "Alkyl halide"],
  // ---- Halides / epoxides / alkylating agents
  ["chloroethane|ethyl chloride", "CCCl", "Chloroethane", "Alkyl halide"],
  ["1-bromobutane|butyl bromide", "CCCCBr", "1-Bromobutane", "Alkyl halide"],
  ["tert-butyl chloride", "CC(C)(C)Cl", "tert-Butyl chloride", "Alkyl halide"],
  ["allyl bromide", "C=CCBr", "Allyl bromide", "Alkyl halide"],
  ["iodomethane|methyl iodide", "CI", "Iodomethane", "Alkyl halide"],
  ["ethylene oxide|oxirane", "C1CO1", "Ethylene oxide", "Epoxide"],
  ["propylene oxide|methyloxirane", "CC1CO1", "Propylene oxide", "Epoxide"],
  ["epichlorohydrin", "ClCC1CO1", "Epichlorohydrin", "Epoxide"],
  ["methyl tosylate", "COS(=O)(=O)c1ccc(C)cc1", "Methyl p-toluenesulfonate", "Sulfonate ester"],
  ["dimethyl sulfate", "COS(=O)(=O)OC", "Dimethyl sulfate", "Sulfate ester"],
  ["methyl methanesulfonate|mms", "COS(C)(=O)=O", "Methyl methanesulfonate", "Sulfonate ester"],
  ["ethyl methanesulfonate|ems", "CCOS(C)(=O)=O", "Ethyl methanesulfonate", "Sulfonate ester"],
  // ---- Sulfur / phosphorus
  ["methanethiol|methyl mercaptan", "CS", "Methanethiol", "Thiol"],
  ["ethanethiol", "CCS", "Ethanethiol", "Thiol"],
  ["thiophenol", "Sc1ccccc1", "Thiophenol", "Thiol"],
  ["dithiothreitol|dtt", "OC(CS)C(O)CS", "Dithiothreitol", "Thiol"],
  ["dimethyl sulfide", "CSC", "Dimethyl sulfide", "Thioether"],
  ["dimethyl disulfide", "CSSC", "Dimethyl disulfide", "Disulfide"],
  ["diphenyl disulfide", "c1ccc(cc1)SSc1ccccc1", "Diphenyl disulfide", "Disulfide"],
  ["benzenesulfonamide", "NS(=O)(=O)c1ccccc1", "Benzenesulfonamide", "Sulfonamide"],
  ["p-toluenesulfonamide|tosylamide", "Cc1ccc(cc1)S(N)(=O)=O", "p-Toluenesulfonamide", "Sulfonamide"],
  ["methanesulfonic acid|msa", "CS(O)(=O)=O", "Methanesulfonic acid", "Sulfonic acid"],
  ["p-toluenesulfonic acid|tsoh|tosic acid", "Cc1ccc(cc1)S(O)(=O)=O", "p-Toluenesulfonic acid", "Sulfonic acid"],
  ["sodium lauryl sulfate|sodium dodecyl sulfate|sls|sds", "CCCCCCCCCCCCOS([O-])(=O)=O.[Na+]", "Sodium dodecyl sulfate", "Alkyl sulfate"],
  ["triphenylphosphine|pph3", "c1ccc(cc1)P(c1ccccc1)c1ccccc1", "Triphenylphosphine", "Phosphine"],
  ["trimethyl phosphate", "COP(=O)(OC)OC", "Trimethyl phosphate", "Phosphate ester"],
  // ---- Mineral acids / bases / salts / redox reagents
  ["hydrochloric acid|hcl|hydrogen chloride", "Cl", "Hydrochloric acid", "Strong acid"],
  ["hydrobromic acid|hbr", "Br", "Hydrobromic acid", "Strong acid"],
  ["sulfuric acid|h2so4", "OS(O)(=O)=O", "Sulfuric acid", "Strong acid"],
  ["nitric acid|hno3", "O[N+]([O-])=O", "Nitric acid", "Strong acid / oxidant"],
  ["phosphoric acid|h3po4", "OP(O)(O)=O", "Phosphoric acid", "Acid"],
  ["sodium hydroxide|naoh|caustic soda", "[Na+].[OH-]", "Sodium hydroxide", "Strong base"],
  ["potassium hydroxide|koh", "[K+].[OH-]", "Potassium hydroxide", "Strong base"],
  ["lithium hydroxide|lioh", "[Li+].[OH-]", "Lithium hydroxide", "Strong base"],
  ["calcium hydroxide|slaked lime", "[Ca+2].[OH-].[OH-]", "Calcium hydroxide", "Base"],
  ["magnesium hydroxide", "[Mg+2].[OH-].[OH-]", "Magnesium hydroxide", "Base"],
  ["ammonium hydroxide|ammonia solution", "[NH4+].[OH-]", "Ammonium hydroxide", "Base"],
  ["magnesium oxide|magnesia|mgo", "[Mg+2].[O-2]", "Magnesium oxide", "Metal oxide (basic)"],
  ["calcium oxide|quicklime|cao", "[Ca+2].[O-2]", "Calcium oxide", "Metal oxide (basic)"],
  ["zinc oxide|zno", "[Zn+2].[O-2]", "Zinc oxide", "Metal oxide"],
  ["silicon dioxide|silica|colloidal silicon dioxide|sio2", "O=[Si]=O", "Silicon dioxide", "Inert oxide"],
  ["titanium dioxide|tio2", "O=[Ti]=O", "Titanium dioxide", "Metal oxide (photocatalyst)"],
  ["sodium bicarbonate|sodium hydrogen carbonate|baking soda", "[Na+].OC([O-])=O", "Sodium bicarbonate", "Carbonate / base"],
  ["sodium carbonate|soda ash", "[Na+].[Na+].[O-]C([O-])=O", "Sodium carbonate", "Carbonate / base"],
  ["potassium carbonate", "[K+].[K+].[O-]C([O-])=O", "Potassium carbonate", "Carbonate / base"],
  ["calcium carbonate|chalk", "[Ca+2].[O-]C([O-])=O", "Calcium carbonate", "Carbonate / base"],
  ["magnesium carbonate", "[Mg+2].[O-]C([O-])=O", "Magnesium carbonate", "Carbonate / base"],
  ["sodium methoxide", "C[O-].[Na+]", "Sodium methoxide", "Strong base"],
  ["sodium ethoxide", "CC[O-].[Na+]", "Sodium ethoxide", "Strong base"],
  ["potassium tert-butoxide|kotbu", "CC(C)(C)[O-].[K+]", "Potassium tert-butoxide", "Strong base"],
  ["sodium chloride|nacl|salt", "[Na+].[Cl-]", "Sodium chloride", "Salt"],
  ["potassium chloride|kcl", "[K+].[Cl-]", "Potassium chloride", "Salt"],
  ["magnesium chloride", "[Mg+2].[Cl-].[Cl-]", "Magnesium chloride", "Salt"],
  ["calcium chloride", "[Ca+2].[Cl-].[Cl-]", "Calcium chloride", "Salt"],
  ["zinc chloride", "[Zn+2].[Cl-].[Cl-]", "Zinc chloride", "Lewis acid"],
  ["aluminium chloride|aluminum chloride|alcl3", "[Al+3].[Cl-].[Cl-].[Cl-]", "Aluminium chloride", "Lewis acid"],
  ["iron(iii) chloride|ferric chloride|fecl3", "[Fe+3].[Cl-].[Cl-].[Cl-]", "Iron(III) chloride", "Lewis acid / oxidant"],
  ["iron(ii) sulfate|ferrous sulfate|feso4", "[Fe+2].[O-]S([O-])(=O)=O", "Iron(II) sulfate", "Transition metal salt"],
  ["copper(ii) sulfate|cupric sulfate|cuso4", "[Cu+2].[O-]S([O-])(=O)=O", "Copper(II) sulfate", "Transition metal salt"],
  ["sodium sulfate", "[Na+].[Na+].[O-]S([O-])(=O)=O", "Sodium sulfate", "Salt"],
  ["magnesium sulfate|mgso4", "[Mg+2].[O-]S([O-])(=O)=O", "Magnesium sulfate", "Salt"],
  ["calcium sulfate|gypsum", "[Ca+2].[O-]S([O-])(=O)=O", "Calcium sulfate", "Salt"],
  ["sodium phosphate|disodium hydrogen phosphate", "[Na+].[Na+].OP([O-])([O-])=O", "Disodium hydrogen phosphate", "Phosphate salt"],
  ["calcium phosphate|dibasic calcium phosphate|calcium hydrogen phosphate", "[Ca+2].OP([O-])([O-])=O", "Calcium hydrogen phosphate", "Phosphate salt"],
  ["ammonium chloride", "[NH4+].[Cl-]", "Ammonium chloride", "Salt"],
  ["potassium iodide|ki", "[K+].[I-]", "Potassium iodide", "Halide salt"],
  ["sodium bromide", "[Na+].[Br-]", "Sodium bromide", "Halide salt"],
  ["sodium fluoride", "[Na+].[F-]", "Sodium fluoride", "Halide salt"],
  ["sodium azide", "[Na+].[N-]=[N+]=[N-]", "Sodium azide", "Azide"],
  ["hydrogen peroxide|h2o2", "OO", "Hydrogen peroxide", "Oxidant"],
  ["tert-butyl hydroperoxide|tbhp", "CC(C)(C)OO", "tert-Butyl hydroperoxide", "Oxidant"],
  ["cumene hydroperoxide", "CC(C)(OO)c1ccccc1", "Cumene hydroperoxide", "Oxidant"],
  ["peracetic acid", "CC(=O)OO", "Peracetic acid", "Oxidant"],
  ["m-chloroperbenzoic acid|mcpba", "OOC(=O)c1cccc(Cl)c1", "m-Chloroperbenzoic acid", "Oxidant"],
  ["sodium hypochlorite|bleach", "[Na+].[O-]Cl", "Sodium hypochlorite", "Oxidant"],
  ["potassium permanganate", "[K+].[O-][Mn](=O)(=O)=O", "Potassium permanganate", "Oxidant"],
  ["oxygen|dioxygen|o2", "O=O", "Molecular oxygen", "Oxidant"],
  ["ozone|o3", "[O-][O+]=O", "Ozone", "Oxidant"],
  ["sodium nitrite|nano2", "[Na+].[O-]N=O", "Sodium nitrite", "Nitrosating agent"],
  ["sodium nitrate", "[Na+].[O-][N+]([O-])=O", "Sodium nitrate", "Nitrate salt"],
  ["sodium borohydride|nabh4", "[Na+].[BH4-]", "Sodium borohydride", "Hydride reductant"],
  ["sodium bisulfite|sodium hydrogen sulfite", "[Na+].OS([O-])=O", "Sodium bisulfite", "Reductant"],
  ["sodium sulfite", "[Na+].[Na+].[O-]S([O-])=O", "Sodium sulfite", "Reductant"],
  ["sodium metabisulfite", "[Na+].[Na+].[O-]S(=O)OS([O-])(=O)=O", "Sodium metabisulfite", "Reductant"],
];

interface LibEntry { name: string; smiles: string; category: string; canonical: string; }
const LIBRARY = new Map<string, LibEntry>();
const CANON_NAMES = new Map<string, string>();

function normName(s: string): string {
  return (s || "").toLowerCase().trim()
    .replace(/[_,;]+/g, " ").replace(/\s*[-–]\s*/g, "-").replace(/\s+/g, " ")
    .replace(/\b(anhydrous|monohydrate|dihydrate|trihydrate|hexahydrate|hydrate|usp|bp|ip|ph\.? ?eur\.?|powder|solution|aqueous|reagent grade|pure)\b/g, "")
    .replace(/\s+/g, " ").trim();
}

/** Register (or override) a compound so names typed by the user resolve to a SMILES. */
export function registerCompound(names: string | string[], smiles: string, displayName?: string, category = "User-defined"): boolean {
  const canon = canonicalSmiles(smiles);
  if (!canon) return false;
  const list = (Array.isArray(names) ? names : [names]).map(normName).filter(Boolean);
  if (!list.length) return false;
  const entry: LibEntry = { name: displayName || list[0], smiles, category, canonical: canon };
  list.forEach((n) => LIBRARY.set(n, entry));
  const prevName = CANON_NAMES.get(canon);
  if (!prevName) CANON_NAMES.set(canon, entry.name);
  else if (!prevName.split(" / ").includes(entry.name)) CANON_NAMES.set(canon, `${prevName} / ${entry.name}`); // stereo-blind: keep both names
  return true;
}

for (const [aliases, smiles, display, cat] of LIB_ROWS) {
  registerCompound(aliases.split("|"), smiles, display, cat);
}

/**
 * Backward-compatible export (older code imported PHARMA_COMPOUNDS). It is now a
 * neutral general-chemistry table; the key names of the old table still resolve.
 * @deprecated use COMPOUND_LIBRARY / lookupCompoundSmiles
 */
export const COMPOUND_LIBRARY: Record<string, { smiles: string; name: string; category?: string }> = {};
LIBRARY.forEach((e, k) => { COMPOUND_LIBRARY[k] = { smiles: e.smiles, name: e.name, category: e.category }; });
export const PHARMA_COMPOUNDS = COMPOUND_LIBRARY;

/**
 * Resolve a chemical name / alias to a SMILES string. Exact (normalised) match first,
 * then whole-word containment (longest alias wins). Never a loose substring match.
 */
export function lookupCompoundSmiles(name: string): string | null {
  if (!name) return null;
  const clean = normName(name);
  if (!clean) return null;
  const direct = LIBRARY.get(clean) || LIBRARY.get(clean.replace(/[\s-]/g, "")) ;
  if (direct) return direct.smiles;
  const noSpace = clean.replace(/[\s()-]/g, "");
  for (const [k, e] of LIBRARY) if (k.replace(/[\s()-]/g, "") === noSpace) return e.smiles;
  const padded = ` ${clean} `;
  let best: LibEntry | null = null; let bestLen = 0;
  for (const [k, e] of LIBRARY) {
    if (k.length < 4 || k.length <= bestLen) continue;
    if (padded.includes(` ${k} `)) { best = e; bestLen = k.length; }
  }
  return best ? best.smiles : null;
}

/** Friendly name for a SMILES if it matches a library compound. */
function libraryNameFor(smiles: string): string | null {
  const c = canonicalSmiles(smiles);
  const n = c ? CANON_NAMES.get(c) : undefined;
  if (!n) return null;
  return n.includes(" / ") ? `${n.split(" / ")[0]} (or a stereoisomer)` : n;
}

/* ------------------------------------------------------------------ */
/*  Functional-group reactivity profiles (data)                        */
/* ------------------------------------------------------------------ */
const CR: Vuln = "Critical", HI: Vuln = "High", MO: Vuln = "Moderate", LO: Vuln = "Low", RS: Vuln = "Resistant";
type T = [Vuln, string];
interface Prof {
  name: string; cat: string; frag: string; site: string;
  acid: T; base: T; hyd: T; photo: T; therm: T; ox: T;
  sec: [Vuln, string, string];
  partners: string[];   // reactive classes of a partner molecule this group reacts with
  provides: string[];   // reactive classes this group contributes to its own molecule
}
const P = (name: string, cat: string, frag: string, site: string, acid: T, base: T, hyd: T, photo: T, therm: T, ox: T,
  sec: [Vuln, string, string], partners: string[] = [], provides: string[] = []): Prof =>
  ({ name, cat, frag, site, acid, base, hyd, photo, therm, ox, sec, partners, provides });

const PROFILES: Record<string, Prof> = {
  // ---------------- carboxylic acid family ----------------
  acid: P("Carboxylic Acid", "Carboxylic Acid", "-C(=O)OH", "Carboxyl proton and carbonyl carbon",
    [LO, "Stays un-ionised; only Fischer esterification with alcohols under strong acid / heat."],
    [HI, "Immediate deprotonation to the carboxylate salt (pKa ~4-5); no covalent breakdown."],
    [RS, "Hydrolytically inert terminus."],
    [LO, "Weak UV absorber unless conjugated; alpha-aryl acids can photodecarboxylate."],
    [LO, "Simple acids resist decarboxylation; beta-keto, malonic, alpha-EWG and electron-rich aryl acids lose CO2 on heating."],
    [LO, "Stable to O2; radical decarboxylation only with strong oxidants / photocatalysis."],
    [HI, "Bases, amines, alcohols, metal oxides", "Acid-base salt formation with bases/amines; esterification with alcohols and amidation with amines on heating."],
    ["base", "base_strong", "base_weak", "N_nuc", "O_nuc", "metal", "alkylating"], ["acid_weak"]),
  carboxylate: P("Carboxylate Salt", "Carboxylic Acid", "-C(=O)O-", "Carboxylate anion / ion pair",
    [HI, "Protonated to the free acid by stronger acids (may precipitate)."],
    [RS, "Stable as the anion."], [RS, "Hydrolytically inert."], [LO, "Photostable unless conjugated."],
    [LO, "Stable to moderate heat; decarboxylates only at high temperature."], [LO, "Resistant to auto-oxidation."],
    [MO, "Multivalent metal ions, acids", "Forms insoluble metal carboxylates (soaps) with Ca2+/Mg2+/Zn2+/Al3+; protonated by acids."],
    ["acid_weak", "acid_strong", "metal", "alkylating"], ["base"]),
  ester: P("Carboxylic Ester", "Carbonyl", "-C(=O)O-", "Ester carbonyl carbon and acyl-oxygen bond",
    [HI, "Acid-catalysed acyl-oxygen cleavage (A_AC2) via protonated carbonyl to acid + alcohol; reversible."],
    [CR, "Saponification (B_AC2): hydroxide attacks the carbonyl; irreversible via carboxylate formation."],
    [LO, "Slow neutral hydrolysis, accelerated by humidity, heat and autocatalysis by liberated acid."],
    [LO, "Weak n->pi* absorber; simple alkyl esters are photostable."],
    [MO, "Thermal transesterification; beta-H syn-elimination (ester pyrolysis) only at high temperature."],
    [LO, "Alkoxy alpha C-H oxidation only under radical conditions."],
    [HI, "Amines, alcohols, water", "Aminolysis to amides and transesterification with nucleophilic co-components."],
    ["N_nuc", "O_nuc", "base", "base_strong", "water", "acid_strong"], ["ester"]),
  ester_aryl: P("Aryl / Vinyl Ester (Activated)", "Carbonyl", "Ar-O-C(=O)-", "Ester carbonyl carbon (good phenoxide/enolate leaving group)",
    [HI, "Acid-catalysed hydrolysis to acid + phenol/enol."],
    [CR, "Rapid saponification because the phenoxide/enolate is an excellent leaving group."],
    [HI, "Activated ester: measurable neutral hydrolysis and buffer catalysis under humidity."],
    [HI, "Photo-Fries rearrangement (aryl esters) or acyl-O homolysis on UV excitation."],
    [MO, "Thermal acyl transfer to nucleophiles; Fries-type rearrangement with Lewis acids."],
    [LO, "Resistant to ambient oxidation."],
    [CR, "Amines, alcohols, bases", "Facile acyl transfer (aminolysis / transesterification) and base-catalysed cleavage."],
    ["N_nuc", "O_nuc", "S_nuc", "base", "base_strong", "water"], ["ester", "acylating"]),
  lactone: P("Cyclic Lactone", "Heterocycle / Ester", "C1OC(=O)CC1", "Cyclic ester carbonyl",
    [MO, "Acid-catalysed ring opening to hydroxy-acid; equilibrium often favours the closed 5-/6-ring lactone."],
    [CR, "Alkaline ring opening to the hydroxy-carboxylate salt."],
    [MO, "Ring-chain equilibrium in water; strained (4-ring) lactones hydrolyse rapidly."],
    [MO, "Photodecarboxylation / photorearrangement possible for conjugated lactones."],
    [MO, "Thermal ring opening/polymerisation (ROP) and dehydration."],
    [LO, "Resistant to direct oxidation."],
    [HI, "Amines, alcohols, bases", "Aminolysis/alcoholysis opens the ring to hydroxy-amide / hydroxy-ester."],
    ["N_nuc", "O_nuc", "base", "base_strong", "water"], ["ester", "acylating"]),
  // ---------------- amides and relatives ----------------
  amide: P("Amide Bond", "Carbonyl / Nitrogen", "-C(=O)N<", "Amide carbonyl carbon and C-N bond",
    [MO, "Acid-catalysed hydrolysis (A_AC2) to carboxylic acid + ammonium; needs strong acid / heat."],
    [MO, "Base-promoted acyl substitution; slowed by amide resonance."],
    [LO, "Very slow neutral hydrolysis (half-lives of years); faster at extreme pH."],
    [MO, "UV-induced C-N cleavage; anilides can undergo photo-Fries rearrangement."],
    [LO, "Thermally robust; deamidation or cyclisation only at high temperature."],
    [LO, "Radical alpha-C-H abstraction at N-alkyl groups (N-dealkylation) under harsh conditions."],
    [LO, "Polar co-components", "Hydrogen-bond donor/acceptor interactions; no covalent reactivity under mild conditions."],
    ["acid_strong", "base_strong", "water"], []),
  lactam: P("Lactam (Cyclic Amide)", "Heterocycle / Amide", "C1NC(=O)CC1", "Cyclic amide carbonyl",
    [MO, "Ring opening to amino acid under strong acid / heat (6- and 7-ring faster than 5-ring)."],
    [MO, "Base-catalysed ring opening to amino-carboxylate; slow for 5-ring."],
    [LO, "Stable under ambient humidity."], [LO, "Photostable unless conjugated."],
    [MO, "Ring-opening polymerisation of medium-ring lactams at elevated temperature."],
    [LO, "Resistant; N-alkyl alpha C-H oxidation under radical stress."],
    [LO, "Strong nucleophiles", "Ring opening only with strong nucleophiles or catalysis."], ["base_strong", "acid_strong"], []),
  beta_lactam: P("Beta-Lactam Core", "Heterocycle / Strained Amide", "N1C(=O)CC1", "Strained four-membered lactam carbonyl",
    [CR, "Protonation then rapid nucleophilic ring opening (ring strain ~26 kcal/mol)."],
    [CR, "Direct hydroxide attack on the strained carbonyl; irreversible ring cleavage."],
    [CR, "Spontaneous neutral solvolysis driven by ring strain."],
    [MO, "Photolytic fragmentation of the four-membered ring."],
    [HI, "Thermally accelerated ring rupture."],
    [MO, "Oxidation at adjacent heteroatoms / activated C-H."],
    [CR, "Amines, alcohols, thiols", "Rapid aminolysis / alcoholysis opening the strained lactam."],
    ["N_nuc", "O_nuc", "S_nuc", "base", "base_strong", "water", "acid_strong"], ["acylating"]),
  imide: P("Imide", "Carbonyl / Nitrogen", "-C(=O)NC(=O)-", "Imide carbonyls (electrophilic, acidic N-H)",
    [MO, "Acid-catalysed ring / C-N cleavage to amic acid, then diacid."],
    [CR, "Hydroxide opens the imide to the amic acid within minutes-hours; acidic N-H forms salts."],
    [MO, "Neutral hydrolysis to amic acid on prolonged exposure to moisture."],
    [MO, "Norrish-type chemistry for aryl imides."], [LO, "Thermally stable."], [LO, "Resistant."],
    [HI, "Amines, bases", "Aminolysis opens the imide ring; N-H deprotonation by bases."], ["N_nuc", "base", "base_strong", "water"], ["acylating"]),
  urea: P("Urea", "Carbonyl / Nitrogen", "-NC(=O)N-", "Urea carbonyl carbon",
    [MO, "Slow acid hydrolysis to amine + CO2 (via carbamic acid)."], [LO, "Slow base hydrolysis; N-H urea can eliminate to isocyanate at high pH/heat."],
    [LO, "Hydrolytically stable at ambient conditions."], [LO, "Photostable unless aryl-conjugated."],
    [MO, "Thermal dissociation to isocyanate + amine above ~130 C."], [LO, "Resistant."],
    [LO, "Aldehydes", "Condensation with formaldehyde / aldehydes (methylol ureas)."], ["carbonyl"], []),
  carbamate: P("Carbamate", "Carbonyl / Nitrogen", "-NC(=O)O-", "Carbamate carbonyl carbon",
    [MO, "Acid hydrolysis to amine + CO2 + alcohol; tert-alkyl carbamates cleave rapidly (A_AL1)."],
    [MO, "Base hydrolysis, faster for aryl carbamates via E1cB isocyanate."], [LO, "Stable at neutral pH."],
    [MO, "Photocleavage of aryl carbamates."], [MO, "Thermolysis to isocyanate + alcohol (>150 C); tert-butyl types lose isobutene."],
    [LO, "Resistant."], [MO, "Amines, alcohols", "Transcarbamoylation with strong nucleophiles."], ["N_nuc", "O_nuc", "base_strong"], []),
  carbonate: P("Carbonate Ester", "Carbonyl", "-OC(=O)O-", "Carbonate carbonyl carbon",
    [MO, "Acid hydrolysis to alcohols + CO2."], [CR, "Rapid base hydrolysis to alcohols + carbonate."], [MO, "Slow neutral hydrolysis."],
    [LO, "Photostable (aliphatic)."], [MO, "Thermal decarboxylation / transesterification."], [LO, "Resistant."],
    [HI, "Amines, alcohols", "Aminolysis to carbamates; transesterification."], ["N_nuc", "O_nuc", "base", "base_strong", "water"], ["acylating"]),
  anhydride: P("Acid Anhydride", "Carbonyl", "-C(=O)OC(=O)-", "Anhydride carbonyl carbons",
    [CR, "Acid-catalysed hydrolysis to two carboxylic acids."], [CR, "Instant saponification."], [CR, "Spontaneous hydrolysis by ambient moisture."],
    [LO, "Photostable (aliphatic)."], [MO, "Thermal disproportionation / decarboxylative fragmentation."], [LO, "Resistant."],
    [CR, "Amines, alcohols, water", "Vigorous acylation of amines (amides) and alcohols (esters)."], ["N_nuc", "O_nuc", "S_nuc", "water", "base"], ["acylating"]),
  acyl_halide: P("Acyl Halide", "Carbonyl", "-C(=O)X", "Acyl carbon (strong electrophile)",
    [CR, "Immediate hydrolysis to the acid + HX."], [CR, "Instant saponification."], [CR, "Reacts violently with atmospheric moisture."],
    [MO, "Photolytic C-X cleavage."], [HI, "Thermally labile; ketene formation with alpha-H."], [LO, "Resistant."],
    [CR, "Amines, alcohols, water", "Fast acylation of any nucleophile."], ["N_nuc", "O_nuc", "S_nuc", "water", "base"], ["acylating", "acid_strong"]),
  thioester: P("Thioester", "Carbonyl / Sulfur", "-C(=O)S-", "Thioester carbonyl carbon and sulfur",
    [MO, "Acid hydrolysis to acid + thiol (slower than oxoesters in acid)."], [HI, "Base hydrolysis / thiolate exchange."], [MO, "Slow neutral hydrolysis."],
    [MO, "Photo-cleavage of the C-S bond."], [LO, "Thermally stable."], [HI, "Sulfur oxidation to sulfoxide-type species."],
    [HI, "Amines, thiols", "Aminolysis and thiol-thioester exchange."], ["N_nuc", "S_nuc", "base"], ["acylating"]),
  // ---------------- carbonyls ----------------
  aldehyde: P("Aldehyde", "Carbonyl", "-CHO", "Formyl carbon (oxidation / condensation centre)",
    [HI, "Acid-catalysed hydration and acetal formation with alcohols."],
    [HI, "Aldol self-condensation (alpha-H present) or Cannizzaro disproportionation (no alpha-H)."],
    [MO, "Reversible hydration to the gem-diol in water."], [HI, "Norrish type I/II cleavage and photo-decarbonylation."],
    [MO, "Thermal oligomerisation / decarbonylation."], [CR, "Radical auto-oxidation to the carboxylic acid via peracid."],
    [CR, "Amines, alcohols", "Imine (Schiff base) formation with amines; hemiacetal/acetal with alcohols."],
    ["N_nuc", "N_nuc_prim", "O_nuc", "base", "base_strong", "oxidant"], ["carbonyl"]),
  ketone: P("Ketone", "Carbonyl", "-C(=O)-", "Carbonyl carbon and enolisable alpha-carbons",
    [MO, "Acid-catalysed enolisation, aldol condensation and ketal formation."], [MO, "Enolate formation; aldol condensation; alpha-epimerisation."],
    [RS, "Resistant to neutral hydrolysis."], [HI, "n->pi* absorption (280-320 nm) initiates Norrish type I/II photocleavage."],
    [LO, "Thermally stable."], [MO, "Baeyer-Villiger oxidation by peroxides/peracids to esters."],
    [MO, "Primary amines, alcohols", "Ketimine/enamine formation; ketalisation."], ["N_nuc_prim", "N_nuc", "O_nuc", "oxidant"], ["carbonyl"]),
  michael: P("Alpha,Beta-Unsaturated Carbonyl (Michael Acceptor)", "Carbonyl / Unsaturation", "C=C-C=O", "Electrophilic beta-carbon and carbonyl",
    [MO, "Acid-catalysed hydration / conjugate addition of water."], [HI, "Hydroxide conjugate addition; base-catalysed polymerisation."],
    [LO, "Slow conjugate hydration in water."], [HI, "[2+2] photocycloaddition and E/Z photo-isomerisation."],
    [MO, "Thermal radical polymerisation / Diels-Alder dimerisation."], [HI, "Nucleophilic epoxidation (Weitz-Scheffer) with peroxides."],
    [CR, "Amines, thiols", "Aza-/thia-Michael addition to the beta-carbon."], ["N_nuc", "S_nuc", "base", "oxidant"], ["michael"]),
  quinone: P("Quinone / Quinone Imine", "Carbonyl / Redox-active", "O=C1C=CC(=O)C=C1", "Electrophilic ring carbons and redox-active carbonyls",
    [LO, "Stable in mild acid; acid-catalysed conjugate addition of water."], [HI, "Hydroxide addition / base-promoted decomposition and polymerisation."],
    [LO, "Slow conjugate hydration."], [HI, "Photoreduction and photo-addition to alkenes."], [MO, "Thermal dimerisation / polymerisation."],
    [LO, "Already in the oxidised state."], [CR, "Thiols, amines, reductants", "Michael addition of thiols/amines; reduced to hydroquinone by reductants."],
    ["S_nuc", "N_nuc", "reductant"], ["oxidant", "michael"]),
  thiocarbonyl: P("Thiocarbonyl (C=S)", "Sulfur / Carbonyl", "C=S", "Thiocarbonyl carbon and sulfur",
    [MO, "Acid hydrolysis to the carbonyl analogue with H2S loss."], [MO, "Base hydrolysis to the carbonyl analogue."], [LO, "Slow hydrolysis."],
    [HI, "Photo-oxidation / desulfurisation."], [MO, "Thermal rearrangement (thione-thiol)."], [CR, "S-oxidation and oxidative desulfurisation to C=O."],
    [MO, "Oxidants, metals", "Oxidative desulfurisation; strong metal coordination via sulfur."], ["oxidant", "metal"], ["reductant"]),
  // ---------------- nitrogen ----------------
  amine: P("Aliphatic Amine", "Amine", "-NH2 / -NHR", "Basic, nucleophilic nitrogen lone pair",
    [HI, "Protonation to the ammonium salt (reversible); no covalent breakdown."], [LO, "Stays as the free base (nucleophilic form)."], [RS, "Hydrolytically inert."],
    [MO, "Photosensitised alpha-C-H abstraction / deamination."], [MO, "Thermal deamination, condensation or oxidative discolouration."],
    [HI, "Air / peroxide oxidation to hydroxylamine, nitrone/imine or N-dealkylation."],
    [CR, "Carbonyls, reducing sugars, esters, acylating agents, nitrite", "Schiff-base / Maillard chemistry with carbonyls; acylation by esters/anhydrides; N-alkylation; N-nitrosation of secondary amines."],
    ["carbonyl", "reducing_sugar", "acylating", "ester", "alkylating", "michael", "nitrosating", "acid_weak", "acid_strong", "oxidant"], ["N_nuc", "N_nuc_prim", "base"]),
  amine_tert: P("Tertiary Amine", "Amine", "-NR3", "Tertiary nitrogen lone pair and N-alkyl alpha C-H",
    [HI, "Protonation to the tertiary ammonium salt."], [RS, "Not ionisable by bases."], [RS, "Hydrolytically stable."],
    [LO, "Weak direct UV absorption (unless aryl-conjugated)."], [MO, "Hofmann-type / Cope elimination pathways at elevated temperature."],
    [CR, "N-oxidation to the N-oxide by peroxides / O2; oxidative N-dealkylation."],
    [MO, "Peroxides, alkyl halides, nitrite", "Peroxide N-oxidation; quaternisation by alkylating agents; nitrosative dealkylation to nitrosamines."],
    ["oxidant", "peroxide_former", "alkylating", "nitrosating", "acid_weak", "acid_strong"], ["base"]),
  arylamine: P("Aromatic Amine (Aniline)", "Amine", "Ar-NH2 / Ar-NHR", "Aniline nitrogen and activated ortho/para ring carbons",
    [MO, "Protonation to anilinium (pKa ~4-5); diazotisation with nitrite/HONO."], [LO, "Stable as the free base."], [RS, "Hydrolytically stable."],
    [CR, "Photo-oxidation giving coloured azo / quinone-imine products."], [MO, "Thermal oxidative coupling / dimerisation."],
    [CR, "Multi-electron oxidation to hydroxylamine, nitroso, nitro, azo and quinone-imine polymers."],
    [HI, "Aldehydes, reducing sugars, nitrite, acylating agents", "Schiff-base condensation with carbonyls; diazotisation / N-nitrosation with nitrite; acylation."],
    ["carbonyl", "reducing_sugar", "acylating", "nitrosating", "oxidant", "alkylating"], ["N_nuc", "N_nuc_prim", "base"]),
  quat: P("Quaternary Ammonium", "Amine", "-N+R4", "Quaternary nitrogen and beta C-H (Hofmann)",
    [RS, "Stable in acid."], [MO, "Hofmann elimination / hydroxide displacement at elevated temperature."], [RS, "Hydrolytically stable."],
    [LO, "Photostable."], [MO, "Hofmann elimination and dealkylation on heating."], [LO, "Resistant to oxidation."],
    [LO, "Anionic species", "Ion-pairing with anionic compounds (surfactant-type complexation)."], ["acid_weak"], []),
  ammonium: P("Ammonium / Protonated Amine", "Amine salt", "-NH3+", "Acidic N-H+",
    [RS, "Stable in acid."], [CR, "Deprotonated to the free (nucleophilic) amine by bases."], [RS, "Hydrolytically inert."],
    [LO, "Photostable."], [MO, "Loss of HX / amine on heating (dissociation of the salt)."], [LO, "Resistant until free-based."],
    [MO, "Bases", "Liberation of the free amine, which then reacts as an amine."], ["base", "base_strong"], ["acid_weak"]),
  amidine: P("Amidine / Guanidine", "Nitrogen", "-C(=N)N<", "Amidinium carbon and basic nitrogen",
    [HI, "Protonation to a highly stabilised amidinium/guanidinium; slow hydrolysis to amide/urea."], [MO, "Base hydrolysis to amide / urea."], [MO, "Slow hydrolysis to amide / urea."],
    [LO, "Photostable."], [LO, "Thermally stable."], [LO, "Resistant."], [MO, "Acids, carbonyls", "Strong base; salt formation; condensation with 1,3-dicarbonyls."], ["acid_weak", "acid_strong", "carbonyl"], ["base"]),
  imine: P("Imine (Schiff Base)", "Nitrogen / C=N", "C=N", "Electrophilic imine carbon",
    [CR, "Hydrolysis to carbonyl + amine (protonated iminium)."], [MO, "Base-catalysed hydrolysis / tautomerisation to enamine."], [HI, "Reversible hydrolysis in water."],
    [MO, "E/Z photo-isomerisation and photo-cleavage."], [MO, "Thermal tautomerisation / oligomerisation."], [MO, "Oxidation to amide or oxaziridine."],
    [HI, "Amines, carbonyls, reductants", "Transimination with amines; reduction to amines."], ["N_nuc_prim", "N_nuc", "reductant", "water"], []),
  oxime: P("Oxime / Nitrone", "Nitrogen / C=N-O", "C=N-OH", "Oxime carbon, N-O bond",
    [MO, "Hydrolysis to carbonyl + hydroxylamine; Beckmann rearrangement with strong acid."], [LO, "Stable; oximate salt formation."], [LO, "Slow hydrolysis."],
    [MO, "Photo-isomerisation (E/Z) and N-O cleavage."], [MO, "Thermal Beckmann / dehydration to nitrile."], [MO, "Oxidative cleavage to nitro / carbonyl."],
    [LO, "Acids", "Acid-promoted rearrangement."], ["acid_strong"], []),
  hydrazone: P("Hydrazone", "Nitrogen / C=N-N", "C=N-N", "Hydrazone carbon, N-N unit",
    [HI, "Hydrolysis to carbonyl + hydrazine."], [LO, "Stable; Wolff-Kishner-type reduction requires strong base/heat."], [MO, "Slow hydrolysis."],
    [MO, "Photo-isomerisation / N-N cleavage."], [MO, "Thermal decomposition."], [HI, "Oxidation to azo/diazo compounds."],
    [MO, "Carbonyls", "Exchange with carbonyl compounds."], ["carbonyl"], []),
  hydrazine: P("Hydrazine (N-N)", "Nitrogen", "N-N", "Nucleophilic N-N nitrogens",
    [MO, "Protonation to hydrazinium salts."], [LO, "Stable."], [RS, "Hydrolytically stable."], [MO, "Photo-oxidation / N-N cleavage."],
    [MO, "Thermal decomposition to N2/NH3."], [CR, "Oxidation to diazene / azo compounds with N2 loss."],
    [CR, "Carbonyls, oxidants", "Condensation with carbonyls to hydrazones; reduces oxidants and metal ions."], ["carbonyl", "reducing_sugar", "oxidant", "acylating"], ["N_nuc", "N_nuc_prim", "reductant"]),
  hydroxylamine: P("Hydroxylamine (N-OH)", "Nitrogen / Oxygen", "N-OH", "Nucleophilic N and N-O bond",
    [MO, "Protonation to hydroxylammonium salts."], [MO, "Free-base is unstable; disproportionation."], [RS, "Hydrolytically stable."],
    [MO, "Photolytic N-O cleavage."], [HI, "Thermal disproportionation (can be energetic)."], [CR, "Oxidation to nitroso / nitrone / nitroxide."],
    [CR, "Carbonyls, oxidants", "Oxime formation with carbonyls; reducing agent for oxidants / metal ions."], ["carbonyl", "reducing_sugar", "oxidant", "acylating"], ["N_nuc", "reductant"]),
  nitrile: P("Nitrile", "Nitrogen", "-C#N", "Electrophilic nitrile carbon",
    [MO, "Hydration to amide then carboxylic acid (needs strong acid / heat)."], [MO, "Base hydrolysis to amide, then carboxylate + NH3."], [RS, "Stable at neutral pH."],
    [LO, "Photostable (aliphatic)."], [LO, "Thermally stable."], [LO, "Resistant."], [LO, "Alcohols + acid", "Pinner reaction with alcohols under HCl."], ["O_nuc", "acid_strong"], []),
  nitro: P("Nitro Group", "Nitrogen / Oxygen", "-NO2", "Nitro nitrogen (strong electron-withdrawing group)",
    [RS, "Stable in acid."], [LO, "Stable; aliphatic nitro compounds form nitronate salts and aromatic nitro activates SNAr / Meisenheimer chemistry."], [RS, "Hydrolytically stable."],
    [HI, "Photoreduction, nitro->nitrite rearrangement and ortho-nitrobenzyl photochemistry."], [MO, "Exothermic thermal decomposition at high temperature."],
    [RS, "Already fully oxidised."], [MO, "Reductants, nucleophiles", "Reduced to nitroso / hydroxylamine / amine; activates SNAr displacement of ortho/para leaving groups."], ["reductant", "N_nuc", "base_strong"], ["oxidant"]),
  nitrosamine: P("N-Nitroso (Nitrosamine)", "Nitrogen / Nitroso", "N-N=O", "N-nitroso group (N-N=O)",
    [HI, "Acid-mediated denitrosation / transnitrosation."], [LO, "Stable in base (alpha-deprotonation only with strong base)."], [LO, "Hydrolytically stable."],
    [CR, "UV (230-350 nm) cleaves the N-N bond giving aminyl radical + NO."], [MO, "Thermally stable at moderate temperature."],
    [MO, "Alpha-hydroxylation (P450-like / radical) leading to alkyl-diazonium ions."], [MO, "Thiols, acids", "Transnitrosation to thiols/amines; denitrosation with strong acid."], ["S_nuc", "N_nuc", "acid_strong"], ["nitrosating"]),
  nitroso: P("C-Nitroso", "Nitrogen / Nitroso", "C-N=O", "Nitroso nitrogen and oxygen",
    [MO, "Acid-catalysed rearrangement / condensation."], [MO, "Base-promoted condensation."], [LO, "Hydrolytically stable."],
    [HI, "Photo-dimerisation / photolysis to radicals."], [HI, "Thermal dimerisation-monomerisation equilibrium and decomposition."], [HI, "Oxidation to nitro."],
    [CR, "Amines, thiols", "Condenses with amines (azo) and thiols; ene reactions."], ["N_nuc", "S_nuc", "reductant"], []),
  azo: P("Azo (N=N)", "Nitrogen", "-N=N-", "Azo nitrogens (chromophore)",
    [LO, "Stable (azo-hydrazone tautomerism for hydroxy-azo dyes)."], [LO, "Stable."], [RS, "Hydrolytically stable."],
    [CR, "E/Z photo-isomerisation and photo-cleavage with N2 loss."], [MO, "Thermal cis-trans relaxation; aliphatic azo compounds lose N2 (radical source)."],
    [MO, "Oxidative cleavage by strong oxidants."], [MO, "Reductants", "Reduced to hydrazo compounds and then amines."], ["reductant"], []),
  azide: P("Azide", "Nitrogen", "-N3", "Azide nitrogens (energetic)",
    [MO, "Protonation gives volatile toxic hydrazoic acid (HN3)."], [LO, "Stable."], [RS, "Hydrolytically stable."],
    [CR, "Photolysis extrudes N2 to give a nitrene."], [HI, "Thermal N2 loss (Curtius-like); explosive hazard."], [LO, "Resistant."],
    [HI, "Phosphines, alkynes, reductants", "Staudinger reduction, azide-alkyne cycloaddition, hydrogenation to amines."], ["reductant", "N_nuc"], []),
  isocyanate: P("Isocyanate", "Nitrogen / Carbonyl", "-N=C=O", "Cumulated N=C=O carbon",
    [CR, "Hydrolysis to amine + CO2."], [CR, "Base-catalysed hydrolysis / trimerisation."], [CR, "Reacts with moisture to give amine + CO2 (then urea)."],
    [MO, "Photolysis to nitrene / CO."], [MO, "Dimerisation / trimerisation on heating."], [LO, "Resistant."],
    [CR, "Amines, alcohols", "Addition to amines gives ureas; to alcohols gives carbamates."], ["N_nuc", "O_nuc", "S_nuc", "water"], ["acylating"]),
  isothiocyanate: P("Isothiocyanate", "Nitrogen / Sulfur", "-N=C=S", "Cumulated N=C=S carbon",
    [MO, "Hydrolysis to amine + COS."], [HI, "Base-promoted hydrolysis."], [MO, "Slow hydrolysis in water."],
    [MO, "Photolysis."], [MO, "Thermal rearrangement."], [MO, "Oxidative desulfurisation."],
    [CR, "Amines, thiols", "Addition of amines gives thioureas; thiols give dithiocarbamates."], ["N_nuc", "S_nuc"], []),
  n_oxide: P("N-Oxide", "Nitrogen / Oxygen", "N+-O-", "N-oxide oxygen",
    [MO, "Protonation to N-hydroxy ammonium."], [LO, "Stable."], [RS, "Hydrolytically stable."],
    [HI, "Photo-deoxygenation / rearrangement."], [HI, "Cope / Meisenheimer / Polonovski rearrangements on heating."], [LO, "Already oxidised."],
    [MO, "Reductants, acylating agents", "Deoxygenation by reductants; Polonovski activation by anhydrides."], ["reductant", "acylating"], ["oxidant"]),
  // ---------------- oxygen ----------------
  alcohol: P("Aliphatic Alcohol", "Hydroxyl", "-CH(OH)-", "Carbinol carbon (dehydration / oxidation centre)",
    [MO, "Acid-promoted dehydration, etherification or esterification."], [RS, "Resistant to base (forms alkoxide only with strong base)."], [RS, "Hydrolytically stable."],
    [LO, "Transparent in the solar UV."], [MO, "Thermal dehydration / dehydrogenation."], [HI, "Oxidation to aldehyde/ketone (primary/secondary carbinols)."],
    [MO, "Acids, esters, carbonyls", "Esterification / transesterification; hemiacetal and acetal formation."], ["acid_weak", "acid_strong", "ester", "acylating", "carbonyl", "alkylating"], ["O_nuc"]),
  alcohol_tert: P("Tertiary Alcohol", "Hydroxyl", "-C(OH)<", "Tertiary carbinol carbon",
    [HI, "E1 dehydration to alkene via tertiary carbocation."], [RS, "Resistant."], [RS, "Hydrolytically stable."],
    [LO, "Transparent in the solar UV."], [MO, "Thermal / acid-catalysed dehydration."], [RS, "Resists oxidation without C-C cleavage."],
    [LO, "Acylating agents", "Sterically hindered esterification."], ["acylating", "acid_strong"], ["O_nuc"]),
  phenol: P("Phenol (Ar-OH)", "Hydroxyl", "Ar-OH", "Phenolic oxygen and activated ortho/para ring carbons",
    [RS, "Resistant to acid solvolysis."], [HI, "Deprotonation to phenolate, which oxidises much faster."], [RS, "Hydrolytically stable."],
    [HI, "UV excitation gives phenoxyl radicals; photo-coupling to dimers."], [MO, "Thermally accelerated oxidative coupling."],
    [CR, "Single-electron oxidation to phenoxyl radical, then C-C/C-O coupling or quinone formation."],
    [MO, "Bases, oxidants, carbonyls, metal ions", "Salt formation with bases; peroxide oxidation; Fe/Cu-catalysed oxidation; Mannich/aldehyde condensation."], ["base", "base_strong", "oxidant", "redox_metal", "carbonyl", "acylating"], ["O_nuc", "acid_weak", "reductant"]),
  enol: P("Enol", "Hydroxyl / Unsaturation", "C=C-OH", "Enolic carbon-oxygen system",
    [MO, "Tautomerises to the carbonyl form."], [HI, "Enolate formation."], [LO, "Tautomerisation in water."],
    [MO, "Photo-tautomerisation."], [MO, "Thermal tautomerisation / oxidation."], [HI, "Air oxidation (e.g. ene-diols such as ascorbate)."],
    [MO, "Oxidants, metals", "Oxidised readily; chelates metal ions."], ["oxidant", "redox_metal", "metal"], ["reductant"]),
  ether: P("Ether", "Ether", "C-O-C", "Ether oxygen and alpha C-H bonds",
    [LO, "Cleavage only with HBr/HI or strong Lewis acids."], [RS, "Stable to base."], [RS, "Hydrolytically inert."],
    [LO, "Weakly UV active."], [LO, "Thermally stable."], [HI, "Autoxidation at alpha C-H to hydroperoxides (peroxide formation)."],
    [MO, "Oxidisable co-components", "Peroxides formed on ageing oxidise co-existing S, N and phenolic groups."], ["oxidant"], ["peroxide_former"]),
  aryl_ether: P("Aryl Ether", "Ether", "Ar-O-C", "Ether oxygen and O-alkyl alpha C-H",
    [LO, "Stable except with HBr/HI or BBr3."], [RS, "Stable to base."], [RS, "Hydrolytically inert."],
    [MO, "UV-induced C-O homolysis for some aryl ethers."], [LO, "Thermally stable."], [MO, "O-dealkylation via alpha C-H abstraction (radical / enzymatic-like)."],
    [LO, "Strong Lewis acids", "Demethylation with BBr3/HBr."], ["acid_strong", "lewis_acid"], []),
  epoxide: P("Epoxide (Oxirane)", "Ether / Strained ring", "C1OC1", "Strained C-O bonds (ring-opening centre)",
    [CR, "Acid-catalysed ring opening to the 1,2-diol (or halohydrin with HX)."], [HI, "SN2 ring opening by hydroxide to the diol."], [MO, "Slow neutral hydrolysis to diol."],
    [LO, "Photostable (aliphatic)."], [MO, "Thermal rearrangement to carbonyl (Meinwald) / polymerisation."], [LO, "Resistant."],
    [CR, "Amines, alcohols, acids, thiols", "Ring opening by nucleophiles giving beta-substituted alcohols."], ["N_nuc", "O_nuc", "S_nuc", "acid_weak", "acid_strong", "base_strong", "water"], ["alkylating"]),
  acetal: P("Acetal / Ketal", "Ether / Carbonyl derivative", "C(OR)2", "Acetal carbon (oxocarbenium precursor)",
    [CR, "Rapid acid hydrolysis via oxocarbenium to carbonyl + alcohols (glycosides likewise)."], [RS, "Stable to base."], [LO, "Stable at neutral pH."],
    [LO, "Photostable."], [LO, "Thermally stable in the absence of acid."], [MO, "Oxidation of the acetal C-H to esters."],
    [MO, "Alcohols + acid", "Transacetalisation with alcohols under acid catalysis."], ["acid_weak", "acid_strong", "O_nuc", "water"], []),
  hemiacetal: P("Hemiacetal / Reducing End (Masked Aldehyde)", "Carbonyl derivative", "C(OH)(OR)", "Anomeric carbon (ring-chain tautomerism)",
    [MO, "Mutarotation and glycoside formation under acid catalysis."], [HI, "Ring opening, enolisation, epimerisation and alkaline degradation (Lobry de Bruyn)."], [MO, "Ring-chain tautomerism (mutarotation) in water."],
    [LO, "Photostable."], [HI, "Dehydration / caramelisation on heating."], [HI, "Oxidation of the aldehyde form to aldonic acid."],
    [CR, "Amines, amino acids", "Reducing end condenses with amines (glycosylamine / Schiff base, Maillard browning)."], ["N_nuc", "N_nuc_prim", "oxidant", "base"], ["reducing_sugar", "carbonyl", "reductant"]),
  peroxide: P("Peroxide / Hydroperoxide", "Peroxide", "-O-O-", "Weak O-O bond (oxidant / radical source)",
    [MO, "Acid-catalysed rearrangement (Hock/Criegee) or heterolysis."], [MO, "Base-promoted decomposition / perhydrolysis."], [LO, "Stable in neutral water (aqueous H2O2 disproportionates slowly)."],
    [CR, "UV homolyses O-O to alkoxy / hydroxyl radicals."], [CR, "Thermal homolysis of O-O initiates radical chains."], [RS, "Already an oxidant."],
    [CR, "Sulfides, amines, phenols, alkenes, metals", "Oxidises S, N, phenol and C=C groups; transition metals catalyse radical (Fenton) decomposition."], ["reductant", "N_nuc", "S_nuc", "redox_metal", "ene"], ["oxidant"]),
  dioxygen: P("Molecular Oxygen / Ozone", "Oxidant", "O=O", "Diradical / electrophilic oxygen",
    [RS, "Inert to acid."], [RS, "Inert to base."], [RS, "Does not hydrolyse."], [MO, "Photosensitised singlet oxygen formation."], [LO, "Stable."], [RS, "Oxidant, not oxidised."],
    [HI, "Autoxidisable groups", "Radical autoxidation of aldehydes, thiols, ethers, enols, phenols, amines, benzylic/allylic C-H."], ["ene", "S_nuc", "N_nuc", "reductant", "peroxide_former"], ["oxidant"]),
  strong_base: P("Hydroxide / Alkoxide / Metal Oxide (Strong Base)", "Base", "OH- / RO- / O2-", "Basic oxygen anion",
    [CR, "Neutralised by acids."], [RS, "Already basic."], [RS, "Does not hydrolyse."], [LO, "Photostable."], [LO, "Thermally stable."], [RS, "Not oxidisable."],
    [CR, "Acids, esters, amides, halides", "Raises microenvironmental pH: catalyses saponification, aldol/Cannizzaro chemistry, E2/SN2; neutralises acids."],
    ["acid_weak", "acid_strong", "ester", "acylating", "alkylating", "carbonyl"], ["base_strong", "base"]),
  water: P("Water", "Solvent / Reagent", "H2O", "Nucleophile / proton donor-acceptor",
    [RS, "Amphoteric spectator."], [RS, "Amphoteric spectator."], [RS, "Reagent for hydrolysis of susceptible groups."], [LO, "Transparent in the solar UV."], [RS, "Stable."], [RS, "Stable."],
    [MO, "Electrophiles", "Hydrolysis reagent and hydrogen-bonding medium; mobilises ions (moisture-mediated solid-state reactions)."], ["acylating", "ester", "alkylating"], ["water"]),
  // ---------------- sulfur / phosphorus ----------------
  thiol: P("Thiol (R-SH)", "Sulfur", "-SH", "Thiol S-H (nucleophile / redox centre)",
    [LO, "Stable in acid."], [HI, "Thiolate formation (pKa 8-10) makes air oxidation very fast."], [RS, "Hydrolytically inert."],
    [MO, "S-H homolysis giving thiyl radicals."], [LO, "Thermally stable."], [CR, "Oxidation to disulfide, then sulfinic / sulfonic acids."],
    [CR, "Oxidants, electrophiles, metals", "Disulfide formation; Michael addition; thiolate-metal binding."], ["oxidant", "michael", "alkylating", "acylating", "redox_metal", "metal"], ["S_nuc", "reductant"]),
  thioether: P("Thioether (Sulfide)", "Sulfur", "-C-S-C-", "Divalent sulfur lone pair (S-oxidation centre)",
    [LO, "Resistant to acid cleavage."], [LO, "Resistant to base."], [RS, "Hydrolytically inert."],
    [MO, "Singlet-oxygen sensitised photo-oxidation; C-S homolysis."], [MO, "Thermal C-S bond homolysis at high temperature."],
    [CR, "Rapid oxidation by air / peroxides to sulfoxide and then sulfone."],
    [HI, "Peroxides, alkylating agents", "Oxidised by peroxides (or peroxide-forming ethers); S-alkylation to sulfonium salts."], ["oxidant", "peroxide_former", "alkylating"], ["S_nuc"]),
  disulfide: P("Disulfide (S-S)", "Sulfur", "-S-S-", "Weak S-S bond",
    [LO, "Stable in acid."], [MO, "Thiolate-disulfide exchange and beta-elimination in base."], [RS, "Hydrolytically stable."],
    [HI, "S-S homolysis by UV to thiyl radicals."], [MO, "Thermal scrambling of S-S bonds."], [MO, "Oxidation to thiosulfinate, then sulfonic acids."],
    [HI, "Thiols, reductants", "Thiol-disulfide exchange; reduction to thiols."], ["S_nuc", "reductant"], []),
  sulfoxide: P("Sulfoxide", "Sulfur", "-S(=O)-", "Sulfinyl sulfur / oxygen",
    [MO, "Pummerer rearrangement under acidic activation."], [LO, "Stable in base."], [RS, "Hydrolytically stable."],
    [MO, "Photo-deoxygenation to sulfide or photo-oxidation."], [HI, "Thermal syn-elimination (sulfoxide pyrolysis) giving alkene + sulfenic acid."],
    [HI, "Further oxidation to the sulfone."], [MO, "Reductants, anhydrides", "Deoxygenation to sulfide; Pummerer activation by anhydrides."], ["reductant", "acylating", "oxidant"], []),
  sulfone: P("Sulfone", "Sulfur", "-SO2-", "Sulfonyl sulfur (very stable)",
    [RS, "Resistant."], [LO, "Alpha-deprotonation only with strong bases."], [RS, "Hydrolytically inert."], [LO, "Photostable."], [LO, "Thermally stable (SO2 extrusion only at very high T)."], [RS, "Fully oxidised."],
    [RS, "None", "Chemically inert under formulation conditions."], [], []),
  sulfonamide: P("Sulfonamide", "Sulfur / Nitrogen", "-SO2NH-", "Sulfonamide N and S-N bond",
    [LO, "Cleaved only by boiling strong acid."], [MO, "N-H deprotonation (pKa ~6-10) gives soluble salts."], [RS, "Hydrolytically stable at ambient conditions."],
    [HI, "Strong UV chromophore; S-N cleavage and SO2 extrusion on photolysis."], [LO, "High thermal stability."], [LO, "Resistant to oxidation."],
    [MO, "Bases, metal cations", "Salt formation with bases; chelation to divalent metals."], ["base", "base_strong", "metal"], ["acid_weak"]),
  sulfonic: P("Sulfonic Acid / Sulfonate", "Sulfur", "-SO3H / -SO3-", "Sulfonate group (strong acid)",
    [RS, "Stable."], [MO, "Fully ionised; forms salts."], [RS, "Hydrolytically stable (aryl sulfonic acids desulfonate only in hot aqueous acid)."],
    [LO, "Photostable (aliphatic)."], [LO, "Thermally stable."], [RS, "Resistant."], [HI, "Bases, amines", "Salt formation with amines/metal ions."], ["base", "N_nuc", "metal", "alkylating"], ["acid_strong"]),
  sulfonate_ester: P("Sulfonate / Sulfate Ester (Alkylating)", "Sulfur / Ester", "-SO2-O-R", "Electrophilic alkyl carbon and sulfonyl sulfur",
    [MO, "Acid hydrolysis to sulfonic acid + alcohol."], [CR, "Hydroxide attacks carbon (SN2) or sulfur, releasing sulfonate."], [MO, "Solvolysis in water / alcohols (alkyl sulfonates are alkylating agents)."],
    [LO, "Photostable."], [MO, "Thermal alkylation / elimination."], [RS, "Resistant."],
    [CR, "Amines, thiols, water", "Alkylates nucleophiles (N-, S-, O-alkylation)."], ["N_nuc", "S_nuc", "O_nuc", "water", "base_strong", "acid_weak"], ["alkylating"]),
  sulfonyl_halide: P("Sulfonyl Halide", "Sulfur", "-SO2X", "Electrophilic sulfonyl sulfur",
    [HI, "Hydrolysis to sulfonic acid + HX."], [CR, "Rapid hydrolysis."], [HI, "Moisture-sensitive; hydrolyses in water."], [MO, "Photolytic S-X cleavage."], [MO, "Thermal SO2 extrusion."], [LO, "Resistant."],
    [CR, "Amines, alcohols", "Sulfonylation of amines (sulfonamides) and alcohols (sulfonate esters)."], ["N_nuc", "O_nuc", "water", "base"], ["acylating", "acid_strong"]),
  sulfite: P("Sulfite / Bisulfite (Reductant)", "Sulfur / Inorganic", "SO3(2-) / HSO3-", "Sulfite sulfur (nucleophile / reductant)",
    [HI, "Releases SO2 in acid."], [LO, "Stable as sulfite."], [RS, "Stable."], [MO, "Photo-oxidation."], [MO, "Loses SO2 on heating (bisulfite/metabisulfite)."],
    [CR, "Rapid oxidation to sulfate by O2 / peroxides."], [HI, "Carbonyls, oxidants, quinones", "Adds to aldehydes/ketones (bisulfite adducts); reduces peroxides, quinones and diazonium compounds."], ["carbonyl", "oxidant", "michael"], ["reductant", "S_nuc"]),
  sulfate: P("Sulfate / Sulfate Ester", "Sulfur / Inorganic", "-OSO3-", "Sulfate oxygen / S-O-C linkage",
    [MO, "Alkyl sulfate esters hydrolyse in acid to alcohol + hydrogen sulfate; inorganic sulfate is inert."], [LO, "Stable."], [LO, "Alkyl sulfates hydrolyse slowly; inorganic sulfate is inert."],
    [LO, "Photostable."], [LO, "Thermally stable."], [RS, "Fully oxidised."], [LO, "Metal ions", "Ion pairing / precipitation with Ca2+, Ba2+, Pb2+."], ["metal"], []),
  phosphate: P("Phosphate / Phosphonate", "Phosphorus", "-P(=O)(O)O-", "Phosphoryl phosphorus and P-O-C ester bonds",
    [MO, "Acid-catalysed P-O-C hydrolysis (esters) - slow."], [MO, "Base-catalysed ester hydrolysis (triesters faster)."], [LO, "Slow neutral hydrolysis."],
    [LO, "Photostable."], [LO, "Thermally stable; condensation to pyrophosphates on dehydration."], [RS, "Fully oxidised."],
    [MO, "Divalent metal ions", "Chelation / precipitation with Ca2+, Mg2+, Al3+, Fe3+."], ["metal", "base_strong"], []),
  phosphine: P("Phosphine", "Phosphorus", "-PR3", "Phosphorus lone pair",
    [MO, "Protonation to phosphonium in strong acid."], [LO, "Stable."], [RS, "Hydrolytically stable."], [MO, "Photo-oxidation."], [LO, "Thermally stable."], [CR, "Rapid oxidation to phosphine oxide."],
    [CR, "Azides, oxidants, alkyl halides", "Staudinger reaction with azides; quaternisation with alkyl halides."], ["oxidant", "alkylating", "reductant"], ["reductant"]),
  // ---------------- halogens ----------------
  alkyl_halide: P("Alkyl Halide", "Halide", "-C-X", "Electrophilic carbon and C-X bond",
    [LO, "Stable; tertiary/benzylic halides solvolyse (SN1)."], [HI, "SN2 substitution or E2 dehydrohalogenation by hydroxide/bases."], [MO, "Solvolysis in water, fast for tertiary/benzylic/allylic and alpha-halo ethers."],
    [MO, "C-X homolysis under UV (C-I > C-Br > C-Cl)."], [MO, "Thermal dehydrohalogenation."], [LO, "Resistant."],
    [HI, "Amines, thiols, carboxylates", "SN2 alkylation of nucleophiles (N-, S-, O-alkylation, quaternisation)."], ["N_nuc", "S_nuc", "O_nuc", "base_strong", "base", "water"], ["alkylating"]),
  aryl_halide: P("Aryl / Vinyl Halide", "Halide", "Ar-X / C=C-X", "sp2 C-X bond",
    [RS, "Resistant."], [LO, "Stable unless activated by ortho/para electron-withdrawing groups (SNAr)."], [RS, "Hydrolytically stable."],
    [MO, "Photodehalogenation (C-I > C-Br > C-Cl >> C-F)."], [LO, "Thermally stable."], [LO, "Resistant."],
    [LO, "Strong nucleophiles, metals", "SNAr with activated rings; metal-catalysed coupling."], ["N_nuc", "S_nuc", "O_nuc", "base_strong", "metal"], []),
  cf3: P("Trifluoromethyl / Perfluoroalkyl", "Halide", "-CF3", "C-F bonds (very strong)",
    [RS, "Resistant."], [LO, "Stable (haloform-type loss only when alpha to carbonyl)."], [RS, "Stable."], [LO, "Photostable."], [RS, "Stable."], [RS, "Resistant."], [RS, "None", "Chemically inert."], [], []),
  halide_ion: P("Halide Anion", "Inorganic", "X-", "Halide ion", [RS, "Spectator."], [RS, "Spectator."], [RS, "Spectator."], [LO, "Photostable (I- can photo-oxidise)."], [RS, "Stable."],
    [MO, "Oxidised to X2 by strong oxidants (I- and Br- readily)."], [MO, "Oxidants, electrophiles", "Nucleophilic catalysis; oxidised by peroxides/metal oxidants."], ["oxidant", "alkylating"], ["halide"]),
  hydrogen_halide: P("Hydrogen Halide (Strong Acid)", "Inorganic acid", "H-X", "Proton (strong acid)", [RS, "Acid itself."], [CR, "Neutralised by bases."], [RS, "Ionises in water."],
    [LO, "Photostable."], [LO, "Stable."], [MO, "Oxidised by strong oxidants to X2."], [CR, "Bases, acid-labile groups", "Protonates amines/bases; catalyses ester/acetal/amide hydrolysis."], ["base", "base_strong", "N_nuc", "ester", "acetal"], ["acid_strong", "halide"]),
  hypohalite: P("Hypohalite / Active Halogen (Oxidant)", "Oxidant", "-O-X / N-X", "Electrophilic halogen (oxidant)", [CR, "Releases X2 / HOX in acid."], [MO, "Disproportionates."], [MO, "Decomposes in water."],
    [CR, "Photolysis to radicals."], [HI, "Thermal decomposition."], [RS, "Oxidant."], [CR, "Amines, sulfides, alkenes", "Oxidises/chlorinates S, N, alkenes and aromatic rings."], ["N_nuc", "S_nuc", "ene", "reductant"], ["oxidant"]),
  // ---------------- carbon frameworks ----------------
  alkene: P("Alkene (Olefin)", "Unsaturation", "-C=C-", "Pi bond and allylic C-H",
    [MO, "Acid-catalysed Markovnikov hydration or carbocation oligomerisation."], [LO, "Resistant unless conjugated with electron-withdrawing groups."], [RS, "Resistant."],
    [HI, "cis-trans photo-isomerisation and [2+2] photodimerisation."], [MO, "Thermal isomerisation / Diels-Alder dimerisation (dienes)."],
    [CR, "Epoxidation by peroxides; allylic autoxidation to hydroperoxides / enones."], [MO, "Peroxides, radicals, electrophiles", "Epoxidation by peroxides; radical addition / crosslinking."], ["oxidant", "peroxide_former", "ene"], ["ene"]),
  alkyne: P("Alkyne", "Unsaturation", "-C#C-", "Pi system (and terminal C-H)",
    [MO, "Acid/metal-catalysed hydration to ketone."], [MO, "Terminal alkynes deprotonate to acetylides."], [RS, "Resistant."], [MO, "Photo-cycloaddition / polymerisation."],
    [MO, "Thermal oligomerisation."], [MO, "Oxidative cleavage to acids / diketones."], [MO, "Cu+/Ag+ ions, azides", "Metal acetylide formation; azide-alkyne cycloaddition."], ["metal", "N_nuc"], []),
  aromatic: P("Aromatic Ring", "Aromatic", "c1ccccc1", "Aromatic pi-system and ring C-H",
    [RS, "Resistant to acid solvolysis."], [RS, "Resistant to base."], [RS, "Hydrolytically inert."], [HI, "Strong UV absorption (254-280 nm) leading to triplet-state chemistry."],
    [RS, "Aromatic stabilisation."], [MO, "Hydroxyl-radical attack forms hydroxylated (phenolic) products; electron-rich rings oxidise faster."],
    [MO, "Electrophiles, radicals", "Electrophilic substitution (nitration, halogenation) with activated rings; pi-stacking / charge-transfer complexes."], ["oxidant", "acid_strong"], []),
  hetero_rich: P("Heteroaromatic Ring (Electron-Rich: Furan / Thiophene / Pyrrole / Azole)", "Aromatic heterocycle", "c1ccoc1", "Electron-rich ring C-H and heteroatom",
    [HI, "Acid-catalysed ring opening / polymerisation (furan, pyrrole)."], [LO, "Stable; N-H azoles deprotonate with strong base."], [LO, "Stable."],
    [HI, "Photo-oxidation with singlet oxygen ([4+2]) and photoisomerisation."], [LO, "Thermally stable."], [HI, "Oxidative ring opening / S-oxidation / polymerisation."],
    [MO, "Oxidants, electrophiles", "Oxidised/halogenated readily; azole N coordinates metal ions."], ["oxidant", "metal", "alkylating"], []),
  hetero_azine: P("Heteroaromatic Ring (Pi-Deficient: Pyridine / Azine / Azole N)", "Aromatic heterocycle", "c1ccncc1", "Ring nitrogen lone pair",
    [MO, "Protonation to the azinium salt (pyridine pKa ~5)."], [RS, "Resistant."], [RS, "Hydrolytically stable."], [MO, "Photo-isomerisation and photo-oxidation."], [LO, "Thermally stable."],
    [MO, "N-oxidation by peroxides/peracids."], [MO, "Alkylating agents, metal ions, acids", "N-alkylation (quaternisation); metal coordination; salt formation."], ["alkylating", "metal", "oxidant", "acid_weak", "acid_strong"], ["base"]),
  benzylic_ch: P("Benzylic / Allylic C-H", "Hydrocarbon", "Ar-CH< / C=C-CH<", "Weak benzylic / allylic C-H bonds (BDE ~85-90 kcal/mol)",
    [RS, "Resistant."], [RS, "Resistant (unless acidified by conjugation)."], [RS, "Resistant."], [MO, "Photo-initiated hydrogen abstraction."], [MO, "Thermal autoxidation initiated by trace radicals."],
    [HI, "Autoxidation to hydroperoxide, then alcohol / ketone / acid."], [MO, "Radical initiators, peroxides, metals", "Radical chain oxidation accelerated by peroxides and redox metals."], ["oxidant", "peroxide_former", "redox_metal"], []),
  // ---------------- inorganic / ionic ----------------
  metal: P("Metal Cation (Main-Group)", "Inorganic ion", "M(n+)", "Lewis-acidic cation",
    [RS, "Spectator."], [MO, "Precipitates as hydroxide/oxide at high pH."], [RS, "Spectator."], [LO, "Photostable."], [RS, "Stable."], [RS, "Not oxidisable (fixed oxidation state)."],
    [HI, "Carboxylates, phenolates, phosphates, carbonyls", "Chelation / salt formation with anionic donors; Lewis-acid activation of carbonyls and esters (Mg2+, Ca2+, Zn2+, Al3+)."], ["acid_weak", "acid_strong", "N_nuc", "carbonyl", "ester"], ["metal", "lewis_acid"]),
  metal_transition: P("Transition-Metal Ion (Redox-Active)", "Inorganic ion", "M(n+)", "Redox-active metal centre",
    [RS, "Spectator."], [MO, "Precipitates as hydroxide / oxide at high pH."], [RS, "Spectator."], [MO, "Ligand-to-metal charge-transfer photochemistry."], [LO, "Stable."], [HI, "Cycles between oxidation states (Fenton-type radical generation)."],
    [CR, "Peroxides, thiols, phenols, chelators", "Catalyses autoxidation and peroxide decomposition (radicals); coordinates N/O/S donors."], ["oxidant", "peroxide_former", "S_nuc", "N_nuc", "acid_weak", "enol"], ["metal", "redox_metal", "lewis_acid"]),
  metal_oxide: P("Inorganic Oxide / Oxo-Metal Species", "Inorganic", "M=O", "Oxide surface / oxo-metal centre",
    [MO, "Basic oxides dissolve in acid; SiO2/TiO2 are inert."], [MO, "Amphoteric / basic surface sites."], [LO, "Insoluble; basic oxides hydrate to hydroxides."],
    [MO, "TiO2 / ZnO are photocatalytic (generate radicals under UV); SiO2 is inert."], [RS, "Stable."], [MO, "High-valent oxo-metals (MnO4-, CrO4 2-) are strong oxidants."],
    [MO, "Acids, adsorbable species", "Surface adsorption; basic oxides neutralise acids; photocatalytic oxidation (TiO2/ZnO)."], ["acid_weak", "acid_strong"], []),
  hydride: P("Hydride Donor (Borohydride)", "Reductant", "BH4-", "Hydridic B-H", [CR, "Evolves H2 in acid."], [LO, "Stable in base."], [HI, "Slowly hydrolyses to H2 + borate in water."],
    [LO, "Photostable."], [MO, "Decomposes at high temperature."], [CR, "Strong reductant (readily oxidised)."], [CR, "Carbonyls, nitro, disulfides", "Reduces aldehydes/ketones, imines, disulfides (and nitro with catalysts)."], ["carbonyl", "acylating"], ["reductant"]),
  carbonate_ion: P("Carbonate / Bicarbonate", "Inorganic base", "CO3(2-) / HCO3-", "Basic carbonate oxygen", [CR, "Neutralised by acid with CO2 evolution."], [RS, "Stable."], [LO, "Buffers water (pH ~8-11)."],
    [LO, "Photostable."], [MO, "Bicarbonate decomposes to carbonate + CO2 + H2O on heating."], [RS, "Not oxidisable."],
    [HI, "Acids, esters", "Neutralises acids (salt formation, CO2); raises pH and catalyses base hydrolysis."], ["acid_weak", "acid_strong", "ester", "acylating"], ["base", "base_strong"]),
  nitrite: P("Nitrite / Nitrosating Agent", "Inorganic", "NO2-", "Nitrite nitrogen / N=O", [CR, "Forms HONO / N2O3 / NO+ (nitrosating species) in acid."], [LO, "Stable as nitrite."], [LO, "Stable."],
    [HI, "UV releases NO / NO2."], [MO, "Disproportionates on heating."], [MO, "Oxidised to nitrate."], [CR, "Secondary/tertiary amines, aromatic amines", "N-nitrosation of amines (nitrosamines); diazotisation of aromatic amines."], ["N_nuc", "N_nuc_sec", "S_nuc", "acid_weak", "acid_strong"], ["nitrosating", "oxidant"]),
  nitrate: P("Nitrate", "Inorganic", "NO3-", "Nitrate oxygen", [RS, "Spectator (HNO3 is a strong oxidising acid)."], [RS, "Spectator."], [RS, "Spectator."], [MO, "Photolysis to nitrite / NO2."], [LO, "Stable."], [RS, "Fully oxidised."],
    [LO, "Reductants", "Weak oxidant in neutral media; strong oxidant as HNO3."], ["reductant"], ["oxidant"]),
  hydrocarbon: P("Aliphatic / Polar Heteroatom Framework", "Other", "C-C / Heteroatom", "Aliphatic C-H and polar heteroatom centres",
    [LO, "No acid-labile functional group."], [LO, "No base-labile functional group."], [LO, "No hydrolysable group."], [LO, "UV absorption depends on chromophores."], [LO, "Thermal bond homolysis at high temperature."],
    [MO, "Radical autoxidation across C-H positions (tertiary > secondary > primary)."], [LO, "Co-reactants", "Non-covalent physical interactions only."], [], []),
};

const VSCORE: Record<Vuln, number> = { Critical: 0.93, High: 0.80, Moderate: 0.58, Low: 0.30, Resistant: 0.10 };
const VORDER: Vuln[] = ["Resistant", "Low", "Moderate", "High", "Critical"];
const stepVuln = (v: Vuln, d: number): Vuln => VORDER[Math.max(0, Math.min(4, VORDER.indexOf(v) + d))];

type CondKey = "acid" | "base" | "hyd" | "photo" | "therm" | "ox";

/** Context-dependent adjustments of a base profile (flags come from the structural site finder). */
function tweakProfile(key: string, f: Record<string, any>): Prof {
  const base = PROFILES[key];
  const p: Prof = { ...base };
  const set = (k: CondKey, v: Vuln, t?: string): void => { p[k] = [v, t ?? p[k][1]] as T; };
  switch (key) {
    case "ester":
      if (f.tert) { set("acid", CR, "Acid-labile tert-alkyl ester: A_AL1/E1 cleavage releases the acid + alkene (e.g. isobutene)."); set("base", LO, "Sterically hindered; resists saponification."); set("therm", HI, "Thermal E1/syn-elimination expels the alkene (isobutene) to give the acid."); }
      if (f.formate) { set("hyd", MO, "Formate esters hydrolyse faster than other alkyl esters."); set("base", CR); }
      if (f.alphaEWG) { set("base", CR, "Alpha electron-withdrawing group accelerates saponification."); set("hyd", MO, "Alpha-EWG activation raises neutral hydrolysis rate."); }
      if (f.hinder) set("base", HI, "Steric hindrance slows saponification relative to unhindered esters.");
      break;
    case "lactone":
      if (f.ring === 4) { set("acid", CR); set("hyd", CR, "Strained beta-lactone hydrolyses rapidly with ring-strain relief (~23 kcal/mol)."); set("therm", HI, "Thermal decarboxylation to alkene / polymerisation."); }
      if (f.alphaEnol) { set("base", LO, "Vinylogous acid (ene-diol lactone): deprotonated to the enolate, which resists hydroxide attack."); set("acid", LO); set("hyd", LO); }
      if (f.arom) { set("base", HI, "Aromatic (coumarin-type) lactone opens in base to the (Z)-hydroxycinnamate."); set("acid", LO); set("photo", HI, "Coumarins undergo [2+2] photodimerisation and photo-oxidation."); }
      break;
    case "amide":
      if (f.arom) { set("acid", LO); set("base", LO); set("hyd", RS, "Aromatic (pyridone-type) amide is hydrolytically inert."); }
      if (f.anilide) { set("photo", HI, "Anilides undergo photo-Fries rearrangement to amino-aryl ketones."); }
      if (f.tertiary) { set("base", LO, "Tertiary amides resist base hydrolysis (no N-H, hindered)."); }
      if (f.formamide) { set("hyd", MO, "Formamides hydrolyse faster than other amides."); set("acid", HI); }
      if (f.hydrazide) { set("acid", HI, "Hydrazides hydrolyse more readily than amides."); }
      break;
    case "carbamate":
      if (f.tert) set("acid", CR, "tert-Butyl-type carbamate (Boc-like): A_AL1 / E1 cleavage loses the alkene, then CO2, giving the free amine.");
      break;
    case "acid":
      if (f.decarb) set("therm", HI, "Activated acid (" + (f.decarbKind || "beta-keto / malonic / alpha-EWG / electron-rich aryl") + "): thermal decarboxylation through a cyclic / zwitterionic transition state releases CO2.");
      break;
    case "arylamine":
      if (f.benz === false) { set("ox", MO, "Heteroaryl amine (amidine-like, delocalised lone pair): much less easily oxidised than an aniline."); set("photo", MO); set("acid", LO, "Weakly basic; protonation occurs on the ring nitrogen."); }
      if (f.tertiary) { set("ox", HI, "Tertiary aryl amine: N-oxidation and oxidative N-dealkylation."); set("photo", HI); }
      break;
    case "amine":
      if (f.ammonia) set("ox", LO, "Ammonia resists mild oxidation.");
      if (f.primary) set("ox", HI, "Oxidises to hydroxylamine / imine / nitroso derivatives.");
      break;
    case "aldehyde":
      if (!f.alphaH) set("base", HI, "No alpha-H: Cannizzaro disproportionation to alcohol + carboxylate (no aldol pathway).");
      else set("base", HI, "Alpha-H present: base-catalysed aldol self-condensation.");
      if (f.aryl) set("ox", CR, "Aromatic aldehydes autoxidise readily to the benzoic acid analogue.");
      break;
    case "ketone":
      if (f.aryl) set("photo", HI, "Aryl ketones populate reactive triplet states (Norrish II if a gamma-H exists; photoreduction).");
      if (!f.alphaH) { set("base", LO, "No enolisable alpha-H."); set("acid", LO); }
      break;
    case "alcohol":
      if (f.benzylic) set("ox", CR, "Benzylic / allylic carbinols oxidise very readily to carbonyls.");
      if (f.methanol) set("ox", HI);
      break;
    case "phenol":
      if (f.hydroquinone) { set("ox", CR, "Hydroquinone/catechol/aminophenol motif: two-electron oxidation to quinone / quinone-imine."); }
      if (f.hindered) { set("ox", HI, "Hindered phenol scavenges radicals (antioxidant) and is converted to phenoxyl / quinone-methide products."); }
      if (f.ewg) { set("ox", MO, "Electron-withdrawing substituents raise the oxidation potential."); set("base", CR, "Acidified phenol (pKa < 8) is fully ionised in mild base."); }
      break;
    case "alkyl_halide":
      if (f.tertiary || f.benzylic || f.allylic) { set("hyd", HI, "Tertiary / benzylic / allylic halides solvolyse via SN1 in water."); set("acid", MO); }
      if (f.primary && !f.benzylic && !f.allylic) set("hyd", LO, "Primary halide: slow neutral solvolysis; SN2 with hydroxide dominates.");
      if (f.gemPoly) { set("base", HI, "Gem-polyhalides undergo base-induced alpha-elimination (carbene) / haloform-type cleavage."); }
      if (f.el === "F") { set("base", LO); set("hyd", RS); set("photo", LO); set("therm", LO); }
      if (f.el === "I") set("photo", CR, "C-I bond homolysis under near-UV / visible light.");
      if (f.el === "Br") set("photo", HI, "C-Br homolysis under UV.");
      break;
    case "aryl_halide":
      if (f.snar) { set("base", HI, "Ortho/para electron-withdrawing group activates SNAr displacement of the halide by hydroxide/amines."); set("hyd", LO); }
      if (f.el === "I") set("photo", HI, "C-I homolysis / photodehalogenation.");
      if (f.el === "Br") set("photo", MO, "C-Br photodehalogenation under UV.");
      if (f.el === "F") set("photo", LO, "C-F is photostable.");
      break;
    case "thioether":
      if (f.aryl) set("ox", HI, "Aryl sulfides oxidise somewhat more slowly than dialkyl sulfides.");
      break;
    case "nitro":
      if (f.aliphatic) { set("base", HI, "Alpha-H nitro compounds ionise to nitronate salts (pKa ~10); nitronates are prone to Nef-type chemistry."); }
      break;
    case "ether":
      if (f.cyclic) set("ox", CR, "Cyclic ethers (THF, dioxane) autoxidise rapidly to hydroperoxides.");
      if (!f.alphaH) set("ox", LO, "No alpha C-H: resists autoxidation.");
      break;
    case "alkene":
      if (f.styrene) { set("ox", CR); set("photo", HI, "Styrene-type alkenes photodimerise and undergo E/Z isomerisation."); set("therm", HI, "Styrenic alkenes polymerise thermally / radically."); }
      if (f.diene) { set("therm", HI, "Conjugated dienes undergo Diels-Alder dimerisation on heating."); set("ox", CR); }
      if (f.electronPoor) set("ox", MO);
      break;
    case "aromatic":
      if (f.activated) set("ox", HI, "Electron-rich (donor-substituted) ring is readily hydroxylated / oxidised by radicals.");
      if (f.deactivated) { set("ox", LO, "Electron-poor ring resists radical hydroxylation."); }
      break;
    case "hetero_rich":
      if (f.kind === "furan") { set("acid", CR, "Furans ring-open under aqueous acid to 1,4-dicarbonyls."); set("ox", CR, "Furans are oxidised (endoperoxide / ring opening) very readily."); }
      if (f.kind === "pyrrole") { set("acid", HI, "Pyrroles polymerise in acid (pyrrole red)."); set("ox", CR, "Pyrroles darken by air / peroxide oxidation."); }
      if (f.kind === "thiophene") { set("acid", LO); set("ox", MO, "S-oxidation to thiophene S-oxide requires strong oxidants; ring is fairly stable."); }
      if (f.kind === "azole") { set("acid", MO, "Protonation of the pyridine-type N (imidazolium-type salts)."); set("ox", MO); }
      if (f.kind === "indole") { set("acid", HI, "Indoles dimerise / oligomerise under acid (C3 protonation)."); set("ox", CR, "Electron-rich indole C2=C3 is oxidised to oxindole / cleavage products."); }
      break;
    case "phosphine":
      if (f.oxide) { set("ox", RS, "Phosphine oxide is already fully oxidised."); set("acid", LO); }
      break;
    default: break;
  }
  return p;
}

/* ------------------------------------------------------------------ */
/*  Structural site finder (graph based)                               */
/* ------------------------------------------------------------------ */
interface Site { g: string; atoms: number[]; f: Record<string, any>; }
interface GroupEntry { fg: FunctionalGroupReactivity; prof: Prof; sites: Site[]; }
interface Analysis { mol: Mol; rings: number[][]; sites: Site[]; entries: GroupEntry[]; classes: Set<string>; }

const REDOX_METALS = new Set(["Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Ag", "Au", "Pd", "Pt", "Mo", "W", "Ce"]);

function findSites(m: Mol, rings: number[][]): Site[] {
  const sites: Site[] = [];
  const el = (i: number): string => m.atoms[i].el;
  const H = (i: number): number => m.atoms[i].h;
  const q = (i: number): number => m.atoms[i].q;
  const arom = (i: number): boolean => m.atoms[i].arom;
  const add = (g: string, atoms: number[], f: Record<string, any> = {}): void => { sites.push({ g, atoms, f }); };
  const hasDbl = (i: number, e: string): boolean => m.nb(i).some((j) => el(j) === e && m.order(i, j) === 2);
  const isAcyl = (x: number, excl = -1): boolean =>
    el(x) === "C" && m.nb(x).some((j) => j !== excl && (el(j) === "O" || el(j) === "S") && m.order(x, j) === 2);
  const carbonylLike = (x: number): boolean =>
    el(x) === "C" && m.nb(x).some((j) => m.order(x, j) === 2 && (el(j) === "O" || el(j) === "S" || el(j) === "N"));
  const vinyl = (x: number): boolean => el(x) === "C" && !arom(x) && m.nb(x).some((j) => m.order(x, j) === 2 && el(j) === "C");
  const sp3 = (x: number): boolean => el(x) === "C" && !arom(x) && q(x) === 0 && m.nb(x).every((j) => m.order(x, j) === 1);
  const acetalOs = (c: number): number[] =>
    !sp3(c) ? [] : m.nb(c).filter((j) => el(j) === "O" && m.order(c, j) === 1 && !arom(j) && (H(j) > 0 || m.deg(j) === 2) &&
      !m.nbx(j, c).some((k) => isAcyl(k, j)));
  const isEWG = (x: number): boolean => {
    if (el(x) === "N") return q(x) === 1 && m.nb(x).some((j) => el(j) === "O");
    if (el(x) === "C") return carbonylLike(x) || m.nb(x).some((j) => m.order(x, j) === 3) || (sp3(x) && m.nb(x).filter((j) => el(j) === "F").length >= 3);
    if (el(x) === "S") return m.nb(x).some((j) => el(j) === "O" && m.order(x, j) === 2);
    return false;
  };
  const isDonor = (x: number): boolean => {
    if (el(x) === "O") return m.nb(x).every((j) => m.order(x, j) === 1) && !m.nb(x).some((j) => isAcyl(j, x));
    if (el(x) === "N") return q(x) === 0 && m.nb(x).every((j) => m.order(x, j) === 1) && !m.nb(x).some((j) => carbonylLike(j) || (el(j) === "S" && hasDbl(j, "O")));
    return false;
  };

  // ---- quinones / quinone imines (ring based)
  const quinoneAtoms = new Set<number>();
  const quinoneRing = new Set<number>();
  for (const r of rings) {
    if (r.length !== 6 || r.some((i) => arom(i))) continue;
    const exo = r.filter((i) => el(i) === "C" && m.nb(i).some((j) => !r.includes(j) && m.order(i, j) === 2 && (el(j) === "O" || el(j) === "N")));
    if (exo.length !== 2) continue;
    let dc = 0;
    for (let p = 0; p < 6; p++) if (m.order(r[p], r[(p + 1) % 6]) === 2) dc++;
    if (dc !== 2) continue;
    const pi = r.indexOf(exo[0]), pj = r.indexOf(exo[1]);
    const dist = Math.min(Math.abs(pi - pj), 6 - Math.abs(pi - pj));
    if (dist !== 3 && dist !== 1) continue;
    add("quinone", r, { para: dist === 3, exo });
    exo.forEach((i) => quinoneAtoms.add(i));
    r.forEach((i) => quinoneRing.add(i));
  }

  // ---- azides (pre-pass: mark all three nitrogens)
  const azide = new Set<number>();
  for (const n of m.alive()) {
    if (el(n) !== "N" || q(n) !== 1) continue;
    const d = m.nb(n).filter((j) => m.order(n, j) === 2 && el(j) === "N");
    if (d.length === 2) { add("azide", [n, ...d]); azide.add(n); d.forEach((j) => azide.add(j)); }
  }

  const enonePairs = new Set<string>();
  const pk = (a: number, b: number): string => (a < b ? `${a}-${b}` : `${b}-${a}`);

  // ================= carbon centres =================
  for (const c of m.alive()) {
    if (el(c) !== "C") continue;
    const nbs = m.nb(c);
    const dO = nbs.filter((j) => el(j) === "O" && m.order(c, j) === 2);
    const dS = nbs.filter((j) => el(j) === "S" && m.order(c, j) === 2);
    const dN = nbs.filter((j) => el(j) === "N" && m.order(c, j) === 2);
    const trip = nbs.find((j) => m.order(c, j) === 3);
    if (trip !== undefined) {
      if (el(trip) === "N") add("nitrile", [c, trip], { aromatic: nbs.some((j) => j !== trip && arom(j)) });
      else if (el(trip) === "C" && c < trip) add("alkyne", [c, trip], { terminal: H(c) > 0 || H(trip) > 0 });
      continue;
    }
    if (quinoneAtoms.has(c)) continue;
    if (dO.length >= 2) continue; // CO2
    if (dN.length && dO.length) { add("isocyanate", [c, dN[0], dO[0]]); continue; }
    if (dN.length && dS.length) { add("isothiocyanate", [c, dN[0], dS[0]]); continue; }

    if (dO.length === 1) {
      const o = dO[0];
      const subs = nbs.filter((j) => j !== o);
      const hs = subs.filter((j) => ["O", "N", "S"].includes(el(j)) || HALOGENS.has(el(j)));
      const Ox = hs.filter((j) => el(j) === "O"), Nx = hs.filter((j) => el(j) === "N"), Sx = hs.filter((j) => el(j) === "S");
      const Xx = hs.filter((j) => HALOGENS.has(el(j)));
      const carbs = subs.filter((j) => el(j) === "C");
      if (Xx.length) add("acyl_halide", [c, o, Xx[0]]);
      else if (Ox.length === 2 && !Nx.length && !Sx.length) {
        const esterified = Ox.filter((x) => m.nbx(x, c).length > 0);
        if (esterified.length === 2) add("carbonate", [c, o, ...Ox]);
        else add("carbonate_ion", [c, o, ...Ox], { acidic: Ox.some((x) => H(x) > 0) });
      } else if (Ox.length === 1 && Nx.length === 1) {
        const ox = Ox[0];
        add("carbamate", [c, o, ox, Nx[0]], { tert: m.nbx(ox, c).some((x) => sp3(x) && H(x) === 0 && m.nb(x).filter((k) => el(k) === "C").length === 3) });
      } else if (Nx.length === 2) add("urea", [c, o, ...Nx], { arom: arom(c) });
      else if (Ox.length === 1 && !Nx.length && !Sx.length) {
        const ox = Ox[0]; const other = m.nbx(ox, c);
        if (H(ox) > 0 && !other.length) {
          const a = carbs[0];
          let beta = false, malonic = false, ewgA = false, arylAct = false, pyr = false;
          if (a !== undefined) {
            beta = m.nbx(a, c).some((x) => el(x) === "C" && hasDbl(x, "O"));
            malonic = m.nbx(a, c).some((x) => el(x) === "C" && hasDbl(x, "O") && m.nb(x).some((y) => el(y) === "O" && H(y) > 0));
            ewgA = sp3(a) && (m.nb(a).filter((k) => HALOGENS.has(el(k))).length >= 2 || m.nb(a).some((k) => el(k) === "N" && q(k) === 1));
            if (arom(a)) {
              const rg = rings.find((r) => r.length === 6 && r.includes(a) && r.every(arom));
              if (rg) {
                const pos = rg.indexOf(a);
                arylAct = rg.some((k, idx) => { const d = Math.min(Math.abs(idx - pos), 6 - Math.abs(idx - pos)); return (d === 1 || d === 3) && m.nb(k).some((y) => !rg.includes(y) && (el(y) === "O" || el(y) === "N") && H(y) > 0); });
              }
              pyr = m.nb(a).some((y) => arom(y) && el(y) === "N");
            }
          }
          add("acid", [c, o, ox], { aryl: carbs.some(arom), decarb: beta || malonic || ewgA || arylAct || pyr, decarbKind: beta ? "beta" : malonic ? "malonic" : ewgA ? "alpha-EWG" : arylAct ? "aryl" : pyr ? "picolinic" : "" });
        }
        else if (q(ox) === -1 && !other.length) add("carboxylate", [c, o, ox]);
        else if (other.length === 1) {
          const x = other[0];
          if (el(x) === "C") {
            if (isAcyl(x, ox)) { if (c < x) add("anhydride", [c, o, ox, x]); }
            else {
              const ringSize = ringSizeOfBond(rings, c, ox);
              const f = {
                aryl: arom(x), vinyl: vinyl(x), formate: H(c) > 0, ring: ringSize, arom: arom(c),
                tert: sp3(x) && H(x) === 0 && m.nb(x).filter((k) => el(k) === "C").length === 3,
                alphaEWG: carbs.some((a) => sp3(a) && (m.nb(a).filter((k) => HALOGENS.has(el(k))).length >= 2 || m.nb(a).some((k) => isEWG(k)))),
                hinder: carbs.some((a) => sp3(a) && H(a) === 0), alphaEnol: false,
              };
              f.alphaEnol = carbs.some((a) => vinyl(a) && m.nb(a).some((k) => el(k) === "O" && H(k) > 0));
              if (ringSize > 0 && ringSize <= 8) add("lactone", [c, o, ox, x], f);
              else if (arom(x) || vinyl(x)) add("ester_aryl", [c, o, ox, x], f);
              else add("ester", [c, o, ox, x], f);
            }
          } else if (el(x) !== "O") add("ester_aryl", [c, o, ox, x], { activated: true });   // O-O (peracid / diacyl peroxide) is reported as a peroxide
        }
      } else if (Sx.length === 1 && !Ox.length && !Nx.length) add("thioester", [c, o, Sx[0]]);
      else if (Nx.length === 1 && !Ox.length && !Sx.length) {
        const n = Nx[0];
        const otherAcyl = m.nbx(n, c).filter((x) => isAcyl(x, n));
        if (otherAcyl.length) { if (c < otherAcyl[0]) add("imide", [c, o, n, otherAcyl[0]]); }
        else {
          const ring = ringSizeOfBond(rings, c, n);
          const f = {
            arom: arom(c), anilide: m.nbx(n, c).some(arom), tertiary: H(n) === 0, formamide: H(c) > 0, hydrazide: m.nbx(n, c).some((x) => el(x) === "N"),
          };
          if (ring === 4) add("beta_lactam", [c, o, n], f);
          else if (ring > 0 && ring <= 8) add("lactam", [c, o, n], f);
          else add("amide", [c, o, n], f);
        }
      } else if (hs.length === 0) {
        if (H(c) >= 1 && carbs.length <= 1) add("aldehyde", [c, o], { aryl: carbs.some(arom), alphaH: carbs.some((a) => sp3(a) && H(a) > 0) });
        else if (carbs.length >= 1) add("ketone", [c, o, ...carbs], { aryl: carbs.some(arom), alphaH: carbs.some((a) => !arom(a) && H(a) > 0), ring: ringSizeOfAtom(rings, c) });
      }
      // Michael acceptor (any carbonyl with alpha,beta C=C)
      for (const a of carbs) {
        if (!vinyl(a)) continue;
        const beta = m.nb(a).find((j) => j !== c && el(j) === "C" && m.order(a, j) === 2);
        if (beta !== undefined) { add("michael", [c, a, beta]); enonePairs.add(pk(a, beta)); }
      }
    } else if (dS.length) {
      add("thiocarbonyl", [c, dS[0]], { thioamide: nbs.some((j) => el(j) === "N") });
    } else if (dN.length) {
      const n = dN[0];
      const oN = m.nb(n).find((j) => el(j) === "O" && m.order(n, j) === 1);
      const nN = m.nb(n).find((j) => el(j) === "N" && m.order(n, j) === 1);
      const csub = nbs.filter((j) => j !== n && (el(j) === "N" || el(j) === "O"));
      if (oN !== undefined) add("oxime", [c, n, oN], { nitrone: q(n) === 1 });
      else if (nN !== undefined) add("hydrazone", [c, n, nN]);
      else if (csub.length) add("amidine", [c, n, ...csub], { guanidine: csub.filter((x) => el(x) === "N").length >= 1 && csub.length >= 2 });
      else add("imine", [c, n], { aryl: nbs.some(arom) });
    }
  }

  // ---- alkenes / enols / enamines
  for (const c of m.alive()) {
    if (el(c) !== "C" || arom(c)) continue;
    for (const j of m.nb(c)) {
      if (j < c || el(j) !== "C" || arom(j) || m.order(c, j) !== 2) continue;
      if (quinoneRing.has(c) && quinoneRing.has(j)) continue;
      const subs = [...m.nbx(c, j), ...m.nbx(j, c)];
      const enolO = subs.find((x) => el(x) === "O" && H(x) > 0);
      if (enolO !== undefined) { add("enol", [c, j, enolO], { enediol: subs.filter((x) => el(x) === "O" && H(x) > 0).length >= 2 }); continue; }
      if (enonePairs.has(pk(c, j))) continue;
      add("alkene", [c, j], {
        styrene: subs.some(arom), diene: subs.some((x) => vinyl(x)), electronPoor: subs.some((x) => isEWG(x)),
        enolEther: subs.some((x) => el(x) === "O"), enamine: subs.some((x) => el(x) === "N" && q(x) === 0),
      });
    }
  }

  // ---- acetals / hemiacetals
  for (const c of m.alive()) {
    const os = acetalOs(c);
    if (os.length >= 2) {
      const oh = os.filter((x) => H(x) > 0).length;
      if (oh >= 1) add("hemiacetal", [c, ...os], { hydrate: oh >= 2, ring: ringSizeOfAtom(rings, c) });
      else add("acetal", [c, ...os], { ketal: H(c) === 0 });
    }
  }

  // ================= oxygen centres =================
  const peroxKeys = new Set<string>();
  for (const o of m.alive()) {
    if (el(o) !== "O" || arom(o)) continue;
    const nbs = m.nb(o);
    const dbl = nbs.filter((j) => m.order(o, j) === 2);
    if (dbl.length) { if (dbl.length === 1 && el(dbl[0]) === "O" && o < dbl[0]) add("dioxygen", [o, dbl[0]]); continue; }
    if (nbs.length === 0) {
      if (H(o) === 2 && q(o) === 0) add("water", [o]);
      else if (q(o) < 0 || (H(o) === 0 && q(o) === 0)) add("strong_base", [o], { kind: H(o) === 1 ? "hydroxide" : "oxide" });
      continue;
    }
    if (nbs.length === 1) {
      const x = nbs[0];
      if (H(o) > 0) {
        if (el(x) === "C") {
          if (m.nb(x).some((j) => el(j) === "O" && m.order(x, j) === 2)) continue;   // acid (C-centre)
          if (arom(x)) {
            const ring = rings.find((r) => r.includes(x) && r.length === 6 && r.every(arom));
            let hq = false;
            if (ring) {
              const pos = ring.indexOf(x);
              for (const k of ring) {
                const d = Math.min(Math.abs(ring.indexOf(k) - pos), 6 - Math.abs(ring.indexOf(k) - pos));
                if ((d === 1 || d === 3) && m.nbx(k, -1).some((y) => !ring.includes(y) && (el(y) === "O" || (el(y) === "N" && isDonor(y))) && (H(y) > 0 || d === 3))) hq = true;
              }
            }
            const ewg = ring ? ring.some((k) => m.nb(k).some((y) => !ring.includes(y) && isEWG(y))) : false;
            const hindered = m.nb(x).some((k) => arom(k) && m.nb(k).some((y) => sp3(y) && H(y) === 0 && m.nb(y).filter((z) => el(z) === "C").length >= 3));
            add("phenol", [x, o], { hydroquinone: hq, ewg, hindered });
          } else if (vinyl(x)) continue;
          else if (acetalOs(x).length >= 2) continue;
          else if (sp3(x)) {
            const cn = m.nbx(x, o).filter((j) => el(j) === "C").length;
            add(cn === 3 ? "alcohol_tert" : "alcohol", [x, o], { benzylic: m.nbx(x, o).some((j) => arom(j) || vinyl(j)), methanol: H(x) === 3, primary: cn <= 1 });
          }
        } else if (el(x) === "O") {
          const k = pk(o, x);
          if (!peroxKeys.has(k)) { peroxKeys.add(k); add("peroxide", [o, x], { hydroperoxide: true }); }
        }
        continue;
      }
      if (q(o) === -1) {
        if (el(x) === "C" && !carbonylLike(x)) { if (arom(x)) add("phenol", [x, o], { phenolate: true }); else add("strong_base", [o], { kind: "alkoxide" }); }
        else if (["Cl", "Br", "I"].includes(el(x))) add("hypohalite", [o, x]);
      }
      continue;
    }
    if (nbs.length === 2) {
      const [a, b] = nbs;
      if (el(a) === "O" || el(b) === "O") {
        const other = el(a) === "O" ? a : b;
        const k = pk(o, other);
        if (!peroxKeys.has(k)) { peroxKeys.add(k); add("peroxide", [o, other], { hydroperoxide: false }); }
        continue;
      }
      if (el(a) === "C" && el(b) === "C") {
        if (isAcyl(a, o) || isAcyl(b, o)) continue;
        if (acetalOs(a).length >= 2 || acetalOs(b).length >= 2) continue;
        const ringSz = ringSizeOfAtom(rings, o);
        const alphaH = [a, b].some((x) => sp3(x) && H(x) > 0);
        if (ringSz === 3) add("epoxide", [o, a, b]);
        else if (arom(a) || arom(b)) add("aryl_ether", [o, a, b], { alphaH });
        else add("ether", [o, a, b], { alphaH, cyclic: ringSz > 0 });
      }
    }
  }

  // ================= nitrogen centres =================
  for (const n of m.alive()) {
    if (el(n) !== "N" || arom(n) || azide.has(n)) continue;
    const nbs = m.nb(n); const a = m.atoms[n];
    const Os = nbs.filter((j) => el(j) === "O");
    const dblN = nbs.filter((j) => m.order(n, j) === 2);
    if (a.q === 1 && Os.length >= 2 && Os.some((o) => m.order(n, o) === 2 || q(o) === -1)) {
      const others = nbs.filter((j) => el(j) !== "O");
      if (others.length === 1 && el(others[0]) === "C") add("nitro", [n, ...Os], { aromatic: arom(others[0]), aliphatic: !arom(others[0]) && !vinyl(others[0]), alphaH: H(others[0]) > 0 });
      else if (others.length === 1) add("nitro", [n, ...Os], { aliphatic: true });
      else if (others.length === 0) add("nitrate", [n]);
      continue;
    }
    const dO = nbs.find((j) => el(j) === "O" && m.order(n, j) === 2);
    if (dO !== undefined && a.q === 0) {
      const partner = nbs.find((j) => j !== dO);
      if (partner !== undefined) {
        if (el(partner) === "N") add("nitrosamine", [n, dO, partner]);
        else if (el(partner) === "C") add("nitroso", [n, dO, partner]);
        else if (el(partner) === "O") add("nitrite", [n], { ester: q(partner) !== -1 });
      } else add("nitrite", [n], { ester: false });
      continue;
    }
    if (a.q === 1 && Os.length === 1 && m.order(n, Os[0]) === 1 && q(Os[0]) === -1) { if (!dblN.length) add("n_oxide", [n, Os[0]]); continue; }
    if (dblN.length) {
      if (dblN.length === 1 && el(dblN[0]) === "N" && a.q === 0 && q(dblN[0]) === 0 && n < dblN[0]) add("azo", [n, dblN[0]], { aryl: nbs.some(arom) });
      continue;
    }
    if (nbs.some((j) => m.order(n, j) === 3)) continue;
    if (a.q === 1) {
      if (nbs.length === 4 && a.h === 0) add("quat", [n]);
      else if (a.h > 0) add("ammonium", [n]);
      continue;
    }
    if (a.q !== 0) continue;
    if (nbs.some((x) => carbonylLike(x) || ((el(x) === "S" || el(x) === "P") && hasDbl(x, "O")))) continue;
    if (nbs.some((x) => el(x) === "N" && m.nb(x).some((j) => m.order(x, j) === 2))) continue;
    if (nbs.some((x) => el(x) === "O")) { add("hydroxylamine", [n, nbs.find((x) => el(x) === "O")!]); continue; }
    const nn = nbs.find((x) => el(x) === "N");
    if (nn !== undefined) { if (n < nn) add("hydrazine", [n, nn]); continue; }
    const primary = a.h === 2, secondary = a.h === 1, ammonia = nbs.length === 0 && a.h === 3;
    if (nbs.some(arom)) {
      const benz = nbs.some((x) => arom(x) && rings.some((r) => r.length === 6 && r.includes(x) && r.every((i) => arom(i) && el(i) === "C")));
      add("arylamine", [n], { primary, secondary, tertiary: a.h === 0, benz });
    }
    else add(a.h === 0 && nbs.length >= 3 ? "amine_tert" : "amine", [n], { primary: primary || ammonia, secondary, ammonia, enamine: nbs.some((x) => vinyl(x)) });
  }

  // ================= sulfur / phosphorus =================
  for (const s of m.alive()) {
    if (el(s) !== "S" || arom(s)) continue;
    const nbs = m.nb(s);
    const dO = nbs.filter((j) => el(j) === "O" && m.order(s, j) === 2);
    const singles = nbs.filter((j) => m.order(s, j) === 1);
    if (dO.length === 0) {
      if (q(s) === -1 && nbs.length === 1) add("thiol", [s], { anion: true });
      else if (nbs.length === 0 && H(s) === 2) add("thiol", [s], { h2s: true });
      else if (nbs.length === 1 && H(s) > 0) add("thiol", [s]);
      else if (nbs.length === 2) {
        const [a, b] = nbs;
        if (el(a) === "C" && el(b) === "C") add("thioether", [s, a, b], { aryl: arom(a) || arom(b) });
        else if ((el(a) === "S" || el(b) === "S")) { const other = el(a) === "S" ? a : b; if (s < other) add("disulfide", [s, other]); }
      }
    } else if (dO.length === 1) {
      const Os = singles.filter((x) => el(x) === "O");
      if (singles.length === 2 && singles.every((x) => el(x) === "C")) add("sulfoxide", [s, dO[0], ...singles]);
      else if (Os.length >= 1) add("sulfite", [s]);
      else add("sulfoxide", [s, dO[0], ...singles]);
    } else if (dO.length === 2) {
      const Ns = singles.filter((x) => el(x) === "N"), Os = singles.filter((x) => el(x) === "O"), Cs = singles.filter((x) => el(x) === "C");
      const Xs = singles.filter((x) => HALOGENS.has(el(x)));
      if (Xs.length) add("sulfonyl_halide", [s, Xs[0]]);
      else if (Ns.length) add("sulfonamide", [s, Ns[0]], { primary: H(Ns[0]) >= 2 });
      else if (Cs.length && Os.length) { if (m.nbx(Os[0], s).length > 0) add("sulfonate_ester", [s, Os[0], ...Cs]); else add("sulfonic", [s, Os[0]], { acid: H(Os[0]) > 0 }); }
      else if (Cs.length >= 2) add("sulfone", [s]);
      else if (Os.length >= 2) {
        const est = Os.filter((x) => m.nbx(x, s).length > 0);
        if (est.length === 2) add("sulfonate_ester", [s, ...Os], { sulfate: true });
        else if (est.length === 1) add("sulfate", [s], { ester: true });
        else add("sulfate", [s], { acid: Os.some((x) => H(x) > 0) });
      }
    }
  }
  for (const p of m.alive()) {
    if (el(p) !== "P" || arom(p)) continue;
    const nbs = m.nb(p);
    const dO = nbs.filter((j) => el(j) === "O" && m.order(p, j) === 2);
    const Os = nbs.filter((j) => el(j) === "O" && m.order(p, j) === 1);
    const Cs = nbs.filter((j) => el(j) === "C");
    if (dO.length >= 1 && Os.length >= 2) add("phosphate", [p], { ester: Os.some((o) => m.nbx(o, p).some((x) => el(x) === "C")), phosphonate: Cs.length > 0, acid: Os.some((o) => H(o) > 0) });
    else if (dO.length >= 1) add("phosphine", [p], { oxide: true });
    else add("phosphine", [p]);
  }

  // ================= halogens / metals / misc inorganic =================
  const ringOf6 = (y: number): number[] | undefined => rings.find((r) => r.length === 6 && r.includes(y) && r.every(arom));
  const snarActivated = (y: number): boolean => {
    const r = ringOf6(y); if (!r) return false;
    const pos = r.indexOf(y);
    for (let k = 0; k < 6; k++) {
      const d = Math.min(Math.abs(k - pos), 6 - Math.abs(k - pos));
      if (d !== 1 && d !== 3) continue;
      const at = r[k];
      if (el(at) === "N") return true;
      if (m.nb(at).some((z) => !r.includes(z) && isEWG(z))) return true;
    }
    return false;
  };
  for (const x of m.alive()) {
    if (HALOGENS.has(el(x))) {
      const nbs = m.nb(x);
      if (!nbs.length) { if (q(x) === -1) add("halide_ion", [x], { el: el(x) }); else if (H(x) === 1) add("hydrogen_halide", [x], { el: el(x) }); continue; }
      if (nbs.length !== 1) continue;
      const y = nbs[0];
      if (el(y) === "C") {
        if (m.nb(y).some((j) => el(j) === "O" && m.order(y, j) === 2)) continue;
        if (arom(y) || vinyl(y)) add("aryl_halide", [x, y], { el: el(x), aromatic: arom(y), snar: arom(y) && snarActivated(y) });
        else if (sp3(y)) {
          const hal = m.nb(y).filter((j) => HALOGENS.has(el(j)));
          const fs = hal.filter((j) => el(j) === "F");
          if (el(x) === "F" && fs.length >= 2) { if (x === Math.min(...fs)) add("cf3", [y]); continue; }
          const cn = m.nbx(y, x).filter((j) => el(j) === "C").length;
          add("alkyl_halide", [x, y], { el: el(x), primary: cn <= 1, secondary: cn === 2, tertiary: cn === 3, benzylic: m.nbx(y, x).some(arom), allylic: m.nbx(y, x).some(vinyl), gemPoly: hal.length >= 2 });
        }
      } else if ((el(y) === "O" || el(y) === "N") && el(x) !== "F") add("hypohalite", [x, y]);
      continue;
    }
    if (METALS.has(el(x))) {
      const oxo = m.nb(x).filter((j) => el(j) === "O" && m.order(x, j) === 2);
      if (oxo.length >= 2 || (oxo.length >= 1 && m.nb(x).length >= 3)) add("metal_oxide", [x], { el: el(x), oxidant: ["Mn", "Cr", "Mo", "V", "W", "Ce"].includes(el(x)) && el(x) !== "Ti", photocat: el(x) === "Ti" });
      else if (REDOX_METALS.has(el(x))) add("metal_transition", [x], { el: el(x), q: q(x), oxidant: q(x) >= 3 || (el(x) === "Cu" && q(x) >= 2) });
      else add("metal", [x], { el: el(x), q: q(x) });
      continue;
    }
    if (el(x) === "Si") { if (m.nb(x).some((j) => el(j) === "O")) add("metal_oxide", [x], { el: "Si", oxidant: false, photocat: false }); continue; }
    if (el(x) === "B" && m.atoms[x].h >= 1 && q(x) < 0) add("hydride", [x]);
  }

  // ================= carbon frameworks =================
  const seenRing = new Set<string>();
  for (const r of rings) {
    if (r.length !== 5 && r.length !== 6) continue;
    if (!r.every(arom)) continue;
    const key = [...r].sort((a, b) => a - b).join(",");
    if (seenRing.has(key)) continue; seenRing.add(key);
    const het = r.filter((i) => el(i) !== "C");
    if (!het.length) {
      const subs = r.flatMap((i) => m.nb(i).filter((j) => !r.includes(j)));
      const donors = subs.filter(isDonor).length, ewgs = subs.filter(isEWG).length;
      add("aromatic", r, { activated: donors > 0 && ewgs === 0, deactivated: ewgs > 0 && donors === 0 });
      continue;
    }
    if (r.length === 6) { if (het.some((i) => el(i) === "N")) add("hetero_azine", r, { kind: "azine" }); continue; }
    const pyridineN = het.some((i) => el(i) === "N" && H(i) === 0 && m.deg(i) === 2);
    const pyrroleN = het.some((i) => el(i) === "N" && (H(i) > 0 || m.deg(i) === 3));
    const fused = r.some((i) => rings.some((rr) => rr !== r && rr.length === 6 && rr.every(arom) && rr.includes(i)));
    let kind = "azole";
    if (het.some((i) => el(i) === "O")) kind = "furan";
    else if (het.some((i) => el(i) === "S")) kind = "thiophene";
    else if (pyrroleN && !pyridineN) kind = fused ? "indole" : "pyrrole";
    add("hetero_rich", r, { kind });
    if (pyridineN) add("hetero_azine", r, { kind: "azole-N" });
  }
  for (const c of m.alive()) {
    if (!sp3(c) || H(c) === 0) continue;
    if (m.nb(c).some((j) => ["O", "N", "S"].includes(el(j)) || HALOGENS.has(el(j)))) continue;
    const conj = m.nb(c).find((j) => arom(j) || vinyl(j));
    if (conj !== undefined) add("benzylic_ch", [c, conj], { h: H(c), benzylic: arom(conj) });
  }
  if (!sites.length) add("hydrocarbon", m.alive().slice(0, 1), {});
  return sites;
}

/* ------------------------------------------------------------------ */
/*  Group assembly, reactive classes, analysis object                  */
/* ------------------------------------------------------------------ */
function deriveClasses(sites: Site[]): Set<string> {
  const c = new Set<string>();
  const add = (...t: string[]): void => t.forEach((x) => c.add(x));
  for (const s of sites) {
    const f = s.f;
    switch (s.g) {
      case "acid": add("acid_weak"); break;
      case "carboxylate": add("base_weak"); break;
      case "ester": case "lactone": add("ester"); break;
      case "ester_aryl": add("ester", "acylating"); break;
      case "beta_lactam": case "anhydride": case "acyl_halide": case "isocyanate": case "isothiocyanate": case "carbonate": case "sulfonyl_halide": case "thioester": case "imide": add("acylating"); break;
      case "aldehyde": case "ketone": add("carbonyl"); break;
      case "hemiacetal": add("reducing_sugar", "carbonyl", "reductant"); break;
      case "michael": add("michael"); break;
      case "quinone": add("michael", "oxidant"); break;
      case "amine": add("N_nuc", "base"); if (f.primary) add("N_nuc_prim"); if (f.secondary) add("N_nuc_sec"); break;
      case "arylamine": add("N_nuc"); if (f.primary) add("N_nuc_prim"); if (f.secondary) add("N_nuc_sec"); break;
      case "amine_tert": case "amidine": case "hetero_azine": add("base"); break;
      case "ammonium": add("acid_weak"); break;
      case "hydrazine": case "hydroxylamine": add("N_nuc", "N_nuc_prim", "reductant"); break;
      case "alcohol": case "alcohol_tert": add("O_nuc"); break;
      case "phenol": add("O_nuc", "acid_weak"); break;
      case "enol": add("O_nuc", "reductant"); break;
      case "thiol": add("S_nuc", "reductant", "acid_weak"); break;
      case "epoxide": case "alkyl_halide": add("alkylating"); break;
      case "sulfonate_ester": add("alkylating"); break;
      case "peroxide": add("oxidant", "peroxide"); break;
      case "dioxygen": case "hypohalite": case "nitrate": case "n_oxide": add("oxidant"); break;
      case "ether": if (f.alphaH) add("peroxide_former"); break;
      case "aryl_ether": break;
      case "benzylic_ch": add("peroxide_former"); break;
      case "strong_base": add("base", "base_strong"); break;
      case "carbonate_ion": add("base", "base_strong"); break;
      case "water": add("water"); break;
      case "metal": add("metal", "lewis_acid"); break;
      case "metal_transition": add("metal", "redox_metal", "lewis_acid"); if (f.oxidant) add("oxidant"); break;
      case "metal_oxide": if (f.oxidant) add("oxidant"); break;
      case "nitrite": add("nitrosating"); break;
      case "hydride": add("reductant"); break;
      case "sulfite": add("reductant", "S_nuc"); break;
      case "halide_ion": add("halide"); break;
      case "hydrogen_halide": add("acid_strong", "halide"); break;
      case "sulfonic": if (f.acid) add("acid_strong"); break;
      case "sulfate": if (f.acid) add("acid_strong"); break;
      case "phosphate": if (f.acid) add("acid_weak"); break;
      case "alkene": case "alkyne": add("ene"); break;
      case "nitrosamine": add("nitrosating"); break;
      default: break;
    }
  }
  return c;
}

const VRANK = (v: Vuln): number => VORDER.indexOf(v);

function buildEntries(sites: Site[]): GroupEntry[] {
  const byName = new Map<string, { prof: Prof; sites: Site[]; g: string }>();
  for (const s of sites) {
    const prof = tweakProfile(s.g, s.f);
    const cur = byName.get(prof.name);
    if (!cur) { byName.set(prof.name, { prof, sites: [s], g: s.g }); continue; }
    cur.sites.push(s);
    // keep the most reactive variant per condition (worst-case reporting)
    const merged: Prof = { ...cur.prof };
    (["acid", "base", "hyd", "photo", "therm", "ox"] as CondKey[]).forEach((k) => { if (VRANK(prof[k][0]) > VRANK(merged[k][0])) merged[k] = prof[k]; });
    if (VRANK(prof.sec[0]) > VRANK(merged.sec[0])) merged.sec = prof.sec;
    cur.prof = merged;
  }
  const out: GroupEntry[] = [];
  byName.forEach(({ prof, sites: ss }) => {
    const n = ss.length;
    const fg: FunctionalGroupReactivity = {
      groupName: prof.name, category: prof.cat, smilesFragment: prof.frag,
      reactiveSite: prof.site + (n > 1 ? ` (x${n})` : ""),
      acidic: { vulnerability: prof.acid[0], mechanism: prof.acid[1] },
      basic: { vulnerability: prof.base[0], mechanism: prof.base[1] },
      hydrolysis: { vulnerability: prof.hyd[0], mechanism: prof.hyd[1] },
      photolytic: { vulnerability: prof.photo[0], mechanism: prof.photo[1] },
      thermal: { vulnerability: prof.therm[0], mechanism: prof.therm[1] },
      oxidative: { vulnerability: prof.ox[0], mechanism: prof.ox[1] },
      secondaryInteraction: { vulnerability: prof.sec[0], partnerGroup: prof.sec[1], mechanism: prof.sec[2] },
    };
    out.push({ fg, prof, sites: ss });
  });
  return out;
}

function analyze(m: Mol): Analysis {
  const rings = findRings(m, 8);
  const sites = findSites(m, rings);
  return { mol: m, rings, sites, entries: buildEntries(sites), classes: deriveClasses(sites) };
}

const CLASS_TEXT: Record<string, string> = {
  N_nuc: "nucleophilic amine / hydrazine", N_nuc_prim: "primary amine", N_nuc_sec: "secondary amine", O_nuc: "hydroxyl nucleophile (alcohol/phenol)", S_nuc: "thiol / sulfur nucleophile",
  carbonyl: "aldehyde / ketone", reducing_sugar: "reducing sugar (masked aldehyde)", acylating: "acylating agent", ester: "ester", acid_weak: "Bronsted acid (carboxylic acid / phenol)",
  acid_strong: "strong acid", base: "base", base_strong: "strong base / alkaline microenvironment", alkylating: "alkylating agent", michael: "Michael acceptor", oxidant: "oxidant",
  peroxide_former: "peroxide-forming / autoxidisable group", reductant: "reductant", metal: "metal ion", redox_metal: "redox-active metal", lewis_acid: "Lewis acid", nitrosating: "nitrosating agent",
  halide: "halide", water: "water", base_weak: "weakly basic carboxylate salt", ene: "C=C / C#C unsaturation", peroxide: "peroxide",
};

/**
 * Re-evaluate the "secondary compound" column of each group against the ACTUAL co-reactant(s):
 * complementary reactive class present -> keep vulnerability and name the partner; otherwise
 * downgrade (no covalent partner) so the table does not promise chemistry that cannot occur.
 */
function refineSecondary(entries: GroupEntry[], coName: string, coClasses: Set<string>): FunctionalGroupReactivity[] {
  return entries.map(({ fg, prof }) => {
    const hit = prof.partners.filter((c) => coClasses.has(c));
    if (hit.length) {
      return { ...fg, secondaryInteraction: { vulnerability: prof.sec[0], partnerGroup: `${coName}: ${[...new Set(hit.map((h) => CLASS_TEXT[h] || h))].slice(0, 3).join(", ")}`, mechanism: prof.sec[2] } };
    }
    return { ...fg, secondaryInteraction: { vulnerability: stepVuln(prof.sec[0], -2), partnerGroup: `${coName} (no complementary reactive group)`, mechanism: "No matching reactive partner class in the co-reactant; only non-covalent (hydrogen-bond, ionic, dipolar) association expected." } };
  });
}

/**
 * Systematic identification of functional groups and reactive centres from the molecular GRAPH
 * (not text patterns), so the result is independent of SMILES writing order / Kekule vs aromatic form.
 * Reactivity is graded per condition: acid, base, neutral hydrolysis, photolysis, thermal, oxidation,
 * plus the cross-interaction column for secondary compounds.
 */
export function detectFunctionalGroupsDetailed(smiles: string): { features: string[]; sites: string[]; functionalGroups: FunctionalGroupReactivity[] } {
  const m = parseSmiles(smiles);
  if (!m) {
    const fb = buildEntries([{ g: "hydrocarbon", atoms: [], f: {} }]);
    return { features: fb.map((e) => e.fg.groupName), sites: fb.map((e) => e.fg.reactiveSite), functionalGroups: fb.map((e) => e.fg) };
  }
  const an = analyze(m);
  return { features: an.entries.map((e) => e.fg.groupName), sites: an.entries.map((e) => e.fg.reactiveSite), functionalGroups: an.entries.map((e) => e.fg) };
}

/** Backward-compatible wrapper for feature and site extraction */
export function detectFunctionalGroups(smiles: string): { features: string[]; sites: string[] } {
  const { features, sites } = detectFunctionalGroupsDetailed(smiles);
  return { features, sites };
}

/* ------------------------------------------------------------------ */
/*  Reaction templates: graph-editing helpers                          */
/* ------------------------------------------------------------------ */
type Source = "Stress degradation" | "Interaction with other compound";
interface Rx {
  cond: Cond; key: string; label: string; mol: Mol | null; mech: string; dG: number;
  ck: CondKey; vuln: Vuln; mult?: number; source?: Source; co?: string; physical?: boolean;
  mergedParent?: boolean; // true: product derives from the primary+co-reactant mixture (spectators = both compounds)
}

function addO(p: Mol, at: number, order = 1, h = order === 1 ? 1 : 0, q = 0): number {
  const o = p.add("O", { h, q });
  p.bond(at, o, order);
  return o;
}
function addOOH(p: Mol, at: number): void {
  const o1 = p.add("O", { h: 0 }); const o2 = p.add("O", { h: 1 });
  p.bond(at, o1, 1); p.bond(o1, o2, 1);
}
function tautomerizeEnols(p: Mol): void {
  for (const o of p.alive()) {
    if (p.el(o) !== "O" || p.atoms[o].h < 1 || p.nb(o).length !== 1) continue;
    const c = p.nb(o)[0];
    if (p.el(c) !== "C" || p.atoms[c].arom) continue;
    const j = p.nb(c).find((x) => p.el(x) === "C" && p.order(c, x) === 2 && !p.atoms[x].arom);
    if (j === undefined) continue;
    p.atoms[o].h -= 1; p.bond(c, o, 2); p.bond(c, j, 1); p.atoms[j].h += 1;
  }
}
/** hydrolytic cleavage of the C(acyl)-X bond: acid on the acyl side, X-H on the leaving side */
function cleaveAcyl(m: Mol, c: number, lg: number): Mol {
  const p = m.clone();
  p.unbond(c, lg); addO(p, c, 1); p.atoms[lg].h += 1;
  tautomerizeEnols(p);
  return p;
}
/** remove a carbamic / carbonic acid group C(=O)(OH)-X (spontaneous CO2 loss); `keep` receives an H */
function dropCO2(p: Mol, c: number, keep: number): void {
  p.unbond(c, keep); p.atoms[keep].h += 1;
  for (const o of p.nb(c)) p.remove(o);
  p.remove(c);
}
/** E1 / syn elimination on carbon x (loses lg), forms alkene towards the most substituted beta-carbon bearing H */
function eliminateAlkene(p: Mol, x: number, lg: number): boolean {
  const beta = p.nbx(x, lg).filter((b) => p.el(b) === "C" && p.atoms[b].h > 0 && !p.atoms[b].arom && p.order(x, b) === 1);
  if (!beta.length) return false;
  beta.sort((a, b) => p.atoms[a].h - p.atoms[b].h);
  const b = beta[0];
  p.unbond(x, lg); p.bond(x, b, 2); p.atoms[b].h -= 1;
  return true;
}
function shortestPath(m: Mol, a: number, b: number, banned: Set<number>): number {
  const dist = new Map<number, number>([[a, 0]]); const q = [a];
  while (q.length) {
    const u = q.shift()!;
    if (u === b) return dist.get(u)!;
    for (const v of m.nb(u)) { if (banned.has(v) || dist.has(v)) continue; dist.set(v, dist.get(u)! + 1); q.push(v); }
  }
  return -1;
}
function isDonorAtom(m: Mol, x: number): boolean {
  const el = m.el(x);
  const acylN = (k: number): boolean => m.nb(k).some((j) => m.el(j) === "C" && m.nb(j).some((y) => m.order(j, y) === 2 && (m.el(y) === "O" || m.el(y) === "S")));
  if (el === "O") return m.nb(x).every((j) => m.order(x, j) === 1) && !m.nb(x).some((j) => m.el(j) === "C" && m.nb(j).some((y) => y !== x && m.order(j, y) === 2 && m.el(y) === "O"));
  if (el === "N") return m.atoms[x].q === 0 && m.nb(x).every((j) => m.order(x, j) === 1) && !acylN(x);
  return false;
}
function isEwgAtom(m: Mol, x: number): boolean {
  const el = m.el(x);
  if (el === "N") return m.atoms[x].q === 1 && m.nb(x).some((j) => m.el(j) === "O");
  if (el === "C") return m.nb(x).some((j) => m.order(x, j) === 2 && ["O", "S", "N"].includes(m.el(j))) || m.nb(x).some((j) => m.order(x, j) === 3) || m.nb(x).filter((j) => m.el(j) === "F").length >= 3;
  if (el === "S") return m.nb(x).some((j) => m.el(j) === "O" && m.order(x, j) === 2);
  return false;
}
/** most electron-rich free aromatic C-H of a benzenoid ring (radical / electrophilic hydroxylation site) */
function bestArylCH(m: Mol, ring: number[]): { atom: number; score: number } | null {
  if (ring.length !== 6) return null;
  let best: { atom: number; score: number } | null = null;
  ring.forEach((a, ia) => {
    if (m.atoms[a].h < 1 || m.el(a) !== "C") return;
    let sc = 0;
    ring.forEach((k, ik) => {
      const d = Math.min(Math.abs(ia - ik), 6 - Math.abs(ia - ik));
      if (d === 0) return;
      for (const s of m.nb(k)) {
        if (ring.includes(s)) continue;
        const o = d === 1 || d === 3;
        if (isDonorAtom(m, s)) sc += o ? 2 : 0.3;
        else if (m.el(s) === "C" && !m.atoms[s].arom && m.nb(s).every((y) => m.order(s, y) === 1) && !isEwgAtom(m, s)) sc += o ? 1 : 0.2;
        else if (HALOGENS.has(m.el(s))) sc += o ? 0.4 : -0.2;
        else if (isEwgAtom(m, s)) sc += o ? -2 : 0.3;
      }
    });
    if (!best || sc > best.score) best = { atom: a, score: sc };
  });
  return best;
}
function toQuinone(m: Mol, x: number, rings: number[][]): Mol | null {
  const r = rings.find((rr) => rr.length === 6 && rr.includes(x) && rr.every((i) => m.atoms[i].arom && m.el(i) === "C"));
  if (!r) return null;
  if (rings.some((rr) => rr !== r && rr.every((i) => m.atoms[i].arom) && rr.some((i) => r.includes(i)))) return null;
  const s0 = r.indexOf(x);
  const ring = [...r.slice(s0), ...r.slice(0, s0)];
  const ox = m.nb(x).find((j) => m.el(j) === "O" && m.atoms[j].h > 0);
  if (ox === undefined) return null;
  const hetOn = (k: number): number | undefined =>
    m.nb(k).find((j) => !ring.includes(j) && ((m.el(j) === "O" && m.atoms[j].h > 0) || (m.el(j) === "N" && m.atoms[j].h > 0 && m.atoms[j].q === 0)));
  let idx = -1;
  if (hetOn(ring[3]) !== undefined) idx = 3;
  else for (const k of [1, 5]) if (hetOn(ring[k]) !== undefined) { idx = k; break; }
  if (idx < 0) return null;
  const partner = hetOn(ring[idx])!;
  const p = m.clone();
  ring.forEach((i) => { p.atoms[i].arom = false; });
  for (let i = 0; i < 6; i++) p.bond(ring[i], ring[(i + 1) % 6], 1);
  if (idx === 3) { p.bond(ring[1], ring[2], 2); p.bond(ring[4], ring[5], 2); }
  else if (idx === 1) { p.bond(ring[2], ring[3], 2); p.bond(ring[4], ring[5], 2); }
  else { p.bond(ring[1], ring[2], 2); p.bond(ring[3], ring[4], 2); }
  p.atoms[ox].h = 0; p.bond(x, ox, 2);
  if (p.el(partner) === "O") { p.atoms[partner].h = 0; } else { p.atoms[partner].h -= 1; }
  p.bond(ring[idx], partner, 2);
  return p;
}


/** Re-establish a consistent Kekule/aromatic description after a ring was opened (returns false if impossible). */
function sanitizeProduct(p: Mol): boolean {
  const aromRings = findRings(p, 8).filter((r) => r.every((i) => p.atoms[i].arom));
  const inArom = new Set<number>(); aromRings.forEach((r) => r.forEach((i) => inArom.add(i)));
  const resid = p.alive().filter((i) => p.atoms[i].arom && !inArom.has(i));
  if (!resid.length) return true;
  resid.forEach((a) => { p.atoms[a].arom = false; });
  for (const a of resid) for (const b of p.nb(a)) if (p.order(a, b) === 1.5) p.bond(a, b, 1);
  const deficit = (i: number): number => {
    const at = p.atoms[i]; const vals = VALENCE[at.el]; if (!vals) return 0;
    let sum = at.h; for (const j of p.nb(i)) { const o = p.order(i, j); sum += o === 1.5 ? 1 : o; }
    const allowed = vals.map((v) => (at.el === "C" ? v - Math.abs(at.q) : at.el === "B" ? v - at.q : v + at.q)).filter((v) => v >= 0).sort((x, y) => x - y);
    const t = allowed.find((v) => v >= sum);
    return t === undefined ? 0 : t - sum;
  };
  for (let pass = 0; pass < 6; pass++) {
    let changed = false;
    for (const a of resid) {
      if (deficit(a) < 1) continue;
      const b = p.nb(a).find((j) => resid.includes(j) && p.order(a, j) === 1 && deficit(j) >= 1);
      if (b !== undefined) { p.bond(a, b, 2); changed = true; }
    }
    if (!changed) break;
  }
  return resid.every((i) => valenceOK(p, i));
}

/* ------------------------------------------------------------------ */
/*  Product finalisation, naming, descriptors                          */
/* ------------------------------------------------------------------ */
function subMol(m: Mol, comp: number[]): Mol {
  const out = new Mol(); const map = new Map<number, number>();
  comp.forEach((i) => map.set(i, out.add(m.atoms[i].el, { ...m.atoms[i] })));
  comp.forEach((i) => m.adj[i].forEach((o, j) => { if (map.has(j) && i < j) out.bond(map.get(i)!, map.get(j)!, o); }));
  return out;
}
function isTrivialComponent(m: Mol, comp: number[]): boolean {
  const els = comp.map((i) => m.el(i)).sort().join("");
  const hs = comp.reduce((s, i) => s + m.atoms[i].h, 0);
  const nH = comp.length;
  if (nH === 1) {
    const e = m.el(comp[0]);
    if (HALOGENS.has(e)) return true;
    if (e === "O" && hs === 2) return true;      // water
    if (e === "N" && hs === 0) return false;
    if (e === "O" && hs === 0 && m.atoms[comp[0]].q === 0) return true;
  }
  if (els === "CO" + "O" && hs === 0) return true;   // CO2
  if (els === "CO" && hs === 0) return true;          // CO
  if (els === "NN" && hs === 0) return true;          // N2
  if (els === "OOS" && hs === 0) return true;         // SO2
  if (els === "NO" && hs === 0 && nH === 2) return true; // NO
  return false;
}
/** drop trivial co-products (H2O, CO2, CO, N2, SO2, HX, halide ions) and return remaining components */
function pruneProducts(m: Mol): { mol: Mol; dropped: string[] } | null {
  const dropped: string[] = []; const p = m.clone();
  const comps = components(p);
  let kept = 0;
  for (const c of comps) {
    if (isTrivialComponent(p, c)) {
      const e = c.map((i) => p.el(i)).sort().join("");
      dropped.push(e === "O" && c.length === 1 ? "H2O" : e === "COO" ? "CO2" : e === "CO" ? "CO" : e === "NN" ? "N2" : e === "OOS" ? "SO2" : c.length === 1 && HALOGENS.has(p.el(c[0])) ? "H" + p.el(c[0]) : e);
      c.forEach((i) => p.remove(i));
    } else kept++;
  }
  return kept ? { mol: p, dropped } : null;
}
function heavyCount(m: Mol, comp: number[]): number { return comp.filter((i) => m.el(i) !== "H").length; }

function fragmentName(m: Mol, comp: number[]): string {
  const sub = subMol(m, comp);
  const lib = libraryNameFor(writeSmiles(sub));
  if (lib) return lib;
  const { formula } = formulaOf(sub);
  const an = analyze(sub);
  const names = an.entries.map((e) => e.fg.groupName).filter((n) => !/^(Aromatic Ring|Benzylic|Aliphatic \/ Polar|Aliphatic Alcohol)/.test(n));
  const main = names.length ? names[0] : an.entries[0]?.fg.groupName ?? "species";
  return `${formula} ${main.split(" / ")[0].toLowerCase().replace(/\s*\(.*\)/, "")}`;
}

interface Cand {
  key: string; cond: Cond; source: Source; label: string; smiles: string; mol: Mol; mech: string; dG: number; like: number;
  co?: string; physical?: boolean; alsoCond: Cond[]; names: string; dropped: string[]; mainFormula: string; mainMass: number; deltaMass: number;
}

function finalizeRx(rx: Rx, parentCanon: string, spectators: Set<string>, refMass: number): Cand | null {
  if (!rx.mol) return null;
  const work = rx.mol.clone();
  if (!sanitizeProduct(work)) return null;
  if (!molValid(work)) return null;
  const pr = pruneProducts(work);
  if (!pr) return null;
  const prod = pr.mol;
  // spectator components (unchanged parts of a salt / mixture) are not products
  if (spectators.size) {
    for (const c of components(prod)) if (spectators.has(writeIsolated(prod, c))) c.forEach((i) => prod.remove(i));
    if (!prod.alive().length) return null;
  }
  if (!molValid(prod)) return null;
  const comps = components(prod).sort((a, b) => heavyCount(prod, b) - heavyCount(prod, a) || a[0] - b[0]);
  const strs = comps.map((c) => ({ c, s: writeIsolated(prod, c) }));
  const smiles = strs.map((x) => x.s).join(".");
  const canon = writeSmiles(prod);
  if (canon === parentCanon) return null;
  const main = comps[0];
  const { formula, mono } = formulaOf(prod, main);
  const names = strs.map((x) => fragmentName(prod, x.c)).join(" + ");
  const mult = rx.mult ?? 1;
  const like = Math.max(0.01, Math.min(0.99, VSCORE[rx.vuln] * mult));
  return {
    key: rx.key, cond: rx.cond, source: rx.source ?? "Stress degradation", label: rx.label, smiles, mol: prod, mech: rx.mech, dG: rx.dG, like,
    co: rx.co, physical: rx.physical, alsoCond: [], names, dropped: pr.dropped, mainFormula: formula, mainMass: mono, deltaMass: mono - refMass,
  };
}

/* ------------------------------------------------------------------ */
/*  Stress-degradation templates                                       */
/* ------------------------------------------------------------------ */
function isSp3C(m: Mol, x: number): boolean {
  return m.el(x) === "C" && !m.atoms[x].arom && m.atoms[x].q === 0 && m.nb(x).every((j) => m.order(x, j) === 1);
}
function isCarbonylC(m: Mol, x: number): boolean {
  return m.el(x) === "C" && m.nb(x).some((j) => m.order(x, j) === 2 && (m.el(j) === "O" || m.el(j) === "S" || m.el(j) === "N"));
}
/** alkyl carbon on N/O that can be hydroxylated (prefers methyl) */
function pickAlkyl(m: Mol, n: number): number | undefined {
  const cands = m.nb(n).filter((x) => isSp3C(m, x) && m.atoms[x].h >= 1 && !isCarbonylC(m, x));
  cands.sort((a, b) => m.atoms[b].h - m.atoms[a].h);
  return cands[0];
}
function dimer22(m: Mol, a: number, b: number): Mol {
  const { mol: p, off } = mergeMol(m, m);
  p.bond(a, b, 1); p.bond(a + off, b + off, 1); p.bond(a, a + off, 1); p.bond(b, b + off, 1);
  return p;
}
function ringOfAromatic(m: Mol, rings: number[][], at: number): number[] | undefined {
  return rings.find((r) => r.length === 6 && r.includes(at) && r.every((i) => m.atoms[i].arom));
}

function genStress(an: Analysis): Rx[] {
  const m = an.mol; const out: Rx[] = [];
  const el = (i: number): string => m.el(i);
  const H = (i: number): number => m.atoms[i].h;
  const V = (s: Site, ck: CondKey): Vuln => tweakProfile(s.g, s.f)[ck][0];
  const R = (s: Site, cond: Cond, ck: CondKey, key: string, label: string, mol: Mol | null, mech: string, dG: number, mult = 1, vuln?: Vuln): void => {
    if (mol) out.push({ cond, key, label, mol, mech, dG, ck, vuln: vuln ?? V(s, ck), mult });
  };
  const hyd3 = (s: Site, key: string, label: string, mol: Mol | null, tx: [string, string, string], dG: [number, number, number], mult: [number, number, number] = [1, 1, 1]): void => {
    R(s, "Acidic Hydrolysis", "acid", key, label, mol, tx[0], dG[0], mult[0]);
    R(s, "Basic Hydrolysis", "base", key, label, mol, tx[1], dG[1], mult[1]);
    R(s, "Hydrolysis", "hyd", key, label, mol, tx[2], dG[2], mult[2]);
  };
  const cyc: Site[] = [];

  for (const s of an.sites) {
    const f = s.f;
    switch (s.g) {
      // ---------------------------------------------------------------- hydrolysis of acyl derivatives
      case "ester": case "ester_aryl": case "lactone": {
        const [c, , ox, x] = s.atoms;
        const ring = f.ring || 0;
        const kind = s.g === "ester" ? "ester" : s.g === "ester_aryl" ? "aryl/vinyl ester" : ring === 4 ? "strained beta-lactone" : "lactone";
        let dG: [number, number, number] = [-2.0, -6.0, -1.5];
        if (s.g === "ester_aryl") dG = [-3.0, -8.0, -3.5];
        if (s.g === "lactone") dG = ring === 4 ? [-8, -9, -7.5] : f.arom ? [-3, -7, -2] : ring === 5 ? [0.5, -4.0, 0.2] : ring === 6 ? [-0.5, -4.5, -0.3] : [-2, -6, -1.5];
        const p = cleaveAcyl(m, c, ox);
        hyd3(s, "acyl-O-cleavage", `${kind} hydrolysis`, p, [
          `Acid-catalysed acyl-oxygen cleavage (A_AC2): the protonated carbonyl is attacked by water and the ${kind} splits into carboxylic acid + alcohol/phenol${s.g === "lactone" ? " (ring opened to the hydroxy-acid)" : ""}.`,
          `Saponification (B_AC2): hydroxide adds to the carbonyl and expels the alkoxide/phenoxide; irreversible once the carboxylate forms.`,
          `Water-mediated hydrolysis of the ${kind} under humidity / neutral pH${s.g === "ester_aryl" ? " (the activated leaving group accelerates it)" : ""}.`,
        ], dG);
        if (f.tert) {
          const p2 = m.clone();
          if (eliminateAlkene(p2, x, ox)) {
            p2.atoms[ox].h += 1;
            R(s, "Acidic Hydrolysis", "acid", "tbu-cleavage", "tert-alkyl ester cleavage", p2, "A_AL1 / E1: protonation and loss of the tertiary carbocation as an alkene (e.g. isobutene) leaves the free carboxylic acid; no water is consumed.", -4.0);
            R(s, "Thermal Degradation", "therm", "tbu-cleavage", "tert-alkyl ester thermolysis", p2, "Thermal E1 / syn-elimination expels the alkene and leaves the carboxylic acid.", -4.0);
          }
        }
        if (s.g !== "lactone") cyc.push(s);
        break;
      }
      case "amide": case "lactam": case "beta_lactam": {
        const [c, , n] = s.atoms;
        if (f.arom) break;
        const ring = ringSizeOfBond(an.rings, c, n);
        const dG: [number, number, number] = s.g === "beta_lactam" ? [-8, -9, -7.5] : s.g === "lactam" ? (ring === 5 ? [-1.0, -2.0, -0.3] : [-2.0, -3.0, -0.8]) : [-2.5, -3.0, -0.5];
        const p = cleaveAcyl(m, c, n);
        const what = s.g === "beta_lactam" ? "strained beta-lactam" : s.g === "lactam" ? "lactam" : "amide";
        hyd3(s, "amide-cleavage", `${what} hydrolysis`, p, [
          `Acid-catalysed hydrolysis (A_AC2): O-protonation, water attack and C-N cleavage give the carboxylic acid + ammonium${s.g === "beta_lactam" ? "; ring strain (~26 kcal/mol) makes this fast" : s.g === "lactam" ? " (ring opened to the amino-acid)" : ""}.`,
          `Hydroxide attacks the ${what} carbonyl; C-N cleavage gives the carboxylate + amine${s.g === "beta_lactam" ? " (irreversible ring opening)" : ""}.`,
          `Neutral hydrolysis of the ${what}${s.g === "beta_lactam" ? " driven by ring strain" : " is very slow (half-lives of years) at ambient pH"}.`,
        ], dG);
        break;
      }
      case "imide": {
        const [c, , n] = s.atoms;
        const p = cleaveAcyl(m, c, n);
        hyd3(s, "imide-opening", "imide hydrolysis (amic acid formation)", p, [
          "Acid-catalysed C-N cleavage of the imide to the amic acid (mono-amide / mono-acid).",
          "Hydroxide opens the imide ring at one carbonyl to give the amic acid; fast because the imide carbonyls are strongly electrophilic.",
          "Slow neutral ring opening of the imide to the amic acid on moisture exposure.",
        ], [-2.0, -6.0, -2.5]);
        break;
      }
      case "urea": {
        if (f.arom) break;
        const [c, , n1, n2] = s.atoms;
        const p = cleaveAcyl(m, c, n1); dropCO2(p, c, n2);
        hyd3(s, "urea-cleavage", "urea hydrolysis (via carbamic acid, CO2 loss)", p, [
          "Acid hydrolysis of the urea carbonyl gives a carbamic acid that loses CO2 to two amines.",
          "Base hydrolysis (slow) gives the carbamate / carbamic acid, which decarboxylates to the amines.",
          "Very slow neutral hydrolysis of the urea to two amines + CO2.",
        ], [-3.0, -2.0, -1.0]);
        break;
      }
      case "carbamate": {
        const [c, , ox, n] = s.atoms;
        const p = cleaveAcyl(m, c, ox); dropCO2(p, c, n);
        hyd3(s, "carbamate-cleavage", "carbamate hydrolysis (via carbamic acid, CO2 loss)", p, [
          "Acid hydrolysis of the carbamate gives the alcohol + carbamic acid, which loses CO2 to the free amine.",
          "Base hydrolysis (BAc2; E1cB via isocyanate for N-H aryl carbamates) releases the alcohol and, after CO2 loss, the amine.",
          "Slow neutral hydrolysis of the carbamate to amine + alcohol + CO2.",
        ], [-3.0, -3.5, -1.5]);
        if (f.tert) {
          const x = m.nbx(ox, c)[0]; const p2 = m.clone();
          if (x !== undefined && eliminateAlkene(p2, x, ox)) {
            p2.atoms[ox].h += 1; dropCO2(p2, c, n);
            R(s, "Acidic Hydrolysis", "acid", "boc-cleavage", "tert-butyl carbamate (Boc-type) cleavage", p2, "A_AL1 / E1 acidolysis: loss of the tert-alkyl cation as an alkene and CO2 releases the free amine.", -4.5);
            R(s, "Thermal Degradation", "therm", "boc-cleavage", "tert-butyl carbamate thermolysis", p2, "Thermolysis (>150 C) expels the alkene and CO2 to give the free amine.", -4.0);
          }
        }
        if (H(n) > 0) {
          const p3 = m.clone(); p3.unbond(c, ox); p3.atoms[ox].h += 1; p3.atoms[n].h -= 1; p3.bond(c, n, 2);
          R(s, "Thermal Degradation", "therm", "carbamate-to-isocyanate", "carbamate thermolysis to isocyanate", p3, "Reversible thermal dissociation of an N-H carbamate into isocyanate + alcohol (above ~150 C).", 4.0);
        }
        break;
      }
      case "carbonate": {
        const [c, , ox1, ox2] = s.atoms;
        const p = cleaveAcyl(m, c, ox1); dropCO2(p, c, ox2);
        hyd3(s, "carbonate-cleavage", "carbonate hydrolysis (CO2 loss)", p, [
          "Acid-catalysed hydrolysis to two alcohols + CO2.", "Hydroxide attacks the carbonate carbonyl; two alkoxides + carbonate result.", "Slow neutral hydrolysis to two alcohols + CO2.",
        ], [-5.5, -6.5, -4.0]);
        break;
      }
      case "anhydride": {
        const [c, , ox] = s.atoms;
        const p = cleaveAcyl(m, c, ox);
        hyd3(s, "anhydride-hydrolysis", "anhydride hydrolysis", p, [
          "Acid-catalysed hydrolysis of the anhydride to two carboxylic acids.", "Instant saponification to two carboxylates.", "Spontaneous hydrolysis by ambient moisture to two carboxylic acids.",
        ], [-8.5, -10, -9]);
        break;
      }
      case "acyl_halide": {
        const [c, , x] = s.atoms;
        const p = cleaveAcyl(m, c, x);
        hyd3(s, "acyl-halide-hydrolysis", "acyl halide hydrolysis", p, [
          "Immediate hydrolysis to the carboxylic acid + HX.", "Instant saponification.", "Reacts violently with atmospheric moisture to the carboxylic acid + HX.",
        ], [-10, -11, -10]);
        break;
      }
      case "thioester": {
        const [c, , sx] = s.atoms;
        const p = cleaveAcyl(m, c, sx);
        hyd3(s, "thioester-hydrolysis", "thioester hydrolysis", p, [
          "Acid hydrolysis of the thioester to carboxylic acid + thiol.", "Base hydrolysis / thiolate exchange releases the carboxylate + thiol.", "Slow neutral hydrolysis to carboxylic acid + thiol.",
        ], [-2.0, -4.5, -2.5]);
        break;
      }
      case "sulfonate_ester": {
        const sAt = s.atoms[0];
        const oxs = m.nb(sAt).filter((j) => el(j) === "O" && m.order(sAt, j) === 1 && m.nbx(j, sAt).length > 0);
        if (!oxs.length) break;
        const oo = oxs[0]; const r = m.nbx(oo, sAt)[0];
        if (el(r) !== "C") break;
        const p = m.clone();
        if (m.atoms[r].arom) { p.unbond(sAt, oo); p.atoms[oo].h += 1; addO(p, sAt, 1); }
        else { p.unbond(oo, r); p.atoms[oo].h += 1; addO(p, r, 1); }
        hyd3(s, "sulfonate-ester-hydrolysis", "sulfonate/sulfate ester hydrolysis (alkylating agent destroyed)", p, [
          "Acid hydrolysis to the sulfonic acid + alcohol.", "Hydroxide attacks carbon (SN2) or sulfur, releasing the sulfonate + alcohol.", "Solvolysis in water/alcohol: the alkyl sulfonate alkylates water to give the alcohol + sulfonic acid.",
        ], [-5.0, -7.0, -5.0]);
        break;
      }
      case "sulfonyl_halide": {
        const [sAt, x] = s.atoms;
        const p = m.clone(); p.unbond(sAt, x); addO(p, sAt, 1); p.atoms[x].h += 1;
        hyd3(s, "sulfonyl-halide-hydrolysis", "sulfonyl halide hydrolysis", p, [
          "Acid-mediated hydrolysis to the sulfonic acid + HX.", "Rapid base hydrolysis to the sulfonate.", "Moisture hydrolyses the sulfonyl halide to the sulfonic acid + HX.",
        ], [-6, -8, -7]);
        break;
      }
      case "acetal": {
        const c = s.atoms[0]; const os = s.atoms.slice(1);
        const ringO = os.filter((o) => ringSizeOfBond(an.rings, c, o) > 0);
        const exo = os.filter((o) => ringSizeOfBond(an.rings, c, o) === 0);
        let p: Mol; let lab: string; let dG = -1.5;
        if (ringO.length === 1 && exo.length >= 1) { p = m.clone(); p.unbond(c, exo[0]); p.atoms[exo[0]].h += 1; addO(p, c, 1); lab = "glycoside / cyclic acetal cleavage (aglycone released)"; dG = -4.0; }
        else { p = m.clone(); p.unbond(c, os[0]); p.unbond(c, os[1]); p.atoms[os[0]].h += 1; p.atoms[os[1]].h += 1; addO(p, c, 2, 0); lab = f.ketal ? "ketal hydrolysis to ketone + alcohols" : "acetal hydrolysis to aldehyde + alcohols"; }
        R(s, "Acidic Hydrolysis", "acid", "acetal-hydrolysis", lab, p, "Protonation of an acetal oxygen, loss of alcohol to the oxocarbenium ion, then water capture (A1 mechanism); very fast even at mild acidity.", dG);
        R(s, "Hydrolysis", "hyd", "acetal-hydrolysis", lab, p, "Neutral hydrolysis is slow, but trace acid autocatalysis in the presence of moisture can cleave the acetal.", dG, 0.5);
        break;
      }
      case "imine": case "oxime": case "hydrazone": {
        const [c, n] = s.atoms;
        if (m.atoms[n].q !== 0) break;
        const p = m.clone(); p.unbond(c, n); addO(p, c, 2, 0); p.atoms[n].h += 2;
        const dG: [number, number, number] = s.g === "imine" ? [-2, -1, -1] : s.g === "oxime" ? [-0.5, -0.5, -0.3] : [-1, -0.5, -0.8];
        hyd3(s, "c=n-hydrolysis", `${s.g === "imine" ? "imine (Schiff base)" : s.g} hydrolysis to the carbonyl compound`, p, [
          "Protonation gives an iminium that is attacked by water; the carbinolamine collapses to the carbonyl + amine/hydroxylamine/hydrazine.",
          "Base-catalysed hydrolysis of the C=N bond to the carbonyl compound + amine.", "Reversible hydrolysis of the C=N bond in water (equilibrium shifts with dilution).",
        ], dG);
        break;
      }
      case "nitrile": {
        const [c, n] = s.atoms;
        const p = m.clone(); p.bond(c, n, 1); p.atoms[n].h = 2; addO(p, c, 2, 0);
        R(s, "Acidic Hydrolysis", "acid", "nitrile-hydration", "nitrile hydration to the primary amide", p, "Acid-catalysed hydration (Ritter-type) of the C#N to the primary amide; further hydrolysis to the acid + NH4+ needs forcing conditions.", -4.0);
        R(s, "Basic Hydrolysis", "base", "nitrile-hydration", "nitrile hydration to the primary amide", p, "Hydroxide adds to the nitrile carbon; tautomerisation gives the primary amide (acid + NH3 on forcing).", -4.5);
        break;
      }
      case "epoxide": {
        const [o, , b] = s.atoms;
        const p = m.clone(); p.unbond(o, b); p.atoms[o].h = 1; addO(p, b, 1);
        hyd3(s, "epoxide-opening", "epoxide ring opening to the 1,2-diol", p, [
          "Protonated epoxide is opened by water at the more substituted carbon (SN1-like) to the vicinal diol.", "Hydroxide opens the strained ring by SN2 at the less hindered carbon to the diol.", "Slow neutral hydrolysis of the strained oxirane to the diol.",
        ], [-9, -8, -7]);
        break;
      }
      case "alkyl_halide": {
        const [x, y] = s.atoms;
        if (f.el === "F") break;
        const sn1 = f.tertiary || f.benzylic || f.allylic;
        const sub = m.clone(); sub.unbond(x, y); sub.atoms[x].q = -1; sub.atoms[x].h = 0; addO(sub, y, 1);
        R(s, "Basic Hydrolysis", "base", "halide-substitution", "alkyl halide substitution by hydroxide", sub, "SN2 (primary/secondary) or SN1 (tertiary/benzylic/allylic) displacement of the halide by hydroxide gives the alcohol.", -4.0, f.primary && !sn1 ? 0.8 : 1);
        R(s, "Hydrolysis", "hyd", "halide-substitution", "alkyl halide solvolysis", sub, "Solvolysis: ionisation (SN1) of tertiary/benzylic/allylic halides, or slow SN2 by water, gives the alcohol + HX.", -4.0, sn1 ? 0.9 : 0.4);
        if (sn1) R(s, "Acidic Hydrolysis", "acid", "halide-substitution", "alkyl halide solvolysis", sub, "Acid-assisted SN1 solvolysis of the labile halide to the alcohol + HX.", -4.0, 0.8);
        const beta = m.nbx(y, x).filter((b) => el(b) === "C" && H(b) > 0 && !m.atoms[b].arom && m.order(y, b) === 1);
        if (beta.length) {
          beta.sort((a, b) => H(a) - H(b)); const b = beta[0];
          const p = m.clone(); p.unbond(x, y); p.atoms[x].q = -1; p.atoms[x].h = 0; p.bond(y, b, 2); p.atoms[b].h -= 1;
          R(s, "Basic Hydrolysis", "base", "dehydrohalogenation", "dehydrohalogenation (E2/E1) to the alkene", p, "Base-induced elimination of HX (Zaitsev alkene); dominant for tertiary halides, minor for primary.", -3.0, f.tertiary ? 1 : f.secondary ? 0.6 : 0.25);
        }
        break;
      }
      case "aryl_halide": {
        const [x, y] = s.atoms;
        if (!f.snar) break;
        const p = m.clone(); p.unbond(x, y); p.atoms[x].q = -1; p.atoms[x].h = 0; addO(p, y, 1);
        R(s, "Basic Hydrolysis", "base", "snar-hydroxide", "SNAr displacement of the halide by hydroxide (phenol formation)", p, "Hydroxide adds to the ring carbon activated by an ortho/para electron-withdrawing group (Meisenheimer complex) and expels halide.", -6.0);
        R(s, "Hydrolysis", "hyd", "snar-hydroxide", "SNAr displacement of the halide by water", p, "Activated aryl halides are slowly hydrolysed to the phenol by water / weak nucleophiles.", -5.0, 0.5);
        break;
      }
      case "alkene": {
        const [c, j] = s.atoms;
        if (f.enolEther || f.enamine) {
          for (const [x, y] of [[c, j], [j, c]]) {
            const hx = m.nbx(x, y).find((k) => el(k) === "O" || (el(k) === "N" && m.atoms[k].q === 0));
            if (hx === undefined || (el(hx) === "O" && m.nb(hx).some((z) => z !== x && isCarbonylC(m, z)))) continue;
            const p = m.clone(); p.unbond(x, hx); p.atoms[hx].h += 1; p.bond(x, y, 1); p.atoms[y].h += 1; addO(p, x, 2, 0);
            R(s, "Acidic Hydrolysis", "acid", "enol-ether-hydrolysis", f.enamine ? "enamine / enamide hydrolysis to the carbonyl compound" : "enol ether hydrolysis to the carbonyl compound", p, "C-protonation of the electron-rich alkene gives an oxocarbenium/iminium ion that water converts to the carbonyl compound + alcohol/amine.", -2.5, 1, HI);
            R(s, "Hydrolysis", "hyd", "enol-ether-hydrolysis", "enol ether / enamine hydrolysis", p, "Moisture-driven hydrolysis of the heteroatom-substituted alkene to the carbonyl compound.", -1.0, 1, MO);
            break;
          }
          break;
        }
        if (f.electronPoor) break;
        const hi = H(c) <= H(j) ? c : j; const lo = hi === c ? j : c;
        const p = m.clone(); p.bond(c, j, 1); addO(p, hi, 1); p.atoms[lo].h += 1;
        R(s, "Acidic Hydrolysis", "acid", "alkene-hydration", "Markovnikov hydration of the C=C bond", p, "Protonation forms the more stable carbocation, trapped by water (Markovnikov hydration to the alcohol).", -1.8);
        break;
      }
      case "michael": {
        const [, a, b] = s.atoms;
        const p = m.clone(); p.bond(a, b, 1); addO(p, b, 1); p.atoms[a].h += 1;
        hyd3(s, "conjugate-hydration", "conjugate (oxa-Michael) hydration to the beta-hydroxy carbonyl", p, [
          "Acid-catalysed conjugate addition of water to the beta-carbon (reversible).", "Hydroxide adds to the beta-carbon (oxa-Michael); reversible retro-Michael / retro-aldol follows.", "Slow reversible conjugate hydration in water.",
        ], [0.8, 0.5, 1.2], [0.5, 0.5, 0.5]);
        break;
      }

      case "isocyanate": {
        const [c, n, o] = s.atoms;
        const p = m.clone(); p.remove(o); p.remove(c); p.atoms[n].h += 2;
        hyd3(s, "isocyanate-hydrolysis", "isocyanate hydrolysis to the amine (CO2 loss)", p, [
          "Water adds to N=C=O to give the carbamic acid, which loses CO2 to the amine.", "Hydroxide adds to the isocyanate carbon; carbamate decarboxylates to the amine.", "Reacts with moisture: carbamic acid then amine + CO2 (the amine can add to remaining isocyanate to form ureas).",
        ], [-8, -9, -8]);
        break;
      }
      case "quat": {
        const n = s.atoms[0];
        let done = false;
        for (const a of m.nb(n).filter((k) => isSp3C(m, k))) {
          const b = m.nbx(a, n).find((k) => isSp3C(m, k) && H(k) >= 1);
          if (b === undefined) continue;
          const p = m.clone(); p.unbond(n, a); p.atoms[n].q = 0; p.bond(a, b, 2); p.atoms[b].h -= 1;
          R(s, "Basic Hydrolysis", "base", "hofmann", "Hofmann elimination to the alkene + tertiary amine", p, "Hydroxide removes a beta-hydrogen (E2, Hofmann rule) and expels the neutral tertiary amine; needs heat.", -3.0, 0.5);
          done = true; break;
        }
        if (!done) {
          const a = m.nb(n).find((k) => isSp3C(m, k) && H(k) === 3);
          if (a !== undefined) {
            const p = m.clone(); p.unbond(n, a); p.atoms[n].q = 0; addO(p, a, 1);
            R(s, "Basic Hydrolysis", "base", "quat-dealkylation", "SN2 dealkylation of the quaternary ammonium by hydroxide", p, "Hydroxide displaces the neutral tertiary amine from a methyl group (SN2; high temperature).", -3.0, 0.4);
          }
        }
        break;
      }
      case "sulfate": {
        if (!f.ester) break;
        const sAt = s.atoms[0];
        const oc = m.nb(sAt).find((k) => el(k) === "O" && m.nbx(k, sAt).some((z) => el(z) === "C"));
        if (oc === undefined) break;
        const r = m.nbx(oc, sAt)[0];
        const p = m.clone(); p.unbond(oc, r); p.atoms[oc].h += 1; addO(p, r, 1);
        R(s, "Acidic Hydrolysis", "acid", "sulfate-ester-hydrolysis", "alkyl sulfate hydrolysis to the alcohol + hydrogen sulfate", p, "Protonation of the sulfate ester and SN2/SN1 attack of water at carbon releases the alcohol + hydrogen sulfate (autocatalytic as acid builds up).", -4.0);
        R(s, "Hydrolysis", "hyd", "sulfate-ester-hydrolysis", "alkyl sulfate hydrolysis to the alcohol + hydrogen sulfate", p, "Slow neutral hydrolysis of the alkyl sulfate; accelerated by heat and by the acid it liberates.", -4.0, 0.7);
        break;
      }
      case "amidine": {
        const [c, n] = s.atoms;
        if (m.atoms[n].q !== 0) break;
        const p = m.clone(); p.unbond(c, n); addO(p, c, 2, 0); p.atoms[n].h += 2;
        hyd3(s, "amidine-hydrolysis", f.guanidine ? "guanidine / biguanide hydrolysis to the urea + amine" : "amidine hydrolysis to the amide + amine", p, [
          "Protonated amidinium is attacked by water; the tetrahedral intermediate expels ammonia/amine to give the amide (urea for guanidines).", "Hydroxide adds to the amidine carbon and expels the amine to give the amide / urea (needs heat).", "Slow neutral hydrolysis of the amidine / guanidine to the amide / urea + amine.",
        ], [-2.0, -2.5, -1.5]);
        break;
      }
      case "phosphate": {
        const pAt = s.atoms[0];
        const ob = m.nb(pAt).find((k) => el(k) === "O" && m.order(pAt, k) === 1 && m.nbx(k, pAt).some((z) => el(z) === "P"));
        if (ob !== undefined) {
          const p = m.clone(); p.unbond(pAt, ob); p.atoms[ob].h += 1; addO(p, pAt, 1);
          R(s, "Acidic Hydrolysis", "acid", "phosphoanhydride-hydrolysis", "phosphoanhydride (P-O-P) hydrolysis", p, "Protonation of the bridging oxygen and water attack at phosphorus cleave the high-energy P-O-P bond to two phosphates.", -7, 1, HI);
          R(s, "Hydrolysis", "hyd", "phosphoanhydride-hydrolysis", "phosphoanhydride (P-O-P) hydrolysis", p, "Water attacks phosphorus of the pyrophosphate-type linkage (metal-ion and pH dependent).", -7, 1, MO);
          R(s, "Basic Hydrolysis", "base", "phosphoanhydride-hydrolysis", "phosphoanhydride (P-O-P) hydrolysis", p, "Hydroxide attack at phosphorus (slowed by charge repulsion of the polyanion).", -7, 0.8, MO);
        }
        if (!f.ester) break;
        const oc = m.nb(pAt).find((k) => el(k) === "O" && m.order(pAt, k) === 1 && m.nbx(k, pAt).some((z) => el(z) === "C"));
        if (oc === undefined) break;
        const p = m.clone(); p.unbond(pAt, oc); p.atoms[oc].h += 1; addO(p, pAt, 1);
        hyd3(s, "phosphate-ester-hydrolysis", "phosphate ester hydrolysis (P-O cleavage)", p, [
          "Acid-catalysed hydrolysis of the P-O-C ester to the phosphate diester/monoester + alcohol (slow).", "Hydroxide attacks phosphorus (SN2@P) and expels the alkoxide (triesters fast, diesters very slow).", "Slow neutral hydrolysis of the phosphate ester.",
        ], [-3.0, -3.5, -2.5]);
        break;
      }
      case "alkyne": {
        const [c, j] = s.atoms;
        const a = H(c) <= H(j) ? c : j; const b = a === c ? j : c;
        const p = m.clone(); p.bond(a, b, 1); addO(p, a, 2, 0); p.atoms[b].h += 2;
        R(s, "Acidic Hydrolysis", "acid", "alkyne-hydration", "alkyne hydration to the ketone / aldehyde (Markovnikov)", p, "Acid (usually Hg2+/Au+ assisted) hydration via the enol gives the methyl ketone (aldehyde for ethyne).", -5.0, 0.6);
        break;
      }
      case "hetero_rich": {
        if (f.kind !== "furan") break;
        const r = s.atoms;
        const oIdx = r.findIndex((i) => el(i) === "O"); if (oIdx < 0) break;
        const ring = [...r.slice(oIdx), ...r.slice(0, oIdx)];       // O, a1, b1, b2, a2
        const [oR, a1, b1, b2, a2] = ring;
        if (ring.length !== 5) break;
        const p = m.clone();
        ring.forEach((i) => { p.atoms[i].arom = false; });
        p.unbond(oR, a2); p.bond(oR, a1, 2); p.atoms[oR].h = 0;
        p.bond(a1, b1, 1); p.bond(b1, b2, 1); p.bond(b2, a2, 1);
        p.atoms[b1].h += 1; p.atoms[b2].h += 1; addO(p, a2, 2, 0);
        R(s, "Acidic Hydrolysis", "acid", "furan-ring-opening", "acid-catalysed furan ring opening to the 1,4-dicarbonyl", p, "C2-protonation of the furan gives an oxocarbenium ion; water addition and ring opening give the saturated 1,4-dicarbonyl (reverse Paal-Knorr).", -3.0);
        break;
      }

      // ---------------------------------------------------------------- oxidation
      case "thioether": {
        const sAt = s.atoms[0];
        let p = m.clone(); addO(p, sAt, 2, 0);
        R(s, "Oxidation", "ox", "s-oxidation", "sulfide to sulfoxide (S-oxidation)", p, "Electrophilic oxygen transfer from peroxide / O2 to the nucleophilic sulfur lone pair gives the sulfoxide (-S(=O)-).", -7.5);
        p = m.clone(); addO(p, sAt, 2, 0); addO(p, sAt, 2, 0);
        R(s, "Oxidation", "ox", "s-oxidation-sulfone", "sulfide to sulfone (over-oxidation)", p, "Excess oxidant converts the sulfoxide onward to the sulfone (-SO2-).", -6.5, 0.5);
        break;
      }
      case "sulfoxide": {
        const p = m.clone(); addO(p, s.atoms[0], 2, 0);
        R(s, "Oxidation", "ox", "sulfoxide-to-sulfone", "sulfoxide to sulfone", p, "Further oxygen transfer to the sulfoxide sulfur gives the sulfone.", -6.0);
        // thermal syn-elimination (needs beta-H) and photo-deoxygenation
        const sAt = s.atoms[0]; const o = s.atoms[1];
        for (const a of m.nb(sAt).filter((k) => el(k) === "C" && isSp3C(m, k))) {
          const b = m.nbx(a, sAt).find((k) => isSp3C(m, k) && H(k) >= 1);
          if (b === undefined) continue;
          const q = m.clone(); q.unbond(sAt, a); q.bond(a, b, 2); q.atoms[b].h -= 1; q.bond(sAt, o, 1); q.atoms[o].h = 1;
          R(s, "Thermal Degradation", "therm", "sulfoxide-syn-elimination", "sulfoxide pyrolysis (syn-elimination) to alkene + sulfenic acid", q, "Concerted five-membered cyclic (Ei) syn-elimination expels the sulfenic acid and forms the alkene.", -0.5);
          break;
        }
        const d = m.clone(); d.remove(o);
        R(s, "Photodegradation", "photo", "sulfoxide-deoxygenation", "sulfoxide photo-deoxygenation to the sulfide", d, "UV excitation cleaves the S=O bond (photo-deoxygenation) to give the sulfide.", 2.5, 0.7);
        break;
      }
      case "thiol": {
        const sAt = s.atoms[0];
        if (H(sAt) < 1) break;
        const { mol: p, off } = mergeMol(m, m); p.atoms[sAt].h -= 1; p.atoms[sAt + off].h -= 1; p.bond(sAt, sAt + off, 1);
        R(s, "Oxidation", "ox", "disulfide-formation", "oxidative dimerisation to the disulfide", p, "Two thiols are oxidised (thiolate / thiyl radical coupling; O2, peroxides, trace metals) to the S-S disulfide.", -5.0);
        const q = m.clone(); q.atoms[sAt].h = 0; addO(q, sAt, 2, 0); addO(q, sAt, 2, 0); addO(q, sAt, 1);
        R(s, "Oxidation", "ox", "sulfonic-acid", "over-oxidation to the sulfonic acid", q, "Strong oxidants take the thiol through sulfenic and sulfinic acids to the sulfonic acid.", -4.5, 0.45);
        break;
      }
      case "amine": case "amine_tert": case "arylamine": {
        const n = s.atoms[0];
        if (m.nb(n).length === 0 && !f.ammonia) break;
        if (s.g === "amine_tert" || (s.g === "arylamine" && f.tertiary && m.nb(n).filter((k) => m.atoms[k].arom).length <= 1)) {
          const p = m.clone(); p.atoms[n].q = 1; addO(p, n, 1, 0, -1);
          R(s, "Oxidation", "ox", "n-oxide", "tertiary amine N-oxidation", p, "Electrophilic oxygen transfer from peroxides / O2 to the nucleophilic nitrogen lone pair gives the polar N-oxide (R3N+-O-).", -4.5);
        }
        if (s.g !== "arylamine" || f.tertiary || f.secondary) {
          const x = pickAlkyl(m, n);
          if (x !== undefined && !f.ammonia) {
            const p = m.clone(); p.unbond(n, x); p.atoms[n].h += 1; p.atoms[x].h -= 1; addO(p, x, 2, 0);
            const prim = H(n) >= 2;
            R(s, "Oxidation", "ox", prim ? "oxidative-deamination" : "n-dealkylation", prim ? "oxidative deamination to the carbonyl compound + NH3" : "oxidative N-dealkylation (carbinolamine collapse)", p, prim ? "Hydrogen abstraction at the alpha C-H, oxygen rebound and collapse of the carbinolamine release the carbonyl compound + ammonia." : "Hydrogen abstraction at the N-alpha C-H, oxygen rebound and collapse of the carbinolamine release the dealkylated amine + aldehyde.", -5.0, prim ? 0.5 : H(x) === 3 ? 0.9 : 0.7);
          }
        }
        if (f.primary || f.secondary) {
          const p = m.clone(); p.atoms[n].h -= 1; addO(p, n, 1);
          R(s, "Oxidation", "ox", "n-hydroxylation", "N-hydroxylation to the hydroxylamine", p, "Two-electron N-oxidation of the amine N-H gives the hydroxylamine (further oxidation gives nitroso / nitrone species).", -2.5, f.ammonia ? 0.4 : 1);
        }
        if (s.g === "arylamine" && f.primary && f.benz) {
          const p = m.clone(); p.atoms[n].h = 0; addO(p, n, 2, 0);
          R(s, "Oxidation", "ox", "nitroso", "aromatic amine oxidation to the nitroso compound", p, "Stepwise oxidation (hydroxylamine, then nitroso) of the aniline nitrogen gives Ar-N=O.", -3.0, 0.8);
          const { mol: q, off } = mergeMol(m, m); q.atoms[n].h = 0; q.atoms[n + off].h = 0; q.bond(n, n + off, 2);
          R(s, "Oxidation", "ox", "azo-dimer", "oxidative coupling to the azo dimer (Ar-N=N-Ar)", q, "One-electron oxidation gives anilino radicals that couple N-N to hydrazo and then azo compounds (coloured impurities).", -6.0, 0.5);
        }
        break;
      }
      case "hetero_azine": {
        const n = s.atoms.find((i) => el(i) === "N" && H(i) === 0 && m.deg(i) === 2 && m.atoms[i].q === 0);
        if (n === undefined) break;
        const p = m.clone(); p.atoms[n].q = 1; addO(p, n, 1, 0, -1);
        R(s, "Oxidation", "ox", "aza-n-oxide", "heteroaromatic N-oxidation", p, "Peroxides / peracids oxidise the pyridine-type ring nitrogen to the N-oxide.", -3.0, f.kind === "azole-N" ? 0.5 : 0.9);
        break;
      }
      case "hemiacetal": {
        const c = s.atoms[0]; const oh = s.atoms.slice(1).find((o) => H(o) > 0);
        if (oh === undefined || H(c) < 1 || f.hydrate) break;
        const p = m.clone(); p.atoms[oh].h = 0; p.bond(c, oh, 2); p.atoms[c].h -= 1;
        R(s, "Oxidation", "ox", "aldose-oxidation", "oxidation of the reducing end to the lactone / aldonic acid", p, "The masked aldehyde (hemiacetal) is oxidised at the anomeric carbon to the lactone (aldonolactone) / aldonic acid.", -8.0);
        break;
      }
      case "aldehyde": {
        const c = s.atoms[0];
        const p = m.clone(); p.atoms[c].h -= 1; addO(p, c, 1);
        R(s, "Oxidation", "ox", "aldehyde-autoxidation", "aldehyde autoxidation to the carboxylic acid", p, "Radical chain autoxidation: the acyl radical adds O2 to give a peracid, which oxidises a second aldehyde; net formation of the carboxylic acid.", -8.5);
        if (m.nb(c).some((k) => el(k) === "C")) {
          const r = m.nb(c).find((k) => el(k) === "C")!;
          const d = m.clone(); d.unbond(c, r); d.remove(s.atoms[1]); d.remove(c); d.atoms[r].h += 1;
          R(s, "Photodegradation", "photo", "photodecarbonylation", "photodecarbonylation (Norrish type I) to the alkane / arene", d, "n->pi* excitation cleaves the acyl C-C bond (Norrish I); loss of CO and H-atom transfer gives the decarbonylated hydrocarbon.", 1.0, f.aryl ? 0.4 : 0.7);
        }
        break;
      }
      case "ketone": {
        const [c, , ...carbs] = s.atoms;
        // Baeyer-Villiger
        const score = (x: number): number => (m.atoms[x].arom ? 3 : isSp3C(m, x) ? (H(x) === 0 ? 4 : H(x) === 1 ? 3 : H(x) === 2 ? 2 : 1) : 2.5);
        const mg = [...carbs].sort((a, b) => score(b) - score(a))[0];
        if (mg !== undefined) {
          const p = m.clone(); p.unbond(c, mg); const o = p.add("O"); p.bond(c, o, 1); p.bond(o, mg, 1);
          R(s, "Oxidation", "ox", "baeyer-villiger", "Baeyer-Villiger oxidation to the ester / lactone", p, "Peroxide/peracid adds to the carbonyl (Criegee intermediate); the higher-aptitude group migrates to oxygen, inserting O between carbonyl and that group.", -7.0, 0.6);
        }
        // Norrish II
        for (const a of carbs) {
          if (!isSp3C(m, a)) continue;
          let done = false;
          for (const b of m.nbx(a, c)) {
            if (!isSp3C(m, b)) continue;
            const g = m.nbx(b, a).find((k) => isSp3C(m, k) && H(k) >= 1);
            if (g === undefined) continue;
            const p = m.clone(); p.unbond(a, b); p.atoms[a].h += 1; p.bond(b, g, 2); p.atoms[g].h -= 1;
            R(s, "Photodegradation", "photo", "norrish-II", "Norrish type II photocleavage (methyl ketone + alkene)", p, "n->pi* excitation, intramolecular gamma-hydrogen abstraction (1,4-biradical) and beta-scission give a shorter ketone (via the enol) + an alkene.", 2.0);
            done = true; break;
          }
          if (done) break;
        }
        break;
      }
      case "enol": {
        if (!f.enediol) break;
        const [c, j] = s.atoms;
        const o1 = m.nb(c).find((k) => el(k) === "O" && H(k) > 0); const o2 = m.nb(j).find((k) => el(k) === "O" && H(k) > 0);
        if (o1 === undefined || o2 === undefined) break;
        const p = m.clone(); p.atoms[o1].h = 0; p.atoms[o2].h = 0; p.bond(c, o1, 2); p.bond(j, o2, 2); p.bond(c, j, 1);
        R(s, "Oxidation", "ox", "enediol-oxidation", "ene-diol oxidation to the 1,2-dicarbonyl (dehydro form)", p, "Two-electron, two-proton oxidation of the ene-diol (via the radical anion; metal ions catalyse) gives the 1,2-dicarbonyl.", -4.0);
        break;
      }
      case "ether": case "aryl_ether": {
        const [o, a, b] = s.atoms;
        const alphas = [a, b].filter((x) => isSp3C(m, x) && H(x) >= 1).sort((x, y) => H(x) - H(y));
        if (!alphas.length) break;
        if (s.g === "ether") {
          const p = m.clone(); p.atoms[alphas[0]].h -= 1; addOOH(p, alphas[0]);
          R(s, "Oxidation", "ox", "ether-hydroperoxide", "ether autoxidation to the alpha-hydroperoxide", p, "Radical autoxidation: H-abstraction at the alpha C-H next to oxygen, O2 capture and chain transfer give the alpha-hydroperoxide (peroxide formation on storage).", -3.0);
        } else {
          const x = alphas.sort((p1, p2) => H(p2) - H(p1))[0];
          const p = m.clone(); p.unbond(o, x); p.atoms[o].h += 1; p.atoms[x].h -= 1; addO(p, x, 2, 0);
          R(s, "Oxidation", "ox", "o-dealkylation", "oxidative O-dealkylation to the phenol + aldehyde", p, "Alpha C-H hydroxylation of the O-alkyl group gives a hemiacetal that collapses to the phenol + aldehyde.", -4.0, 0.7);
        }
        break;
      }
      case "benzylic_ch": {
        const [c] = s.atoms;
        const h = H(c);
        if (h === 1) {
          const p = m.clone(); p.atoms[c].h -= 1; addOOH(p, c);
          R(s, "Oxidation", "ox", "benzylic-hydroperoxide", "benzylic/allylic autoxidation to the hydroperoxide", p, "Radical chain: H-abstraction at the weak benzylic/allylic C-H and O2 capture give the hydroperoxide (e.g. cumene hydroperoxide).", -3.0);
        } else {
          const p = m.clone(); p.atoms[c].h -= 1; addO(p, c, 1);
          R(s, "Oxidation", "ox", "benzylic-hydroxylation", "benzylic/allylic hydroxylation to the alcohol", p, "Autoxidation via the hydroperoxide followed by reduction/decomposition gives the benzylic/allylic alcohol.", -4.0);
          if (h === 2) {
            const q = m.clone(); q.atoms[c].h = 0; addO(q, c, 2, 0);
            R(s, "Oxidation", "ox", "benzylic-oxidation", "benzylic oxidation to the ketone", q, "Further oxidation of the secondary hydroperoxide / alcohol gives the ketone.", -4.5, 0.6);
          }
        }
        break;
      }
      case "hydrazine": {
        const [n1, n2] = s.atoms;
        if (H(n1) < 1 || H(n2) < 1) break;
        const p = m.clone(); p.atoms[n1].h -= 1; p.atoms[n2].h -= 1; p.bond(n1, n2, 2);
        R(s, "Oxidation", "ox", "hydrazine-oxidation", "hydrazine oxidation to the diazene", p, "Two-electron oxidation of the N-N unit gives the diazene (N=N), which can lose N2.", -5.0);
        break;
      }
      case "hydroxylamine": {
        const [n, o] = s.atoms;
        if (H(n) < 1 || H(o) < 1) break;
        const p = m.clone(); p.atoms[n].h -= 1; p.atoms[o].h = 0; p.bond(n, o, 2);
        R(s, "Oxidation", "ox", "hydroxylamine-oxidation", "hydroxylamine oxidation to the nitroso compound", p, "Oxidation of the N-hydroxy-amine removes two hydrogens to give the C-nitroso compound.", -3.0);
        break;
      }
      case "nitroso": {
        const [n] = s.atoms;
        const p = m.clone(); p.atoms[n].q = 1; addO(p, n, 1, 0, -1);
        R(s, "Oxidation", "ox", "nitroso-to-nitro", "nitroso oxidation to the nitro compound", p, "Oxidation of the nitroso nitrogen (peroxide / O2) gives the nitro group.", -8.0);
        break;
      }
      case "thiocarbonyl": {
        const [c, sAt] = s.atoms;
        const p = m.clone(); p.atoms[sAt].el = "O"; p.atoms[sAt].h = 0;
        R(s, "Oxidation", "ox", "oxidative-desulfurisation", "oxidative desulfurisation (C=S to C=O)", p, "S-oxidation to the sulfine / S-oxide is followed by loss of sulfur, converting the thiocarbonyl to the carbonyl.", -5.0);
        void c; break;
      }
      case "phosphine": {
        if (f.oxide) break;
        const p = m.clone(); addO(p, s.atoms[0], 2, 0);
        R(s, "Oxidation", "ox", "p-oxidation", "phosphine oxidation to the phosphine oxide", p, "Air / peroxide oxidises the P(III) lone pair to the strong P=O bond.", -9.0);
        break;
      }
      // ---------------------------------------------------------------- photolysis / thermal
      case "acid": {
        const [c, o, ox] = s.atoms;
        const a = m.nb(c).find((k) => el(k) === "C");
        if (a === undefined) break;
        const p = m.clone(); p.unbond(c, a); p.remove(o); p.remove(ox); p.remove(c); p.atoms[a].h += 1;
        if (f.decarb) {
          const kind = f.decarbKind || "activated";
          R(s, "Thermal Degradation", "therm", "decarboxylation", `decarboxylation (${kind} acid)`, p, kind === "beta" ? "Cyclic six-membered transition state (beta-keto acid) releases CO2 and gives the enol, which tautomerises to the ketone." : kind === "malonic" ? "Malonic-type acids lose CO2 through a cyclic transition state to give the mono-acid." : kind === "aryl" ? "Ipso-protonation of the electron-rich aryl acid (protodecarboxylation) releases CO2." : "The electron-withdrawing alpha-substituent stabilises the carbanion-like transition state, releasing CO2.", kind === "beta" ? -4.0 : kind === "malonic" ? -3.5 : -1.0);
        }
        if (isSp3C(m, a) && m.nb(a).some((k) => m.atoms[k].arom)) {
          R(s, "Photodegradation", "photo", "photodecarboxylation", "photodecarboxylation of the alpha-aryl acid", p, "UV excitation of the aryl chromophore drives electron transfer / decarboxylation, releasing CO2 and a benzylic radical that abstracts H.", 0.5, 0.7);
        }
        cyc.push(s);
        break;
      }
      case "peroxide": {
        const [o1, o2] = s.atoms;
        const term = [o1, o2].find((k) => m.nb(k).length === 1 && H(k) > 0);
        let p: Mol;
        if (term !== undefined) { const other = term === o1 ? o2 : o1; p = m.clone(); p.remove(term); p.atoms[other].h += 1; }
        else { p = m.clone(); p.unbond(o1, o2); p.atoms[o1].h += 1; p.atoms[o2].h += 1; }
        R(s, "Photodegradation", "photo", "peroxide-homolysis", "peroxide O-O homolysis (reduction to alcohol/acid)", p, "UV homolysis of the weak O-O bond gives alkoxy/hydroxyl radicals that abstract hydrogen to give the alcohol (acid for peracids).", -2.0);
        R(s, "Thermal Degradation", "therm", "peroxide-homolysis", "peroxide thermolysis (O-O homolysis)", p, "Thermal O-O homolysis initiates radical chains and gives the alcohol/acid after H-abstraction.", -2.0);
        break;
      }
      default: break;
    }
    // ----- photochemistry of selected groups (not in the switch above to keep hydrolysis grouped)
    if (s.g === "alkene") {
      const [c, j] = s.atoms;
      const p = m.clone(); p.bond(c, j, 1); const o = p.add("O"); p.bond(o, c, 1); p.bond(o, j, 1);
      R(s, "Oxidation", "ox", "epoxidation", "alkene epoxidation", p, "Electrophilic oxygen transfer from peroxide / peracid (or radical addition-cyclisation) converts the C=C bond to the epoxide.", -4.5, f.enolEther || f.enamine ? 0.4 : f.electronPoor ? 0.6 : 1);
      if (f.styrene || f.electronPoor) {
        const d = dimer22(m, c, j);
        R(s, "Photodegradation", "photo", "[2+2]-dimer", "[2+2] photodimerisation to the cyclobutane", d, "Excited-state alkene adds to a ground-state alkene ([2+2] cycloaddition) giving the cyclobutane dimer.", 2.5, 0.7);
      }
    }
    if (s.g === "aryl_halide" || s.g === "alkyl_halide") {
      const [x, y] = s.atoms;
      const p = m.clone(); p.unbond(x, y); p.atoms[x].q = -1; p.atoms[x].h = 0; p.atoms[y].h += 1;
      const dG = s.g === "aryl_halide" ? ({ I: 1.5, Br: 2.2, Cl: 3.0, F: 5.0 } as Record<string, number>)[f.el] ?? 3.0 : 2.0;
      R(s, "Photodegradation", "photo", "photodehalogenation", s.g === "aryl_halide" ? "photodehalogenation (Ar-X to Ar-H)" : "photo-hydrodehalogenation (R-X to R-H)", p, "UV homolysis of the C-X bond (weakest for I, then Br, Cl) gives an aryl/alkyl radical that abstracts hydrogen from the medium.", dG);
    }
    if (s.g === "ester_aryl" && f.aryl) {
      const [c, , ox, x] = s.atoms;
      const r = ringOfAromatic(m, an.rings, x);
      if (r) {
        const pos = r.indexOf(x);
        const targets: [number, string][] = [[r[(pos + 1) % 6], "ortho"], [r[(pos + 5) % 6], "ortho"], [r[(pos + 3) % 6], "para"]];
        const seen = new Set<string>();
        for (const [k, kind] of targets) {
          if (H(k) < 1 || seen.has(kind)) continue;
          seen.add(kind);
          const p = m.clone(); p.unbond(c, ox); p.bond(c, k, 1); p.atoms[k].h -= 1; p.atoms[ox].h += 1;
          R(s, "Photodegradation", "photo", `photo-fries-${kind}`, `photo-Fries rearrangement (${kind}-hydroxyaryl ketone)`, p, "UV homolysis of the aryl ester C(acyl)-O bond gives a radical pair in the solvent cage; recombination at the " + kind + " ring position gives the hydroxyaryl ketone.", 1.5);
        }
      }
    }
    if (s.g === "amide" && f.anilide) {
      const [c, , n] = s.atoms;
      const xa = m.nb(n).find((k) => m.atoms[k].arom);
      const r = xa !== undefined ? ringOfAromatic(m, an.rings, xa) : undefined;
      if (r && xa !== undefined) {
        const pos = r.indexOf(xa);
        for (const [k, kind] of [[r[(pos + 1) % 6], "ortho"], [r[(pos + 5) % 6], "ortho"], [r[(pos + 3) % 6], "para"]] as [number, string][]) {
          if (H(k) < 1) continue;
          const p = m.clone(); p.unbond(c, n); p.bond(c, k, 1); p.atoms[k].h -= 1; p.atoms[n].h += 1;
          R(s, "Photodegradation", "photo", `photo-fries-anilide-${kind}`, `photo-Fries rearrangement of the anilide (${kind}-aminoaryl ketone)`, p, "UV homolysis of the anilide C(acyl)-N bond and radical-cage recombination at the " + kind + " ring position gives the amino-aryl ketone.", 2.0, 0.8);
          break;
        }
      }
    }
    if (s.g === "nitro" && f.aromatic) {
      const [n] = s.atoms;
      const om = m.nb(n).find((k) => el(k) === "O" && m.atoms[k].q === -1);
      if (om !== undefined) {
        const p = m.clone(); p.remove(om); p.atoms[n].q = 0;
        R(s, "Photodegradation", "photo", "nitro-photoreduction", "nitroarene photoreduction to the nitroso compound", p, "Excited nitroarene (n,pi* triplet) abstracts hydrogen / transfers an O atom, giving the nitroso compound (further reduction gives hydroxylamine/amine).", 3.5);
      }
    }
    if (s.g === "michael") {
      const [, a, b] = s.atoms;
      const d = dimer22(m, a, b);
      R(s, "Photodegradation", "photo", "[2+2]-dimer", "[2+2] photodimerisation of the enone to the cyclobutane", d, "Triplet enone adds to a ground-state alkene ([2+2] photocycloaddition) forming the cyclobutane dimer.", 2.5, 0.7);
    }
    if (s.g === "sulfonamide") {
      const [sAt, n] = s.atoms;
      const ca = m.nb(sAt).find((k) => el(k) === "C" && m.atoms[k].arom);
      if (ca !== undefined) {
        const p = m.clone(); p.unbond(sAt, ca); p.unbond(sAt, n); p.bond(ca, n, 1);
        for (const o of m.nb(sAt).filter((k) => el(k) === "O")) p.remove(o);
        p.remove(sAt);
        R(s, "Photodegradation", "photo", "so2-extrusion", "photo-desulfonylation (SO2 extrusion, Ar-N bond formation)", p, "UV cleavage of the S-N / S-C bonds with SO2 extrusion (Smiles-type rearrangement) joins the aryl carbon to nitrogen.", 2.0, 0.8);
      }
    }
    if (s.g === "azide") {
      const nc = s.atoms.find((k) => m.nb(k).some((z) => el(z) === "C"));
      if (nc !== undefined) {
        const p = m.clone(); for (const k of s.atoms) if (k !== nc) p.remove(k); p.atoms[nc].h = 2;
        for (const z of p.nb(nc)) if (p.el(z) === "N") p.unbond(nc, z);
        R(s, "Photodegradation", "photo", "azide-photolysis", "azide photolysis (nitrene) to the amine", p, "UV photolysis extrudes N2 to the nitrene, which abstracts hydrogen to give the primary amine (or inserts / rearranges).", -1.0);
      }
    }
    if (s.g === "nitrosamine") {
      const [n, dO, partner] = s.atoms;
      const p = m.clone(); p.remove(dO); p.remove(n); p.atoms[partner].h += 1;
      R(s, "Photodegradation", "photo", "denitrosation", "photolytic N-N cleavage (denitrosation) to the secondary amine", p, "UV (230-350 nm) cleaves the N-N bond giving an aminyl radical + NO; H-abstraction gives the parent amine.", 2.0);
      R(s, "Acidic Hydrolysis", "acid", "denitrosation", "acid-mediated denitrosation to the secondary amine", p, "Protonated nitrosamine undergoes nucleophile-assisted (Br-/SCN-/thiourea) transnitrosation, releasing the amine + NO+.", 1.0, 0.7);
    }
    if (s.g === "disulfide") {
      const [s1, s2] = s.atoms;
      const p = m.clone(); p.unbond(s1, s2); p.atoms[s1].h += 1; p.atoms[s2].h += 1;
      R(s, "Photodegradation", "photo", "disulfide-homolysis", "disulfide S-S photolysis to the thiols", p, "UV homolysis of the S-S bond gives thiyl radicals that abstract hydrogen to give two thiols (or scramble / recombine).", 3.0);
    }
    if (s.g === "n_oxide") {
      const [n, o] = s.atoms;
      const p = m.clone(); p.remove(o); p.atoms[n].q = 0;
      R(s, "Photodegradation", "photo", "n-oxide-deoxygenation", "N-oxide photo-deoxygenation to the amine", p, "Photolysis of the N-O bond releases atomic oxygen / O(3P) and regenerates the tertiary amine / pyridine (competes with rearrangement).", 2.0);
    }
    if (s.g === "alcohol" || s.g === "alcohol_tert") {
      // thermal dehydration and beta-hydroxy carbonyl (aldol-type) dehydration
      const [x, oh] = s.atoms;
      for (const a of m.nb(x).filter((k) => isSp3C(m, k) && H(k) >= 1)) {
        const cc = m.nbx(a, x).find((k) => isCarbonylC(m, k));
        if (cc === undefined) continue;
        const p = m.clone(); p.unbond(x, oh); p.atoms[oh].h += 1; p.bond(x, a, 2); p.atoms[a].h -= 1;
        R(s, "Thermal Degradation", "therm", "aldol-dehydration", "dehydration of the beta-hydroxy carbonyl to the enone (E1cB)", p, "Enolisation and E1cB loss of water from the aldol-type beta-hydroxy carbonyl gives the conjugated alpha,beta-unsaturated carbonyl.", -1.0, 0.9);
        break;
      }
      if (s.g === "alcohol_tert" || f.benzylic) {
        const p = m.clone();
        if (eliminateAlkene(p, x, oh)) {
          p.atoms[oh].h += 1;
          R(s, "Thermal Degradation", "therm", "dehydration", "acid-/heat-catalysed dehydration to the alkene", p, "E1 dehydration through the tertiary/benzylic carbocation (trace acid, heat) gives the alkene + water.", 0.8, s.g === "alcohol_tert" ? 1 : 0.7);
        }
      }
    }
    if (s.g === "urea" && !f.arom) {
      const [c, , n1, n2] = s.atoms;
      const nh = H(n1) > 0 ? n1 : H(n2) > 0 ? n2 : undefined;
      if (nh !== undefined) {
        const other = nh === n1 ? n2 : n1;
        const p = m.clone(); p.unbond(c, other); p.atoms[other].h += 1; p.atoms[nh].h -= 1; p.bond(c, nh, 2);
        R(s, "Thermal Degradation", "therm", "urea-dissociation", "thermal dissociation of the urea to isocyanate + amine", p, "N-H ureas dissociate reversibly (>130 C) into an isocyanate and an amine.", 3.5);
      }
    }
  }

  // ----- alcohol oxidation (cap: two most relevant carbinols)
  const alcs = an.sites.filter((x) => x.g === "alcohol").sort((a, b) => (b.f.benzylic ? 2 : 0) + (b.f.primary ? 1 : 0) - ((a.f.benzylic ? 2 : 0) + (a.f.primary ? 1 : 0)));
  alcs.slice(0, 2).forEach((s) => {
    const [x, o] = s.atoms;
    const p = m.clone(); p.atoms[o].h = 0; p.bond(x, o, 2); p.atoms[x].h -= 1;
    R(s, "Oxidation", "ox", "alcohol-oxidation", H(x) >= 2 ? "alcohol oxidation to the aldehyde" : "alcohol oxidation to the ketone", p, "Hydride / hydrogen-atom abstraction from the carbinol C-H (peroxide, O2, metal catalysis) gives the carbonyl compound.", H(x) >= 2 ? -3.5 : -4.0);
  });

  // ----- phenols: quinone, C-C coupling
  for (const s of an.sites.filter((x) => x.g === "phenol" && !x.f.phenolate)) {
    const [x] = s.atoms;
    const q = toQuinone(m, x, an.rings);
    if (q) {
      R(s, "Oxidation", "ox", "quinone", "oxidation to the quinone / quinone-imine", q, "Two-electron, two-proton oxidation (or SET then disproportionation) of the hydroquinone/aminophenol-type ring gives the quinone / quinone-imine.", -2.5);
      break;
    }
  }
  const ph = an.sites.find((x) => x.g === "phenol" && !x.f.phenolate && !x.f.hindered && !x.f.hydroquinone);
  if (ph) {
    const r = ringOfAromatic(m, an.rings, ph.atoms[0]);
    if (r) {
      const pos = r.indexOf(ph.atoms[0]);
      const k = [r[(pos + 1) % 6], r[(pos + 5) % 6]].find((z) => H(z) >= 1);
      if (k !== undefined) {
        const { mol: p, off } = mergeMol(m, m); p.bond(k, k + off, 1); p.atoms[k].h -= 1; p.atoms[k + off].h -= 1;
        R(ph, "Oxidation", "ox", "phenol-coupling", "oxidative ortho-ortho C-C coupling to the biaryl-diol", p, "Phenoxyl radicals (from SET / H-abstraction) couple at the ortho carbons; tautomerisation restores aromaticity to give the 2,2'-biaryl diol.", -3.0, 0.5);
      }
    }
  }

  // ----- aromatic hydroxylation (best free C-H over all benzenoid rings)
  let bestAr: { s: Site; atom: number; score: number } | null = null;
  for (const s of an.sites.filter((x) => x.g === "aromatic")) {
    const b = bestArylCH(m, s.atoms);
    if (b && (!bestAr || b.score > bestAr.score)) bestAr = { s, atom: b.atom, score: b.score };
  }
  if (bestAr) {
    const p = m.clone(); p.atoms[bestAr.atom].h -= 1; addO(p, bestAr.atom, 1);
    R(bestAr.s, "Oxidation", "ox", "aromatic-hydroxylation", "aromatic hydroxylation by hydroxyl radical / peroxide (phenol formation)", p, "Electrophilic HO. / peroxide attack at the most electron-rich free ring C-H (ortho/para to donors, away from electron-withdrawing groups) gives the phenol.", -2.0, bestAr.s.f.activated ? 0.9 : bestAr.s.f.deactivated ? 0.3 : 0.6);
  }

  // ----- saturated hydrocarbon fallback (autoxidation of the weakest C-H)
  if (an.sites.length === 1 && an.sites[0].g === "hydrocarbon") {
    const cs = m.alive().filter((i) => isSp3C(m, i) && H(i) >= 1).sort((a, b) => H(a) - H(b));
    if (cs.length) {
      const p = m.clone(); p.atoms[cs[0]].h -= 1; addO(p, cs[0], 1);
      R(an.sites[0], "Oxidation", "ox", "alkane-autoxidation", "C-H autoxidation to the alcohol (tertiary > secondary > primary)", p, "Radical chain autoxidation abstracts the weakest C-H; O2 capture and decomposition of the hydroperoxide give the alcohol.", -3.0, 0.7);
    }
  }

  // ----- intramolecular cyclisation (lactam / lactone formation), thermal
  const nuSites = an.sites.filter((x) => x.g === "alcohol" || x.g === "alcohol_tert" || (x.g === "amine" && (x.f.primary || x.f.secondary) && !x.f.ammonia));
  let cycCount = 0;
  for (const e of cyc) {
    const c = e.atoms[0]; const o = e.atoms[1]; const ox = e.atoms[2];
    for (const nsite of nuSites) {
      if (cycCount >= 3) break;
      const nu = nsite.g === "amine" ? nsite.atoms[0] : nsite.atoms[1];
      if (nu === ox || nu === o) continue;
      const d = shortestPath(m, c, nu, new Set([o, ox]));
      if (d !== 4 && d !== 5) continue;
      if (H(nu) < 1) continue;
      const p = m.clone(); p.unbond(c, ox);
      if (e.g === "acid") p.remove(ox); else p.atoms[ox].h += 1;
      p.bond(c, nu, 1); p.atoms[nu].h -= 1;
      const isN = nsite.g === "amine";
      out.push({ cond: "Thermal Degradation", key: "cyclisation", label: isN ? "intramolecular lactamisation" : "intramolecular lactonisation", mol: p, ck: "therm",
        mech: `Intramolecular nucleophilic acyl substitution: the ${isN ? "amine" : "hydroxyl"} attacks the ${e.g === "acid" ? "carboxylic acid" : "ester"} carbonyl through a ${d + 1}-membered transition state, closing the ${d + 1}-membered ${isN ? "lactam" : "lactone"} and expelling ${e.g === "acid" ? "water" : "the alcohol"} (heat / acid catalysis).`,
        dG: isN ? -2.5 : 0.5, vuln: isN ? HI : MO, mult: 1 });
      cycCount++;
    }
  }

  // ----- decarboxylation-independent thermal: nothing else
  return out;
}

/* ------------------------------------------------------------------ */
/*  Co-reactant (secondary compound) interaction templates             */
/* ------------------------------------------------------------------ */
const XSRC: Source = "Interaction with other compound";

function shiftSites(sites: Site[], off: number): Site[] {
  return sites.map((s) => ({ g: s.g, atoms: s.atoms.map((a) => a + off), f: { ...s.f, ...(s.f.exo ? { exo: (s.f.exo as number[]).map((x) => x + off) } : {}) } }));
}
interface Nu { at: number; kind: "N" | "O" | "S"; w: number; site: Site; }
function nucleophilesOf(M: Mol, sites: Site[]): Nu[] {
  const out: Nu[] = []; const h = (i: number): number => M.atoms[i].h;
  for (const s of sites) {
    switch (s.g) {
      case "amine": if (h(s.atoms[0]) >= 1) out.push({ at: s.atoms[0], kind: "N", w: 1, site: s }); break;
      case "arylamine": if (h(s.atoms[0]) >= 1) out.push({ at: s.atoms[0], kind: "N", w: 0.6, site: s }); break;
      case "hydrazine": { const [a, b] = s.atoms; const n = h(a) >= h(b) ? a : b; if (h(n) >= 1) out.push({ at: n, kind: "N", w: 1, site: s }); break; }
      case "hydroxylamine": if (h(s.atoms[0]) >= 1) out.push({ at: s.atoms[0], kind: "N", w: 0.8, site: s }); break;
      case "alcohol": out.push({ at: s.atoms[1], kind: "O", w: 0.7, site: s }); break;
      case "alcohol_tert": out.push({ at: s.atoms[1], kind: "O", w: 0.3, site: s }); break;
      case "phenol": if (!s.f.phenolate) out.push({ at: s.atoms[1], kind: "O", w: 0.4, site: s }); break;
      case "thiol": if (h(s.atoms[0]) >= 1) out.push({ at: s.atoms[0], kind: "S", w: 1, site: s }); break;
      default: break;
    }
  }
  return out;
}
const secVuln = (s: Site): Vuln => tweakProfile(s.g, s.f).sec[0];
const minVuln = (a: Vuln, b: Vuln): Vuln => (VRANK(a) <= VRANK(b) ? a : b);

/** every covalent / ionic pathway in which A (site list `ea`) acts as electrophile/acid and B (`nb`) as nucleophile/base */
function crossOriented(M: Mol, ea: Site[], eb: Site[], nameE: string, nameN: string): Rx[] {
  const out: Rx[] = [];
  const H = (i: number): number => M.atoms[i].h;
  const nus = nucleophilesOf(M, eb);
  const X = (E: Site, N: Site | null, cond: Cond, key: string, label: string, mol: Mol | null, mech: string, dG: number, mult = 1, vuln?: Vuln, physical = false): void => {
    if (!mol) return;
    const v = vuln ?? (N ? minVuln(secVuln(E), secVuln(N)) : secVuln(E));
    out.push({ cond, key, label, mol, mech, dG, ck: "acid", vuln: v, mult, source: XSRC, co: nameN, physical, mergedParent: true });
  };
  const nuText: Record<string, string> = { N: "amine nitrogen", O: "hydroxyl oxygen", S: "thiol sulfur" };

  for (const E of ea) {
    const f = E.f;
    // ------------------------------------------------ acyl transfer to nucleophiles
    const acyl = ((): { c: number; lg: number; kind: string; dN: number; dO: number; w: number } | null => {
      switch (E.g) {
        case "ester": return { c: E.atoms[0], lg: E.atoms[2], kind: "ester", dN: -3.5, dO: 0.4, w: 0.9 };
        case "ester_aryl": return { c: E.atoms[0], lg: E.atoms[2], kind: "aryl/vinyl ester", dN: -5.0, dO: -2.0, w: 1 };
        case "lactone": return { c: E.atoms[0], lg: E.atoms[2], kind: "lactone", dN: -3.0, dO: 0.4, w: 0.9 };
        case "anhydride": return { c: E.atoms[0], lg: E.atoms[2], kind: "anhydride", dN: -9, dO: -8, w: 1 };
        case "acyl_halide": return { c: E.atoms[0], lg: E.atoms[2], kind: "acyl halide", dN: -9, dO: -9, w: 1 };
        case "thioester": return { c: E.atoms[0], lg: E.atoms[2], kind: "thioester", dN: -4, dO: -2, w: 0.9 };
        case "carbonate": return { c: E.atoms[0], lg: E.atoms[2], kind: "carbonate", dN: -6, dO: -3, w: 0.9 };
        case "beta_lactam": return { c: E.atoms[0], lg: E.atoms[2], kind: "beta-lactam", dN: -9, dO: -8, w: 1 };
        default: return null;
      }
    })();
    if (acyl) {
      let nN = 0, nO = 0;
      for (const N of nus) {
        if (N.kind === "S") continue;
        if (N.kind === "N" ? nN >= 2 : nO >= 2) continue;
        const p = M.clone(); p.unbond(acyl.c, acyl.lg); p.atoms[acyl.lg].h += 1; p.bond(acyl.c, N.at, 1); p.atoms[N.at].h -= 1;
        tautomerizeEnols(p);
        if (N.kind === "N") {
          nN++;
          X(E, N.site, "Thermal Degradation", "aminolysis", `aminolysis: ${acyl.kind} + amine to amide`, p,
            `Nucleophilic addition-elimination: the ${nuText.N} of ${nameN} attacks the electrophilic ${acyl.kind} carbonyl of ${nameE}, expels the leaving group and forms a covalent amide adduct.`, acyl.dN, acyl.w * N.w);
        } else {
          nO++;
          X(E, N.site, "Thermal Degradation", "transesterification", `acyl transfer: ${acyl.kind} + alcohol/phenol (transesterification)`, p,
            `Intermolecular acyl transfer: the ${nuText.O} of ${nameN} attacks the ${acyl.kind} carbonyl of ${nameE} (accelerated by heat, moisture and trace acid/base), giving a new ester and the released alcohol.`, acyl.dO, acyl.w * N.w * 0.8);
        }
      }
      for (const N of nus.filter((x) => x.kind === "S").slice(0, 1)) {
        const p = M.clone(); p.unbond(acyl.c, acyl.lg); p.atoms[acyl.lg].h += 1; p.bond(acyl.c, N.at, 1); p.atoms[N.at].h -= 1;
        X(E, N.site, "Thermal Degradation", "thioester-formation", `acyl transfer: ${acyl.kind} + thiol (thioester formation)`, p, `Thiolate/thiol of ${nameN} attacks the ${acyl.kind} carbonyl of ${nameE} giving a thioester.`, -1.5, acyl.w * 0.7);
      }
    }
    if (E.g === "isocyanate") {
      const [c, n] = E.atoms;
      for (const N of nus.slice(0, 3)) {
        const p = M.clone(); p.bond(c, N.at, 1); p.atoms[N.at].h -= 1; p.atoms[n].h += 1;
        p.bond(c, n, 1);
        X(E, N.site, "Thermal Degradation", "isocyanate-addition", N.kind === "N" ? "isocyanate + amine to urea" : N.kind === "O" ? "isocyanate + alcohol to carbamate" : "isocyanate + thiol to thiocarbamate", p,
          `Nucleophilic addition of the ${nuText[N.kind]} of ${nameN} to the cumulated N=C=O carbon of ${nameE} gives the ${N.kind === "N" ? "urea" : N.kind === "O" ? "carbamate" : "thiocarbamate"} (very fast).`, -9, N.kind === "N" ? 1 : 0.85);
      }
    }
    // ------------------------------------------------ carboxylic acid + alcohol / amine (esterification / amidation)
    if (E.g === "acid") {
      const [c, , ox] = E.atoms;
      for (const N of nus.filter((x) => x.kind !== "S").slice(0, 3)) {
        const p = M.clone(); p.unbond(c, ox); p.remove(ox); p.bond(c, N.at, 1); p.atoms[N.at].h -= 1;
        if (N.kind === "O") X(E, N.site, "Thermal Degradation", "esterification", "Fischer esterification (acid + alcohol)", p, `Acid-catalysed condensation of the carboxylic acid of ${nameE} with the hydroxyl of ${nameN} gives an ester + water (equilibrium; driven by heat / water loss).`, 0.8, 0.6 * N.w, MO);
        else X(E, N.site, "Thermal Degradation", "amidation", "thermal amidation (acid + amine)", p, `Heating the carboxylic acid of ${nameE} with the amine of ${nameN} (via the ammonium carboxylate salt) dehydrates to the amide.`, -1.0, 0.6 * N.w, MO);
      }
    }
    // ------------------------------------------------ carbonyl + amine: imine / N-glycosylamine (Maillard)
    if (E.g === "aldehyde" || E.g === "ketone") {
      const [c, o] = E.atoms;
      for (const N of nus.filter((x) => x.kind === "N" && H(x.at) >= 2).slice(0, 2)) {
        const p = M.clone(); p.remove(o); p.bond(c, N.at, 2); p.atoms[N.at].h -= 2;
        const nm = N.site.g === "hydrazine" ? "hydrazone" : N.site.g === "hydroxylamine" ? "oxime" : "imine (Schiff base)";
        X(E, N.site, "Thermal Degradation", "condensation", `${E.g} + primary amine condensation to the ${nm}`, p, `Nucleophilic addition of the ${nuText.N} of ${nameN} to the ${E.g} carbonyl of ${nameE} and dehydration gives the ${nm} (reversible; water removal drives it).`, E.g === "aldehyde" ? 1.5 : 3.0, 1, E.g === "ketone" ? MO : undefined);
      }
    }
    if (E.g === "hemiacetal" && !f.hydrate) {
      const c = E.atoms[0]; const oh = E.atoms.slice(1).find((o) => H(o) > 0);
      if (oh !== undefined) {
        for (const N of nus.filter((x) => x.kind === "N").slice(0, 2)) {
          const p = M.clone(); p.remove(oh); p.bond(c, N.at, 1); p.atoms[N.at].h -= 1;
          X(E, N.site, "Thermal Degradation", "maillard", "Maillard condensation: reducing sugar + amine to the N-glycosylamine", p, `The masked aldehyde (anomeric hemiacetal) of ${nameE} condenses with the unprotonated amine of ${nameN} to the N-glycosylamine / Schiff base; the irreversible Amadori rearrangement and browning (melanoidins) that follow on heating / moisture pull the equilibrium forward.`, -3.0);
        }
      }
    }
    // ------------------------------------------------ epoxide ring opening
    if (E.g === "epoxide") {
      const [o, , b] = E.atoms;
      for (const N of nus.slice(0, 3)) {
        const p = M.clone(); p.unbond(o, b); p.atoms[o].h = 1; p.bond(b, N.at, 1); p.atoms[N.at].h -= 1;
        X(E, N.site, "Thermal Degradation", "epoxide-opening", `epoxide ring opening by ${N.kind === "N" ? "amine" : N.kind === "O" ? "alcohol" : "thiol"}`, p, `SN2 opening of the strained epoxide of ${nameE} by the ${nuText[N.kind]} of ${nameN} gives the beta-substituted alcohol.`, -8, N.kind === "O" ? 0.5 : 1);
      }
    }
    // ------------------------------------------------ alkylation
    if (E.g === "alkyl_halide" && f.el !== "F") {
      const [x, y] = E.atoms;
      for (const N of nus.slice(0, 3)) {
        const p = M.clone(); p.unbond(x, y); p.atoms[x].h += 1; p.bond(y, N.at, 1); p.atoms[N.at].h -= 1;
        X(E, N.site, "Thermal Degradation", "alkylation", `N/O/S-alkylation by the alkyl halide (${N.kind === "N" ? "amine" : N.kind === "O" ? "alcohol" : "thiol"})`, p, `SN2 displacement of halide from ${nameE} by the ${nuText[N.kind]} of ${nameN} gives the alkylated product + HX (alkyl halides are genotoxic-type alkylators).`, -6, (f.tertiary ? 0.5 : 1) * (N.kind === "O" ? 0.4 : 1));
      }
      for (const s2 of eb.filter((z) => z.g === "amine_tert").slice(0, 1)) {
        const n = s2.atoms[0];
        const p = M.clone(); p.unbond(x, y); p.atoms[x].q = -1; p.atoms[x].h = 0; p.bond(y, n, 1); p.atoms[n].q = 1;
        X(E, s2, "Thermal Degradation", "quaternisation", "quaternisation of the tertiary amine (Menshutkin)", p, `Menshutkin reaction: the tertiary amine of ${nameN} displaces halide at ${nameE}, giving a quaternary ammonium salt.`, -6, 0.8);
      }
    }
    if (E.g === "sulfonate_ester") {
      const sAt = E.atoms[0];
      const oo = M.nb(sAt).find((k) => M.el(k) === "O" && M.order(sAt, k) === 1 && M.nbx(k, sAt).some((z) => M.el(z) === "C" && isSp3C(M, z)));
      if (oo !== undefined) {
        const r = M.nbx(oo, sAt).find((z) => isSp3C(M, z))!;
        for (const N of nus.slice(0, 3)) {
          const p = M.clone(); p.unbond(oo, r); p.atoms[oo].h += 1; p.bond(r, N.at, 1); p.atoms[N.at].h -= 1;
          X(E, N.site, "Thermal Degradation", "alkylation-sulfonate", `alkylation by the sulfonate/sulfate ester (${N.kind === "N" ? "N" : N.kind === "O" ? "O" : "S"}-alkylation)`, p, `${nameE} is an alkylating agent: SN2 attack by the ${nuText[N.kind]} of ${nameN} on the alkyl carbon releases the sulfonic acid and gives the alkylated adduct.`, -7, N.kind === "O" ? 0.4 : 1);
        }
      }
    }
    // ------------------------------------------------ Michael addition
    if (E.g === "michael") {
      const [, a, b] = E.atoms;
      for (const N of nus.filter((z) => z.kind !== "O").slice(0, 3)) {
        const p = M.clone(); p.bond(a, b, 1); p.bond(b, N.at, 1); p.atoms[a].h += 1; p.atoms[N.at].h -= 1;
        X(E, N.site, "Thermal Degradation", "michael-addition", `${N.kind === "S" ? "thia" : "aza"}-Michael addition to the alpha,beta-unsaturated carbonyl`, p, `Conjugate addition of the ${nuText[N.kind]} of ${nameN} to the electrophilic beta-carbon of ${nameE} gives the beta-substituted carbonyl adduct.`, -5, 1);
      }
    }
  }
  return out;
}

/** interactions where B supplies a NITROSATING agent (nitrite) and A supplies amines */
function nitrosation(M: Mol, aSites: Site[], bSites: Site[], nameA: string, nameB: string, acidPresent: boolean): Rx[] {
  const out: Rx[] = [];
  const nit = bSites.find((s) => s.g === "nitrite" && !s.f.ester);
  if (!nit) return out;
  const nN = nit.atoms[0];
  const drop = (p: Mol): void => { for (const o of p.nb(nN)) p.remove(o); p.remove(nN); };
  const acidMult = acidPresent ? 1 : 0.55;
  const push = (cond: Cond, key: string, label: string, mol: Mol, mech: string, dG: number, mult: number, vuln: Vuln): void => {
    out.push({ cond, key, label, mol, mech, dG, ck: "acid", vuln, mult, source: XSRC, co: nameB, mergedParent: true });
  };
  for (const s of aSites) {
    const n = s.atoms[0];
    if ((s.g === "amine" && s.f.secondary) || (s.g === "arylamine" && s.f.secondary)) {
      const p = M.clone(); p.atoms[n].h -= 1; const nn = p.add("N"); const on = p.add("O"); p.bond(n, nn, 1); p.bond(nn, on, 2); drop(p);
      push("Acidic Hydrolysis", "n-nitrosation", "N-nitrosation of the secondary amine (nitrosamine formation)", p, `Under acidic conditions nitrite forms HONO/N2O3/NO+, which nitrosates the secondary amine of ${nameA} to the N-nitrosamine (a mutagenic impurity class).`, -8.5, acidMult, CR);
    } else if (s.g === "amine_tert" || (s.g === "arylamine" && s.f.tertiary)) {
      const x = pickAlkyl(M, n);
      if (x !== undefined) {
        const p = M.clone(); p.unbond(n, x); p.atoms[x].h -= 1; addO(p, x, 2, 0); const nn = p.add("N"); const on = p.add("O"); p.bond(n, nn, 1); p.bond(nn, on, 2); drop(p);
        push("Acidic Hydrolysis", "nitrosative-dealkylation", "nitrosative dealkylation of the tertiary amine to the N-nitrosamine", p, `Nitrosation of the tertiary amine of ${nameA} forms an unstable N-nitrosammonium species that loses an aldehyde (dealkylation), leaving the secondary N-nitrosamine.`, -5, 0.4 * acidMult, MO);
      }
    } else if (s.g === "arylamine" && s.f.primary) {
      const p = M.clone(); p.atoms[n].h = 0; p.atoms[n].q = 1; const nn = p.add("N"); p.bond(n, nn, 3); drop(p);
      push("Acidic Hydrolysis", "diazotisation", "diazotisation of the primary aromatic amine (aryl diazonium ion)", p, `In acid, nitrite of ${nameB} diazotises the primary aromatic amine of ${nameA} to the aryl diazonium ion (which couples, hydrolyses to phenols, or loses N2).`, -6, acidMult, HI);
    }
  }
  return out;
}

/** proton-transfer salts (physical/ionic interactions): acid donor + basic acceptor, or stronger acid displacing a weaker carboxylate */
function saltFormation(M: Mol, dSites: Site[], aSites: Site[], nameD: string, nameB: string, dIsPrimary: boolean): Rx[] {
  const out: Rx[] = [];
  const H = (i: number): number => M.atoms[i].h;
  const donors: { at: number; strong: boolean; site: Site }[] = [];
  for (const s of dSites) {
    if (s.g === "acid") donors.push({ at: s.atoms[2], strong: false, site: s });
    else if (s.g === "sulfonic" && s.f.acid) donors.push({ at: s.atoms[1], strong: true, site: s });
    else if (s.g === "hydrogen_halide") donors.push({ at: s.atoms[0], strong: true, site: s });
  }
  const accs: { at: number; w: number; strongOnly: boolean; site: Site; carbox?: boolean }[] = [];
  for (const s of aSites) {
    if (s.g === "amine" || s.g === "amine_tert") accs.push({ at: s.atoms[0], w: 1, strongOnly: false, site: s });
    else if (s.g === "amidine") accs.push({ at: s.atoms[1], w: 1, strongOnly: false, site: s });
    else if (s.g === "arylamine") accs.push({ at: s.atoms[0], w: 0.4, strongOnly: true, site: s });
    else if (s.g === "hetero_azine") { const n = s.atoms.find((i) => M.el(i) === "N" && H(i) === 0 && M.deg(i) === 2 && M.atoms[i].q === 0); if (n !== undefined) accs.push({ at: n, w: 0.5, strongOnly: true, site: s }); }
    else if (s.g === "carboxylate") accs.push({ at: s.atoms[2], w: 0.9, strongOnly: false, site: s, carbox: true });
  }
  let count = 0;
  for (const d of donors.slice(0, 2)) for (const a of accs.slice(0, 2)) {
    if (a.strongOnly && !d.strong) continue;
    if (a.carbox) {
      const cAt = a.site.atoms[0];
      const accAryl = M.nb(cAt).some((k) => M.atoms[k].arom);
      const dAcidic = d.strong || d.site.f.aryl || d.site.f.decarbKind === "alpha-EWG";
      if (!dAcidic || (accAryl && !d.strong)) continue;
    }
    if (count++ >= 3) break;
    const p = M.clone();
    p.atoms[d.at].h -= 1; p.atoms[d.at].q = -1;
    if (a.carbox) { p.atoms[a.at].q = 0; p.atoms[a.at].h += 1; } else { p.atoms[a.at].h += 1; p.atoms[a.at].q = 1; }
    const donorName = dIsPrimary ? nameD : nameB, accName = dIsPrimary ? nameB : nameD;
    out.push({ cond: dIsPrimary ? "Basic Hydrolysis" : "Acidic Hydrolysis", key: a.carbox ? "acid-base-exchange" : "proton-transfer-salt",
      label: a.carbox ? "acid-base exchange (stronger acid displaces the weaker carboxylate, forming its salt)" : "acid-base salt formation (ionic pair, proton transfer)", mol: p, ck: "acid",
      mech: a.carbox
        ? `${accName} is the salt of a weak acid (a mild base); the stronger ${d.strong ? "acid" : "aromatic/activated carboxylic acid"} of ${donorName} protonates the carboxylate, liberating the free weaker acid and forming the salt of ${donorName} with the metal cation (alters solubility and raises the local pH).`
        : `Proton transfer from the ${d.strong ? "strong" : "carboxylic"} acid of ${donorName} to the basic nitrogen of ${accName} gives an ammonium-anion salt (changes solubility, hygroscopicity and the pH of the microenvironment; may precipitate).`,
      dG: d.strong ? -6.5 : a.carbox ? -3 : -4.5, vuln: d.strong ? CR : HI, mult: a.w, source: XSRC, co: nameB, physical: true, mergedParent: true });
  }
  return out;
}

/** acid + basic metal species (oxide / hydroxide / carbonate) => metal carboxylate (soap) salt */
function metalSalt(A: Analysis, B: Analysis, nameA: string, nameB: string): Rx[] {
  const out: Rx[] = [];
  const baseSite = B.sites.find((s) => s.g === "strong_base" || s.g === "carbonate_ion");
  const metal = B.sites.find((s) => s.g === "metal" || s.g === "metal_transition");
  if (!baseSite || !metal) return out;
  const z = Math.max(1, Math.abs(B.mol.atoms[metal.atoms[0]].q));
  const acids = A.sites.filter((s) => s.g === "acid").slice(0, 1);
  for (const s of acids) {
    const ox = s.atoms[2];
    let p = A.mol.clone(); p.atoms[ox].h -= 1; p.atoms[ox].q = -1;
    for (let k = 1; k < z && k < 3; k++) { const { mol, off } = mergeMol(p, A.mol); mol.atoms[ox + off].h -= 1; mol.atoms[ox + off].q = -1; p = mol; }
    const mAt = p.add(B.mol.el(metal.atoms[0]), { q: z });
    void mAt;
    out.push({ cond: "Basic Hydrolysis", key: "metal-carboxylate", label: `metal carboxylate (soap) formation with ${B.mol.el(metal.atoms[0])}${z}+`, mol: p, ck: "base",
      mech: `The basic ${nameB} (oxide / hydroxide / carbonate) neutralises the carboxylic acid of ${nameA} stoichiometrically, forming the ${B.mol.el(metal.atoms[0])}(${z}+) carboxylate salt (ionic; Lewis-acid cation coordination; often insoluble).`,
      dG: -6, vuln: CR, mult: 1, source: XSRC, co: nameB, physical: true, mergedParent: false });
  }
  return out;
}

/** reductive interactions: B = reductant (hydride / thiol / sulfite / hydrazine / enediol) acting on reducible groups of A */
function reductions(A: Analysis, B: Analysis, nameA: string, nameB: string): Rx[] {
  const out: Rx[] = [];
  const hydride = B.sites.some((s) => s.g === "hydride");
  const strongRed = hydride || B.sites.some((s) => s.g === "hydrazine" || s.g === "sulfite");
  const push = (key: string, label: string, mol: Mol, mech: string, dG: number, mult: number, vuln: Vuln): void => {
    out.push({ cond: "Oxidation", key, label, mol, mech, dG, ck: "ox", vuln, mult, source: XSRC, co: nameB, mergedParent: false });
  };
  const m = A.mol;
  for (const s of A.sites) {
    if ((s.g === "aldehyde" || s.g === "ketone") && hydride) {
      const [c, o] = s.atoms; const p = m.clone(); p.bond(c, o, 1); p.atoms[o].h = 1; p.atoms[c].h += 1;
      push("hydride-reduction", `reduction of the ${s.g} to the alcohol by ${nameB}`, p, `Hydride transfer from ${nameB} to the carbonyl carbon of ${nameA}, then protonation, gives the alcohol (redox interaction).`, -8, 1, CR);
    }
    if (s.g === "quinone") {
      const p = m.clone(); const ring = s.atoms; ring.forEach((i) => { p.atoms[i].arom = true; });
      for (let i = 0; i < ring.length; i++) p.bond(ring[i], ring[(i + 1) % ring.length], 1.5);
      for (const c of (s.f.exo as number[])) {
        const ex = p.nb(c).find((k) => !ring.includes(k) && p.order(c, k) === 2);
        if (ex !== undefined) { p.bond(c, ex, 1); p.atoms[ex].h += 1; }
      }
      push("quinone-reduction", `reduction of the quinone to the hydroquinone by ${nameB}`, p, `Two-electron, two-proton reduction of the quinone/quinone-imine of ${nameA} by ${nameB} gives the aromatic hydroquinone / aminophenol.`, -6, strongRed ? 1 : 0.8, HI);
    }
    if (s.g === "disulfide") {
      const [s1, s2] = s.atoms; const p = m.clone(); p.unbond(s1, s2); p.atoms[s1].h += 1; p.atoms[s2].h += 1;
      push("disulfide-reduction", `disulfide reduction to thiols by ${nameB}`, p, `Thiol-disulfide exchange / hydride reduction by ${nameB} cleaves the S-S bond of ${nameA} to two thiols.`, -4, strongRed || B.sites.some((z) => z.g === "thiol") ? 1 : 0.6, HI);
    }
    if (s.g === "nitro" && s.f.aromatic && strongRed) {
      const n = s.atoms[0]; const p = m.clone(); for (const o of m.nb(n).filter((k) => m.el(k) === "O")) p.remove(o); p.atoms[n].q = 0; p.atoms[n].h = 2;
      push("nitro-reduction", `nitro group reduction to the amine by ${nameB}`, p, `${nameB} reduces the aromatic nitro group of ${nameA} through nitroso and hydroxylamine to the aniline (usually needs a catalyst / forcing conditions).`, -9, 0.6, MO);
    }
    if (s.g === "peroxide" && strongRed) {
      const [o1, o2] = s.atoms; const term = [o1, o2].find((k) => m.nb(k).length === 1 && m.atoms[k].h > 0);
      if (term !== undefined) { const other = term === o1 ? o2 : o1; const p = m.clone(); p.remove(term); p.atoms[other].h += 1;
        push("peroxide-reduction", `hydroperoxide reduction to the alcohol by ${nameB}`, p, `${nameB} reduces the O-O bond of the hydroperoxide of ${nameA} to the alcohol.`, -12, 1, CR); }
    }
    if (s.g === "azide" && strongRed) {
      const nc = s.atoms.find((k) => m.nb(k).some((z) => m.el(z) === "C"));
      if (nc !== undefined) { const p = m.clone(); for (const k of s.atoms) if (k !== nc) p.remove(k); p.atoms[nc].h = 2; for (const z of p.nb(nc)) if (p.el(z) === "N") p.unbond(nc, z);
        push("azide-reduction", `azide reduction to the amine by ${nameB}`, p, `${nameB} reduces the azide of ${nameA} to the primary amine with N2 loss.`, -10, 0.8, HI); }
    }
  }
  return out;
}

/**
 * Environment effects: copies of the primary compound's own stress pathways, re-graded because the
 * co-reactant changes the microenvironment (alkaline / acidic / moist / oxidising / photocatalytic).
 */
function environmentBoost(stress: Rx[], B: Analysis, nameB: string): Rx[] {
  const out: Rx[] = []; const cl = B.classes;
  const strongAcid = cl.has("acid_strong");
  const photocat = B.sites.some((s) => s.g === "metal_oxide" && s.f.photocat);
  const push = (r: Rx, step: number, why: string, mult = 1): void => {
    out.push({ ...r, vuln: stepVuln(r.vuln, step), mult: (r.mult ?? 1) * mult, source: XSRC, co: nameB, mergedParent: false,
      mech: `${why} ${r.mech}` });
  };
  for (const r of stress) {
    if (r.cond === "Basic Hydrolysis" && cl.has("base_strong")) push(r, 1, `Alkaline microenvironment created by ${nameB} (local pH >> 9) accelerates this pathway.`);
    else if (r.cond === "Basic Hydrolysis" && cl.has("base_weak") && cl.has("metal")) push(r, 1, `${nameB} is the metal salt of a weak acid: dissolving in the moisture film it raises the local pH (mildly alkaline microenvironment) and its cation acts as a Lewis acid, accelerating base-catalysed cleavage.`);
    else if (r.cond === "Acidic Hydrolysis" && strongAcid) push(r, 1, `Acidic microenvironment created by ${nameB} accelerates this pathway.`);
    else if (r.cond === "Hydrolysis" && cl.has("water")) push(r, 1, `Water / moisture supplied by ${nameB} drives this hydrolysis.`);
    else if (r.cond === "Oxidation" && cl.has("oxidant")) push(r, 1, `${nameB} is an oxidant / contains peroxide and supplies the oxidising equivalents.`);
    else if (r.cond === "Oxidation" && cl.has("peroxide_former")) push(r, 1, `Autoxidation of ${nameB} generates trace hydroperoxides / radicals that oxidise susceptible groups.`);
    else if (r.cond === "Oxidation" && cl.has("redox_metal")) push(r, 1, `Redox-active metal ions in ${nameB} catalyse radical (Fenton-type) oxidation.`);
    else if (r.cond === "Photodegradation" && photocat) push(r, 1, `${nameB} (TiO2/ZnO-type semiconductor) is a photocatalyst that generates reactive oxygen species under UV.`);
  }
  return out;
}

function genCross(A: Analysis, B: Analysis, nameA: string, nameB: string, stress: Rx[], ambientAcid = false): Rx[] {
  const { mol: M, off } = mergeMol(A.mol, B.mol);
  const sa = A.sites; const sb = shiftSites(B.sites, off);
  const out: Rx[] = [];
  out.push(...crossOriented(M, sa, sb, nameA, nameB));
  out.push(...crossOriented(M, sb, sa, nameB, nameA));
  const acidPresent = ambientAcid || A.classes.has("acid_weak") || A.classes.has("acid_strong") || B.classes.has("acid_weak") || B.classes.has("acid_strong");
  out.push(...nitrosation(M, sa, sb, nameA, nameB, acidPresent));
  out.push(...nitrosation(M, sb, sa, nameB, nameA, acidPresent));
  out.push(...saltFormation(M, sa, sb, nameA, nameB, true));
  out.push(...saltFormation(M, sb, sa, nameB, nameA, false));
  out.push(...metalSalt(A, B, nameA, nameB));
  out.push(...reductions(A, B, nameA, nameB));
  out.push(...environmentBoost(stress, B, nameB));
  out.forEach((r) => { r.co = nameB; });   // "co" is always the co-reactant (never the primary compound)
  return out;
}

/* ------------------------------------------------------------------ */
/*  Public API: input resolution, ranking, reporting                   */
/* ------------------------------------------------------------------ */
interface Resolved { name: string; smiles: string; mol: Mol | null; typed: string; }

function resolveInput(inp: { type: string; value: string; originalName?: string }, fallback: string): Resolved {
  const value = (inp?.value ?? "").trim();
  const typed = value;
  let smiles = ""; let mol: Mol | null = null;
  if (value) {
    if ((inp.type || "").toUpperCase() === "SMILES") {
      mol = parseSmiles(value); if (mol) smiles = value;
    }
    if (!mol) {
      const lk = lookupCompoundSmiles(value);
      if (lk) { smiles = lk; mol = parseSmiles(lk); }
    }
    if (!mol && !/\s/.test(value)) { mol = parseSmiles(value); if (mol) smiles = value; }   // bare SMILES typed as a name
  }
  const canon = mol ? writeSmiles(mol) : "";
  const libName = canon ? libraryNameFor(canon) : null;
  const name = inp.originalName || ((inp.type || "").toUpperCase() !== "SMILES" && value ? value : libName || fallback);
  return { name, smiles, mol, typed };
}

function fmtMass(x: number): string { return (x >= 0 ? "+" : "") + x.toFixed(4); }
function cap(s: string): string { return s.charAt(0).toUpperCase() + s.slice(1); }

const CORE_CONDS: [CondKey, string][] = [
  ["acid", "Acidic stress"], ["base", "Basic stress"], ["hyd", "Neutral hydrolysis"], ["photo", "Photolytic stress"], ["therm", "Thermal stress"], ["ox", "Oxidative stress"],
];

function unresolvedResult(primary: Resolved, others: Resolved[]): PredictionResult {
  const fb = buildEntries([{ g: "hydrocarbon", atoms: [], f: {} }]);
  return {
    chainOfThought: `[Systematic Functional Group Reactivity & Computational Degradation Assessment]

INPUT RESOLUTION FAILED:
   - The primary input "${primary.typed || "(empty)"}" could not be interpreted as a valid SMILES string or as a compound name known to the local library.
   - No structure was assumed and no impurity was fabricated. Provide a valid SMILES (e.g. CC(=O)OC1=CC=CC=C1C(=O)O) or add the name with registerCompound(name, smiles).`,
    compounds: [{ name: primary.name, smiles: primary.typed, features: ["Unresolved input"], interactionSites: [] },
      ...others.map((o) => ({ name: o.name, smiles: o.smiles || o.typed, features: o.mol ? analyze(o.mol).entries.map((e) => e.fg.groupName) : ["Unresolved input"], interactionSites: [] as string[] }))],
    interactionType: "None",
    mechanism: "The primary compound could not be resolved to a structure; no mechanism can be derived.",
    functionalGroupAnalysis: fb.map((e) => e.fg),
    degradationImpurities: [],
  };
}

/**
 * Functional-group driven degradation and interaction engine (domain-agnostic).
 *  1. Graph perception of functional groups for the primary compound and each co-reactant.
 *  2. Condition-specific reactivity (acid / base / neutral hydrolysis / photolysis / thermal / oxidation).
 *  3. Graph-edit reaction templates -> valid, atom-balanced product structures.
 *  4. Cross-reactivity between complementary functional groups of the co-reactant(s).
 *  5. Ranking: heuristic likelihood from the vulnerability grades + Boltzmann weights (T = 298.15 K).
 */
export function generateComputationalPrediction(
  inputs: { type: string; value: string; originalName?: string }[],
  method: PredictionMethod = "Both"
): PredictionResult {
  const primary = resolveInput(inputs[0] ?? { type: "Name", value: "" }, "Primary Compound");
  const coRes: Resolved[] = inputs.slice(1).filter((c) => (c?.value ?? "").trim()).map((c, i) => resolveInput(c, `Co-reactant ${i + 1}`));
  if (!primary.mol) return unresolvedResult(primary, coRes);
  try { return runEngine(primary, coRes, method); }
  catch (err) {
    const fb = buildEntries([{ g: "hydrocarbon", atoms: [], f: {} }]);
    return {
      chainOfThought: `[Systematic Functional Group Reactivity & Computational Degradation Assessment]\n\nThe engine raised an unexpected error while analysing "${primary.name}": ${String((err as Error)?.message ?? err)}.`,
      compounds: [{ name: primary.name, smiles: primary.smiles, features: [], interactionSites: [] }], interactionType: "None",
      mechanism: "Analysis aborted due to an internal error.", functionalGroupAnalysis: fb.map((e) => e.fg), degradationImpurities: [],
    };
  }
}

function runEngine(primary: Resolved, coRes: Resolved[], method: PredictionMethod): PredictionResult {
  const pm = primary.mol as Mol;
  const an = analyze(pm);
  const pComps = components(pm);
  const pCanon = writeSmiles(pm);
  const pSpect = pComps.length > 1 ? new Set(pComps.map((c) => writeIsolated(pm, c))) : new Set<string>();
  const big = [...pComps].sort((a, b) => heavyCount(pm, b) - heavyCount(pm, a))[0];
  const refMass = formulaOf(pm, big).mono;
  const pFormula = formulaOf(pm, big).formula;

  const cands: Cand[] = [];
  const add = (rx: Rx[], parentCanon: string, spect: Set<string>): void => {
    for (const r of rx) {
      try { const c = finalizeRx(r, parentCanon, spect, refMass); if (c) cands.push(c); } catch { /* skip malformed template output */ }
    }
  };

  // ---- intrinsic stress pathways of the primary compound
  const stress = genStress(an).slice(0, 600);
  add(stress, pCanon, pSpect);

  // ---- co-reactant analysis and cross pathways
  const coAn: { r: Resolved; a: Analysis }[] = [];
  for (const r of coRes) if (r.mol) coAn.push({ r, a: analyze(r.mol) });
  const ambientAcid = coAn.some(({ a }) => a.classes.has("acid_strong") || a.classes.has("acid_weak")) || an.classes.has("acid_strong") || an.classes.has("acid_weak");
  for (const { r, a } of coAn) {
    const rx = genCross(an, a, primary.name, r.name, stress, ambientAcid);
    const merged = rx.filter((x) => x.mergedParent);
    const other = rx.filter((x) => !x.mergedParent);
    const { mol: M } = mergeMol(pm, a.mol);
    const mSpect = new Set<string>([...components(pm).map((c) => writeIsolated(pm, c)), ...components(a.mol).map((c) => writeIsolated(a.mol, c))]);
    add(merged, writeSmiles(M), mSpect);
    add(other, pCanon, pSpect);
  }

  // ---- de-duplicate identical products (keep the most likely origin, remember other conditions)
  const byKey = new Map<string, Cand>();
  for (const c of cands) {
    const k = writeSmiles(c.mol) + (c.physical ? "|ionic" : "");
    const cur = byKey.get(k);
    if (!cur) { byKey.set(k, c); continue; }
    const crossC = c.source === "Interaction with other compound", crossCur = cur.source === "Interaction with other compound";
    if (c.like > cur.like + 1e-9 || (Math.abs(c.like - cur.like) <= 1e-9 && crossC && !crossCur)) { c.alsoCond = [...new Set([...cur.alsoCond, cur.cond].filter((x) => x !== c.cond))]; byKey.set(k, c); }
    else if (c.cond !== cur.cond && !cur.alsoCond.includes(c.cond)) cur.alsoCond.push(c.cond);
  }
  const unique = [...byKey.values()];

  // ---- Boltzmann distribution (T = 298.15 K): P_i = exp(-dG_i / RT) / sum_j exp(-dG_j / RT)
  const R = 0.0019872; // kcal / (mol * K)
  const T = 298.15;
  const RT = R * T;
  const expTerms = unique.map((c) => Math.exp(-c.dG / RT));
  const sumExp = expTerms.reduce((acc, v) => acc + v, 0) || 1;
  const coNamesAll = coAn.map((c) => c.r.name);

  const calculated = unique.map((cand, idx) => {
    const pBoltzmann = Math.min(0.99, Math.max(0.01, Number((expTerms[idx] / sumExp).toFixed(4))));
    const pHeuristic = Math.min(0.99, Math.max(0.01, Number(cand.like.toFixed(4))));
    let probability = pHeuristic;
    if (method === "Boltzmann") probability = pBoltzmann;
    else if (method === "Both") probability = Number(((pHeuristic + pBoltzmann) / 2).toFixed(4));
    const cross = cand.source === "Interaction with other compound";
    const frags = cand.names;
    const sub = cand.smiles.includes(".") ? ` Fragments: ${frags}.` : "";
    const trivial = cand.dropped.length ? ` Co-product(s) not listed: ${[...new Set(cand.dropped)].join(", ")}.` : "";
    const desc = cand.physical
      ? `Ionic / non-covalent association product (${cand.mainFormula}; no covalent bond change).${sub}`
      : `${cap(cand.label)}. Main species ${cand.mainFormula}, monoisotopic mass ${cand.mainMass.toFixed(4)} Da (${fmtMass(cand.deltaMass)} Da vs ${pFormula} parent).${sub}${trivial}`;
    const also = cand.alsoCond.length ? ` Also formed under: ${cand.alsoCond.join(", ")}.` : "";
    return {
      iupacName: `${cap(cand.label)}: ${frags}`,
      smiles: cand.smiles,
      structureDescription: desc,
      origin: cross ? `${primary.name} + ${cand.co ?? coNamesAll[0] ?? "co-reactant"}` : primary.name,
      condition: cand.cond,
      source: cand.source,
      mechanismExplanation: cand.mech + also,
      relativeEnergy: Number(cand.dG.toFixed(2)),
      probability,
      probabilityHeuristic: pHeuristic,
      probabilityBoltzmann: pBoltzmann,
      _cand: cand,
    };
  });

  // Sort strictly by formation probability descending
  calculated.sort((a, b) => b.probability - a.probability);
  let top = calculated.slice(0, 5);
  // when co-reactants are supplied, make sure the most probable interaction product is reported
  const crossPool = calculated.filter((c) => c._cand.source === "Interaction with other compound" && c.probabilityHeuristic >= 0.5);
  if (coAn.length && crossPool.length && !top.some((t) => crossPool.includes(t))) top = [...calculated.slice(0, 4), crossPool[0]];
  const degradationImpurities = top.map(({ _cand, ...rest }) => { void _cand; return rest; });

  // ---- interaction classification
  const crossAll = cands.filter((c) => c.source === "Interaction with other compound" && c.like >= 0.5).sort((a, b) => b.like - a.like)
    .map((c) => ({ _cand: c }));
  const chemical = crossAll.filter((c) => !c._cand.physical);
  const interactionType: "Physical" | "Chemical" | "None" = !coAn.length ? "None" : chemical.length ? "Chemical" : crossAll.length ? "Physical" : "None";
  const bestCross = (chemical[0] ?? crossAll[0])?._cand;
  const noPathways = unique.length === 0;

  // ---- functional-group table (secondary column evaluated against actual co-reactants)
  let fgTable = an.entries.map((e) => e.fg);
  if (coAn.length) {
    const union = new Set<string>(); coAn.forEach(({ a }) => a.classes.forEach((c) => union.add(c)));
    fgTable = refineSecondary(an.entries, coNamesAll.join(" + "), union);
  }

  // ---- compounds list
  const compoundsList = [{ name: primary.name, smiles: primary.smiles, features: an.entries.map((e) => e.fg.groupName), interactionSites: an.entries.map((e) => e.fg.reactiveSite) }];
  for (const r of coRes) {
    const a = coAn.find((x) => x.r === r)?.a;
    compoundsList.push({ name: r.name, smiles: r.smiles || r.typed, features: a ? a.entries.map((e) => e.fg.groupName) : ["Unresolved input (no structure available)"], interactionSites: a ? a.entries.map((e) => e.fg.reactiveSite) : [] });
  }

  // ---- narrative
  const grade = (e: GroupEntry, k: CondKey): Vuln => e.prof[k][0];
  const stressLines = CORE_CONDS.map(([k, lab]) => {
    const hits = an.entries.filter((e) => VRANK(grade(e, k)) >= 2).sort((a, b) => VRANK(grade(b, k)) - VRANK(grade(a, k))).slice(0, 3);
    return `   - ${lab}: ` + (hits.length ? hits.map((e) => `${e.fg.groupName} [${grade(e, k)}]`).join("; ") : "no group above Moderate susceptibility");
  }).join("\n");
  const coBlock = coAn.length ? `CROSS-FUNCTIONAL INTERACTIONS WITH CO-REACTANTS:
${coAn.map(({ r, a }) => `   - Co-reactant evaluated: ${r.name} (SMILES: ${r.smiles})
   - Co-reactant functional groups: ${a.entries.map((e) => e.fg.groupName).join(", ") || "none detected"}
   - Reactive classes detected: ${[...a.classes].map((c) => CLASS_TEXT[c] || c).join("; ") || "none"}`).join("\n")}
   - Templates evaluated: acyl transfer (aminolysis / transesterification), Schiff-base / Maillard condensation, esterification / amidation, epoxide & alkyl-halide alkylation, Michael addition, isocyanate addition, N-nitrosation / diazotisation, proton-transfer salt formation, metal carboxylate formation, redox (oxidant / reductant) and microenvironment (pH, moisture, peroxide, photocatalyst) effects.
   - Result: ${bestCross ? `${bestCross.label} (${bestCross.physical ? "ionic / physical" : "covalent"})` : "no complementary reactive pair with meaningful likelihood; only non-covalent association expected"}.` : "MULTI-COMPONENT MATRIX: Single-compound intrinsic forced degradation analysis.";
  const unresolvedNote = coRes.filter((r) => !r.mol).length ? `\n   - Unresolved co-reactant input(s) ignored (no structure assumed): ${coRes.filter((r) => !r.mol).map((r) => `"${r.typed}"`).join(", ")}.` : "";

  const chainOfThought = `[Systematic Functional Group Reactivity & Computational Degradation Assessment]

PRIMARY MOLECULAR FUNCTIONAL GROUP INVENTORY:
   - Target Structure: ${primary.name} (SMILES: ${primary.smiles}; ${pFormula}, monoisotopic ${refMass.toFixed(4)} Da)
   - Identified Functional Groups: ${an.entries.map((e) => `${e.fg.groupName} [${e.fg.category}]`).join(", ")}
   - Identified Reactive Centers: ${an.entries.map((e) => e.fg.reactiveSite).join("; ")}

REACTION SUSCEPTIBILITY BY STRESS CONDITION (most vulnerable groups):
${stressLines}

${coBlock}${unresolvedNote}

THERMODYNAMIC BOLTZMANN PARTITION & KINETIC PROBABILITIES (T = 298.15 K):
   - ${cands.length} candidate pathways were generated from graph-edit templates; ${unique.length} distinct products remained after de-duplication (identical products from several conditions are merged).
   - Relative formation free energies (dG, kcal/mol) are reaction-class estimates on a compressed scale (ring strain, bond-energy and ionisation terms), not quantum-chemical values - confirm key pathways with DFT if needed.
   - Boltzmann distribution derived via P_i = exp(-dG_i / RT) / sum exp(-dG_j / RT); heuristic likelihood derives from the vulnerability grade of the reacting group under the same condition (Critical 0.93, High 0.80, Moderate 0.58, Low 0.30, Resistant 0.10) scaled by pathway-specific factors. Method used: ${method}.
   - Ranked top ${top.length} dominant degradation impurities and reaction adducts.`;

  const groupsTxt = an.entries.slice(0, 4).map((e) => e.fg.groupName).join(", ");
  const mechanism = noPathways
    ? `${primary.name} shows no functional-group liability under the evaluated stress conditions (${groupsTxt}); it is chemically inert or ionic under acid, base, neutral hydrolysis, light, heat and oxidative stress.`
    : coAn.length
    ? (bestCross
      ? `Cross-functional interaction (${bestCross.label}): ${bestCross.mech}`
      : `No complementary reactive functional-group pair was found between ${primary.name} and ${coNamesAll.join(" + ")}; interaction is limited to non-covalent association. Intrinsic stress degradation is governed by ${groupsTxt}.`)
    : `Intrinsic stress degradation governed by hydrolytic, oxidative, photolytic and thermal reactivity of the functional groups present in ${primary.name} (${groupsTxt}).`;

  return {
    chainOfThought,
    compounds: compoundsList,
    interactionType,
    mechanism,
    functionalGroupAnalysis: fgTable,
    degradationImpurities,
  };
}
