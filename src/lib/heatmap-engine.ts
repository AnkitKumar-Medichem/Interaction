export interface ReactiveInteractionDetail {
  row: number;
  col: number;
  rowLabel: string;
  colLabel: string;
  score: number;
  percentage: string;
  mechanism: string;
  severity: "Critical" | "High" | "Moderate" | "Low";
  description: string;
}

export interface HeatmapResponse {
  success: boolean;
  compound1Name: string;
  compound2Name: string;
  rowLabels: string[];
  colLabels: string[];
  matrix: number[][];
  details: ReactiveInteractionDetail[];
  topInteractions: ReactiveInteractionDetail[];
  heatmapBase64: string;
  summary: string;
  error?: string;
}

export interface HeatmapRequestCompound {
  name: string;
  smiles?: string;
  features?: string[];
  interactionSites?: string[];
}

export const DEFAULT_CONDITIONS = [
  "Acidic",
  "Basic",
  "Hydrolysis",
  "Photolysis",
  "Thermal",
  "Oxidative",
] as const;

export type ForcedDegradationCondition = (typeof DEFAULT_CONDITIONS)[number];

export interface HeatmapOptions {
  compounds: HeatmapRequestCompound[];
  cmap?: string;
  title?: string;
  conditions?: string[];
}

/**
 * Evaluates chemical site susceptibility under specific forced degradation / environmental stress conditions:
 * Acidic, Basic, Hydrolysis, Photolysis, Thermal, Oxidative.
 * Also accounts for co-formulated excipients and counter-compounds in the formulation matrix.
 */
export function calculateConditionVulnerability(
  site: string,
  condition: string,
  compounds: HeatmapRequestCompound[] = []
): {
  score: number;
  mechanism: string;
  severity: "Critical" | "High" | "Moderate" | "Low";
  description: string;
} {
  const s = site.toLowerCase();
  const c = condition.toLowerCase();

  // Detect co-formulated microenvironment influences
  const allNames = compounds.map((cp) => cp.name.toLowerCase()).join(" ");
  const allSites = compounds
    .flatMap((cp) => [...(cp.interactionSites || []), ...(cp.features || [])])
    .join(" ")
    .toLowerCase();

  const hasAlkalinePartner = [
    "magnesium oxide",
    "mgo",
    "hydroxide",
    "calcium carbonate",
    "basic",
    "alkaline",
    "sodium bicarbonate",
    "nahco3",
    "meglumine",
  ].some((k) => allNames.includes(k) || allSites.includes(k));

  const hasOxidantPartner = [
    "povidone",
    "pvp",
    "crospovidone",
    "polysorbate",
    "peroxide",
    "peg",
    "polyethylene glycol",
  ].some((k) => allNames.includes(k) || allSites.includes(k));

  const hasReducingSugar = [
    "lactose",
    "glucose",
    "dextrose",
    "maltose",
    "reducing sugar",
  ].some((k) => allNames.includes(k) || allSites.includes(k));

  const hasMetalExcipient = [
    "magnesium stearate",
    "mg stearate",
    "calcium",
    "zinc",
    "metal",
    "mg2+",
  ].some((k) => allNames.includes(k) || allSites.includes(k));

  // Site classification
  const isEster = ["ester", "acetyl", "acetoxy", "carbonyl", "acyl"].some((k) => s.includes(k));
  const isPhenol = ["phenol", "phenolic", "hydroxyl", "hydroxy"].some((k) => s.includes(k)) && !s.includes("carboxyl");
  const isCatechol = ["catechol", "ortho-diphenol", "dihydroxyphenyl", "dioh"].some((k) => s.includes(k));
  const isAmine = ["amine", "amino", "methylamino", "biguanide", "guanidine", "nh2", "nh"].some((k) => s.includes(k));
  const isAcid = ["carboxylic acid", "carboxyl", "acid group", "cooh", "salicylate"].some((k) => s.includes(k));
  const isAmide = ["amide", "acetamide", "anilide", "lactam"].some((k) => s.includes(k));
  const isAromatic = ["aromatic", "benzene", "phenyl", "ring"].some((k) => s.includes(k));
  const isSugar = ["aldose", "lactose", "sugar", "hemiacetal", "anomeric"].some((k) => s.includes(k));
  const isAlkalineCenter = ["alkaline", "oxide", "basic surface", "hydroxide"].some((k) => s.includes(k));
  const isStearateCenter = ["stearate", "fatty acid", "aliphatic"].some((k) => s.includes(k));
  const isInterfacial = ["interface", "coupling", "cross-interaction", "contact"].some((k) => s.includes(k));

  // 1. ACIDIC CONDITION
  if (c.includes("acid")) {
    if (isAlkalineCenter) {
      return {
        score: 0.96,
        mechanism: "Acid-Base Stoichiometric Neutralization",
        severity: "Critical",
        description: "Rapid stoichiometric reaction with hydronium ions (H3O+) causing dissolution and counter-cation dissociation.",
      };
    }
    if (isSugar) {
      return {
        score: 0.86,
        mechanism: "Acidic Glycosidic / Acetal Cleavage",
        severity: "Critical",
        description: "Protonation of glycosidic oxygen generates oxocarbenium intermediates, causing hydrolytic chain scission.",
      };
    }
    if (isEster) {
      const boost = hasMetalExcipient ? 0.08 : 0;
      return {
        score: Math.min(0.92, 0.74 + boost),
        mechanism: "Specific Acid-Catalyzed Ester Solvolysis",
        severity: "High",
        description: "Carbonyl oxygen protonation increases electrophilicity toward nucleophilic water attack and cleavage.",
      };
    }
    if (isAmide) {
      return {
        score: 0.44,
        mechanism: "Acid-Promoted Amide Hydrolysis",
        severity: "Moderate",
        description: "Hydronium activation enables slow nucleophilic cleavage of secondary amide bond.",
      };
    }
    if (isAmine) {
      return {
        score: 0.16,
        mechanism: "Ammonium Salt Stabilization",
        severity: "Low",
        description: "Protonation forms a stable conjugate ammonium cation, shielding lone pairs from nucleophilic addition.",
      };
    }
    if (isAcid || isPhenol || isAromatic) {
      return {
        score: 0.12,
        mechanism: "Acidic Protonation Inertness",
        severity: "Low",
        description: "Suppressed ionization maintains unionized covalent state with low degradation potential.",
      };
    }
    return {
      score: 0.20,
      mechanism: "Acidic Environmental Exposure",
      severity: "Low",
      description: "Low-level susceptibility under low pH aqueous conditions.",
    };
  }

  // 2. BASIC CONDITION
  if (c.includes("basic") || c.includes("alkaline")) {
    if (isEster) {
      const boost = hasAlkalinePartner ? 0.05 : 0;
      return {
        score: Math.min(0.98, 0.94 + boost),
        mechanism: "Base-Promoted Ester Saponification (B_Ac2)",
        severity: "Critical",
        description: "Direct hydroxide (OH-) nucleophilic attack at ester carbonyl generates tetrahedral intermediate followed by cleavage.",
      };
    }
    if (isCatechol) {
      return {
        score: 0.97,
        mechanism: "Base-Catalyzed Catecholate Deprotonation",
        severity: "Critical",
        description: "Rapid deprotonation of ortho-hydroxyls creates electron-rich catecholate, accelerating autoxidation into reactive ortho-quinones.",
      };
    }
    if (isPhenol) {
      return {
        score: 0.88,
        mechanism: "Phenolate Ion Formation & Autoxidation Activation",
        severity: "Critical",
        description: "Deprotonation to phenolate anion drastically lowers oxidation potential, triggering oxidative discoloration and dimerization.",
      };
    }
    if (isSugar) {
      return {
        score: 0.85,
        mechanism: "Lobry de Bruyn-Alberda van Ekenstein Rearrangement",
        severity: "High",
        description: "Base-catalyzed enolization and dicarbonyl fragmentation yielding reactive advanced degradation products.",
      };
    }
    if (isAmide) {
      return {
        score: 0.64,
        mechanism: "Hydroxide-Assisted Amide Cleavage",
        severity: "Moderate",
        description: "Hydroxide attacks secondary amide carbonyl at elevated alkaline pH, generating carboxylate and free amine.",
      };
    }
    if (isInterfacial) {
      return {
        score: hasAlkalinePartner ? 0.95 : 0.60,
        mechanism: "Alkaline Interfacial Hydrolytic Destabilization",
        severity: hasAlkalinePartner ? "Critical" : "Moderate",
        description: "Basic microenvironment promotes localized interfacial proton transfer and catalytic degradation.",
      };
    }
    if (isAcid) {
      return {
        score: 0.28,
        mechanism: "Carboxylate Salt Formation",
        severity: "Low",
        description: "Forms stable deprotonated carboxylate salt with alkaline counter-ions.",
      };
    }
    if (isAlkalineCenter) {
      return {
        score: 0.10,
        mechanism: "Basic Phase Stability",
        severity: "Low",
        description: "Compatible with basic microenvironment with minimal chemical change.",
      };
    }
    return {
      score: 0.22,
      mechanism: "Alkaline Surface Interaction",
      severity: "Low",
      description: "Mild electrostatic interaction under alkaline conditions.",
    };
  }

  // 3. HYDROLYSIS (Neutral / Moisture)
  if (c.includes("hydrolysis") || c.includes("aqueous") || c.includes("moisture")) {
    if (isEster) {
      const boost = hasAlkalinePartner ? 0.20 : 0;
      return {
        score: Math.min(0.95, 0.72 + boost),
        mechanism: hasAlkalinePartner ? "Alkaline-Accelerated Moisture Hydrolysis" : "Aqueous Ester Solvolysis",
        severity: hasAlkalinePartner ? "Critical" : "High",
        description: hasAlkalinePartner
          ? "Excipient microenvironment raises local pH, drastically accelerating moisture-driven ester cleavage."
          : "Water nucleophilic attack across ester carbonyl initiates hydrolytic degradation.",
      };
    }
    if (isSugar) {
      return {
        score: 0.48,
        mechanism: "Aqueous Anomeric Mutarotation",
        severity: "Moderate",
        description: "Solvent-assisted equilibration between anomers via transient open-chain carbonyl form.",
      };
    }
    if (isAmide) {
      return {
        score: 0.32,
        mechanism: "Neutral Amide Solvolysis",
        severity: "Low",
        description: "Amide resonance stabilization provides resistance against neutral water-mediated cleavage.",
      };
    }
    if (isAcid || isPhenol) {
      return {
        score: 0.14,
        mechanism: "Aqueous Hydration Equilibrium",
        severity: "Low",
        description: "Hydrogen-bonded solvation shell formation with negligible covalent bond rupture.",
      };
    }
    if (isStearateCenter || isAromatic) {
      return {
        score: 0.08,
        mechanism: "Hydrophobic Barrier Shield",
        severity: "Low",
        description: "Hydrophobic properties retard aqueous ingress and interfacial hydrolysis.",
      };
    }
    return {
      score: 0.18,
      mechanism: "Moisture-Mediated Interaction",
      severity: "Low",
      description: "Low solvolytic susceptibility under neutral humidity conditions.",
    };
  }

  // 4. PHOTOLYSIS (Light / UV-Vis exposure)
  if (c.includes("photo") || c.includes("light") || c.includes("uv")) {
    if (isEster && isPhenol) {
      return {
        score: 0.82,
        mechanism: "Photo-Fries Intramolecular Rearrangement",
        severity: "High",
        description: "UV photon absorption triggers homolytic acyl-oxygen cleavage followed by radical cage recombination at ortho/para ring positions.",
      };
    }
    if (isCatechol || isPhenol) {
      return {
        score: 0.78,
        mechanism: "Photochemical Phenoxyl Radical Generation",
        severity: "High",
        description: "UV irradiation excites phenolic chromophores, generating radical intermediates prone to oxidative dimerization.",
      };
    }
    if (isAromatic) {
      return {
        score: 0.62,
        mechanism: "Pi-Pi* Conjugated Chromophore Excitation",
        severity: "Moderate",
        description: "Absorption of UV wavelengths promotes molecules to reactive triplet excited states.",
      };
    }
    if (isAmine || isAmide) {
      return {
        score: 0.50,
        mechanism: "Photo-Induced Electron Transfer (PET)",
        severity: "Moderate",
        description: "Nitrogen lone pairs undergo photo-sensitized electron transfer yielding reactive aminium radical species.",
      };
    }
    if (isAlkalineCenter || isStearateCenter || isSugar) {
      return {
        score: 0.08,
        mechanism: "UV-Vis Photostability",
        severity: "Low",
        description: "Absence of conjugated chromophores renders center transparent to near-UV radiation.",
      };
    }
    return {
      score: 0.25,
      mechanism: "Photochemical Excitation",
      severity: "Low",
      description: "Low photochemical cross-section under ambient light.",
    };
  }

  // 5. THERMAL STRESS (High temperature / Thermolysis)
  if (c.includes("thermal") || c.includes("heat") || c.includes("thermolysis")) {
    if (isAmine && hasReducingSugar) {
      return {
        score: 0.95,
        mechanism: "Thermally-Driven Maillard Condensation",
        severity: "Critical",
        description: "Elevated temperature drives nucleophilic condensation of amine with reducing sugar carbonyl into glycosylamine and Amadori intermediates.",
      };
    }
    if (isEster) {
      return {
        score: 0.78,
        mechanism: "Thermally-Activated Transesterification",
        severity: "High",
        description: "Thermal kinetic energy overcomes activation barrier (ΔH‡) for intermolecular acyl transfer and oligomerization.",
      };
    }
    if (isAcid) {
      return {
        score: 0.72,
        mechanism: "Thermal Decarboxylation & Anhydride Condensation",
        severity: "High",
        description: "Thermal stress drives concerted decarboxylation or intermolecular condensation into anhydride dimers.",
      };
    }
    if (isCatechol || isPhenol) {
      return {
        score: 0.68,
        mechanism: "Thermal Radical Coupling",
        severity: "Moderate",
        description: "Accelerated radical generation facilitates oxidative polymerization and colored adduct formation.",
      };
    }
    if (isStearateCenter) {
      return {
        score: 0.55,
        mechanism: "Thermal Phase Transition & Softening",
        severity: "Moderate",
        description: "Excipient softening increases molecular mobility and solid-state reaction kinetics.",
      };
    }
    if (isAlkalineCenter) {
      return {
        score: 0.15,
        mechanism: "High Lattice Thermal Stability",
        severity: "Low",
        description: "Inorganic lattice structure maintains stability under pharmaceutical thermal stress temperatures.",
      };
    }
    return {
      score: 0.35,
      mechanism: "Thermally-Induced Kinetic Acceleration",
      severity: "Low",
      description: "Thermal acceleration of baseline reaction pathways according to Arrhenius kinetics.",
    };
  }

  // 6. OXIDATIVE STRESS (Peroxides / Radicals / Auto-oxidation)
  if (c.includes("oxid") || c.includes("peroxide") || c.includes("radical")) {
    if (isCatechol) {
      return {
        score: 0.98,
        mechanism: "Rapid Ortho-Quinone Auto-Oxidation",
        severity: "Critical",
        description: "Electrophilic radical abstraction of catechol protons yields reactive ortho-quinones prone to rapid polymerization.",
      };
    }
    if (isPhenol) {
      const boost = hasOxidantPartner ? 0.08 : 0;
      return {
        score: Math.min(0.96, 0.88 + boost),
        mechanism: "Radical Phenoxyl Auto-Oxidation",
        severity: "Critical",
        description: "Peroxide and oxygen radicals abstract phenolic hydrogen, initiating quinonoid and dimer formation.",
      };
    }
    if (isAmine) {
      const boost = hasOxidantPartner ? 0.10 : 0;
      return {
        score: Math.min(0.90, 0.76 + boost),
        mechanism: "Peroxide-Mediated N-Oxidation",
        severity: "High",
        description: "Electrophilic oxygen transfer from peroxides onto nitrogen lone pair yielding N-oxide and hydroxylamine species.",
      };
    }
    if (isAromatic) {
      return {
        score: 0.65,
        mechanism: "Electrophilic Aromatic Ring Hydroxylation",
        severity: "Moderate",
        description: "Hydroxyl radical attack at activated aromatic positions generates dihydroxybenzoic and phenolic derivatives.",
      };
    }
    if (isEster) {
      return {
        score: 0.36,
        mechanism: "Secondary Oxidative Cleavage",
        severity: "Low",
        description: "Ester group displays moderate resistance; cleavage occurs secondary to aromatic radical intermediates.",
      };
    }
    if (isAcid || isAlkalineCenter || isStearateCenter) {
      return {
        score: 0.12,
        mechanism: "Oxidative Inertness",
        severity: "Low",
        description: "Fully oxidized or saturated functional center exhibits high resistance to oxygen transfer.",
      };
    }
    return {
      score: 0.28,
      mechanism: "Peroxide Surface Exposure",
      severity: "Low",
      description: "Moderate susceptibility to ambient radical oxidation.",
    };
  }

  // Fallback
  return {
    score: 0.20,
    mechanism: "Environmental Stress Interaction",
    severity: "Low",
    description: "Standard baseline degradation potential under experimental conditions.",
  };
}

// ---------------------------------------------------------------------------
// Divergent WarmCool Palette & Color Interpolation
// ---------------------------------------------------------------------------
const PALETTES: Record<string, number[][]> = {
  warmcool: [
    [59, 76, 192],    // 0.00: Deep cool blue
    [98, 130, 234],   // 0.12: Soft blue
    [141, 176, 254],  // 0.25: Light sky blue
    [184, 208, 249],  // 0.37: Soft ice blue
    [221, 221, 221],  // 0.50: Neutral slate
    [245, 183, 162],  // 0.62: Soft warm salmon
    [237, 132, 107],  // 0.75: Warm coral
    [211, 74, 62],    // 0.87: Strong crimson
    [180, 4, 38],     // 1.00: Deep warm red
  ],
  coolwarm: [
    [59, 76, 192],
    [98, 130, 234],
    [141, 176, 254],
    [184, 208, 249],
    [221, 221, 221],
    [245, 183, 162],
    [237, 132, 107],
    [211, 74, 62],
    [180, 4, 38],
  ],
};

function interpolateColor(t: number, paletteKey: string): string {
  const stops = PALETTES[paletteKey] || PALETTES["warmcool"];
  const clamped = Math.max(0, Math.min(1, t));
  const pos = clamped * (stops.length - 1);
  const idx = Math.floor(pos);
  const frac = pos - idx;

  const c1 = stops[idx];
  const c2 = stops[Math.min(stops.length - 1, idx + 1)];

  const r = Math.round(c1[0] + frac * (c2[0] - c1[0]));
  const g = Math.round(c1[1] + frac * (c2[1] - c1[1]));
  const b = Math.round(c1[2] + frac * (c2[2] - c1[2]));
  return `rgb(${r}, ${g}, ${b})`;
}

function getBrightness(rgbStr: string): number {
  const match = rgbStr.match(/\d+/g);
  if (!match || match.length < 3) return 128;
  const [r, g, b] = match.map(Number);
  return (r * 299 + g * 587 + b * 114) / 1000;
}

function escapeXml(unsafe: string): string {
  return unsafe
    .replace(/&/g, "&amp;")
    .replace(/</g, "&lt;")
    .replace(/>/g, "&gt;")
    .replace(/"/g, "&quot;")
    .replace(/'/g, "&apos;");
}

// ---------------------------------------------------------------------------
// Main Heatmap Generator
// ---------------------------------------------------------------------------
export function generateInteractionHeatmap(options: HeatmapOptions): HeatmapResponse {
  const compounds = options.compounds || [];
  if (compounds.length === 0) {
    throw new Error("At least one compound is required for interaction heatmap.");
  }

  const comp1 = compounds[0];
  const comp1Name = comp1.name || "Compound 1";
  const comp2 = compounds.length > 1 ? compounds[1] : null;
  const comp2Name = comp2 ? comp2.name || "Compound 2" : "Intramolecular";

  // Helper to extract clean functional group name without & or parentheses
  function cleanFgName(raw: string): string {
    return raw
      .replace(/^\[[^\]]*\]\s*/, "")
      .replace(/\s*\([^)]*\)/g, "")
      .replace(/\s*&.*$/, "")
      .replace(/\s*↔.*$/, "")
      .replace(/Reactive Center/i, "Aliphatic Center")
      .trim();
  }

  // Functional groups placed on X-axis (colLabels)
  const colLabels: string[] = [];

  const c1Sites = (comp1.interactionSites && comp1.interactionSites.length > 0)
    ? comp1.interactionSites
    : (comp1.features && comp1.features.length > 0)
    ? comp1.features
    : ["Aliphatic Center"];

  for (const s of c1Sites.slice(0, 5)) {
    const c = cleanFgName(s);
    if (c && !colLabels.includes(c)) {
      colLabels.push(c);
    }
  }

  if (comp2) {
    const c2Sites = (comp2.interactionSites && comp2.interactionSites.length > 0)
      ? comp2.interactionSites
      : (comp2.features && comp2.features.length > 0)
      ? comp2.features
      : ["Additive Center"];
    for (const s of c2Sites.slice(0, 4)) {
      const c = cleanFgName(s);
      if (c && !colLabels.includes(c)) {
        colLabels.push(c);
      }
    }
  }

  if (colLabels.length === 0) {
    colLabels.push("Aliphatic Framework");
  }

  // Conditions placed on Y-axis (rowLabels)
  const rowLabels: string[] =
    options.conditions && options.conditions.length > 0
      ? options.conditions
      : [...DEFAULT_CONDITIONS];

  const matrix: number[][] = [];
  const details: ReactiveInteractionDetail[] = [];

  // Deterministic hash helper for reproducible micro-variations
  function hashPair(str: string): number {
    let hash = 0;
    for (let i = 0; i < str.length; i++) {
      hash = (hash << 5) - hash + str.charCodeAt(i);
      hash |= 0;
    }
    return Math.abs(hash);
  }

  for (let r = 0; r < rowLabels.length; r++) {
    const rowRow: number[] = [];
    const conditionLabel = rowLabels[r];
    for (let c = 0; c < colLabels.length; c++) {
      const siteLabel = colLabels[c];

      const res = calculateConditionVulnerability(siteLabel, conditionLabel, compounds);
      const hashOffset = ((hashPair(siteLabel + conditionLabel) % 9) - 4) * 0.01;
      const score = Math.max(0.05, Math.min(0.99, Number((res.score + hashOffset).toFixed(2))));

      rowRow.push(score);
      details.push({
        row: r,
        col: c,
        rowLabel: conditionLabel,
        colLabel: siteLabel,
        score,
        percentage: `${Math.round(score * 100)}%`,
        mechanism: res.mechanism,
        severity: res.severity,
        description: res.description,
      });
    }
    matrix.push(rowRow);
  }

  // Rank top interactions / vulnerabilities
  const topInteractions = [...details]
    .sort((a, b) => b.score - a.score)
    .slice(0, 4);

  // SVG Heatmap Construction
  const cmap = options.cmap || "warmcool";
  const numRows = rowLabels.length;
  const numCols = colLabels.length;

  const cellWidth = 112;
  const cellHeight = 48;

  const margin = { top: 96, right: 100, bottom: 85, left: 240 };
  const gridWidth = numCols * cellWidth;
  const gridHeight = numRows * cellHeight;
  const totalWidth = margin.left + gridWidth + margin.right;
  const totalHeight = margin.top + gridHeight + margin.bottom;

  let svgCells = "";
  for (let r = 0; r < numRows; r++) {
    for (let c = 0; c < numCols; c++) {
      const val = matrix[r][c];
      const color = interpolateColor(val, cmap);
      const x = margin.left + c * cellWidth;
      const y = margin.top + r * cellHeight;
      const textColor = getBrightness(color) < 140 ? "#FFFFFF" : "#0F172A";

      svgCells += `
        <g class="cell-group">
          <rect
            x="${x + 2}"
            y="${y + 2}"
            width="${cellWidth - 4}"
            height="${cellHeight - 4}"
            rx="6"
            fill="${color}"
            stroke="#CBD5E1"
            stroke-width="0.75"
          />
        </g>
      `;
    }
  }

  // Row labels (Chemical Centers)
  let svgRowLabels = "";
  for (let r = 0; r < numRows; r++) {
    const y = margin.top + r * cellHeight + cellHeight / 2 + 4;
    const rawLabel = rowLabels[r];
    const label = rawLabel.length > 30 ? rawLabel.substring(0, 28) + "…" : rawLabel;
    svgRowLabels += `
      <text
        x="${margin.left - 14}"
        y="${y}"
        font-family="system-ui, -apple-system, sans-serif"
        font-size="12"
        font-weight="600"
        fill="#334155"
        text-anchor="end"
      >${escapeXml(label)}</text>
    `;
  }

  // Col labels (Conditions)
  let svgColLabels = "";
  for (let c = 0; c < numCols; c++) {
    const x = margin.left + c * cellWidth + cellWidth / 2;
    const y = margin.top - 14;
    const label = colLabels[c];
    svgColLabels += `
      <text
        x="${x}"
        y="${y}"
        font-family="system-ui, -apple-system, sans-serif"
        font-size="13"
        font-weight="700"
        fill="#0F172A"
        text-anchor="middle"
      >${escapeXml(label)}</text>
    `;
  }

  // Vertical color bar legend on right side
  const cbX = margin.left + gridWidth + 30;
  const cbY = margin.top + 8;
  const cbW = 16;
  const cbH = Math.max(120, gridHeight - 16);

  let gradientStops = "";
  const numStops = 10;
  for (let i = 0; i <= numStops; i++) {
    const t = 1 - i / numStops;
    const col = interpolateColor(t, cmap);
    gradientStops += `<stop offset="${(i / numStops) * 100}%" stop-color="${col}" />\n`;
  }

  const cbSvg = `
    <defs>
      <linearGradient id="legendGrad" x1="0" y1="0" x2="0" y2="1">
        ${gradientStops}
      </linearGradient>
    </defs>
    <rect x="${cbX}" y="${cbY}" width="${cbW}" height="${cbH}" rx="4" fill="url(#legendGrad)" stroke="#CBD5E1" stroke-width="0.75"/>
    <text x="${cbX + cbW + 8}" y="${cbY + 10}" font-family="system-ui, sans-serif" font-size="10" font-weight="700" fill="#0F172A">1.00</text>
    <text x="${cbX + cbW + 8}" y="${cbY + cbH / 2 + 3}" font-family="system-ui, sans-serif" font-size="10" font-weight="600" fill="#64748B">0.50</text>
    <text x="${cbX + cbW + 8}" y="${cbY + cbH}" font-family="system-ui, sans-serif" font-size="10" font-weight="700" fill="#0F172A">0.00</text>
    <text x="${cbX - 4}" y="${cbY - 10}" font-family="system-ui, sans-serif" font-size="9" font-weight="700" fill="#64748B">RISK</text>
  `;

  const chartTitle =
    options.title ||
    (comp2
      ? `Stress Degradation & Incompatibility Heatmap: ${comp1Name} & ${comp2Name}`
      : `Stress Degradation Heatmap: ${comp1Name}`);

  const svgFull = `
    <svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 ${totalWidth} ${totalHeight}" width="100%" height="100%">
      <rect width="100%" height="100%" fill="#FFFFFF" rx="16"/>
      <!-- Header -->
      <text x="${totalWidth / 2}" y="36" font-family="'Playfair Display', Georgia, serif" font-size="18" font-weight="700" fill="#0F172A" text-anchor="middle">
        ${escapeXml(chartTitle)}
      </text>
      <text x="${totalWidth / 2}" y="56" font-family="system-ui, sans-serif" font-size="11" font-weight="500" fill="#64748B" text-anchor="middle">
        Vulnerability Assessment across Forced Degradation Stress Conditions (WarmCool Palette)
      </text>

      <!-- Cells -->
      ${svgCells}

      <!-- Labels -->
      ${svgRowLabels}
      ${svgColLabels}

      <!-- Legend -->
      ${cbSvg}

      <!-- Bottom Axis Label -->
      <text x="${margin.left + gridWidth / 2}" y="${margin.top + gridHeight + 35}" font-family="system-ui, sans-serif" font-size="11" font-weight="700" fill="#64748B" text-anchor="middle">
        Forced Degradation &amp; Environmental Stress Conditions (ICH Q1A/Q1B Standards)
      </text>
    </svg>
  `.trim();

  const svgBase64 = `data:image/svg+xml;base64,${Buffer.from(svgFull).toString("base64")}`;

  const maxScore = Math.max(...details.map((d) => d.score));
  const topItem = topInteractions[0];
  const summary = `Evaluated ${numRows} reactive site(s) across ${numCols} stress degradation condition(s) (Acidic, Basic, Hydrolysis, Photolysis, Thermal, Oxidative). Peak susceptibility detected at '${topItem?.rowLabel}' under '${topItem?.colLabel}' (${maxScore.toFixed(2)} - ${topItem?.mechanism || "Degradation Pathway"}).`;

  return {
    success: true,
    compound1Name: comp1Name,
    compound2Name: comp2 ? comp2Name : "Intramolecular",
    rowLabels,
    colLabels,
    matrix,
    details,
    topInteractions,
    heatmapBase64: svgBase64,
    summary,
  };
}
