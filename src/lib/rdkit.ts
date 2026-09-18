interface RDKitModule {
  get_mol: (smiles: string) => RDKitMolecule | null;
}

interface RDKitMolecule {
  is_valid: () => boolean;
  get_smiles: () => string;
  delete: () => void;
  get_descriptors: () => string;
  get_svg: (width?: number, height?: number) => string;
}

export interface MolecularDescriptors {
  MolWt?: number;
  MolLogP?: number;
  TPSA?: number;
  NumRotatableBonds?: number;
  HBD?: number;
  HBA?: number;
  HeavyAtomCount?: number;
  NumAromaticRings?: number;
  NumHeteroatoms?: number;
  FractionCSP3?: number;
}

let rdkitModule: RDKitModule | null = null;
let initializationPromise: Promise<RDKitModule> | null = null;

// High-efficiency in-memory caches to prevent redundant WASM re-evaluations
const svgCache = new Map<string, string>();
const descriptorCache = new Map<string, MolecularDescriptors | null>();

export function sanitizeSmiles(raw: string): string {
  if (!raw) return "";
  let s = raw.trim();
  s = s.replace(/^[`"']+|[`"']+$/g, '');
  s = s.replace(/^(?:canonical\s+)?smiles\s*:\s*/i, '');
  s = s.replace(/\s+\(.*?\)$/, '');
  s = s.replace(/[;,. \t]+$/, '');
  s = s.replace(/\s+/g, '');
  return s;
}

export async function initRDKit(): Promise<RDKitModule> {
  if (rdkitModule) return rdkitModule;
  if (initializationPromise) return initializationPromise;

  initializationPromise = new Promise((resolve, reject) => {
    const locateFile = (filePath: string) => {
      if (filePath.endsWith(".wasm")) {
        return "/RDKit_minimal.wasm";
      }
      return filePath;
    };

    const startInit = () => {
      // @ts-ignore
      if (typeof window.initRDKitModule === "function") {
        // @ts-ignore
        window.initRDKitModule({ locateFile })
          .then((module: RDKitModule) => {
            rdkitModule = module;
            resolve(module);
          })
          .catch((err: unknown) => {
            console.warn("RDKit init with locateFile failed, trying default loader:", err);
            // @ts-ignore
            window.initRDKitModule()
              .then((module: RDKitModule) => {
                rdkitModule = module;
                resolve(module);
              })
              .catch((e: unknown) => {
                initializationPromise = null;
                reject(e);
              });
          });
      } else {
        initializationPromise = null;
        reject(new Error("initRDKitModule function not found on window"));
      }
    };

    // @ts-ignore
    if (typeof window.initRDKitModule === "function") {
      startInit();
    } else {
      // First try loading local script from /RDKit_minimal.js
      const script = document.createElement("script");
      script.src = "/RDKit_minimal.js";
      script.onload = () => startInit();
      script.onerror = () => {
        // Fallback: try to load from CDN if local script fails
        const fallbackScript = document.createElement("script");
        fallbackScript.src = "https://unpkg.com/@rdkit/rdkit/dist/RDKit_minimal.js";
        fallbackScript.onload = () => startInit();
        fallbackScript.onerror = () => {
          initializationPromise = null;
          reject(new Error("Failed to load RDKit script"));
        };
        document.head.appendChild(fallbackScript);
      };
      document.head.appendChild(script);
    }
  });

  return initializationPromise;
}

export async function validateSmiles(smiles: string): Promise<{ isValid: boolean; canonicalSmiles?: string; error?: string }> {
  try {
    const clean = sanitizeSmiles(smiles);
    if (!clean) return { isValid: false, error: "Empty chemical structure" };

    const rdkit = await initRDKit();
    let mol = rdkit.get_mol(clean);
    if (!mol || !mol.is_valid()) {
      if (mol) mol.delete();
      const relaxed = clean.replace(/[@\\/]/g, "").replace(/\(\)/g, "");
      mol = rdkit.get_mol(relaxed);
    }
    if (!mol) {
      return { isValid: false, error: "Invalid chemical structure (RDKit could not parse SMILES)" };
    }
    const isValid = mol.is_valid();
    const canonical = mol.get_smiles();
    mol.delete();
    
    if (!isValid) {
      return { isValid: false, error: "Chemical structure is invalid (Valence or bonding errors)" };
    }
    
    return { isValid: true, canonicalSmiles: canonical };
  } catch (err) {
    console.error("RDKit Validation Error:", err);
    return { isValid: false, error: "Validation engine error" };
  }
}

export async function getMoleculeSvg(smiles: string, width: number = 200, height: number = 200): Promise<string | null> {
  if (!smiles) return null;
  const cleanInputSmiles = sanitizeSmiles(smiles);
  if (!cleanInputSmiles) return null;
  const cacheKey = `${cleanInputSmiles}_${width}x${height}`;
  
  if (svgCache.has(cacheKey)) {
    return svgCache.get(cacheKey)!;
  }

  try {
    const rdkit = await initRDKit();
    let mol = rdkit.get_mol(cleanInputSmiles);
    
    // Fallback: try parsing with relaxed stereochemistry or empty branch pruning
    if (!mol || !mol.is_valid()) {
      if (mol) mol.delete();
      const relaxed = cleanInputSmiles.replace(/[@\\/]/g, "").replace(/\(\)/g, "");
      if (relaxed && relaxed !== cleanInputSmiles) {
        mol = rdkit.get_mol(relaxed);
      }
    }
    
    if (!mol) return null;
    
    const isValid = mol.is_valid();
    if (!isValid) {
      mol.delete();
      return null;
    }
    
    // RDKit minimal get_svg takes width and height parameters directly
    const rawSvg = mol.get_svg(width, height);
    mol.delete();
    
    if (rawSvg) {
      // Strip XML declaration for valid inline SVG rendering in React
      const cleanSvg = rawSvg.replace(/<\?xml[^>]*\?>/i, '').trim();
      svgCache.set(cacheKey, cleanSvg);
      return cleanSvg;
    }
    return null;
  } catch (err) {
    console.error("RDKit SVG Generation Error:", err);
    return null;
  }
}

export async function getMolecularDescriptors(smiles: string): Promise<MolecularDescriptors | null> {
  if (!smiles) return null;
  const clean = smiles.trim();
  if (descriptorCache.has(clean)) {
    return descriptorCache.get(clean)!;
  }

  try {
    const rdkit = await initRDKit();
    let mol = rdkit.get_mol(clean);
    
    // Fallback: try to canonicalize if first attempt fails
    if (!mol || !mol.is_valid()) {
      if (mol) mol.delete();
      descriptorCache.set(clean, null);
      return null;
    }
    
    const descriptorsJson = mol.get_descriptors();
    const raw = JSON.parse(descriptorsJson);
    mol.delete();
    
    // Normalize keys - extracting all requested RDKit descriptors
    const normalized: MolecularDescriptors = {
      MolWt: raw.MolWt != null ? Number(Number(raw.MolWt).toFixed(2)) : (raw.amw != null ? Number(Number(raw.amw).toFixed(2)) : undefined),
      MolLogP: raw.MolLogP != null ? Number(Number(raw.MolLogP).toFixed(2)) : (raw.CrippenClogP != null ? Number(Number(raw.CrippenClogP).toFixed(2)) : undefined),
      TPSA: raw.TPSA != null ? Number(Number(raw.TPSA).toFixed(2)) : undefined,
      NumRotatableBonds: raw.NumRotatableBonds ?? raw.numRotatableBonds,
      HBD: raw.NumHBD ?? raw.lipinskiHBD,
      HBA: raw.NumHBA ?? raw.lipinskiHBA,
      HeavyAtomCount: raw.NumHeavyAtoms ?? raw.HeavyAtomCount,
      NumAromaticRings: raw.NumAromaticRings,
      NumHeteroatoms: raw.NumHeteroatoms,
      FractionCSP3: raw.FractionCSP3 != null ? Number(Number(raw.FractionCSP3).toFixed(3)) : undefined,
    };
    
    descriptorCache.set(clean, normalized);
    return normalized;
  } catch (err) {
    console.error("RDKit Descriptors Error:", err);
    return null;
  }
}

export async function computeStrainEnergy(smiles: string): Promise<number | null> {
  try {
    const rdkit = await initRDKit();
    let mol = rdkit.get_mol(smiles);
    if (!mol || !mol.is_valid()) {
      if (mol) mol.delete();
      return null;
    }
    
    let energy = null;
    if (typeof (mol as any).add_hs === 'function') {
      (mol as any).add_hs();
    }
    
    // Some RDKit WASM builds expose an optimization or force field API
    if (typeof (mol as any).optimize_geometry === 'function') {
       const res = (mol as any).optimize_geometry(); // returns energy
       if (typeof res === 'number') energy = res;
    }

    mol.delete();
    return typeof energy === 'number' ? energy : null;
  } catch(e) {
    return null;
  }
}
