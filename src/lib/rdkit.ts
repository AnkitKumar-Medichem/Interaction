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
}

let rdkitModule: RDKitModule | null = null;
let initializationPromise: Promise<RDKitModule> | null = null;

// High-efficiency in-memory caches to prevent redundant WASM re-evaluations
const svgCache = new Map<string, string>();
const descriptorCache = new Map<string, MolecularDescriptors | null>();

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
              .catch(reject);
          });
      } else {
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
        fallbackScript.onerror = () => reject(new Error("Failed to load RDKit script"));
        document.head.appendChild(fallbackScript);
      };
      document.head.appendChild(script);
    }
  });

  return initializationPromise;
}

export async function validateSmiles(smiles: string): Promise<{ isValid: boolean; canonicalSmiles?: string; error?: string }> {
  try {
    const clean = smiles.trim();
    if (!clean) return { isValid: false, error: "Empty chemical structure" };

    const rdkit = await initRDKit();
    const mol = rdkit.get_mol(clean);
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
  const cleanInputSmiles = smiles.trim().replace(/\s+/g, '');
  const cacheKey = `${cleanInputSmiles}_${width}x${height}`;
  
  if (svgCache.has(cacheKey)) {
    return svgCache.get(cacheKey)!;
  }

  try {
    const rdkit = await initRDKit();
    const mol = rdkit.get_mol(cleanInputSmiles);
    
    if (!mol) return null;
    
    const isValid = mol.is_valid();
    if (!isValid) {
      mol.delete();
      return null;
    }
    
    // RDKit minimal get_svg takes width and height parameters directly
    const svg = mol.get_svg(width, height);
    mol.delete();
    
    if (svg) {
      svgCache.set(cacheKey, svg);
    }
    return svg;
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
    
    // Normalize keys - extracting MolWt, LogP, TPSA, and RotatableBonds
    const normalized: MolecularDescriptors = {
      MolWt: raw.MolWt ?? raw.amw ?? raw.MolWeight ?? raw.mw,
      MolLogP: raw.MolLogP ?? raw.logp ?? raw.CrippenClogP,
      TPSA: raw.TPSA ?? raw.tpsa,
      NumRotatableBonds: raw.NumRotatableBonds ?? raw.numRotatableBonds ?? raw.rotatableBonds,
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
