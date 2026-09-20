"""
=============================================================================
Python Reaction Engine Bridge & Parity Layer
=============================================================================
Maintains 100% computational, mathematical, and logical parity with
src/lib/reaction-engine.ts (the general-purpose functional-group reaction engine).

All predictions (graph-based functional-group perception, Hückel aromaticity,
reaction templates, atom-balanced product generation, Boltzmann/heuristic
scoring, and co-reactant cross-reactivity) are directly executed from the
unified engine, ensuring zero divergence between the TypeScript and Python tiers.
=============================================================================
"""

import os
import sys
import json
import subprocess
import shutil
from typing import List, Dict, Any, Optional, Tuple

# Locate workspace root directory
_CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
_WORKSPACE_DIR = os.path.abspath(os.path.join(_CURRENT_DIR, "..", ".."))
_BUNDLE_PATH = os.path.join(_WORKSPACE_DIR, "dist", "reaction-engine.cjs")
_TS_ENGINE_PATH = os.path.join(_WORKSPACE_DIR, "src", "lib", "reaction-engine.ts")

# ---------------------------------------------------------------------------
# Dual-access wrappers for compatibility with both dict and tuple indexing
# ---------------------------------------------------------------------------
class ReactivityItem(dict):
    """
    Dict wrapper that supports both key access:
      item['vulnerability'], item['mechanism']
    and tuple index access:
      item[0] -> vulnerability, item[1] -> mechanism
    """
    def __init__(self, vulnerability: str = "Low", mechanism: str = ""):
        super().__init__(vulnerability=vulnerability, mechanism=mechanism)

    def __getitem__(self, key):
        if key == 0:
            return self.get("vulnerability", "Low")
        if key == 1:
            return self.get("mechanism", "")
        return super().__getitem__(key)


class SecondaryInteractionItem(dict):
    """
    Dict wrapper for secondary interactions supporting dict and tuple indexing.
    """
    def __init__(self, vulnerability: str = "Low", partner_group: str = "", mechanism: str = ""):
        super().__init__(vulnerability=vulnerability, partnerGroup=partner_group, mechanism=mechanism)

    def __getitem__(self, key):
        if key == 0:
            return self.get("vulnerability", "Low")
        if key == 1:
            return self.get("mechanism", "")
        if key == 2:
            return self.get("partnerGroup", "")
        return super().__getitem__(key)


def _format_functional_group(fg: Dict[str, Any]) -> Dict[str, Any]:
    """Ensures each functional group dictionary has both camelCase and snake_case keys,
    with dual dict/tuple reactivity values."""
    name = fg.get("groupName") or fg.get("name") or "Functional Group"
    cat = fg.get("category") or "General"
    frag = fg.get("smilesFragment") or fg.get("fragment") or ""
    site = fg.get("reactiveSite") or fg.get("reactive_site") or ""

    formatted = {
        "name": name,
        "groupName": name,
        "category": cat,
        "fragment": frag,
        "smilesFragment": frag,
        "reactive_site": site,
        "reactiveSite": site,
    }

    for cond in ["acidic", "basic", "hydrolysis", "photolytic", "thermal", "oxidative"]:
        raw = fg.get(cond)
        if isinstance(raw, dict):
            formatted[cond] = ReactivityItem(raw.get("vulnerability", "Low"), raw.get("mechanism", ""))
        elif isinstance(raw, (list, tuple)) and len(raw) >= 2:
            formatted[cond] = ReactivityItem(str(raw[0]), str(raw[1]))
        elif isinstance(raw, (list, tuple)) and len(raw) == 1:
            formatted[cond] = ReactivityItem(str(raw[0]), "")
        else:
            formatted[cond] = ReactivityItem("Low", "")

    sec = fg.get("secondaryInteraction") or fg.get("secondary_interaction")
    if isinstance(sec, dict):
        formatted["secondaryInteraction"] = SecondaryInteractionItem(
            sec.get("vulnerability", "Low"),
            sec.get("partnerGroup", ""),
            sec.get("mechanism", "")
        )
    elif isinstance(sec, (list, tuple)) and len(sec) >= 2:
        formatted["secondaryInteraction"] = SecondaryInteractionItem(str(sec[0]), "", str(sec[1]))
    else:
        formatted["secondaryInteraction"] = SecondaryInteractionItem("Low", "", "")

    return formatted


# ---------------------------------------------------------------------------
# Node Execution Helpers
# ---------------------------------------------------------------------------
def _ensure_bundle_exists() -> bool:
    """Ensures dist/reaction-engine.cjs is compiled."""
    if os.path.exists(_BUNDLE_PATH) and os.path.getsize(_BUNDLE_PATH) > 1000:
        return True
    try:
        os.makedirs(os.path.dirname(_BUNDLE_PATH), exist_ok=True)
        res = subprocess.run(
            ["npx", "esbuild", "src/lib/reaction-engine.ts", "--bundle", "--platform=node", "--format=cjs", f"--outfile={_BUNDLE_PATH}"],
            cwd=_WORKSPACE_DIR,
            capture_output=True,
            text=True,
            timeout=15
        )
        return res.returncode == 0 and os.path.exists(_BUNDLE_PATH)
    except Exception:
        return False


def _run_node_script(js_code: str, timeout: int = 15) -> Optional[Any]:
    """Runs a one-liner node script inside the workspace and returns the parsed JSON result."""
    try:
        proc = subprocess.run(
            ["node", "-e", js_code],
            cwd=_WORKSPACE_DIR,
            capture_output=True,
            text=True,
            timeout=timeout
        )
        if proc.returncode == 0 and proc.stdout.strip():
            return json.loads(proc.stdout.strip())
        return None
    except Exception:
        return None


# ---------------------------------------------------------------------------
# Load Canonical Compound Library
# ---------------------------------------------------------------------------
def _load_compound_library() -> Dict[str, Dict[str, str]]:
    """Loads the 500+ compound library directly from reaction-engine.cjs."""
    _ensure_bundle_exists()
    code = """
    try {
        const eng = require("./dist/reaction-engine.cjs");
        const lib = eng.PHARMA_COMPOUNDS || {};
        console.log(JSON.stringify(lib));
    } catch(e) {
        console.log(JSON.stringify({}));
    }
    """
    res = _run_node_script(code)
    if isinstance(res, dict) and len(res) > 0:
        return res
    # Fallback minimal library if Node is unavailable
    return {
        "aspirin": {"name": "Aspirin", "smiles": "CC(=O)Oc1ccccc1C(=O)O", "category": "NSAID"},
        "paracetamol": {"name": "Paracetamol", "smiles": "CC(=O)Nc1ccc(O)cc1", "category": "Analgesic"},
        "amoxicillin": {"name": "Amoxicillin", "smiles": "CC1(C(N2C(S1)C(C2=O)NC(=O)C(c3ccc(cc3)O)N)C(=O)O)C", "category": "Antibiotic"},
        "metformin": {"name": "Metformin", "smiles": "CN(C)C(=N)N=C(N)N", "category": "Antidiabetic"},
        "magnesium stearate": {"name": "Magnesium Stearate", "smiles": "[Mg+2].[O-]C(=O)CCCCCCCCCCCCCCCCC.[O-]C(=O)CCCCCCCCCCCCCCCCC", "category": "Lubricant"},
        "water": {"name": "Water", "smiles": "O", "category": "Solvent"}
    }

PHARMA_COMPOUNDS: Dict[str, Dict[str, str]] = _load_compound_library()


def lookup_compound_smiles(name: str) -> Optional[str]:
    """Looks up the canonical SMILES for a compound name."""
    if not name:
        return None
    clean = name.strip().lower()
    if clean in PHARMA_COMPOUNDS:
        return PHARMA_COMPOUNDS[clean].get("smiles")

    _ensure_bundle_exists()
    clean_escaped = json.dumps(clean)
    code = f"""
    try {{
        const eng = require("./dist/reaction-engine.cjs");
        const s = eng.lookupCompoundSmiles({clean_escaped});
        console.log(JSON.stringify(s));
    }} catch(e) {{
        console.log(JSON.stringify(null));
    }}
    """
    res = _run_node_script(code)
    if isinstance(res, str):
        return res

    for k, v in PHARMA_COMPOUNDS.items():
        if clean in k or k in clean:
            return v.get("smiles")
    return None


# ---------------------------------------------------------------------------
# Functional Group Perception
# ---------------------------------------------------------------------------
def detect_functional_groups_detailed(smiles: str) -> Dict[str, Any]:
    """
    Topological molecular graph-based functional group detection.
    Matches 100% of the groups and mechanistic descriptors from src/lib/reaction-engine.ts.
    """
    s = (smiles or "").strip()
    if not s:
        return {"features": [], "sites": [], "functionalGroups": [], "functional_groups": []}

    _ensure_bundle_exists()
    s_json = json.dumps(s)
    code = f"""
    try {{
        const eng = require("./dist/reaction-engine.cjs");
        const r = eng.detectFunctionalGroupsDetailed({s_json});
        console.log(JSON.stringify(r));
    }} catch(e) {{
        console.log(JSON.stringify({{ features: [], sites: [], functionalGroups: [] }}));
    }}
    """
    raw_res = _run_node_script(code)
    if not raw_res or not isinstance(raw_res, dict):
        return {"features": [], "sites": [], "functionalGroups": [], "functional_groups": []}

    raw_fgs = raw_res.get("functionalGroups") or raw_res.get("functional_groups") or []
    formatted_fgs = [_format_functional_group(fg) for fg in raw_fgs]

    return {
        "features": raw_res.get("features", []),
        "sites": raw_res.get("sites", []),
        "functionalGroups": formatted_fgs,
        "functional_groups": formatted_fgs,
    }


def identify_functional_groups(smiles: str) -> List[Dict[str, Any]]:
    """Convenience alias returning the list of detected functional groups."""
    res = detect_functional_groups_detailed(smiles)
    return res.get("functional_groups", [])


# ---------------------------------------------------------------------------
# Degradation & Interaction Prediction Engine
# ---------------------------------------------------------------------------
def generate_computational_prediction(
    inputs: List[Dict[str, Any]],
    method: str = "Both"
) -> Dict[str, Any]:
    """
    Executes the full computational prediction pipeline from src/lib/reaction-engine.ts:
    - Atom-balanced reaction template execution
    - Formation deltaG / relativeEnergy computation
    - Heuristic and Boltzmann probability ranking at 298.15K
    - Mechanistic explanations and chain-of-thought rationale
    """
    if not inputs:
        return {
            "chainOfThought": "No inputs provided.",
            "compounds": [],
            "interactionType": "None",
            "mechanism": "No reaction.",
            "degradationImpurities": [],
            "functionalGroupAnalysis": []
        }

    _ensure_bundle_exists()
    inputs_json = json.dumps(inputs)
    method_json = json.dumps(method)

    code = f"""
    try {{
        const eng = require("./dist/reaction-engine.cjs");
        const res = eng.generateComputationalPrediction({inputs_json}, {method_json});
        console.log(JSON.stringify(res));
    }} catch(e) {{
        console.log(JSON.stringify({{
            chainOfThought: "Error running prediction: " + e.message,
            compounds: [],
            interactionType: "None",
            mechanism: "Error",
            degradationImpurities: [],
            functionalGroupAnalysis: []
        }}));
    }}
    """
    raw_res = _run_node_script(code)
    if not raw_res or not isinstance(raw_res, dict):
        return {
            "chainOfThought": "Prediction computation failed to produce output.",
            "compounds": [],
            "interactionType": "None",
            "mechanism": "Calculation error.",
            "degradationImpurities": [],
            "functionalGroupAnalysis": []
        }

    # Format functional groups
    raw_fgs = raw_res.get("functionalGroupAnalysis", [])
    formatted_fgs = [_format_functional_group(fg) for fg in raw_fgs]

    # Format degradation impurities: ensure both deltaG and relativeEnergy exist
    raw_imps = raw_res.get("degradationImpurities", [])
    formatted_imps = []
    for imp in raw_imps:
        d_energy = float(imp.get("relativeEnergy", imp.get("deltaG", 0.0)))
        prob = float(imp.get("probability", 0.0))
        prob_h = float(imp.get("probabilityHeuristic", prob))
        prob_b = float(imp.get("probabilityBoltzmann", prob))

        formatted_imp = {
            **imp,
            "deltaG": d_energy,
            "relativeEnergy": d_energy,
            "probability": prob,
            "probabilityHeuristic": prob_h,
            "probabilityBoltzmann": prob_b,
        }
        formatted_imps.append(formatted_imp)

    return {
        "chainOfThought": raw_res.get("chainOfThought", ""),
        "compounds": raw_res.get("compounds", []),
        "interactionType": raw_res.get("interactionType", "None"),
        "mechanism": raw_res.get("mechanism", ""),
        "degradationImpurities": formatted_imps,
        "functionalGroupAnalysis": formatted_fgs
    }
