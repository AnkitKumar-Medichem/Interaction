"""
Interaction - Chemical Reaction & Impurity Prediction Platform
Unified Computational Chemistry & Cheminformatics Application.
"""

import os
import re
import math
import json
import io
import base64
import urllib.parse
import datetime
import csv
import html
from typing import List, Dict, Any, Tuple

# Optional 3rd-party dependencies with resilient fallbacks
try:
    import pandas as pd
except ImportError:
    pd = None

try:
    import numpy as np
except ImportError:
    np = None

try:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
except Exception:
    plt = None

try:
    import seaborn as sns
except Exception:
    sns = None

try:
    import streamlit as st
except ImportError:
    raise RuntimeError("Streamlit is required to run this application. Please run: pip install streamlit")

# Reaction engine (RDKit + SMARTS, lives in reaction_engine.py next to this file).
# If it cannot be imported the app still starts and explains what is missing.
try:
    import reaction_engine as _rxe
    _RXE_ERROR = None
except Exception as _rxe_exc:  # pragma: no cover
    _rxe = None
    _RXE_ERROR = _rxe_exc


def safe_rerun():
    """Reruns the Streamlit application safely across different Streamlit versions."""
    if hasattr(st, "rerun"):
        st.rerun()
    elif hasattr(st, "experimental_rerun"):
        st.experimental_rerun()


# ==============================================================================
# Page Configuration & Styling
# ==============================================================================
st.set_page_config(
    page_title="Interaction - Chemical Reaction & Impurity Prediction",
    layout="wide",
    initial_sidebar_state="collapsed",
)

# Custom CSS matching the AI Studio design system
st.markdown("""
<link rel="preconnect" href="https://fonts.googleapis.com">
<link rel="preconnect" href="https://fonts.gstatic.com" crossorigin>
<link href="https://fonts.googleapis.com/css2?family=Inter:wght@400;500;600;700;800&family=Playfair+Display:ital,wght@0,600;0,700;0,800;1,600&family=JetBrains+Mono:wght@400;500;600&display=swap" rel="stylesheet">

<style>
    /* Global Typography & Reset */
    html, body, [class*="css"] {
        font-family: 'Inter', -apple-system, BlinkMacSystemFont, sans-serif;
        color: #0F172A;
    }
    
    /* Top Header */
    .app-header {
        border-bottom: 1px solid #E2E8F0;
        padding-bottom: 1.25rem;
        margin-bottom: 1.75rem;
    }
    .app-title {
        font-family: 'Playfair Display', Georgia, serif;
        font-size: 2.2rem;
        font-weight: 700;
        color: #0F172A;
        letter-spacing: -0.02em;
        line-height: 1.2;
        margin-bottom: 0.35rem;
    }
    .app-subtitle {
        color: #64748B;
        font-size: 0.95rem;
        font-weight: 400;
        line-height: 1.5;
    }

    /* Section Headings */
    .section-title {
        font-family: 'Playfair Display', Georgia, serif;
        font-size: 1.35rem;
        font-weight: 700;
        color: #0F172A;
        margin-bottom: 0.25rem;
    }
    .section-desc {
        color: #64748B;
        font-size: 0.85rem;
        margin-bottom: 1.25rem;
    }

    /* Input Card */
    .input-card {
        background: #FFFFFF;
        border: 1px solid #E2E8F0;
        border-radius: 16px;
        padding: 1.75rem;
        box-shadow: 0 1px 3px rgba(0,0,0,0.02);
        margin-bottom: 1.5rem;
    }

    /* Primary & Secondary Compound Display Cards */
    .ap1-comp-card {
        background: #FFFFFF;
        border: 1px solid #E2E8F0;
        border-radius: 14px;
        overflow: hidden;
        display: flex;
        flex-direction: row;
        margin-bottom: 1rem;
        width: 100%;
        box-shadow: 0 1px 3px rgba(0,0,0,0.02);
    }
    .ap1-comp-mol {
        width: 240px;
        min-width: 240px;
        background: #FFFFFF;
        border-right: 1px solid #F1F5F9;
        display: flex;
        align-items: center;
        justify-content: center;
        padding: 1rem;
        position: relative;
    }
    .ap1-comp-badge {
        position: absolute;
        top: 10px;
        left: 10px;
        background: #EFF6FF;
        color: #2563EB;
        font-size: 0.65rem;
        font-weight: 700;
        padding: 0.2rem 0.55rem;
        border-radius: 4px;
        text-transform: uppercase;
        letter-spacing: 0.05em;
    }
    .ap1-comp-info {
        padding: 1.25rem 1.5rem;
        flex: 1;
        display: flex;
        flex-direction: column;
        justify-content: center;
    }
    .ap1-comp-name {
        font-family: 'Playfair Display', Georgia, serif;
        font-size: 1.35rem;
        font-weight: 700;
        color: #0F172A;
        line-height: 1.2;
    }
    .ap1-comp-role {
        display: inline-block;
        font-size: 0.65rem;
        font-weight: 700;
        text-transform: uppercase;
        letter-spacing: 0.05em;
        padding: 0.2rem 0.55rem;
        border-radius: 4px;
        margin-left: 0.5rem;
    }
    .role-primary {
        background: #EFF6FF;
        color: #1D4ED8;
    }
    .role-secondary {
        background: #F1F5F9;
        color: #475569;
    }
    .ap1-smiles-box {
        font-family: 'JetBrains Mono', monospace;
        font-size: 0.75rem;
        color: #64748B;
        background: #F8FAFC;
        border: 1px solid #F1F5F9;
        padding: 0.35rem 0.65rem;
        border-radius: 6px;
        margin: 0.5rem 0 0.75rem 0;
        word-break: break-all;
    }
    .ap1-tag-group {
        display: flex;
        flex-wrap: wrap;
        gap: 0.35rem;
        margin-top: 0.25rem;
    }
    .ap1-pill {
        font-size: 0.7rem;
        font-weight: 500;
        padding: 0.2rem 0.6rem;
        border-radius: 9999px;
        border: 1px solid #E2E8F0;
        background: #FFFFFF;
        color: #475569;
    }
    .ap1-pill.mw {
        font-family: 'JetBrains Mono', monospace;
        font-weight: 600;
        background: #F1F5F9;
        color: #334155;
        border-color: #E2E8F0;
    }
    .ap1-pill.site {
        background: #EFF6FF;
        border-color: #DBEAFE;
        color: #1D4ED8;
        font-weight: 600;
    }

    /* Impurity Card Presentation */
    .ap1-imp-card {
        background: #FFFFFF;
        border: 1px solid #E2E8F0;
        border-radius: 14px;
        overflow: hidden;
        display: flex;
        flex-direction: row;
        margin-bottom: 1.25rem;
        width: 100%;
        box-shadow: 0 1px 3px rgba(0,0,0,0.02);
    }
    .ap1-imp-svg {
        width: 250px;
        min-width: 250px;
        background: #FFFFFF;
        border-right: 1px solid #F1F5F9;
        display: flex;
        align-items: center;
        justify-content: center;
        padding: 1.25rem;
        position: relative;
    }
    .ap1-imp-idx-badge {
        position: absolute;
        top: 12px;
        left: 12px;
        background: #F1F5F9;
        color: #475569;
        font-size: 0.75rem;
        font-weight: 800;
        padding: 0.15rem 0.5rem;
        border-radius: 4px;
    }
    .ap1-imp-body {
        padding: 1.5rem 1.75rem;
        flex: 1;
    }
    .ap1-imp-header {
        display: flex;
        justify-content: space-between;
        align-items: flex-start;
        gap: 1rem;
        margin-bottom: 0.5rem;
    }
    .ap1-imp-title {
        font-size: 1.15rem;
        font-weight: 700;
        color: #0F172A;
        line-height: 1.25;
    }
    .ap1-imp-prob-val {
        font-size: 1.45rem;
        font-weight: 800;
        color: #B91C1C;
        text-align: right;
        line-height: 1;
    }
    .ap1-imp-prob-sub {
        font-size: 0.68rem;
        font-weight: 600;
        color: #94A3B8;
        text-transform: uppercase;
        letter-spacing: 0.05em;
        text-align: right;
        margin-top: 0.25rem;
    }
    .ap1-prob-bar-bg {
        width: 100%;
        height: 6px;
        background: #F1F5F9;
        border-radius: 9999px;
        overflow: hidden;
        margin: 0.65rem 0 1rem 0;
    }
    .ap1-prob-bar-fill {
        height: 100%;
        background: linear-gradient(90deg, #3B82F6, #EA580C, #B91C1C);
        border-radius: 9999px;
    }
    .ap1-mech-box {
        background: #F8FAFC;
        border: 1px solid #E2E8F0;
        border-radius: 8px;
        padding: 0.85rem 1rem;
        font-size: 0.82rem;
        color: #334155;
        line-height: 1.55;
        margin-bottom: 1rem;
    }
    .ap1-mech-title {
        font-weight: 700;
        color: #0F172A;
        margin-bottom: 0.25rem;
    }
    .ap1-badge-cond {
        font-size: 0.68rem;
        font-weight: 600;
        text-transform: uppercase;
        letter-spacing: 0.05em;
        padding: 0.2rem 0.55rem;
        border-radius: 4px;
    }
    .cond-hydro { background: #EFF6FF; color: #1D4ED8; border: 1px solid #DBEAFE; }
    .cond-oxid { background: #FEF2F2; color: #B91C1C; border: 1px solid #FEE2E2; }
    .cond-therm { background: #FFF7ED; color: #C2410C; border: 1px solid #FFEDD5; }
    .cond-photo { background: #FEFCE8; color: #A16207; border: 1px solid #FEF08A; }
    .cond-react { background: #FAF5FF; color: #7E22CE; border: 1px solid #F3E8FF; }

    /* Vulnerability Badges */
    .vuln-critical { background: #FFF1F2; color: #BE123C; border: 1px solid #FECDD3; font-weight: 600; font-size: 0.72rem; padding: 2px 8px; border-radius: 4px; }
    .vuln-high { background: #FFFBEB; color: #B45309; border: 1px solid #FDE68A; font-weight: 600; font-size: 0.72rem; padding: 2px 8px; border-radius: 4px; }
    .vuln-moderate { background: #FEFCE8; color: #854D0E; border: 1px solid #FEF08A; font-weight: 600; font-size: 0.72rem; padding: 2px 8px; border-radius: 4px; }
    .vuln-low { background: #EFF6FF; color: #1D4ED8; border: 1px solid #DBEAFE; font-weight: 600; font-size: 0.72rem; padding: 2px 8px; border-radius: 4px; }
    .vuln-resistant { background: #F8FAFC; color: #475569; border: 1px solid #E2E8F0; font-weight: 600; font-size: 0.72rem; padding: 2px 8px; border-radius: 4px; }

    /* Condition Grid Card */
    .cond-card {
        background: #FFFFFF;
        border: 1px solid #E2E8F0;
        border-radius: 8px;
        padding: 0.75rem;
        height: 100%;
    }
    .cond-card-title {
        font-size: 0.75rem;
        font-weight: 700;
        color: #0F172A;
        display: flex;
        justify-content: space-between;
        align-items: center;
        margin-bottom: 0.35rem;
    }
    .cond-card-desc {
        font-size: 0.7rem;
        color: #475569;
        line-height: 1.45;
    }

    /* Mechanistic Framework & Disclaimer */
    .ap1-cot-box {
        background: #F8FAFC;
        border: 1px solid #E2E8F0;
        border-radius: 12px;
        padding: 1.25rem 1.5rem;
        font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif;
        font-size: 0.85rem;
        color: #334155;
        line-height: 1.65;
        white-space: pre-wrap;
        margin-bottom: 1rem;
    }
    .ap1-cot-note {
        background: #F8FAFC;
        border: 1px solid #E2E8F0;
        border-radius: 10px;
        padding: 1rem 1.25rem;
        font-size: 0.78rem;
        color: #475569;
        line-height: 1.55;
    }
    .ap1-disclaimer-box {
        background: #FAFAFA;
        border: 1px solid #E2E8F0;
        border-radius: 10px;
        padding: 1rem 1.25rem;
        font-size: 0.78rem;
        color: #64748B;
        font-style: italic;
        line-height: 1.55;
        margin-top: 1.75rem;
    }

    /* Responsive adjustments */
    @media (max-width: 768px) {
        .ap1-comp-card, .ap1-imp-card {
            flex-direction: column !important;
        }
        .ap1-comp-mol, .ap1-imp-svg {
            width: 100% !important;
            min-width: 100% !important;
            border-right: none !important;
            border-bottom: 1px solid #F1F5F9 !important;
        }
    }
</style>
""", unsafe_allow_html=True)

# ==============================================================================
# Clean HTML Render Utility (Prevents Markdown Code Block Indentation Glitches)
# ==============================================================================
def render_html(html_str: str):
    """
    Renders an HTML snippet cleanly in Streamlit without risk of Markdown
    parsers treating indented lines or blank lines as code blocks (<pre><code>).
    """
    clean_lines = [line.strip() for line in html_str.splitlines() if line.strip()]
    st.markdown(" ".join(clean_lines), unsafe_allow_html=True)

# ==============================================================================
# Chemical Structure Rendering & Descriptors Engine
# ==============================================================================
def sanitize_smiles_py(raw: Any) -> str:
    """Sanitizes raw chemical SMILES input with strict type safety."""
    if not raw or not isinstance(raw, str):
        return ""
    s = raw.strip()
    s = re.sub(r'^[`"\'\s]+|[`"\'\s]+$', '', s)
    s = re.sub(r'^(?:canonical\s+)?smiles\s*:\s*', '', s, flags=re.I)
    s = re.sub(r'\s*\(.*?\)$', '', s)
    s = re.sub(r'[;,. \t]+$', '', s)
    s = re.sub(r'\s+', '', s)
    return s

def get_chemical_structure_img(smiles: str, width: int = 240, height: int = 200) -> str:
    """
    Renders a crisp 2D chemical structure image.
    Uses native cheminformatics drawing if available, with resilient remote fallbacks (PubChem then Cactus).
    """
    clean = sanitize_smiles_py(smiles)
    if not clean:
        return ""

    # 1. Try local cheminformatics engine (RDKit vector SVG)
    try:
        from rdkit import Chem
        from rdkit.Chem import rdDepictor
        from rdkit.Chem.Draw import rdMolDraw2D
        mol = Chem.MolFromSmiles(clean)
        if mol is None:
            relaxed = re.sub(r'[@\\/]', '', clean)
            relaxed = re.sub(r'\(\)', '', relaxed)
            if relaxed and relaxed != clean:
                mol = Chem.MolFromSmiles(relaxed)
        if mol is not None:
            try:
                rdDepictor.Compute2DCoords(mol)
            except Exception:
                pass
            drawer = rdMolDraw2D.MolDraw2DSVG(width, height)
            opts = drawer.drawOptions()
            opts.clearBackground = True
            opts.padding = 0.05
            drawer.DrawMolecule(mol)
            drawer.FinishDrawing()
            svg = drawer.GetDrawingText()
            b64 = base64.b64encode(svg.encode("utf-8")).decode("utf-8")
            return f"data:image/svg+xml;base64,{b64}"
    except Exception:
        pass

    # 2. Resilient fallback to PubChem depiction API
    try:
        encoded = urllib.parse.quote(clean)
        return f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/{encoded}/PNG?record_type=2d&image_size={width}x{height}"
    except Exception:
        pass

    # 3. Secondary fallback to Cactus NCI depiction
    try:
        encoded = urllib.parse.quote(clean)
        return f"https://cactus.nci.nih.gov/chemical/structure/{encoded}/image"
    except Exception:
        return ""

def get_molecular_descriptors(smiles: str) -> Dict[str, Any]:
    """
    Calculates molecular descriptors using RDKit:
    logP, TPSA, HBD, HBA, NumRotatableBonds, HeavyAtomCount,
    NumAromaticRings, NumHeteroatoms, FractionCSP3, and MolWt.
    Includes lightweight formula approximation fallback if RDKit is unavailable.
    """
    clean = sanitize_smiles_py(smiles)
    if not clean:
        return {}
    try:
        from rdkit import Chem
        from rdkit.Chem import Descriptors, rdMolDescriptors
        mol = Chem.MolFromSmiles(clean)
        if mol is None:
            relaxed = re.sub(r'[@\\/]', '', clean)
            relaxed = re.sub(r'\(\)', '', relaxed)
            if relaxed and relaxed != clean:
                mol = Chem.MolFromSmiles(relaxed)
        if mol is not None:
            return {
                "mw": round(float(Descriptors.MolWt(mol)), 2),
                "logp": round(float(Descriptors.MolLogP(mol)), 2),
                "tpsa": round(float(Descriptors.TPSA(mol)), 2),
                "hbd": int(rdMolDescriptors.CalcNumHBD(mol)),
                "hba": int(rdMolDescriptors.CalcNumHBA(mol)),
                "rotatable_bonds": int(Descriptors.NumRotatableBonds(mol)),
                "heavy_atom_count": int(Descriptors.HeavyAtomCount(mol)),
                "aromatic_rings": int(Descriptors.NumAromaticRings(mol)),
                "heteroatoms": int(Descriptors.NumHeteroatoms(mol)),
                "fraction_csp3": round(float(Descriptors.FractionCSP3(mol)), 3)
            }
    except Exception:
        pass

    # Lightweight heuristic fallback if RDKit is not installed
    try:
        heavy_atoms = len(re.findall(r'[A-Za-z]', re.sub(r'\[H\]|[Hh]', '', clean)))
        aromatic_count = 1 if ('c1' in clean or 'c2' in clean) else 0
        rot_bonds = max(0, len(re.findall(r'CC|CO|CN|CS', clean)) - 1)
        hba = len(re.findall(r'[O|N|o|n]', clean))
        hbd = len(re.findall(r'O(?![C|c])|N(?![C|c])|OH|NH', clean))
        approx_mw = round(heavy_atoms * 13.5, 2)
        if heavy_atoms > 0:
            return {
                "mw": approx_mw,
                "logp": 1.5,
                "tpsa": round(hba * 18.0, 1),
                "hbd": hbd,
                "hba": hba,
                "rotatable_bonds": rot_bonds,
                "heavy_atom_count": heavy_atoms,
                "aromatic_rings": aromatic_count,
                "heteroatoms": hba,
                "fraction_csp3": 0.35
            }
    except Exception:
        pass
    return {}

def format_descriptor_pills(desc: Dict[str, Any]) -> str:
    """Formats calculated RDKit molecular descriptors into styled HTML pills."""
    if not desc:
        return ""
    pills = []
    if "mw" in desc and desc["mw"] is not None:
        pills.append(f'<span class="ap1-pill mw" title="Molecular Weight (g/mol)">MW: {desc["mw"]:.2f} g/mol</span>')
    if "logp" in desc and desc["logp"] is not None:
        pills.append(f'<span class="ap1-pill" title="Partition Coefficient (LogP)">LogP: {desc["logp"]}</span>')
    if "tpsa" in desc and desc["tpsa"] is not None:
        pills.append(f'<span class="ap1-pill" title="Topological Polar Surface Area (sq A)">TPSA: {desc["tpsa"]} sq A</span>')
    if "hbd" in desc and desc["hbd"] is not None:
        pills.append(f'<span class="ap1-pill" title="Hydrogen Bond Donors">HBD: {desc["hbd"]}</span>')
    if "hba" in desc and desc["hba"] is not None:
        pills.append(f'<span class="ap1-pill" title="Hydrogen Bond Acceptors">HBA: {desc["hba"]}</span>')
    if "rotatable_bonds" in desc and desc["rotatable_bonds"] is not None:
        pills.append(f'<span class="ap1-pill" title="Number of Rotatable Bonds">RotB: {desc["rotatable_bonds"]}</span>')
    if "heavy_atom_count" in desc and desc["heavy_atom_count"] is not None:
        pills.append(f'<span class="ap1-pill" title="Heavy Atom Count">HeavyAtoms: {desc["heavy_atom_count"]}</span>')
    if "aromatic_rings" in desc and desc["aromatic_rings"] is not None:
        pills.append(f'<span class="ap1-pill" title="Number of Aromatic Rings">AromRings: {desc["aromatic_rings"]}</span>')
    if "heteroatoms" in desc and desc["heteroatoms"] is not None:
        pills.append(f'<span class="ap1-pill" title="Number of Heteroatoms">Heteroatoms: {desc["heteroatoms"]}</span>')
    if "fraction_csp3" in desc and desc["fraction_csp3"] is not None:
        pills.append(f'<span class="ap1-pill" title="Fraction of sp3 Carbons (Fsp3)">Fsp3: {desc["fraction_csp3"]}</span>')
    return "".join(pills)

# ==============================================================================
# Persistent CSV Logbook Configuration (Rolling 100 Queries)
# ==============================================================================
LOGBOOK_FILE = "query_logbook.csv"
LOGBOOK_COLUMNS = ["timestamp", "primary_compound_smiles", "secondary_compounds_smiles", "predicted_impurities"]

def load_logbook():
    """Loads the query logbook into a DataFrame or list of dicts with date-only timestamps."""
    if os.path.exists(LOGBOOK_FILE):
        try:
            if pd is not None:
                df = pd.read_csv(LOGBOOK_FILE)
                for col in LOGBOOK_COLUMNS:
                    if col not in df.columns:
                        df[col] = ""
                if "timestamp" in df.columns:
                    df["timestamp"] = df["timestamp"].astype(str).apply(lambda x: x.split(" ")[0].split("T")[0] if x and str(x).lower() != 'nan' else "")
                return df
            else:
                rows = []
                with open(LOGBOOK_FILE, mode='r', encoding='utf-8', errors='ignore') as f:
                    reader = csv.DictReader(f)
                    for r in reader:
                        if "timestamp" in r and r["timestamp"]:
                            r["timestamp"] = str(r["timestamp"]).split(" ")[0].split("T")[0]
                        rows.append(r)
                return rows
        except Exception:
            return pd.DataFrame(columns=LOGBOOK_COLUMNS) if pd is not None else []
    return pd.DataFrame(columns=LOGBOOK_COLUMNS) if pd is not None else []

def append_to_logbook(primary_smiles: str, secondary_smiles_list: List[str], impurities: List[Dict[str, Any]]):
    """Appends a new prediction query entry into the persistent rolling CSV logbook with date-only timestamp."""
    try:
        imp_summary = "; ".join([f"{i.get('iupacName', 'Unknown')} ({i.get('smiles', '')}) [{(i.get('probability', 0)*100):.1f}%]" for i in impurities[:5]])
        new_entry = {
            "timestamp": datetime.date.today().strftime("%Y-%m-%d"),
            "primary_compound_smiles": primary_smiles,
            "secondary_compounds_smiles": "; ".join([s for s in secondary_smiles_list if s.strip()]),
            "predicted_impurities": imp_summary
        }
        if pd is not None:
            df = load_logbook()
            if isinstance(df, pd.DataFrame):
                df = pd.concat([pd.DataFrame([new_entry]), df], ignore_index=True)
                df = df.head(100)
                df.to_csv(LOGBOOK_FILE, index=False)
                return

        # Pure-Python fallback using built-in csv
        existing = []
        if os.path.exists(LOGBOOK_FILE):
            with open(LOGBOOK_FILE, mode='r', encoding='utf-8', errors='ignore') as f:
                reader = csv.DictReader(f)
                for r in reader:
                    existing.append(r)
        all_entries = [new_entry] + existing
        all_entries = all_entries[:100]
        with open(LOGBOOK_FILE, mode='w', encoding='utf-8', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=LOGBOOK_COLUMNS)
            writer.writeheader()
            writer.writerows(all_entries)
    except Exception:
        pass


# ==============================================================================
# Functional Group Identification Engine
# ==============================================================================
DEFAULT_CONDITIONS = ["Acidic", "Basic", "Hydrolysis", "Photolysis", "Thermal", "Oxidative"]

def identify_functional_groups(smiles: str) -> List[Dict[str, Any]]:
    """
    Identifies functional groups and maps mechanistic reactivity across:
    Acidic, Basic, Hydrolysis, Photolytic, Thermal, Oxidative conditions (plus cross-interaction).
    Graph-based SMARTS perception in reaction_engine.py (independent of how a SMILES is written);
    returns [] for input that is not a valid molecule instead of guessing.
    """
    s = sanitize_smiles_py(smiles)
    if not s or _rxe is None:
        return []
    try:
        return _rxe.identify_functional_groups(s)
    except Exception:
        return []

# ==============================================================================
# Heatmap Plotter (Pure Publication-Quality Matrix)
# ==============================================================================
def plot_heatmap(matrix, row_labels: List[str], col_labels: List[str], title: str):
    """
    Generates a publication-quality vulnerability matrix using Matplotlib/Seaborn.
    Returns None if graphics libraries are unavailable or if rendering fails.
    """
    if pd is None or sns is None or plt is None:
        return None

    fig = None
    try:
        display_rows = [r if len(r) <= 35 else r[:32] + "..." for r in row_labels]
        df = pd.DataFrame(matrix, index=display_rows, columns=col_labels)

        n_rows = len(display_rows)
        n_cols = len(col_labels)
        fig_width = max(8.5, n_cols * 1.6 + 2.0)
        fig_height = max(4.6, n_rows * 0.7 + 1.8)

        fig, ax = plt.subplots(figsize=(fig_width, fig_height), dpi=150)
        fig.patch.set_facecolor('#FFFFFF')
        ax.set_facecolor('#F8FAFC')

        sns.heatmap(
            df,
            annot=False,
            cmap="coolwarm",
            vmin=0.0,
            vmax=1.0,
            cbar_kws={'label': 'Degradation / Incompatibility Potential', 'shrink': 0.85},
            linewidths=2.0,
            linecolor='#FFFFFF',
            square=False,
            ax=ax
        )

        ax.set_title(title, fontsize=13, fontweight='bold', pad=18, color='#0F172A')
        ax.set_xticklabels(ax.get_xticklabels(), rotation=15 if n_cols > 3 else 0, ha='center', fontsize=9.5, fontweight='600', color='#334155')
        ax.set_yticklabels(ax.get_yticklabels(), rotation=0, fontsize=9.5, fontweight='600', color='#334155')

        plt.tight_layout()
        return fig
    except Exception:
        if fig is not None:
            try:
                plt.close(fig)
            except Exception:
                pass
        return None

def render_html_heatmap(matrix, row_labels: List[str], col_labels: List[str], title: str):
    """
    Bulletproof pure-HTML heatmap fallback. Never fails, requires 0 external dependencies.
    """
    def get_color_style(val: float) -> Tuple[str, str, str]:
        if val >= 0.85:
            return "#FEE2E2", "#991B1B", "Critical"
        elif val >= 0.65:
            return "#FFEDD5", "#9A3412", "High"
        elif val >= 0.35:
            return "#FEF9C3", "#854D0E", "Moderate"
        elif val >= 0.15:
            return "#E0F2FE", "#075985", "Low"
        return "#F8FAFC", "#475569", "Resistant"

    header_cols = "".join([f'<th style="padding: 10px 14px; background: #F1F5F9; color: #1E293B; font-size: 0.85rem; font-weight: 700; border: 1px solid #CBD5E1; text-align: center;">{c}</th>' for c in col_labels])
    
    rows_html = []
    for r_idx, r_name in enumerate(row_labels):
        cells = []
        for c_idx in range(len(col_labels)):
            val = float(matrix[r_idx][c_idx]) if (r_idx < len(matrix) and c_idx < len(matrix[r_idx])) else 0.15
            bg, text_color, label = get_color_style(val)
            cells.append(f'<td style="padding: 10px 14px; background: {bg}; color: {text_color}; font-size: 0.8rem; font-weight: 600; text-align: center; border: 1px solid #E2E8F0;">{label} ({int(val*100)}%)</td>')
        rows_html.append(f'<tr><td style="padding: 10px 14px; background: #F8FAFC; color: #0F172A; font-weight: 700; font-size: 0.85rem; border: 1px solid #CBD5E1; white-space: nowrap;">{r_name}</td>{"".join(cells)}</tr>')

    table_html = f"""
    <div style="background: #FFFFFF; border: 1px solid #E2E8F0; border-radius: 12px; padding: 1.25rem; margin: 1rem 0; overflow-x: auto; box-shadow: 0 1px 3px rgba(0,0,0,0.02);">
        <div style="font-weight: 700; color: #0F172A; font-size: 1.05rem; margin-bottom: 0.75rem;">{title}</div>
        <table style="width: 100%; border-collapse: collapse; font-family: \'Inter\', -apple-system, sans-serif;">
            <thead>
                <tr>
                    <th style="padding: 10px 14px; background: #F1F5F9; color: #1E293B; font-size: 0.85rem; font-weight: 700; border: 1px solid #CBD5E1; text-align: left;">Stress Condition</th>
                    {header_cols}
                </tr>
            </thead>
            <tbody>
                {"".join(rows_html)}
            </tbody>
        </table>
    </div>
    """
    render_html(table_html)


# ==============================================================================
# Reaction & Degradation Prediction Engine
# ==============================================================================
def predict_degradation_and_reactions(
    primary_smiles: str,
    secondary_smiles_list: List[str],
    method: str = "Both"
) -> Dict[str, Any]:
    """
    Calculates degradation products, free energies (Delta G), Boltzmann & Heuristic probabilities
    from the functional groups of the primary compound and of any co-reactants.
    Chemistry is done by reaction_engine.py (graph-based, RDKit); this function keeps the
    dictionary contract used by the rest of the UI (impurities, heatmap, chain of thought).
    """
    clean_primary = sanitize_smiles_py(primary_smiles)
    clean_secondaries = [sanitize_smiles_py(s) for s in secondary_smiles_list if s and sanitize_smiles_py(s)]

    def _absent_result(message: str, chain: str) -> Dict[str, Any]:
        return {
            "functional_groups": [],
            "impurities": [{
                "iupacName": "Reactive functional group is absent",
                "smiles": "",
                "condition": "Hydrolysis",
                "source": "Stress degradation",
                "mechanismExplanation": message,
                "deltaG": 0.0,
                "kineticLikelihood": 0.0,
                "probability": 0.0,
                "probabilityBoltzmann": 0.0,
                "probabilityHeuristic": 0.0
            }],
            "heatmap_matrix": [[0.15]*1 for _ in range(6)],
            "row_labels": ["Acidic", "Basic", "Hydrolysis", "Photolysis", "Thermal", "Oxidative"],
            "col_labels": ["Aliphatic Framework"],
            "chain_of_thought": chain
        }

    if not clean_primary:
        return _absent_result("Reactive functional group is absent. Please enter a valid molecular SMILES string.",
                              "No primary compound provided.")

    if _rxe is None:
        note = ("The reaction engine module could not be loaded (%s). Place reaction_engine.py next to app.py "
                "and make sure RDKit is installed (pip install rdkit)." % (_RXE_ERROR,))
        return _absent_result(note, note)

    engine = _rxe.predict(clean_primary, clean_secondaries, method)
    if not engine.get("ok"):
        return _absent_result(engine.get("reason", "The primary SMILES could not be parsed."),
                              engine.get("reason", "The primary SMILES could not be parsed."))

    p_groups = engine["functional_groups"]
    top_5 = engine["impurities"]
    if not top_5:
        top_5 = [{
            "iupacName": "Reactive functional group is absent",
            "smiles": "",
            "condition": "Hydrolysis",
            "source": "Stress degradation",
            "mechanismExplanation": "Reactive functional group is absent. " + engine.get("mechanism", ""),
            "deltaG": 0.0,
            "kineticLikelihood": 0.0,
            "probability": 0.0,
            "probabilityBoltzmann": 0.0,
            "probabilityHeuristic": 0.0
        }]


    # Build Heatmap matrix with functional groups on X-axis (col_labels) and conditions on Y-axis (row_labels)
    def clean_fg_name(name_str: str) -> str:
        clean = re.sub(r"\s*\(.*?\)", "", name_str)
        clean = re.sub(r"\s*&.*$", "", clean)
        clean = re.sub(r"^\[.*?\]\s*", "", clean)
        clean = re.sub(r"\s*(?:<->|\u2194).*$", "", clean)
        clean = re.sub(r"Reactive Center", "Aliphatic Center", clean, flags=re.I)
        return clean.strip()

    col_labels = []
    # Primary compound functional groups
    for g in p_groups:
        c_name = clean_fg_name(g.get("name", "Functional Group"))
        if c_name and c_name not in col_labels:
            col_labels.append(c_name)

    # Secondary compound functional groups
    for sec_s in secondary_smiles_list:
        if sec_s and sec_s.strip():
            for sg in identify_functional_groups(sec_s.strip()):
                c_name = clean_fg_name(sg.get("name", ""))
                if c_name and c_name not in col_labels:
                    col_labels.append(c_name)

    if len(col_labels) == 0:
        col_labels = ["Aliphatic Framework"]

    # Conditions on Y-axis
    row_labels = ["Acidic", "Basic", "Hydrolysis", "Photolysis", "Thermal", "Oxidative"]
    vuln_map = {"Critical": 0.92, "High": 0.75, "Moderate": 0.45, "Low": 0.20, "Resistant": 0.05}

    # Map each functional group to vulnerability ratings
    cond_keys = ["acidic", "basic", "hydrolysis", "photolytic", "thermal", "oxidative"]
    fg_dict_lookup = {}
    for g in p_groups:
        c_name = clean_fg_name(g.get("name", ""))
        fg_dict_lookup[c_name] = g
    for sec_s in secondary_smiles_list:
        if sec_s and sec_s.strip():
            for sg in identify_functional_groups(sec_s.strip()):
                c_name = clean_fg_name(sg.get("name", ""))
                if c_name not in fg_dict_lookup:
                    fg_dict_lookup[c_name] = sg

    matrix = []
    for cond_key in cond_keys:
        row_vals = []
        for fg_col in col_labels:
            g_obj = fg_dict_lookup.get(fg_col)
            if g_obj and cond_key in g_obj:
                score = vuln_map.get(g_obj[cond_key][0], 0.20)
            else:
                score = 0.15
            row_vals.append(score)
        matrix.append(row_vals)


    chain_of_thought = engine["chain_of_thought"]

    return {
        "functional_groups": p_groups,
        "impurities": top_5,
        "heatmap_matrix": np.array(matrix) if np is not None else matrix,
        "row_labels": row_labels,
        "col_labels": col_labels,
        "chain_of_thought": chain_of_thought
    }

# ==============================================================================
# Main Application UI Header
# ==============================================================================
render_html("""
<div class="app-header">
    <div class="app-title">Interaction</div>
    <div class="app-subtitle">Computational Chemical Reaction & Impurity Prediction Platform</div>
</div>
""")

# Top Navigation Tabs matching AI Studio
tab_predict, tab_logbook = st.tabs([
    "Analysis & Predictions",
    "Query Logbook (100 Queries)"
])

# ==============================================================================
# TAB 1: Analysis & Predictions
# ==============================================================================
with tab_predict:
    # Input Chemical Data Form
    with st.container():
        st.markdown('<div class="section-title">Input Molecular Structures (SMILES Only)</div>', unsafe_allow_html=True)
        st.markdown('<div class="section-desc">Enter canonical SMILES representations. Calculations evaluate functional group reactivity, stress degradation pathways, and intermolecular incompatibilities.</div>', unsafe_allow_html=True)

        col_input1, col_input2 = st.columns([1.2, 1], gap="large")

        with col_input1:
            primary_smiles = st.text_input(
                "Primary Compound (SMILES) *",
                value="",
                placeholder="Enter canonical SMILES string"
            )

            if "num_secondary" not in st.session_state:
                st.session_state.num_secondary = 1

            sec_smiles_list = []
            st.markdown("<div style='font-size: 0.85rem; font-weight: 700; color: #334155; margin-top: 0.75rem; margin-bottom: 0.35rem;'>Secondary Compounds (Optional Co-reactants & Excipients)</div>", unsafe_allow_html=True)
            for i in range(st.session_state.num_secondary):
                sec_val = st.text_input(
                    f"Secondary Compound {i+1} (SMILES)",
                    key=f"sec_smiles_{i}",
                    placeholder="Enter secondary compound SMILES",
                    help=f"Co-reactant, excipient, or secondary ingredient {i+1} SMILES."
                )
                sec_smiles_list.append(sec_val)

            col_add, col_rem = st.columns(2)
            with col_add:
                if st.button("Add Secondary Compound (SMILES)", disabled=st.session_state.num_secondary >= 4):
                    st.session_state.num_secondary += 1
                    safe_rerun()
            with col_rem:
                if st.button("Remove Secondary Compound", disabled=st.session_state.num_secondary <= 1):
                    st.session_state.num_secondary -= 1
                    safe_rerun()

        with col_input2:
            st.markdown("<div style='font-size: 0.85rem; font-weight: 700; color: #334155; margin-bottom: 0.35rem;'>Prediction Engine & Methodology</div>", unsafe_allow_html=True)
            method_choice = st.radio(
                "Select Prediction Framework",
                [
                    "Dual Engine (Heuristic Kinetic Rules + Boltzmann Thermodynamic Delta G)",
                    "Heuristic (Expert Kinetic Activation & Transition States)",
                    "Boltzmann (Thermodynamic Free Energy Delta G Distribution at 298.15K)"
                ],
                index=0,
                label_visibility="collapsed"
            )

            # Map choice to short method key
            if "Heuristic" in method_choice and "Boltzmann" not in method_choice:
                method_key = "Heuristic"
            elif "Boltzmann" in method_choice and "Heuristic" not in method_choice:
                method_key = "Boltzmann"
            else:
                method_key = "Both"

            st.markdown("<div style='height: 1.5rem;'></div>", unsafe_allow_html=True)
            predict_btn = st.button("Predict Chemical Interactions", type="primary", use_container_width=True)

    # Perform Analysis on Click
    if predict_btn:
        clean_primary = sanitize_smiles_py(primary_smiles)
        if not clean_primary:
            st.error("Please provide a valid Primary Compound SMILES string before starting prediction.")
        else:
            with st.spinner("Analyzing functional groups and calculating condition reactivity..."):
                try:
                    calc_results = predict_degradation_and_reactions(primary_smiles, sec_smiles_list, method_key)
                    append_to_logbook(primary_smiles, sec_smiles_list, calc_results["impurities"])
                    st.session_state.last_results = calc_results
                    st.session_state.last_primary = primary_smiles
                    st.session_state.last_secondary = [s for s in sec_smiles_list if s.strip()]
                except Exception as ex:
                    st.error(f"Prediction calculation error: {ex}")


    # Display Results Dashboard matching AI Studio
    if "last_results" in st.session_state:
        res = st.session_state.last_results
        cur_primary = st.session_state.last_primary
        cur_secondary = st.session_state.get("last_secondary", [])

        st.markdown("<hr style='border: none; border-top: 1px solid #E2E8F0; margin: 2rem 0;'/>", unsafe_allow_html=True)

        # ----------------------------------------------------------------------
        # 1. Input Chemical Data Card (With 2D Chemical Structure!)
        # ----------------------------------------------------------------------
        st.markdown('<div class="section-title">1. Input Chemical Data</div>', unsafe_allow_html=True)
        st.markdown('<div class="section-desc">Calculated molecular descriptors, functional group features, and 2D chemical structure diagrams.</div>', unsafe_allow_html=True)

        # Primary Compound Card
        primary_img = get_chemical_structure_img(cur_primary, width=220, height=200)
        primary_desc = get_molecular_descriptors(cur_primary)
        primary_desc_pills = format_descriptor_pills(primary_desc)

        p_groups = res.get("functional_groups", [])
        reactive_centers = list(set([g["reactive_site"] for g in p_groups if "reactive_site" in g]))

        primary_mol_html = f'<img src="{primary_img}" alt="Primary Compound" style="max-width: 100%; max-height: 180px; object-fit: contain;"/>' if primary_img else '<div style="color: #94A3B8; font-size: 0.75rem;">Structure diagram unavailable</div>'

        sites_html = ""
        if reactive_centers:
            sites_pills = "".join([f'<span class="ap1-pill site">{site}</span>' for site in reactive_centers[:4]])
            sites_html = f"""
            <div style="margin-top: 0.75rem;">
                <div style="font-size: 0.7rem; font-weight: 700; color: #2563EB; text-transform: uppercase; letter-spacing: 0.05em; margin-bottom: 0.25rem;">
                    Reactive Interaction Centers:
                </div>
                <div class="ap1-tag-group">
                    {sites_pills}
                </div>
            </div>
            """

        safe_primary = html.escape(cur_primary)
        render_html(f"""
        <div class="ap1-comp-card">
            <div class="ap1-comp-mol">
                <span class="ap1-comp-badge">C1</span>
                {primary_mol_html}
            </div>
            <div class="ap1-comp-info">
                <div style="display: flex; align-items: center; margin-bottom: 0.25rem;">
                    <span class="ap1-comp-name">Primary Compound</span>
                    <span class="ap1-comp-role role-primary">Primary Active</span>
                </div>
                <div class="ap1-smiles-box" title="{safe_primary}">
                    {safe_primary}
                </div>
                <div class="ap1-tag-group">
                    {primary_desc_pills}
                </div>
                {sites_html}
            </div>
        </div>
        """)

        # Secondary Compound Cards (If provided)
        for s_idx, sec_sm in enumerate(cur_secondary):
            safe_sec = html.escape(sec_sm)
            sec_img = get_chemical_structure_img(sec_sm, width=220, height=200)
            sec_desc = get_molecular_descriptors(sec_sm)
            sec_desc_pills = format_descriptor_pills(sec_desc)
            sec_mol_html = f'<img src="{sec_img}" alt="Secondary Compound {s_idx+1}" style="max-width: 100%; max-height: 180px; object-fit: contain;"/>' if sec_img else '<div style="color: #94A3B8; font-size: 0.75rem;">Structure diagram unavailable</div>'

            render_html(f"""
            <div class="ap1-comp-card">
                <div class="ap1-comp-mol">
                    <span class="ap1-comp-badge">C{s_idx+2}</span>
                    {sec_mol_html}
                </div>
                <div class="ap1-comp-info">
                    <div style="display: flex; align-items: center; margin-bottom: 0.25rem;">
                        <span class="ap1-comp-name">Secondary Co-reactant {s_idx+1}</span>
                        <span class="ap1-comp-role role-secondary">Co-reactant / Excipient</span>
                    </div>
                    <div class="ap1-smiles-box" title="{safe_sec}">
                        {safe_sec}
                    </div>
                    <div class="ap1-tag-group">
                        {sec_desc_pills}
                    </div>
                </div>
            </div>
            """)

        st.markdown("<hr style='border: none; border-top: 1px solid #E2E8F0; margin: 2rem 0;'/>", unsafe_allow_html=True)

        # ----------------------------------------------------------------------
        # Stress Degradation & Incompatibility Heatmap
        # ----------------------------------------------------------------------
        st.markdown('<div class="section-title">Stress Degradation & Incompatibility Heatmap</div>', unsafe_allow_html=True)
        st.markdown('<div class="section-desc">Quantitative stress matrix modeling reactive center vulnerability across Acidic, Basic, Hydrolysis, Photolysis, Thermal, and Oxidative conditions using the WarmCool spectrum.</div>', unsafe_allow_html=True)

        fig = None
        try:
            fig = plot_heatmap(
                res["heatmap_matrix"],
                res["row_labels"],
                res["col_labels"],
                title=f"Stress Incompatibility Profile: {cur_primary}"
            )
        except Exception:
            fig = None

        if fig is not None:
            try:
                st.pyplot(fig)
            except Exception:
                render_html_heatmap(
                    res["heatmap_matrix"],
                    res["row_labels"],
                    res["col_labels"],
                    title=f"Stress Incompatibility Profile: {cur_primary}"
                )
        else:
            render_html_heatmap(
                res["heatmap_matrix"],
                res["row_labels"],
                res["col_labels"],
                title=f"Stress Incompatibility Profile: {cur_primary}"
            )


        st.markdown("<hr style='border: none; border-top: 1px solid #E2E8F0; margin: 2rem 0;'/>", unsafe_allow_html=True)

        # ----------------------------------------------------------------------
        # Mechanistic Framework Evaluation
        # ----------------------------------------------------------------------
        st.markdown('<div class="section-title">Mechanistic Framework Evaluation</div>', unsafe_allow_html=True)
        st.markdown('<div class="section-desc">Comprehensive kinetic pathways, microenvironmental influences, and thermodynamic justification.</div>', unsafe_allow_html=True)

        cot_text = res.get("chain_of_thought", "")
        if cot_text:
            cot_html = cot_text.replace("\n", "<br>")
            render_html(f"""
            <div class="ap1-cot-box">{cot_html}</div>
            """)

        st.markdown("<hr style='border: none; border-top: 1px solid #E2E8F0; margin: 2rem 0;'/>", unsafe_allow_html=True)

        # ----------------------------------------------------------------------
        # Degradation Products and Details (Top 5 Ranked)
        # ----------------------------------------------------------------------
        st.markdown('<div class="section-title">Degradation Products and Details</div>', unsafe_allow_html=True)
        st.markdown('<div class="section-desc">Ranked strictly by formation probability and thermodynamic stability (Top 5 maximum).</div>', unsafe_allow_html=True)

        for idx, imp in enumerate(res["impurities"]):
            prob_pct = imp["probability"] * 100
            cond = imp.get("condition", "Direct Degradation")
            cond_lower = cond.lower()

            if "oxid" in cond_lower:
                cond_class = "cond-oxid"
            elif "therm" in cond_lower:
                cond_class = "cond-therm"
            elif "photo" in cond_lower:
                cond_class = "cond-photo"
            elif "react" in cond_lower or "incomp" in cond_lower:
                cond_class = "cond-react"
            else:
                cond_class = "cond-hydro"

            imp_smiles = imp.get("smiles", "")
            imp_desc = get_molecular_descriptors(imp_smiles) if imp_smiles else {}
            imp_desc_pills = format_descriptor_pills(imp_desc) if imp_desc else ""
            enc_smiles = urllib.parse.quote(sanitize_smiles_py(imp_smiles)) if imp_smiles else ""

            if not imp_smiles:
                imp_mol_html = '<div style="color: #64748B; font-size: 0.8rem; text-align: center; padding: 2.5rem 0.5rem; font-weight: 500; font-family: \'Inter\', sans-serif;">Reactive functional group is absent</div>'
            else:
                imp_img = get_chemical_structure_img(imp_smiles, width=240, height=200)
                imp_mol_html = f'<img src="{imp_img}" alt="Structure of {imp.get("iupacName", "Impurity")}" onerror="if(!this.dataset.fallback){{this.dataset.fallback=\'1\';this.src=\'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/{enc_smiles}/PNG?record_type=2d&image_size=300x300\';}}else{{this.style.display=\'none\';}}" style="max-width: 100%; max-height: 180px; object-fit: contain;"/>' if imp_img else '<div style="color: #94A3B8; font-size: 0.75rem; text-align: center;">Structure diagram unavailable</div>'

            render_html(f"""
            <div class="ap1-imp-card">
                <div class="ap1-imp-svg">
                    <span class="ap1-imp-idx-badge">#{idx + 1}</span>
                    {imp_mol_html}
                </div>
                <div class="ap1-imp-body">
                    <div class="ap1-imp-header">
                        <div>
                            <div class="ap1-imp-title">{html.escape(str(imp.get('smiles') or imp.get('iupacName', 'Impurity')))}</div>
                            <div class="ap1-tag-group" style="margin-top: 0.4rem;">
                                {imp_desc_pills}
                            </div>
                        </div>
                        <div style="text-align: right;">
                            <div class="ap1-imp-prob-val">{prob_pct:.1f}%</div>
                            <div class="ap1-imp-prob-sub">
                                Heuristic: {(imp['probabilityHeuristic']*100):.1f}% | Boltzmann: {(imp['probabilityBoltzmann']*100):.1f}%
                            </div>
                            <div style="font-size: 0.75rem; font-family: 'JetBrains Mono', monospace; color: #64748B; margin-top: 0.25rem;">
                                Delta G: {imp['deltaG']:.2f} kcal/mol
                            </div>
                        </div>
                    </div>

                    <div class="ap1-mech-box">
                        <div class="ap1-mech-title">Chemical Mechanism:</div>
                        <div>{imp['mechanismExplanation']}</div>
                    </div>

                    <div style="display: flex; flex-wrap: wrap; gap: 0.5rem; align-items: center;">
                        <span class="ap1-badge-cond {cond_class}">{cond}</span>
                        <span class="ap1-pill" style="font-weight: 600; color: #1D4ED8; background: #EFF6FF; border-color: #DBEAFE;">
                            Origin: {imp.get('source', 'Stress degradation')}
                        </span>
                    </div>
                </div>
            </div>
            """)

        # ----------------------------------------------------------------------
        # CSV Report Export (Integrated directly in Output Page)
        # ----------------------------------------------------------------------
        report_rows = []
        for idx, imp in enumerate(res["impurities"]):
            prob_pct = round(imp.get("probability", 0.0) * 100, 2)
            i_smiles = imp.get("smiles", "")
            i_d = get_molecular_descriptors(i_smiles)
            report_rows.append({
                "Rank": idx + 1,
                "Byproduct IUPAC Name": imp.get("iupacName", ""),
                "SMILES": i_smiles,
                "MW (g/mol)": i_d.get("mw", ""),
                "LogP": i_d.get("logp", ""),
                "TPSA (sq A)": i_d.get("tpsa", ""),
                "HBD": i_d.get("hbd", ""),
                "HBA": i_d.get("hba", ""),
                "NumRotatableBonds": i_d.get("rotatable_bonds", ""),
                "HeavyAtomCount": i_d.get("heavy_atom_count", ""),
                "NumAromaticRings": i_d.get("aromatic_rings", ""),
                "NumHeteroatoms": i_d.get("heteroatoms", ""),
                "FractionCSP3": i_d.get("fraction_csp3", ""),
                "Degradation Condition": imp.get("condition", ""),
                "Origin": imp.get("source", "Stress degradation"),
                "Combined Formation Probability (%)": prob_pct,
                "Heuristic Probability (%)": round(imp.get("probabilityHeuristic", 0.0) * 100, 2),
                "Boltzmann Probability (%)": round(imp.get("probabilityBoltzmann", 0.0) * 100, 2),
                "Free Energy Delta G (kcal/mol)": imp.get("deltaG", 0.0),
                "Chemical Mechanism": imp.get("mechanismExplanation", ""),
                "Primary Compound": cur_primary,
                "Timestamp": datetime.date.today().strftime("%Y-%m-%d")
            })

        if pd is not None:
            df_report = pd.DataFrame(report_rows)
            report_csv_data = df_report.to_csv(index=False).encode('utf-8')
        else:
            out = io.StringIO()
            if report_rows:
                writer = csv.DictWriter(out, fieldnames=list(report_rows[0].keys()))
                writer.writeheader()
                writer.writerows(report_rows)
            report_csv_data = out.getvalue().encode('utf-8')
            df_report = report_rows

        st.markdown("<hr style='border: none; border-top: 1px solid #E2E8F0; margin: 2rem 0;'/>", unsafe_allow_html=True)
        st.markdown('<div class="section-title">Computational CSV Report</div>', unsafe_allow_html=True)
        st.markdown('<div class="section-desc">Download complete structured analysis data including reaction pathways, thermodynamic free energy values, and kinetic formation probabilities.</div>', unsafe_allow_html=True)

        st.download_button(
            label="Download Degradation Analysis CSV Report",
            data=report_csv_data,
            file_name=f"interaction_prediction_report_{datetime.date.today()}.csv",
            mime="text/csv",
            type="primary",
            use_container_width=True
        )

        with st.expander("View Tabular CSV Report Data", expanded=False):
            st.dataframe(df_report, use_container_width=True)

        # ----------------------------------------------------------------------
        # Disclaimer Card
        # ----------------------------------------------------------------------
        render_html("""
        <div class="ap1-disclaimer-box" style="margin-top: 1.5rem;">
            Disclaimer: INTERACTION is an AI-assisted computational chemistry modeling tool designed for reaction pathway exploration and byproduct screening. Predictions should be verified by experimental analytical assays.
        </div>
        """)

# ==============================================================================
# TAB 2: Query Logbook (100 Queries)
# ==============================================================================
with tab_logbook:
    st.markdown('<div class="section-title">Persistent Reaction Query Logbook</div>', unsafe_allow_html=True)
    st.markdown('<div class="section-desc">Maintaining the 100 most recent calculation queries with SMILES inputs, predicted impurities, and dates.</div>', unsafe_allow_html=True)

    df_log = load_logbook()

    is_empty = False
    if pd is not None and isinstance(df_log, pd.DataFrame):
        is_empty = df_log.empty
    elif isinstance(df_log, list):
        is_empty = len(df_log) == 0
    else:
        is_empty = True

    if is_empty:
        st.info("No queries recorded yet. Run a prediction on the Analysis & Predictions tab to record data.")
    else:
        if pd is not None and isinstance(df_log, pd.DataFrame):
            csv_data = df_log.to_csv(index=False).encode('utf-8')
        else:
            out = io.StringIO()
            if df_log:
                writer = csv.DictWriter(out, fieldnames=LOGBOOK_COLUMNS)
                writer.writeheader()
                writer.writerows(df_log)
            csv_data = out.getvalue().encode('utf-8')

        st.download_button(
            label="Download Logbook CSV",
            data=csv_data,
            file_name=f"query_logbook_{datetime.date.today()}.csv",
            mime="text/csv",
            type="primary"
        )
        st.dataframe(df_log, use_container_width=True, height=500)
