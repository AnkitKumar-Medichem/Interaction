"""
A-Pi1: AI-Powered Drug-Excipient Interaction & Degradation Predictor
Streamlit Application for Pharmaceutical Formulation & Stability Evaluation
"""

import io
import json
import os
import time
from typing import List, Optional
import pandas as pd
from pydantic import BaseModel, Field
import streamlit as st

# Attempt RDKit import for chemical structure rendering
try:
    from rdkit import Chem
    from rdkit.Chem import Descriptors, Draw
    from rdkit.Chem.Draw import rdMolDraw2D
    RDKIT_AVAILABLE = True
except ImportError:
    RDKIT_AVAILABLE = False

# Google GenAI SDK
try:
    from google import genai
    from google.genai import types
    GENAI_AVAILABLE = True
except ImportError:
    GENAI_AVAILABLE = False


# ==============================================================================
# Page Configuration & Styling
# ==============================================================================
st.set_page_config(
    page_title="A-Pi1 | Drug-Excipient Degradation Predictor",
    page_icon="💊",
    layout="wide",
    initial_sidebar_state="expanded",
)

# Modern, clean pharmaceutical laboratory styling
st.markdown("""
<style>
    /* Metric Card Styling */
    .metric-card {
        background: #f8fafc;
        border: 1px solid #e2e8f0;
        border-radius: 10px;
        padding: 16px;
        text-align: center;
    }
    .metric-title {
        font-size: 0.85rem;
        font-weight: 600;
        color: #64748b;
        text-transform: uppercase;
        letter-spacing: 0.05em;
    }
    .metric-value {
        font-size: 1.8rem;
        font-weight: 700;
        color: #0f172a;
        margin-top: 4px;
    }
    
    /* Condition Badges */
    .badge-oxidation { background-color: #fef3c7; color: #92400e; border: 1px solid #fde68a; padding: 2px 8px; border-radius: 9999px; font-size: 0.75rem; font-weight: 600; }
    .badge-acidic { background-color: #e0f2fe; color: #0369a1; border: 1px solid #bae6fd; padding: 2px 8px; border-radius: 9999px; font-size: 0.75rem; font-weight: 600; }
    .badge-basic { background-color: #ecfeff; color: #0e7490; border: 1px solid #a5f3fc; padding: 2px 8px; border-radius: 9999px; font-size: 0.75rem; font-weight: 600; }
    .badge-photodegradation { background-color: #f3e8ff; color: #6b21a8; border: 1px solid #e9d5ff; padding: 2px 8px; border-radius: 9999px; font-size: 0.75rem; font-weight: 600; }
    .badge-thermal { background-color: #ffe4e6; color: #9f1239; border: 1px solid #fecdd3; padding: 2px 8px; border-radius: 9999px; font-size: 0.75rem; font-weight: 600; }
    
    /* Impurity Card */
    .impurity-card {
        background: #ffffff;
        border: 1px solid #e2e8f0;
        border-radius: 12px;
        padding: 18px;
        box-shadow: 0 1px 3px rgba(0,0,0,0.05);
        margin-bottom: 16px;
    }
</style>
""", unsafe_allow_html=True)


# ==============================================================================
# Pydantic Output Schemas
# ==============================================================================
class CompoundItem(BaseModel):
    name: str = Field(description="Strict original or identified compound name")
    smiles: str = Field(description="Canonical valid SMILES string")
    features: List[str] = Field(description="Key structural moieties and functional groups")
    interactionSites: List[str] = Field(description="Specific sites likely to interact or undergo degradation")

class DegradantItem(BaseModel):
    iupacName: str = Field(description="IUPAC or common name of the newly generated degradant/impurity")
    smiles: str = Field(description="Valid, canonical SMILES string of the impurity with strictly correct valences")
    structureDescription: str = Field(description="Key structural difference or transformation from parent")
    origin: str = Field(description="Originating compound(s), e.g. Compound 1 or Compound 1 + Compound 2")
    probability: float = Field(description="Primary ranking probability between 0.01 and 0.99")
    probabilityHeuristic: Optional[float] = Field(default=None, description="Heuristic reasoning probability")
    probabilityBoltzmann: Optional[float] = Field(default=None, description="Boltzmann thermodynamic probability at 298.15K")
    relativeEnergy: Optional[float] = Field(default=None, description="Estimated relative formation energy ΔG in kcal/mol")
    condition: str = Field(description="One of: Oxidation, Acidic Hydrolysis, Basic Hydrolysis, Photodegradation, Thermal Degradation")
    source: str = Field(description="One of: Stress degradation, Interaction with other compound")
    mechanismExplanation: str = Field(description="Chemical mechanism including pH shift, transesterification, hydrolysis, or complexation")

class AnalysisReport(BaseModel):
    chainOfThought: str = Field(description="Detailed step-by-step chemical reasoning before prediction")
    compounds: List[CompoundItem] = Field(description="Profile of starting compounds")
    interactionType: str = Field(description="Physical, Chemical, or None")
    mechanism: str = Field(description="Overall interaction narrative and formulation implications")
    degradationImpurities: List[DegradantItem] = Field(description="Exactly 10 predicted degradation impurities")


# ==============================================================================
# Preset Chemical Formulations
# ==============================================================================
PRESETS = {
    "Aspirin + Lactose Monohydrate (Maillard/Transesterification)": [
        {"name": "Aspirin", "smiles": "CC(=O)Oc1ccccc1C(=O)O", "type": "SMILES"},
        {"name": "Lactose Monohydrate", "smiles": "C([C@@H]1[C@@H]([C@@H]([C@H]([C@H](O1)O[C@@H]2[C@H](O[C@H]([C@@H]([C@H]2O)O)O)CO)O)O)O", "type": "SMILES"}
    ],
    "Paracetamol + Povidone (K-30) (Oxidation & Complexation)": [
        {"name": "Paracetamol", "smiles": "CC(=O)Nc1ccc(cc1)O", "type": "SMILES"},
        {"name": "Povidone (Monomer Unit)", "smiles": "C=CN1CCCC1=O", "type": "SMILES"}
    ],
    "Metformin HCl + Magnesium Stearate (Physical/Chemical Shift)": [
        {"name": "Metformin", "smiles": "CN(C)C(=N)NC(=N)N", "type": "SMILES"},
        {"name": "Stearic Acid", "smiles": "CCCCCCCCCCCCCCCCCC(=O)O", "type": "SMILES"}
    ],
    "Ciprofloxacin + Calcium Carbonate (Chelation Complex)": [
        {"name": "Ciprofloxacin", "smiles": "C1CC1N2C=C(C(=O)C3=CC(=C(C=C32)N4CCNCC4)F)C(=O)O", "type": "SMILES"},
        {"name": "Calcium Carbonate", "smiles": "[Ca+2].[O-]C(=O)[O-]", "type": "SMILES"}
    ]
}


# ==============================================================================
# Chemical Utilities (RDKit)
# ==============================================================================
def render_mol_svg(smiles: str, width: int = 300, height: int = 200) -> Optional[str]:
    """Render a 2D molecule as an SVG string using RDKit."""
    if not RDKIT_AVAILABLE or not smiles:
        return None
    try:
        mol = Chem.MolFromSmiles(smiles.strip())
        if not mol:
            return None
        drawer = rdMolDraw2D.MolDraw2DSVG(width, height)
        opts = drawer.drawOptions()
        opts.clearBackground = True
        drawer.DrawMolecule(mol)
        drawer.FinishDrawing()
        return drawer.GetDrawingText()
    except Exception:
        return None

def get_molecular_descriptors(smiles: str) -> dict:
    """Calculate basic physicochemical properties using RDKit."""
    if not RDKIT_AVAILABLE or not smiles:
        return {}
    try:
        mol = Chem.MolFromSmiles(smiles.strip())
        if not mol:
            return {}
        return {
            "MolWt": round(Descriptors.MolWt(mol), 2),
            "LogP": round(Descriptors.MolLogP(mol), 2),
            "TPSA": round(Descriptors.TPSA(mol), 2),
            "RotatableBonds": Descriptors.NumRotatableBonds(mol),
            "HBD": Descriptors.NumHDonors(mol),
            "HBA": Descriptors.NumHAcceptors(mol),
        }
    except Exception:
        return {}


# ==============================================================================
# Gemini Prediction Engine
# ==============================================================================
def call_gemini_prediction(
    api_key: str,
    compounds: List[dict],
    method: str,
    model_name: str = "gemini-3.1-flash-lite"
) -> AnalysisReport:
    """Invokes the Google Gemini Generative AI API with structured Pydantic output."""
    if not GENAI_AVAILABLE:
        raise RuntimeError("google-genai Python package is not installed. Please run: pip install google-genai")

    client = genai.Client(api_key=api_key)

    # Format compound details
    compounds_info = []
    for i, comp in enumerate(compounds):
        name = comp.get("name", f"Compound {i+1}").strip() or f"Compound {i+1}"
        val = comp.get("value", "").strip()
        comp_type = comp.get("type", "Name")
        
        desc_str = ""
        if comp_type == "SMILES" and val:
            desc = get_molecular_descriptors(val)
            if desc:
                desc_str = f" [Calculated: MW={desc['MolWt']}, LogP={desc['LogP']}, TPSA={desc['TPSA']}]"

        spec = f"Compound {i+1}: Name='{name}', Format={comp_type}, Data='{val}'.{desc_str}"
        compounds_info.append(spec)

    compounds_text = "\n".join(compounds_info)

    probability_instruction = ""
    if method == "Boltzmann":
        probability_instruction = "Calculate probability of formation based on Boltzmann distribution at 298.15K. You MUST provide relative formation energy (relativeEnergy) in kcal/mol."
    elif method == "Heuristic":
        probability_instruction = "Calculate probability of formation based on expert chemical reasoning, reactive site vulnerability, and reaction kinetics."
    else:  # Both
        probability_instruction = (
            "Provide BOTH 'probabilityHeuristic' (expert chemical reasoning) and 'probabilityBoltzmann' "
            "(thermodynamic stability from relative ΔG at 298.15K). Also provide 'relativeEnergy' in kcal/mol. "
            "The main 'probability' field must reflect the Boltzmann probability for ranking."
        )

    system_instruction = f"""You are a senior pharmaceutical degradation evaluator and preformulation stability scientist.
Your task is to predict the degradation profile of Compound 1 (Active Pharmaceutical Ingredient) in the provided mixture using the {method} framework.

Analytical Mandates:
1. Identify exact chemical structures, functional groups, and reactive interaction sites for all input compounds.
2. Evaluate degradation arising from:
   - Stress degradation of Compound 1 under accelerated stability testing conditions.
   - Chemical and physical interactions between Compound 1 and excipients/co-compounds.
3. Test under all 5 mandatory stress conditions:
   - Oxidation
   - Acidic Hydrolysis
   - Basic Hydrolysis
   - Photodegradation
   - Thermal Degradation
4. Predict exactly 10 distinct degradation impurities / interaction products.
5. For each impurity:
   - Provide a unique, descriptive IUPAC or common name for the NEW degradant. Do NOT repeat the starting material name.
   - Output a strictly valid, canonical SMILES string. Obey standard valences (Carbon max 4, Nitrogen 3/4, Oxygen 2).
   - Detail the formation mechanism (e.g. transesterification, Maillard reaction, oxidation of amine/thioether, ester hydrolysis).
   - {probability_instruction}
6. Classify the overall interaction type as Physical, Chemical, or None.
"""

    prompt = f"""Predict and evaluate the pharmaceutical degradation of Compound 1 in the following mixture using the {method} analytical method:

INPUT COMPOUNDS:
{compounds_text}

Provide exactly 10 degradation impurities ranked by probability."""

    # Models to attempt in order of reliability/quota
    models_to_try = [model_name, "gemini-3.1-flash-lite", "gemini-flash-latest", "gemini-3.8-flash"]
    # Deduplicate while preserving order
    seen = set()
    candidate_models = [m for m in models_to_try if not (m in seen or seen.add(m))]

    last_error = None
    for candidate in candidate_models:
        for attempt in range(2):
            try:
                response = client.models.generate_content(
                    model=candidate,
                    contents=prompt,
                    config=types.GenerateContentConfig(
                        system_instruction=system_instruction,
                        temperature=0.1,
                        response_mime_type="application/json",
                        response_schema=AnalysisReport,
                    ),
                )
                
                if not response.text:
                    raise ValueError("Model returned an empty response.")
                
                # Parse structured output
                data = json.loads(response.text)
                return AnalysisReport.model_validate(data)
                
            except Exception as e:
                last_error = e
                err_str = str(e)
                # If rate limit or temporary server spike, retry with small backoff
                if "429" in err_str or "quota" in err_str.lower() or "503" in err_str or "overloaded" in err_str.lower():
                    time.sleep(2 + attempt * 2)
                    continue
                # If model not found or fatal, jump to next candidate model
                break

    raise RuntimeError(f"Failed to generate prediction with Gemini API: {last_error}")


# ==============================================================================
# Excel Export Generator
# ==============================================================================
def generate_excel_report(report: AnalysisReport) -> bytes:
    """Generate a multi-tab Excel (.xlsx) workbook for pharmaceutical reporting."""
    output = io.BytesIO()
    with pd.ExcelWriter(output, engine="openpyxl") as writer:
        # Sheet 1: Executive Summary
        summary_data = {
            "Parameter": [
                "Interaction Type",
                "Primary Degradation Mechanism",
                "Total Predicted Impurities",
                "Evaluation Timestamp",
            ],
            "Value": [
                report.interactionType,
                report.mechanism,
                len(report.degradationImpurities),
                time.strftime("%Y-%m-%d %H:%M:%S UTC", time.gmtime()),
            ]
        }
        pd.DataFrame(summary_data).to_excel(writer, sheet_name="Executive Summary", index=False)

        # Sheet 2: Starting Materials
        compounds_data = []
        for c in report.compounds:
            desc = get_molecular_descriptors(c.smiles)
            compounds_data.append({
                "Compound Name": c.name,
                "Canonical SMILES": c.smiles,
                "Key Features": ", ".join(c.features),
                "Interaction Sites": ", ".join(c.interactionSites),
                "Mol Wt (g/mol)": desc.get("MolWt", "N/A"),
                "LogP": desc.get("LogP", "N/A"),
                "TPSA (Å²)": desc.get("TPSA", "N/A"),
            })
        pd.DataFrame(compounds_data).to_excel(writer, sheet_name="Starting Materials", index=False)

        # Sheet 3: Degradation Impurities
        impurities_data = []
        for imp in report.degradationImpurities:
            desc = get_molecular_descriptors(imp.smiles)
            impurities_data.append({
                "IUPAC / Impurity Name": imp.iupacName,
                "SMILES String": imp.smiles,
                "Condition": imp.condition,
                "Source Type": imp.source,
                "Originating Species": imp.origin,
                "Primary Probability": imp.probability,
                "Heuristic Probability": imp.probabilityHeuristic if imp.probabilityHeuristic is not None else "N/A",
                "Boltzmann Probability": imp.probabilityBoltzmann if imp.probabilityBoltzmann is not None else "N/A",
                "Relative ΔG (kcal/mol)": imp.relativeEnergy if imp.relativeEnergy is not None else "N/A",
                "Mol Wt (g/mol)": desc.get("MolWt", "N/A"),
                "LogP": desc.get("LogP", "N/A"),
                "TPSA (Å²)": desc.get("TPSA", "N/A"),
                "Mechanism Explanation": imp.mechanismExplanation,
            })
        pd.DataFrame(impurities_data).to_excel(writer, sheet_name="Impurities Profile", index=False)

    return output.getvalue()


# ==============================================================================
# UI Application Layout
# ==============================================================================
def main():
    # Header Banner
    st.title("💊 A-Pi1 Degradation Predictor")
    st.caption("AI-Powered Pharmaceutical Drug-Excipient Compatibility & Forced Stress Degradation Evaluator")

    # Sidebar: Configuration
    with st.sidebar:
        st.header("⚙️ Configuration")
        
        # API Key Resolution
        default_api_key = (
            st.secrets.get("GEMINI_API_KEY", "") 
            if hasattr(st, "secrets") else ""
        ) or os.environ.get("GEMINI_API_KEY", "")
        
        api_key = st.text_input(
            "Gemini API Key",
            value=default_api_key,
            type="password",
            help="Provide your Google Gemini API key. On Streamlit Cloud, configure it in App Secrets."
        )

        st.divider()

        # Method Selection
        st.subheader("🔬 Analytical Framework")
        method = st.radio(
            "Evaluation Method",
            ["Both", "Heuristic", "Boltzmann"],
            format_func=lambda x: {
                "Both": "Both (Comparative: Heuristic + Boltzmann)",
                "Heuristic": "Heuristic (Reaction Kinetics & Expert Logic)",
                "Boltzmann": "Boltzmann (ΔG Thermodynamics @ 298.15K)"
            }[x],
            help="Select whether to evaluate degradation via chemical reaction kinetics, Boltzmann thermodynamics, or a dual comparative view."
        )

        st.divider()

        # Model Selector
        model_selection = st.selectbox(
            "AI Core Engine",
            ["gemini-3.1-flash-lite", "gemini-flash-latest", "gemini-3.8-flash"],
            index=0,
            help="Default is gemini-3.1-flash-lite for maximum speed, accuracy, and high quota resilience."
        )

        st.divider()
        st.markdown("""
        **About A-Pi1:**
        - Predicts drug-excipient interactions and stress degradants.
        - Evaluates Oxidation, Hydrolysis (Acid/Base), Photolysis, & Thermolysis.
        - Computes 2D chemical structures via RDKit.
        - Exports full stability profiles to Excel.
        """)

    # Main Interface Tabs
    tab_input, tab_results = st.tabs(["🧪 Formulation Input", "📊 Stability Analysis"])

    # Handle Preset Selection
    with tab_input:
        col_top1, col_top2 = st.columns([2, 1])
        with col_top1:
            st.subheader("Starting Formulation Components")
            st.write("Specify Compound 1 (Active Pharmaceutical Ingredient) and up to 4 excipients or co-compounds.")
        with col_top2:
            preset_choice = st.selectbox(
                "Quick Load Formulation Preset:",
                ["-- Select a Preset --"] + list(PRESETS.keys())
            )
            if preset_choice != "-- Select a Preset --":
                st.session_state["compounds"] = [
                    {"name": c["name"], "value": c["smiles"], "type": c["type"]}
                    for c in PRESETS[preset_choice]
                ]

        # Initialize compounds in session state
        if "compounds" not in st.session_state:
            st.session_state["compounds"] = [
                {"name": "Aspirin", "value": "CC(=O)Oc1ccccc1C(=O)O", "type": "SMILES"},
                {"name": "Lactose Monohydrate", "value": "C([C@@H]1[C@@H]([C@@H]([C@H]([C@H](O1)O[C@@H]2[C@H](O[C@H]([C@@H]([C@H]2O)O)O)CO)O)O)O", "type": "SMILES"},
            ]

        # Render compound input forms
        compounds_list = st.session_state["compounds"]
        
        for idx in range(len(compounds_list)):
            comp = compounds_list[idx]
            is_api = (idx == 0)
            
            with st.container():
                st.markdown(f"**{'💊 Active Drug Substance (Compound 1)' if is_api else f'📦 Excipient / Co-Compound {idx+1}'}**")
                c1, c2, c3, c4 = st.columns([2, 1, 3, 0.5])
                
                with c1:
                    comp["name"] = st.text_input(
                        f"Compound {idx+1} Name",
                        value=comp.get("name", ""),
                        key=f"name_{idx}",
                        placeholder="e.g. Aspirin, Lactose, Povidone"
                    )
                with c2:
                    comp["type"] = st.selectbox(
                        "Input Type",
                        ["SMILES", "Name"],
                        index=0 if comp.get("type", "SMILES") == "SMILES" else 1,
                        key=f"type_{idx}"
                    )
                with c3:
                    comp["value"] = st.text_input(
                        "SMILES String / Identifier",
                        value=comp.get("value", ""),
                        key=f"val_{idx}",
                        placeholder="e.g. CC(=O)Oc1ccccc1C(=O)O" if comp["type"] == "SMILES" else "Chemical name"
                    )
                with c4:
                    if not is_api:
                        if st.button("🗑️", key=f"del_{idx}", help="Remove this component"):
                            compounds_list.pop(idx)
                            st.rerun()

                # Live 2D Structure Preview if SMILES provided
                if comp["type"] == "SMILES" and comp["value"].strip() and RDKIT_AVAILABLE:
                    svg_data = render_mol_svg(comp["value"], width=220, height=130)
                    if svg_data:
                        desc = get_molecular_descriptors(comp["value"])
                        p_col1, p_col2 = st.columns([1, 4])
                        with p_col1:
                            st.image(f"data:image/svg+xml;utf8,{svg_data}", width=180)
                        with p_col2:
                            st.caption(f"**Molecular Weight:** {desc.get('MolWt')} g/mol | **LogP:** {desc.get('LogP')} | **TPSA:** {desc.get('TPSA')} Å² | **Rotatable Bonds:** {desc.get('RotatableBonds')}")

                st.divider()

        # Add Excipient Button
        if len(compounds_list) < 5:
            if st.button("➕ Add Excipient / Additional Component"):
                compounds_list.append({"name": "", "value": "", "type": "Name"})
                st.rerun()

        st.write("")
        run_col1, run_col2 = st.columns([1, 3])
        with run_col1:
            run_btn = st.button("🚀 Predict Degradation & Compatibility", type="primary", use_container_width=True)
        with run_col2:
            st.caption("Predicts stress degradation pathways (Oxidation, Hydrolysis, Photodegradation, Thermal) and chemical incompatibilities.")

        if run_btn:
            if not api_key:
                st.error("⚠️ Gemini API Key is required. Enter it in the sidebar or set GEMINI_API_KEY.")
                return

            if not compounds_list[0]["value"].strip() and not compounds_list[0]["name"].strip():
                st.error("⚠️ Please provide at least Compound 1 (Active Drug Substance).")
                return

            with st.spinner("Analyzing molecular reactivity, stress degradation pathways, and thermodynamics..."):
                try:
                    result = call_gemini_prediction(
                        api_key=api_key,
                        compounds=compounds_list,
                        method=method,
                        model_name=model_selection
                    )
                    st.session_state["prediction_report"] = result
                    st.success("✅ Degradation prediction and compatibility evaluation complete!")
                except Exception as exc:
                    st.error(f"❌ Analysis failed: {exc}")

    # Results Display Tab
    with tab_results:
        report: Optional[AnalysisReport] = st.session_state.get("prediction_report")
        
        if not report:
            st.info("👈 No prediction report available yet. Enter compounds in the 'Formulation Input' tab and click 'Predict Degradation'.")
            return

        # Top Executive Metrics
        st.subheader("📋 Executive Stability Summary")
        m_col1, m_col2, m_col3, m_col4 = st.columns(4)
        
        with m_col1:
            int_color = "#ef4444" if report.interactionType == "Chemical" else ("#f59e0b" if report.interactionType == "Physical" else "#10b981")
            st.markdown(f"""
            <div class="metric-card">
                <div class="metric-title">Interaction Type</div>
                <div class="metric-value" style="color: {int_color};">{report.interactionType}</div>
            </div>
            """, unsafe_allow_html=True)
        with m_col2:
            st.markdown(f"""
            <div class="metric-card">
                <div class="metric-title">Predicted Impurities</div>
                <div class="metric-value">{len(report.degradationImpurities)}</div>
            </div>
            """, unsafe_allow_html=True)
        with m_col3:
            high_prob = sum(1 for imp in report.degradationImpurities if imp.probability >= 0.5)
            st.markdown(f"""
            <div class="metric-card">
                <div class="metric-title">High Probability (>50%)</div>
                <div class="metric-value" style="color: #b91c1c;">{high_prob}</div>
            </div>
            """, unsafe_allow_html=True)
        with m_col4:
            # Count most common condition
            cond_counts = pd.Series([imp.condition for imp in report.degradationImpurities]).value_counts()
            top_cond = cond_counts.index[0] if not cond_counts.empty else "N/A"
            st.markdown(f"""
            <div class="metric-card">
                <div class="metric-title">Most Vulnerable Pathway</div>
                <div class="metric-value" style="font-size: 1.25rem; margin-top: 10px;">{top_cond}</div>
            </div>
            """, unsafe_allow_html=True)

        st.write("")
        st.markdown(f"**Mechanism Narrative:** {report.mechanism}")

        with st.expander("🔍 Analytical Step-by-Step Chemical Reasoning (Chain-of-Thought)"):
            st.write(report.chainOfThought)

        st.divider()

        # Section: Starting Materials Breakdown
        st.subheader("🧬 Starting Components & Reactive Sites")
        c_cols = st.columns(min(len(report.compounds), 4))
        for i, comp in enumerate(report.compounds):
            with c_cols[i % len(c_cols)]:
                st.markdown(f"**{comp.name}**")
                st.caption(f"`{comp.smiles}`")
                
                if RDKIT_AVAILABLE:
                    svg_data = render_mol_svg(comp.smiles, width=240, height=150)
                    if svg_data:
                        st.image(f"data:image/svg+xml;utf8,{svg_data}", width=220)
                
                desc = get_molecular_descriptors(comp.smiles)
                if desc:
                    st.write(f"**MW:** {desc.get('MolWt')} | **LogP:** {desc.get('LogP')}")
                
                st.markdown("**Interaction Sites:**")
                for site in comp.interactionSites:
                    st.markdown(f"- `{site}`")

        st.divider()

        # Section: Degradation Impurities Dashboard
        st.subheader("⚠️ Predicted Degradation Impurities")
        
        # Interactive Filtering & Sorting Controls
        f_col1, f_col2, f_col3 = st.columns(3)
        with f_col1:
            condition_filter = st.selectbox(
                "Filter by Stress Condition:",
                ["All Conditions", "Oxidation", "Acidic Hydrolysis", "Basic Hydrolysis", "Photodegradation", "Thermal Degradation"]
            )
        with f_col2:
            source_filter = st.selectbox(
                "Filter by Degradation Source:",
                ["All Sources", "Stress degradation", "Interaction with other compound"]
            )
        with f_col3:
            sort_by = st.selectbox(
                "Sort Impurities By:",
                ["Primary Probability (High to Low)", "Relative Energy ΔG (Low to High)", "Alphabetical"]
            )

        # Filter impurities
        filtered_impurities = [
            imp for imp in report.degradationImpurities
            if (condition_filter == "All Conditions" or imp.condition == condition_filter)
            and (source_filter == "All Sources" or imp.source == source_filter)
        ]

        # Sort impurities
        if sort_by == "Primary Probability (High to Low)":
            filtered_impurities.sort(key=lambda x: x.probability, reverse=True)
        elif sort_by == "Relative Energy ΔG (Low to High)":
            filtered_impurities.sort(key=lambda x: (x.relativeEnergy if x.relativeEnergy is not None else 999.0))
        elif sort_by == "Alphabetical":
            filtered_impurities.sort(key=lambda x: x.iupacName)

        st.caption(f"Showing {len(filtered_impurities)} of {len(report.degradationImpurities)} predicted degradation products.")

        # Grid of Impurity Cards
        for imp in filtered_impurities:
            badge_class = {
                "Oxidation": "badge-oxidation",
                "Acidic Hydrolysis": "badge-acidic",
                "Basic Hydrolysis": "badge-basic",
                "Photodegradation": "badge-photodegradation",
                "Thermal Degradation": "badge-thermal",
            }.get(imp.condition, "badge-oxidation")

            with st.container():
                card_col1, card_col2 = st.columns([1, 2.5])
                
                with card_col1:
                    if RDKIT_AVAILABLE:
                        svg_data = render_mol_svg(imp.smiles, width=280, height=180)
                        if svg_data:
                            st.image(f"data:image/svg+xml;utf8,{svg_data}", width=250)
                        else:
                            st.caption("2D Structure Unavailable")
                    else:
                        st.caption("RDKit not installed")

                with card_col2:
                    st.markdown(f"""
                    <h4>{imp.iupacName}</h4>
                    <span class="{badge_class}">{imp.condition}</span>
                    <span style="background: #f1f5f9; color: #475569; padding: 2px 8px; border-radius: 9999px; font-size: 0.75rem; font-weight: 600; margin-left: 6px;">{imp.source}</span>
                    """, unsafe_allow_html=True)
                    
                    st.markdown(f"**SMILES:** `{imp.smiles}`")
                    st.markdown(f"**Origin:** {imp.origin} | **Structural Shift:** {imp.structureDescription}")

                    # Probability Bars & Energy Metrics
                    stat_col1, stat_col2, stat_col3 = st.columns(3)
                    with stat_col1:
                        st.metric("Primary Probability", f"{round(imp.probability * 100, 1)}%")
                    with stat_col2:
                        if imp.probabilityHeuristic is not None:
                            st.metric("Heuristic Probability", f"{round(imp.probabilityHeuristic * 100, 1)}%")
                        elif imp.probabilityBoltzmann is not None:
                            st.metric("Boltzmann Probability", f"{round(imp.probabilityBoltzmann * 100, 1)}%")
                        else:
                            st.metric("Probability Score", f"{round(imp.probability * 100, 1)}%")
                    with stat_col3:
                        energy_str = f"{imp.relativeEnergy:.2f} kcal/mol" if imp.relativeEnergy is not None else "N/A"
                        st.metric("Relative ΔG", energy_str)

                    st.markdown(f"**Mechanism:** {imp.mechanismExplanation}")
                
                st.divider()

        # Section: Comparative Chart
        st.subheader("📈 Energy & Probability Distribution")
        chart_data = []
        for imp in report.degradationImpurities:
            chart_data.append({
                "Impurity": imp.iupacName[:25] + "...",
                "Primary Probability": imp.probability,
                "Heuristic": imp.probabilityHeuristic if imp.probabilityHeuristic is not None else imp.probability,
                "Boltzmann": imp.probabilityBoltzmann if imp.probabilityBoltzmann is not None else imp.probability,
                "ΔG (kcal/mol)": imp.relativeEnergy if imp.relativeEnergy is not None else 0.0,
            })
        df_chart = pd.DataFrame(chart_data)
        st.bar_chart(df_chart.set_index("Impurity")[["Primary Probability", "Heuristic", "Boltzmann"]])

        # Section: Export Center
        st.subheader("💾 Export Reports")
        ex_col1, ex_col2 = st.columns(2)
        
        with ex_col1:
            excel_bytes = generate_excel_report(report)
            st.download_button(
                label="📥 Download Comprehensive Excel Report (.xlsx)",
                data=excel_bytes,
                file_name="A-Pi1_Degradation_Report.xlsx",
                mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
                use_container_width=True
            )
        
        with ex_col2:
            json_str = report.model_dump_json(indent=2)
            st.download_button(
                label="📥 Download Raw JSON Stability Report",
                data=json_str,
                file_name="A-Pi1_Degradation_Report.json",
                mime="application/json",
                use_container_width=True
            )


if __name__ == "__main__":
    main()
