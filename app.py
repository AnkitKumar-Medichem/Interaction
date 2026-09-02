"""
A-Pi1: AI-Powered Drug-Excipient Interaction & Degradation Predictor
Streamlit Application Entry Point

This file serves and hosts the exact A-Pi1 pharmaceutical formulation GUI
(React + Vite + Tailwind CSS + Lucide + RDKit engine) inside Streamlit with 
complete fidelity, zero visual compromises, and responsive layout.
"""

import os
import shutil
import subprocess
import streamlit as st
import streamlit.components.v1 as components

# ==============================================================================
# 1. Streamlit Page Configuration
# ==============================================================================
st.set_page_config(
    page_title="A-Pi1 | Drug-Excipient Degradation Predictor",
    page_icon="💊",
    layout="wide",
    initial_sidebar_state="collapsed",
)

# ==============================================================================
# 2. Exact GUI Viewport Styling
# Hides Streamlit's default chrome so the exact A-Pi1 interface fills the screen
# ==============================================================================
st.markdown("""
<style>
    /* Remove default Streamlit padding, header, toolbar, and footer */
    header[data-testid="stHeader"] {
        display: none !important;
        height: 0 !important;
    }
    footer {
        display: none !important;
    }
    #MainMenu {
        display: none !important;
    }
    div[data-testid="stToolbar"] {
        display: none !important;
    }
    div[data-testid="stDecoration"] {
        display: none !important;
    }
    div[data-testid="stStatusWidget"] {
        display: none !important;
    }
    .main .block-container {
        padding: 0 !important;
        margin: 0 !important;
        max-width: 100vw !important;
        width: 100vw !important;
        height: 100vh !important;
        overflow: hidden !important;
    }
    iframe {
        position: fixed !important;
        top: 0 !important;
        left: 0 !important;
        width: 100vw !important;
        height: 100vh !important;
        border: none !important;
        margin: 0 !important;
        padding: 0 !important;
        z-index: 999999 !important;
    }
</style>
""", unsafe_allow_html=True)

# ==============================================================================
# 3. Static Asset Preparation
# Ensures built production files are accessible in the static directory
# ==============================================================================
root_dir = os.path.dirname(os.path.abspath(__file__))
static_dir = os.path.join(root_dir, "static")
dist_dir = os.path.join(root_dir, "dist")

if not os.path.exists(static_dir) or not os.path.exists(os.path.join(static_dir, "index.html")):
    if os.path.exists(dist_dir) and os.path.exists(os.path.join(dist_dir, "index.html")):
        os.makedirs(static_dir, exist_ok=True)
        for item in os.listdir(dist_dir):
            s = os.path.join(dist_dir, item)
            d = os.path.join(static_dir, item)
            if os.path.isdir(s):
                shutil.copytree(s, d, dirs_exist_ok=True)
            else:
                shutil.copy2(s, d)
    else:
        try:
            subprocess.run(["npm", "run", "build"], cwd=root_dir, check=True)
            if os.path.exists(dist_dir):
                os.makedirs(static_dir, exist_ok=True)
                for item in os.listdir(dist_dir):
                    s = os.path.join(dist_dir, item)
                    d = os.path.join(static_dir, item)
                    if os.path.isdir(s):
                        shutil.copytree(s, d, dirs_exist_ok=True)
                    else:
                        shutil.copy2(s, d)
        except Exception:
            pass

# ==============================================================================
# 4. Resolve Gemini API Key from Streamlit Secrets or Environment
# ==============================================================================
api_key = os.environ.get("GEMINI_API_KEY", "")
if not api_key:
    try:
        api_key = st.secrets.get("GEMINI_API_KEY", "")
    except Exception:
        pass
if not api_key:
    try:
        api_key = st.secrets.get("VITE_GEMINI_API_KEY", "")
    except Exception:
        pass

# ==============================================================================
# 5. Serve the Exact GUI inside Streamlit
# ==============================================================================
iframe_src = "app/static/index.html"
if api_key:
    iframe_src += f"?gemini_key={api_key}"

if os.path.exists(os.path.join(static_dir, "index.html")):
    components.iframe(src=iframe_src, height=1200, scrolling=True)
else:
    st.error("Application build files not found. Please run 'npm run build' to generate the production assets.")
