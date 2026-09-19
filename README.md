<div align="center">
<img width="1200" height="475" alt="GHBanner" src="https://github.com/user-attachments/assets/0aa67016-6eaf-458a-adb2-6e31a0763ed6" />
</div>

# INTERACTION: AI-Powered Chemical Interaction & Reaction Product Predictor

INTERACTION predicts chemical cross-interactions, reactivity, and degradation pathways (Oxidation, Acidic Hydrolysis, Basic Hydrolysis, Photolysis, and Thermolysis) between chemical compounds. It evaluates reaction pathways using both **Heuristic Kinetic Logic** and **Boltzmann Thermodynamic Distribution ($\Delta G$ at 298.15 K)**.

---

## 🚀 Streamlit Cloud Deployment Guide

This repository is optimized for one-click deployment on **[Streamlit Community Cloud](https://share.streamlit.io/)**.

### Step 1: Deploy on Streamlit Cloud
1. Fork or push this repository to your GitHub account.
2. Go to **[share.streamlit.io](https://share.streamlit.io/)** and click **"Create app"**.
3. Select your repository and branch (`main`).
4. Set the **Main file path** to:
   ```text
   app.py
   ```
5. Click **"Deploy!"**

Streamlit Cloud will automatically detect:
- `requirements.txt` to install the chemistry & visualization libraries (`rdkit`, `pandas`, `numpy`, `matplotlib`, `seaborn`).
- `packages.txt` to install the system graphic libraries (`libgl1`, `libglib2.0-0`, `libxrender1`, etc.) required for 2D molecular structure rendering.
- `.streamlit/config.toml` for UI theme presets and optimal server settings.

---

## 💻 Running the Streamlit App Locally

### Prerequisites
- Python 3.10+
- A Google Gemini API Key ([Get one from Google AI Studio](https://aistudio.google.com/))

### Steps:
1. **Clone the repository:**
   ```bash
   git clone <your-repo-url>
   cd <repo-folder>
   ```

2. **Create and activate a virtual environment:**
   ```bash
   python3 -m venv venv
   source venv/bin/activate    # On Windows: venv\Scripts\activate
   ```

3. **Install Python dependencies:**
   ```bash
   pip install -r requirements.txt
   ```

4. **Set your environment variable:**
   ```bash
   export GEMINI_API_KEY="your_gemini_api_key_here"
   # On Windows Command Prompt: set GEMINI_API_KEY=your_gemini_api_key_here
   # On Windows PowerShell: $env:GEMINI_API_KEY="your_gemini_api_key_here"
   ```

5. **Launch the Streamlit app:**
   ```bash
   streamlit run app.py
   ```
   The app will open automatically in your browser at `http://localhost:8501`.

---

## 🧪 Key Capabilities in the Streamlit App

- **Reaction Presets**: Instant one-click loading of benchmark chemical reaction mixtures (e.g. *Aspirin + Magnesium Stearate*, *Metformin + Lactose*).
- **Dual Analytical Framework**:
  - **Heuristic Kinetics**: Evaluates reactive site vulnerability and reaction kinetics.
  - **Boltzmann Thermodynamics**: Evaluates relative formation energy $\Delta G$ (kcal/mol) at $298.15\text{ K}$.
  - **Both**: Comparative dual perspective with ranking and distribution charts.
- **2D Chemical Structure Rendering**: Native RDKit vector SVG generation for all input compounds and predicted transformation products.
- **Physicochemical Properties**: Real-time calculation of Molecular Weight, LogP, TPSA, Rotatable Bonds, and Hydrogen Bond Donors/Acceptors.
- **Excel & JSON Export**: Download complete multi-tab `.xlsx` reaction reports (Executive Summary, Starting Materials, Reaction Products Profile) and raw JSON.

---

## 🌐 Running the React/Vite Version (Alternative)

If you wish to run the Node.js/React frontend locally:

1. Install dependencies:
   ```bash
   npm install
   ```
2. Set your `GEMINI_API_KEY` in `.env.local` or `.env`.
3. Run development server:
   ```bash
   npm run dev
   ```
4. Open `http://localhost:3000`.

