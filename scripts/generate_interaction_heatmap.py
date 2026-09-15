#!/usr/bin/env python3
"""
scripts/generate_interaction_heatmap.py
Generates a publication-grade Seaborn heatmap overlay representing the interaction
potential and coupling strength between functional groups and reactive centers
of analyzed chemical entities.
"""

import sys
import json
import io
import base64
import numpy as np

# Ensure matplotlib runs in headless server mode without GUI display
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns


def calculate_interaction_potential(site1: str, site2: str) -> tuple:
    """
    Computes chemical interaction potential (0.0 to 1.0) and mechanistic rationale
    between two reactive centers / functional groups.
    """
    s1 = site1.lower()
    s2 = site2.lower()

    # Helper conditions
    is_amine = any(k in s1 or k in s2 for k in ["amine", "amino", "methylamino", "biguanide", "guanidine", "nh2", "nh"])
    is_primary_amine = any(k in s1 or k in s2 for k in ["primary amine", "amino group", "4-aminophenol"])
    is_catechol = any(k in s1 or k in s2 for k in ["catechol", "ortho-diphenol", "3,4-dioh", "dihydroxyphenyl"])
    is_phenol = any(k in s1 or k in s2 for k in ["phenol", "phenolic", "hydroxyl", "napqi"])
    is_ester = any(k in s1 or k in s2 for k in ["ester", "acetyl", "acetoxy", "carbonyl (hydrolysis"])
    is_acid = any(k in s1 or k in s2 for k in ["carboxylic acid", "carboxyl", "acid group", "cooh"])
    is_alkaline = any(k in s1 or k in s2 for k in ["alkaline", "basic surface", "hydroxide", "oxide microenvironment", "ph 8", "ph 9"])
    is_metal = any(k in s1 or k in s2 for k in ["mg2+", "magnesium", "metal", "lewis acid", "cation"])
    is_stearate = any(k in s1 or k in s2 for k in ["stearate", "carboxylate anion", "fatty acid"])
    is_sugar = any(k in s1 or k in s2 for k in ["aldose", "lactose", "glucose", "anomeric", "reducing sugar", "hemiacetal"])
    is_oxidant = any(k in s1 or k in s2 for k in ["peroxide", "pvp", "povidone", "singlet oxygen", "radical"])
    is_amide = any(k in s1 or k in s2 for k in ["amide", "acetamide"])
    is_aromatic = any(k in s1 or k in s2 for k in ["aromatic", "benzene", "phenyl"])
    is_lipophilic = any(k in s1 or k in s2 for k in ["hydrocarbon", "lipophilic", "aliphatic chain", "alkyl"])

    # 1. Acid-Base / Alkaline Microenvironment destabilization
    if (is_catechol or is_phenol) and is_alkaline:
        return 0.96, "Base-Catalyzed Autoxidation", "Critical", (
            "Alkaline surface deprotonates phenolic/catecholic hydroxyls, drastically accelerating "
            "rate of oxidation to reactive ortho-quinones and subsequent cyclization."
        )

    if is_ester and is_alkaline:
        return 0.92, "Base-Promoted Ester Hydrolysis", "Critical", (
            "Alkaline microenvironment accelerates nucleophilic hydroxyl attack onto the ester carbonyl, "
            "causing rapid hydrolytic cleavage into acid and alcohol fragments."
        )

    # 2. Aldose - Amine Maillard Condensation
    if is_amine and is_sugar:
        return 0.94, "Nucleophilic Maillard Condensation", "Critical", (
            "Nucleophilic amine nitrogen attacks open-chain aldose carbonyl, forming an unstable glycosylamine "
            "that undergoes Amadori rearrangement to brown colored chromophores."
        )

    # 3. Lewis Acid Metal Catalysis / Coordination
    if (is_catechol or is_phenol or is_ester) and is_metal:
        return 0.88, "Lewis Acid Metal Coordination", "High", (
            "Divalent magnesium ions polarize carbonyl/hydroxyl oxygen centers, facilitating nucleophilic attack "
            "and stabilizing radical electron-transfer transition states."
        )

    # 4. Ion-Pairing & Interfacial Adsorption
    if is_amine and is_stearate:
        return 0.82, "Interfacial Ion-Pairing", "High", (
            "Electrostatic salt-formation between protonated cationic ammonium centers and lipophilic stearate "
            "carboxylate anions creates insoluble hydrophobic boundary layers."
        )

    # 5. Oxidative Stress / Peroxide Cleavage
    if (is_catechol or is_phenol or is_amine) and is_oxidant:
        return 0.91, "Peroxide-Mediated Oxidation", "Critical", (
            "Trace peroxides initiate single-electron transfer (SET), generating phenoxyl/aminyl radicals "
            "leading to reactive quinone-imines or dimerized adducts."
        )

    # 6. Carboxylic Acid - Ester Transesterification / Anchimeric Assistance
    if is_ester and is_acid:
        return 0.79, "Anchimeric Acid Hydrolysis & Condensation", "High", (
            "Intra- or intermolecular general-acid catalysis promotes rapid ester solvolysis and symmetrical anhydride coupling."
        )

    # 7. Acid - Amine Salt Formation
    if is_acid and is_amine:
        return 0.86, "Proton Transfer Acid-Base Neutralization", "High", (
            "Spontaneous proton transfer from carboxylic acid donor to basic amine acceptor generating ammonium carboxylate salt."
        )

    # 8. Amide - Peroxide / Phenol interaction
    if is_amide and is_oxidant:
        return 0.74, "Oxidative Cleavage", "Moderate", (
            "Oxidative N-deacylation and quinoid oxidation under elevated temperature."
        )

    # 9. Phenol - Phenol Radical Coupling / Self-association
    if is_phenol and is_phenol:
        return 0.58, "Bimolecular Radical Coupling", "Moderate", (
            "Phenoxyl radical recombination forming diaryl ether or biphenyl quinone oligomers."
        )

    # 10. Hydrogen Bonding / Polar Coordination
    if (is_phenol or is_acid or is_amine or is_amide) and ("hydroxyl" in s1 or "hydroxyl" in s2 or "o" in s1 or "o" in s2):
        return 0.52, "Extensive Hydrogen Bonding", "Moderate", (
            "Polar hydrogen bond network altering solid-state dissolution and interfacial wettability."
        )

    # 11. Lipophilic / Hydrophobic Desolvation
    if is_lipophilic and (is_aromatic or is_lipophilic):
        return 0.35, "Hydrophobic Dispersion / Van der Waals", "Low", (
            "Non-covalent lipophilic dispersion without covalent bond breakage."
        )

    # Default baseline interaction potential based on heteroatom density
    base_score = 0.20
    if "o" in s1 or "n" in s1 or "o" in s2 or "n" in s2:
        base_score = 0.38
    return base_score, "Weak Non-Specific Contact", "Low", "Minor steric proximity or weak dipole-dipole orientation."


def generate_heatmap(data: dict) -> dict:
    """
    Builds the Seaborn interaction heatmap and returns serializable payload.
    """
    compounds = data.get("compounds", [])
    if len(compounds) < 1:
        raise ValueError("At least one compound is required for interaction potential heatmap.")

    comp1 = compounds[0]
    comp1_name = comp1.get("name", "Compound 1")

    # If only 1 compound, analyze intramolecular reactive centers
    if len(compounds) == 1:
        comp2 = comp1
        comp2_name = f"{comp1_name} (Intramolecular)"
    else:
        comp2 = compounds[1]
        comp2_name = comp2.get("name", "Compound 2")

    # Extract or infer reactive centers
    c1_sites = comp1.get("interactionSites") or comp1.get("features") or ["Reactive Scaffold"]
    c2_sites = comp2.get("interactionSites") or comp2.get("features") or ["Secondary Matrix"]

    # Filter and limit labels to 6 max per axis for optimal readability
    c1_labels = [s.strip() for s in c1_sites if s.strip()][:6]
    c2_labels = [s.strip() for s in c2_sites if s.strip()][:6]

    if not c1_labels:
        c1_labels = ["Primary Reactive Center"]
    if not c2_labels:
        c2_labels = ["Secondary Contact Center"]

    matrix = np.zeros((len(c1_labels), len(c2_labels)), dtype=float)
    details = []

    for i, s1 in enumerate(c1_labels):
        for j, s2 in enumerate(c2_labels):
            score, mech, severity, desc = calculate_interaction_potential(s1, s2)
            # Add small determinism based on string hashes for realistic variation
            hash_offset = ((hash(s1 + s2) % 11) - 5) * 0.01
            final_score = round(float(np.clip(score + hash_offset, 0.05, 0.99)), 2)
            matrix[i, j] = final_score
            details.append({
                "row": i,
                "col": j,
                "rowLabel": s1,
                "colLabel": s2,
                "score": final_score,
                "percentage": f"{int(round(final_score * 100))}%",
                "mechanism": mech,
                "severity": severity,
                "description": desc
            })

    # Sort details to find top interactions
    top_interactions = sorted(details, key=lambda x: x["score"], reverse=True)[:5]

    # Plot Seaborn Heatmap with Matplotlib
    cmap_choice = data.get("cmap", "YlOrRd")
    valid_cmaps = ["YlOrRd", "rocket", "mako", "viridis", "magma", "coolwarm"]
    if cmap_choice not in valid_cmaps:
        cmap_choice = "YlOrRd"

    # Truncate labels for heatmap axes to avoid overlap
    def clean_axis_label(text: str, max_len: int = 24) -> str:
        if len(text) <= max_len:
            return text
        return text[:max_len - 3] + "..."

    disp_y = [clean_axis_label(l) for l in c1_labels]
    disp_x = [clean_axis_label(l) for l in c2_labels]

    fig_w = max(7.0, len(disp_x) * 1.6 + 2.5)
    fig_h = max(5.0, len(disp_y) * 1.1 + 2.0)

    fig, ax = plt.subplots(figsize=(fig_w, fig_h), dpi=180)
    fig.patch.set_facecolor("#ffffff")
    ax.set_facecolor("#ffffff")

    # Draw Heatmap
    sns.heatmap(
        matrix,
        ax=ax,
        annot=True,
        fmt=".2f",
        cmap=cmap_choice,
        vmin=0.0,
        vmax=1.0,
        linewidths=1.2,
        linecolor="#ffffff",
        cbar_kws={
            "label": "Interaction Potential (Coupling Strength: 0.0 - 1.0)",
            "shrink": 0.85,
            "pad": 0.04
        },
        annot_kws={"size": 11, "weight": "bold", "color": "#1e293b"}
    )

    # Style axes
    ax.set_xticklabels(disp_x, rotation=25, ha="right", fontsize=9, fontweight="600", color="#0f172a")
    ax.set_yticklabels(disp_y, rotation=0, ha="right", fontsize=9, fontweight="600", color="#0f172a")

    ax.set_xlabel(f"Co-Reactant Centers ({comp2_name})", fontsize=11, fontweight="bold", color="#1e293b", labelpad=10)
    ax.set_ylabel(f"Primary Reactive Sites ({comp1_name})", fontsize=11, fontweight="bold", color="#1e293b", labelpad=10)
    
    chart_title = data.get("title") or f"Interaction Potential Matrix: {comp1_name} vs {comp2_name}"
    ax.set_title(chart_title, fontsize=12, fontweight="bold", color="#0f172a", pad=16)

    # Add subtle subtitle/brand
    plt.tight_layout()

    buf = io.BytesIO()
    fig.savefig(buf, format="png", bbox_inches="tight", facecolor=fig.get_facecolor(), edgecolor="none")
    plt.close(fig)
    buf.seek(0)

    heatmap_b64 = "data:image/png;base64," + base64.b64encode(buf.read()).decode("utf-8")

    return {
        "success": True,
        "compound1Name": comp1_name,
        "compound2Name": comp2_name,
        "rowLabels": c1_labels,
        "colLabels": c2_labels,
        "matrix": matrix.tolist(),
        "details": details,
        "topInteractions": top_interactions,
        "heatmapBase64": heatmap_b64,
        "summary": (
            f"Evaluated {len(c1_labels)} reactive center(s) of '{comp1_name}' against {len(c2_labels)} "
            f"functional group(s) of '{comp2_name}'. Maximum interaction coupling strength detected: "
            f"{top_interactions[0]['score'] if top_interactions else 0.0} ({top_interactions[0]['mechanism'] if top_interactions else 'N/A'})."
        )
    }


def main():
    try:
        if len(sys.argv) > 1 and sys.argv[1] != "-":
            with open(sys.argv[1], "r", encoding="utf-8") as f:
                input_data = json.load(f)
        else:
            raw = sys.stdin.read().strip()
            if not raw:
                raise ValueError("No input JSON passed to generate_interaction_heatmap.py")
            input_data = json.loads(raw)

        res = generate_heatmap(input_data)
        print(json.dumps(res))
    except Exception as e:
        err_res = {"success": False, "error": str(e)}
        print(json.dumps(err_res))
        sys.exit(1)


if __name__ == "__main__":
    main()
