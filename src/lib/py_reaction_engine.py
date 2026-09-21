"""
reaction_engine.py - General-purpose functional-group reaction engine (RDKit based)
====================================================================================

Domain-agnostic: works for any organic / organometallic / inorganic species given as SMILES.
It replaces the regex/keyword engine of the Streamlit app (identify_functional_groups and the
candidate generation inside predict_degradation_and_reactions).

Pipeline
  1. RDKit graph perception of functional groups with SMARTS (independent of how a SMILES is
     written; aromatic and Kekule input give the same answer).
  2. Per-group reactivity profile for 7 conditions (acid, base, neutral hydrolysis, photolysis,
     thermal, oxidation, cross-interaction), with variants (aryl vs alkyl ester, ring strain,
     tert-alkyl, hindered phenol, ...).
  3. Reaction templates (reaction SMARTS + a few RDKit graph edits) generate valid, atom-balanced
     products; every product is sanitised, checked for element conservation and de-duplicated.
  4. Co-reactants are matched group-to-group (acyl transfer, Maillard, alkylation, Michael,
     nitrosation, salts, metal carboxylates, redox, pH / moisture / peroxide / photocatalyst effects).
  5. Ranking: heuristic likelihood from the vulnerability grade + Boltzmann weights at 298.15 K.

Public functions
  identify_functional_groups(smiles)                       -> list of group dicts (app format)
  predict(primary_smiles, secondary_smiles_list, method)   -> dict(impurities, functional_groups,
                                                              chain_of_thought, ...)
Limits: stereochemistry is not propagated to products; free energies are reaction-class estimates
on a compressed kcal/mol scale (not QM); names are descriptive, not systematic IUPAC.
"""
from __future__ import annotations

import math
import re
from typing import Any, Callable, Dict, Iterable, List, Optional, Sequence, Set, Tuple

from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem, Descriptors, rdMolDescriptors

RDLogger.DisableLog("rdApp.*")

# --------------------------------------------------------------------------------------
# Vulnerability grading
# --------------------------------------------------------------------------------------
VORDER = ["Resistant", "Low", "Moderate", "High", "Critical"]
VSCORE = {"Critical": 0.93, "High": 0.80, "Moderate": 0.58, "Low": 0.30, "Resistant": 0.10}
CKEYS = ("acid", "base", "hyd", "photo", "therm", "ox")
COND_LABEL = {
    "acid": "Acidic Hydrolysis", "base": "Basic Hydrolysis", "hyd": "Hydrolysis",
    "photo": "Photodegradation", "therm": "Thermal Degradation", "ox": "Oxidation",
}
COND_TITLE = {
    "acid": "Acidic stress", "base": "Basic stress", "hyd": "Neutral hydrolysis",
    "photo": "Photolytic stress", "therm": "Thermal stress", "ox": "Oxidative stress",
}
SRC_STRESS = "Stress degradation"
SRC_CROSS = "Interaction with other compound"


def vrank(v: str) -> int:
    return VORDER.index(v) if v in VORDER else 0


def step_vuln(v: str, d: int) -> str:
    return VORDER[max(0, min(4, vrank(v) + d))]


def clean_text(t: str) -> str:
    """The app injects text into HTML unescaped: keep it free of markup characters."""
    return (t or "").replace("->", "\u2192").replace("<", "below ").replace(">", "above ").replace("&", "and")


# --------------------------------------------------------------------------------------
# RDKit helpers
# --------------------------------------------------------------------------------------
def parse_smiles(smi: str) -> Optional[Chem.Mol]:
    s = (smi or "").strip()
    if not s:
        return None
    m = Chem.MolFromSmiles(s)
    if m is None:
        relaxed = re.sub(r"[@\\/]", "", s)
        if relaxed and relaxed != s:
            m = Chem.MolFromSmiles(relaxed)
    return m


def canon(m: Chem.Mol) -> str:
    return Chem.MolToSmiles(m, isomericSmiles=False)


def frag_mols(m: Chem.Mol) -> List[Chem.Mol]:
    return list(Chem.GetMolFrags(m, asMols=True, sanitizeFrags=True))


def formula_of(m: Chem.Mol) -> str:
    try:
        return rdMolDescriptors.CalcMolFormula(m)
    except Exception:
        return ""


def mono_mass(m: Chem.Mol) -> float:
    try:
        return float(Descriptors.ExactMolWt(m))
    except Exception:
        return 0.0


def element_counts(m: Chem.Mol) -> Dict[str, int]:
    c: Dict[str, int] = {}
    for a in m.GetAtoms():
        c[a.GetSymbol()] = c.get(a.GetSymbol(), 0) + 1
    return c


# --------------------------------------------------------------------------------------
# Compound library (used to give products readable names)
# --------------------------------------------------------------------------------------
LIB_ROWS = [
    ('O', 'Water', 'Solvent'),
    ('CO', 'Methanol', 'Alcohol'),
    ('CCO', 'Ethanol', 'Alcohol'),
    ('CC(C)O', 'Propan-2-ol', 'Alcohol'),
    ('CCCO', 'Propan-1-ol', 'Alcohol'),
    ('CCCCO', 'Butan-1-ol', 'Alcohol'),
    ('CC(C)(C)O', 'tert-Butanol', 'Alcohol'),
    ('OC1CCCCC1', 'Cyclohexanol', 'Alcohol'),
    ('OCc1ccccc1', 'Benzyl alcohol', 'Alcohol'),
    ('OCCO', 'Ethylene glycol', 'Polyol'),
    ('CC(O)CO', 'Propylene glycol', 'Polyol'),
    ('OCC(O)CO', 'Glycerol', 'Polyol'),
    ('OCCOCCO', 'Diethylene glycol', 'Polyether'),
    ('OCCOCCOCCO', 'Polyethylene glycol (triethylene glycol model)', 'Polyether'),
    ('CC(C)=O', 'Acetone', 'Ketone'),
    ('CCC(C)=O', 'Butan-2-one', 'Ketone'),
    ('O=C1CCCCC1', 'Cyclohexanone', 'Ketone'),
    ('CC(=O)c1ccccc1', 'Acetophenone', 'Ketone'),
    ('O=C(c1ccccc1)c1ccccc1', 'Benzophenone', 'Ketone'),
    ('CC#N', 'Acetonitrile', 'Nitrile'),
    ('N#Cc1ccccc1', 'Benzonitrile', 'Nitrile'),
    ('CN(C)C=O', 'N,N-Dimethylformamide', 'Amide'),
    ('CN(C)C(C)=O', 'N,N-Dimethylacetamide', 'Amide'),
    ('CN1CCCC1=O', 'N-Methyl-2-pyrrolidone', 'Lactam'),
    ('CS(C)=O', 'Dimethyl sulfoxide', 'Sulfoxide'),
    ('CS(C)(=O)=O', 'Dimethyl sulfone', 'Sulfone'),
    ('C1CCOC1', 'Tetrahydrofuran', 'Ether'),
    ('C1COCCO1', '1,4-Dioxane', 'Ether'),
    ('CCOCC', 'Diethyl ether', 'Ether'),
    ('COC(C)(C)C', 'Methyl tert-butyl ether', 'Ether'),
    ('COc1ccccc1', 'Anisole', 'Ether'),
    ('CCOC(C)=O', 'Ethyl acetate', 'Ester'),
    ('COC(C)=O', 'Methyl acetate', 'Ester'),
    ('CCCCOC(C)=O', 'Butyl acetate', 'Ester'),
    ('CC(=O)OC(C)(C)C', 'tert-Butyl acetate', 'Ester'),
    ('CC(=O)Oc1ccccc1', 'Phenyl acetate', 'Ester'),
    ('CC(=O)OC=C', 'Vinyl acetate', 'Ester'),
    ('COC(=O)c1ccccc1', 'Methyl benzoate', 'Ester'),
    ('CCOC=O', 'Ethyl formate', 'Ester'),
    ('C=CC(=O)OCC', 'Ethyl acrylate', 'Michael acceptor'),
    ('COC(=O)C(C)=C', 'Methyl methacrylate', 'Michael acceptor'),
    ('COC(=O)OC', 'Dimethyl carbonate', 'Carbonate'),
    ('CCOC(=O)c1ccccc1C(=O)OCC', 'Diethyl phthalate', 'Ester'),
    ('CC(=O)OCC(COC(C)=O)OC(C)=O', 'Triacetin', 'Ester'),
    ('CC(=O)Oc1ccccc1C(=O)O', '2-Acetoxybenzoic acid (acetylsalicylic acid)', 'Aryl ester / acid'),
    ('ClCCl', 'Dichloromethane', 'Alkyl halide'),
    ('ClC(Cl)Cl', 'Chloroform', 'Alkyl halide'),
    ('ClC(Cl)(Cl)Cl', 'Carbon tetrachloride', 'Alkyl halide'),
    ('ClCCCl', '1,2-Dichloroethane', 'Alkyl halide'),
    ('CCCCCC', 'Hexane', 'Hydrocarbon'),
    ('CCCCCCC', 'Heptane', 'Hydrocarbon'),
    ('C1CCCCC1', 'Cyclohexane', 'Hydrocarbon'),
    ('C', 'Methane', 'Hydrocarbon'),
    ('C=C', 'Ethene', 'Alkene'),
    ('CC=C', 'Propene', 'Alkene'),
    ('CC(C)=C', 'Isobutene', 'Alkene'),
    ('CCCCC=C', '1-Hexene', 'Alkene'),
    ('C1=CCCCC1', 'Cyclohexene', 'Alkene'),
    ('C=Cc1ccccc1', 'Styrene', 'Alkene'),
    ('C=CC=C', '1,3-Butadiene', 'Diene'),
    ('C1=CCC=C1', 'Cyclopentadiene', 'Diene'),
    ('C#C', 'Ethyne', 'Alkyne'),
    ('C#Cc1ccccc1', 'Phenylacetylene', 'Alkyne'),
    ('c1ccccc1', 'Benzene', 'Aromatic'),
    ('Cc1ccccc1', 'Toluene', 'Aromatic'),
    ('Cc1ccccc1C', 'o-Xylene', 'Aromatic'),
    ('CCc1ccccc1', 'Ethylbenzene', 'Aromatic'),
    ('CC(C)c1ccccc1', 'Cumene', 'Aromatic'),
    ('c1ccc2ccccc2c1', 'Naphthalene', 'Aromatic'),
    ('c1ccc(cc1)-c1ccccc1', 'Biphenyl', 'Aromatic'),
    ('Clc1ccccc1', 'Chlorobenzene', 'Aryl halide'),
    ('Brc1ccccc1', 'Bromobenzene', 'Aryl halide'),
    ('Ic1ccccc1', 'Iodobenzene', 'Aryl halide'),
    ('Fc1ccccc1', 'Fluorobenzene', 'Aryl halide'),
    ('FC(F)(F)c1ccccc1', 'Benzotrifluoride', 'Aryl halide'),
    ('Clc1ccc(cc1[N+]([O-])=O)[N+]([O-])=O', '1-Chloro-2,4-dinitrobenzene', 'Aryl halide'),
    ('[O-][N+](=O)c1ccccc1', 'Nitrobenzene', 'Nitro'),
    ('Oc1ccc(cc1)[N+]([O-])=O', '4-Nitrophenol', 'Phenol'),
    ('C[N+]([O-])=O', 'Nitromethane', 'Nitro'),
    ('CC[N+]([O-])=O', 'Nitroethane', 'Nitro'),
    ('CC(=O)O', 'Acetic acid', 'Carboxylic acid'),
    ('OC=O', 'Formic acid', 'Carboxylic acid'),
    ('CCC(=O)O', 'Propanoic acid', 'Carboxylic acid'),
    ('CCCC(=O)O', 'Butanoic acid', 'Carboxylic acid'),
    ('OC(=O)c1ccccc1', 'Benzoic acid', 'Carboxylic acid'),
    ('OC(=O)c1ccccc1O', 'Salicylic acid', 'Carboxylic acid / phenol'),
    ('OC(=O)c1cc(O)ccc1O', 'Gentisic acid', 'Carboxylic acid / phenol'),
    ('OC(=O)c1cc(O)c(O)c(O)c1', 'Gallic acid', 'Carboxylic acid / phenol'),
    ('OC(=O)c1ccccc1C(O)=O', 'Phthalic acid', 'Carboxylic acid'),
    ('OC(=O)C=Cc1ccccc1', 'Cinnamic acid', 'Michael acceptor'),
    ('OC(=O)C(O)=O', 'Oxalic acid', 'Carboxylic acid'),
    ('OC(=O)CC(O)=O', 'Malonic acid', 'Carboxylic acid'),
    ('OC(=O)CCC(O)=O', 'Succinic acid', 'Carboxylic acid'),
    ('OC(=O)CCCCC(O)=O', 'Adipic acid', 'Carboxylic acid'),
    ('OC(=O)C=CC(O)=O', 'Maleic acid', 'Michael acceptor'),
    ('OC(=O)C=CC(O)=O', 'Fumaric acid', 'Michael acceptor'),
    ('CC(=O)CC(O)=O', 'Acetoacetic acid', 'Keto acid'),
    ('CC(=O)C(O)=O', 'Pyruvic acid', 'Keto acid'),
    ('CC(O)C(O)=O', 'Lactic acid', 'Hydroxy acid'),
    ('OCC(O)=O', 'Glycolic acid', 'Hydroxy acid'),
    ('OC(C(O)C(O)=O)C(O)=O', 'Tartaric acid', 'Hydroxy acid'),
    ('OC(=O)CC(O)(CC(O)=O)C(O)=O', 'Citric acid', 'Hydroxy acid'),
    ('OCC(O)C1OC(=O)C(O)=C1O', 'Ascorbic acid', 'Enediol lactone'),
    ('OC(=O)C(F)(F)F', 'Trifluoroacetic acid', 'Carboxylic acid'),
    ('OC(=O)CCl', 'Chloroacetic acid', 'Carboxylic acid'),
    ('OC(=O)C(Cl)(Cl)Cl', 'Trichloroacetic acid', 'Carboxylic acid'),
    ('CCCCCCCCCCCC(O)=O', 'Lauric acid', 'Fatty acid'),
    ('CCCCCCCCCCCCCCCC(O)=O', 'Palmitic acid', 'Fatty acid'),
    ('CCCCCCCCCCCCCCCCCC(O)=O', 'Stearic acid', 'Fatty acid'),
    ('CCCCCCCCC=CCCCCCCCC(O)=O', 'Oleic acid', 'Fatty acid'),
    ('CCCCCC=CCC=CCCCCCCCC(O)=O', 'Linoleic acid', 'Fatty acid'),
    ('CCCCCCCCC=CCCCCCCCC(=O)OC', 'Methyl oleate', 'Ester'),
    ('CC([O-])=O.[Na+]', 'Sodium acetate', 'Carboxylate salt'),
    ('[O-]C(=O)c1ccccc1.[Na+]', 'Sodium benzoate', 'Carboxylate salt'),
    ('OC(CC([O-])=O)(CC([O-])=O)C([O-])=O.[Na+].[Na+].[Na+]', 'Trisodium citrate', 'Carboxylate salt'),
    ('[Mg+2].[O-]C(=O)CCCCCCCCCCCCCCCCC.[O-]C(=O)CCCCCCCCCCCCCCCCC', 'Magnesium stearate', 'Metal carboxylate'),
    ('[Ca+2].[O-]C(=O)CCCCCCCCCCCCCCCCC.[O-]C(=O)CCCCCCCCCCCCCCCCC', 'Calcium stearate', 'Metal carboxylate'),
    ('[Na+].[O-]C(=O)CCCCCCCCCCCCCCCCC', 'Sodium stearate', 'Metal carboxylate'),
    ('CCCCCCCCCCCCCCCCCCOC(=O)C=CC([O-])=O.[Na+]', 'Sodium stearyl fumarate', 'Ester / carboxylate'),
    ('CC(=O)OC(C)=O', 'Acetic anhydride', 'Anhydride'),
    ('O=C1CCC(=O)O1', 'Succinic anhydride', 'Anhydride'),
    ('O=C1C=CC(=O)O1', 'Maleic anhydride', 'Anhydride'),
    ('O=C1OC(=O)c2ccccc12', 'Phthalic anhydride', 'Anhydride'),
    ('CC(Cl)=O', 'Acetyl chloride', 'Acyl halide'),
    ('O=C(Cl)c1ccccc1', 'Benzoyl chloride', 'Acyl halide'),
    ('O=C(OOC(=O)c1ccccc1)c1ccccc1', 'Benzoyl peroxide', 'Peroxide'),
    ('O=C1CCCO1', 'gamma-Butyrolactone', 'Lactone'),
    ('O=C1CCCCO1', 'delta-Valerolactone', 'Lactone'),
    ('O=C1CCO1', 'beta-Propiolactone', 'Lactone'),
    ('O=C1CCCCCO1', 'epsilon-Caprolactone', 'Lactone'),
    ('O=c1ccc2ccccc2o1', 'Coumarin', 'Lactone'),
    ('O=C1CCCCCN1', 'epsilon-Caprolactam', 'Lactam'),
    ('O=C1CCCN1', '2-Pyrrolidone', 'Lactam'),
    ('O=C1CCN1', '2-Azetidinone (beta-lactam)', 'Lactam'),
    ('C=CN1CCCC1=O', 'N-Vinylpyrrolidone (PVP repeat-unit model)', 'Lactam'),
    ('O=C1CCC(=O)N1', 'Succinimide', 'Imide'),
    ('O=C1NC(=O)c2ccccc12', 'Phthalimide', 'Imide'),
    ('O=C1C=CC(=O)N1', 'Maleimide', 'Imide'),
    ('CC(N)=O', 'Acetamide', 'Amide'),
    ('NC=O', 'Formamide', 'Amide'),
    ('NC(=O)c1ccccc1', 'Benzamide', 'Amide'),
    ('CC(=O)Nc1ccccc1', 'Acetanilide', 'Anilide'),
    ('CC(=O)Nc1ccc(O)cc1', 'N-(4-Hydroxyphenyl)acetamide', 'Anilide / phenol'),
    ('NC(N)=O', 'Urea', 'Urea'),
    ('NC(N)=S', 'Thiourea', 'Thiocarbonyl'),
    ('CCOC(N)=O', 'Ethyl carbamate', 'Carbamate'),
    ('CC(C)(C)OC(=O)OC(=O)OC(C)(C)C', 'Di-tert-butyl dicarbonate', 'Anhydride'),
    ('CC(C)(C)OC(N)=O', 'tert-Butyl carbamate', 'Carbamate'),
    ('C=O', 'Formaldehyde', 'Aldehyde'),
    ('CC=O', 'Acetaldehyde', 'Aldehyde'),
    ('CCC=O', 'Propanal', 'Aldehyde'),
    ('O=Cc1ccccc1', 'Benzaldehyde', 'Aldehyde'),
    ('O=CC=O', 'Glyoxal', 'Aldehyde'),
    ('C=CC=O', 'Acrolein', 'Michael acceptor'),
    ('CC(=O)C=C', 'Methyl vinyl ketone', 'Michael acceptor'),
    ('COc1cc(C=O)ccc1O', 'Vanillin', 'Aldehyde / phenol'),
    ('O=C1C=CC(=O)C=C1', '1,4-Benzoquinone', 'Quinone'),
    ('OCC1OC(O)C(O)C(O)C1O', 'Glucose (pyranose)', 'Reducing sugar'),
    ('OCC1(O)OCC(O)C(O)C1O', 'Fructose (pyranose)', 'Reducing sugar'),
    ('OCC1OC(O)C(O)C(O)C1O', 'Galactose (pyranose)', 'Reducing sugar'),
    ('OCC1OC(OC2C(CO)OC(O)C(O)C2O)C(O)C(O)C1O', 'Lactose', 'Reducing disaccharide'),
    ('OCC1OC(OC2C(CO)OC(O)C(O)C2O)C(O)C(O)C1O', 'Maltose', 'Reducing disaccharide'),
    ('OCC1OC(OC2(CO)OC(CO)C(O)C2O)C(O)C(O)C1O', 'Sucrose', 'Non-reducing disaccharide'),
    ('OCC(O)C(O)C(O)C(O)CO', 'Mannitol', 'Polyol'),
    ('OCC(O)C(O)C(O)C(O)CO', 'Sorbitol', 'Polyol'),
    ('OCC(O)C(O)C(O)CO', 'Xylitol', 'Polyol'),
    ('OCC1OC(O)C(O)C(O)C1O', 'Glucose repeat unit (cellulose / starch model)', 'Carbohydrate'),
    ('N', 'Ammonia', 'Amine'),
    ('CN', 'Methylamine', 'Amine'),
    ('CNC', 'Dimethylamine', 'Amine'),
    ('CN(C)C', 'Trimethylamine', 'Amine'),
    ('CCN', 'Ethylamine', 'Amine'),
    ('CCNCC', 'Diethylamine', 'Amine'),
    ('CCN(CC)CC', 'Triethylamine', 'Amine'),
    ('CCN(C(C)C)C(C)C', 'N,N-Diisopropylethylamine', 'Amine'),
    ('NCCN', 'Ethylenediamine', 'Amine'),
    ('NCCO', 'Ethanolamine', 'Amino alcohol'),
    ('OCCNCCO', 'Diethanolamine', 'Amino alcohol'),
    ('OCCN(CCO)CCO', 'Triethanolamine', 'Amino alcohol'),
    ('NC(CO)(CO)CO', 'Tris(hydroxymethyl)aminomethane', 'Amino alcohol'),
    ('CNCC(O)C(O)C(O)C(O)CO', 'N-Methylglucamine (meglumine)', 'Amino polyol'),
    ('C1CCNCC1', 'Piperidine', 'Amine'),
    ('C1COCCN1', 'Morpholine', 'Amine'),
    ('C1CNCCN1', 'Piperazine', 'Amine'),
    ('C1CCNC1', 'Pyrrolidine', 'Amine'),
    ('C1CCC2=NCCCN2CC1', 'DBU', 'Amidine base'),
    ('NCc1ccccc1', 'Benzylamine', 'Amine'),
    ('Nc1ccccc1', 'Aniline', 'Aromatic amine'),
    ('CNc1ccccc1', 'N-Methylaniline', 'Aromatic amine'),
    ('CN(C)c1ccccc1', 'N,N-Dimethylaniline', 'Aromatic amine'),
    ('Nc1ccc(O)cc1', '4-Aminophenol', 'Aromatic amine / phenol'),
    ('Nc1ccc(cc1)[N+]([O-])=O', '4-Nitroaniline', 'Aromatic amine'),
    ('NC(N)=N', 'Guanidine', 'Guanidine'),
    ('NN', 'Hydrazine', 'Hydrazine'),
    ('NNc1ccccc1', 'Phenylhydrazine', 'Hydrazine'),
    ('NO', 'Hydroxylamine', 'Hydroxylamine'),
    ('CC(C)=NO', 'Acetone oxime', 'Oxime'),
    ('ON=Cc1ccccc1', 'Benzaldehyde oxime', 'Oxime'),
    ('CN(C)N=O', 'N-Nitrosodimethylamine', 'N-Nitrosamine'),
    ('CCN(CC)N=O', 'N-Nitrosodiethylamine', 'N-Nitrosamine'),
    ('c1ccc(cc1)N=Nc1ccccc1', 'Azobenzene', 'Azo'),
    ('[N-]=[N+]=Nc1ccccc1', 'Phenyl azide', 'Azide'),
    ('CN=C=O', 'Methyl isocyanate', 'Isocyanate'),
    ('S=C=Nc1ccccc1', 'Phenyl isothiocyanate', 'Isothiocyanate'),
    ('NCC(O)=O', 'Glycine', 'Amino acid'),
    ('CC(N)C(O)=O', 'Alanine', 'Amino acid'),
    ('NCCCCC(N)C(O)=O', 'Lysine', 'Amino acid'),
    ('NC(CCC(O)=O)C(O)=O', 'Glutamic acid', 'Amino acid'),
    ('NC(CS)C(O)=O', 'Cysteine', 'Amino acid'),
    ('CSCCC(N)C(O)=O', 'Methionine', 'Amino acid'),
    ('NC(Cc1ccccc1)C(O)=O', 'Phenylalanine', 'Amino acid'),
    ('NC(Cc1ccc(O)cc1)C(O)=O', 'Tyrosine', 'Amino acid'),
    ('NC(Cc1c[nH]c2ccccc12)C(O)=O', 'Tryptophan', 'Amino acid'),
    ('OC(=O)C1CCCN1', 'Proline', 'Amino acid'),
    ('CCCCCCCCCCCCCCCC[N+](C)(C)C.[Br-]', 'Cetyltrimethylammonium bromide', 'Quaternary ammonium'),
    ('CCCCCCCCCCCC[N+](C)(C)Cc1ccccc1.[Cl-]', 'Benzalkonium chloride (C12)', 'Quaternary ammonium'),
    ('C[N+](C)(C)C.[OH-]', 'Tetramethylammonium hydroxide', 'Strong base'),
    ('Oc1ccccc1', 'Phenol', 'Phenol'),
    ('Oc1ccccc1O', 'Catechol', 'Phenol'),
    ('Oc1cccc(O)c1', 'Resorcinol', 'Phenol'),
    ('Oc1ccc(O)cc1', 'Hydroquinone', 'Phenol'),
    ('Cc1ccc(O)cc1', 'p-Cresol', 'Phenol'),
    ('Cc1cc(c(O)c(c1)C(C)(C)C)C(C)(C)C', 'Butylated hydroxytoluene', 'Hindered phenol'),
    ('c1ccc2[nH]ccc2c1', 'Indole', 'Heteroaromatic'),
    ('c1ccncc1', 'Pyridine', 'Heteroaromatic'),
    ('c1cc[nH]c1', 'Pyrrole', 'Heteroaromatic'),
    ('c1ccoc1', 'Furan', 'Heteroaromatic'),
    ('c1ccsc1', 'Thiophene', 'Heteroaromatic'),
    ('c1c[nH]cn1', 'Imidazole', 'Heteroaromatic'),
    ('Cn1cnc2c1c(=O)n(C)c(=O)n2C', '1,3,7-Trimethylxanthine', 'Heteroaromatic'),
    ('ClCc1ccccc1', 'Benzyl chloride', 'Alkyl halide'),
    ('CCCl', 'Chloroethane', 'Alkyl halide'),
    ('CCCCBr', '1-Bromobutane', 'Alkyl halide'),
    ('CC(C)(C)Cl', 'tert-Butyl chloride', 'Alkyl halide'),
    ('C=CCBr', 'Allyl bromide', 'Alkyl halide'),
    ('CI', 'Iodomethane', 'Alkyl halide'),
    ('C1CO1', 'Ethylene oxide', 'Epoxide'),
    ('CC1CO1', 'Propylene oxide', 'Epoxide'),
    ('ClCC1CO1', 'Epichlorohydrin', 'Epoxide'),
    ('COS(=O)(=O)c1ccc(C)cc1', 'Methyl p-toluenesulfonate', 'Sulfonate ester'),
    ('COS(=O)(=O)OC', 'Dimethyl sulfate', 'Sulfate ester'),
    ('COS(C)(=O)=O', 'Methyl methanesulfonate', 'Sulfonate ester'),
    ('CCOS(C)(=O)=O', 'Ethyl methanesulfonate', 'Sulfonate ester'),
    ('CS', 'Methanethiol', 'Thiol'),
    ('CCS', 'Ethanethiol', 'Thiol'),
    ('Sc1ccccc1', 'Thiophenol', 'Thiol'),
    ('OC(CS)C(O)CS', 'Dithiothreitol', 'Thiol'),
    ('CSC', 'Dimethyl sulfide', 'Thioether'),
    ('CSSC', 'Dimethyl disulfide', 'Disulfide'),
    ('c1ccc(cc1)SSc1ccccc1', 'Diphenyl disulfide', 'Disulfide'),
    ('NS(=O)(=O)c1ccccc1', 'Benzenesulfonamide', 'Sulfonamide'),
    ('Cc1ccc(cc1)S(N)(=O)=O', 'p-Toluenesulfonamide', 'Sulfonamide'),
    ('CS(O)(=O)=O', 'Methanesulfonic acid', 'Sulfonic acid'),
    ('Cc1ccc(cc1)S(O)(=O)=O', 'p-Toluenesulfonic acid', 'Sulfonic acid'),
    ('CCCCCCCCCCCCOS([O-])(=O)=O.[Na+]', 'Sodium dodecyl sulfate', 'Alkyl sulfate'),
    ('c1ccc(cc1)P(c1ccccc1)c1ccccc1', 'Triphenylphosphine', 'Phosphine'),
    ('COP(=O)(OC)OC', 'Trimethyl phosphate', 'Phosphate ester'),
    ('Cl', 'Hydrochloric acid', 'Strong acid'),
    ('Br', 'Hydrobromic acid', 'Strong acid'),
    ('OS(O)(=O)=O', 'Sulfuric acid', 'Strong acid'),
    ('O[N+]([O-])=O', 'Nitric acid', 'Strong acid / oxidant'),
    ('OP(O)(O)=O', 'Phosphoric acid', 'Acid'),
    ('[Na+].[OH-]', 'Sodium hydroxide', 'Strong base'),
    ('[K+].[OH-]', 'Potassium hydroxide', 'Strong base'),
    ('[Li+].[OH-]', 'Lithium hydroxide', 'Strong base'),
    ('[Ca+2].[OH-].[OH-]', 'Calcium hydroxide', 'Base'),
    ('[Mg+2].[OH-].[OH-]', 'Magnesium hydroxide', 'Base'),
    ('[NH4+].[OH-]', 'Ammonium hydroxide', 'Base'),
    ('[Mg+2].[O-2]', 'Magnesium oxide', 'Metal oxide (basic)'),
    ('[Ca+2].[O-2]', 'Calcium oxide', 'Metal oxide (basic)'),
    ('[Zn+2].[O-2]', 'Zinc oxide', 'Metal oxide'),
    ('O=[Si]=O', 'Silicon dioxide', 'Inert oxide'),
    ('O=[Ti]=O', 'Titanium dioxide', 'Metal oxide (photocatalyst)'),
    ('[Na+].OC([O-])=O', 'Sodium bicarbonate', 'Carbonate / base'),
    ('[Na+].[Na+].[O-]C([O-])=O', 'Sodium carbonate', 'Carbonate / base'),
    ('[K+].[K+].[O-]C([O-])=O', 'Potassium carbonate', 'Carbonate / base'),
    ('[Ca+2].[O-]C([O-])=O', 'Calcium carbonate', 'Carbonate / base'),
    ('[Mg+2].[O-]C([O-])=O', 'Magnesium carbonate', 'Carbonate / base'),
    ('C[O-].[Na+]', 'Sodium methoxide', 'Strong base'),
    ('CC[O-].[Na+]', 'Sodium ethoxide', 'Strong base'),
    ('CC(C)(C)[O-].[K+]', 'Potassium tert-butoxide', 'Strong base'),
    ('[Na+].[Cl-]', 'Sodium chloride', 'Salt'),
    ('[K+].[Cl-]', 'Potassium chloride', 'Salt'),
    ('[Mg+2].[Cl-].[Cl-]', 'Magnesium chloride', 'Salt'),
    ('[Ca+2].[Cl-].[Cl-]', 'Calcium chloride', 'Salt'),
    ('[Zn+2].[Cl-].[Cl-]', 'Zinc chloride', 'Lewis acid'),
    ('[Al+3].[Cl-].[Cl-].[Cl-]', 'Aluminium chloride', 'Lewis acid'),
    ('[Fe+3].[Cl-].[Cl-].[Cl-]', 'Iron(III) chloride', 'Lewis acid / oxidant'),
    ('[Fe+2].[O-]S([O-])(=O)=O', 'Iron(II) sulfate', 'Transition metal salt'),
    ('[Cu+2].[O-]S([O-])(=O)=O', 'Copper(II) sulfate', 'Transition metal salt'),
    ('[Na+].[Na+].[O-]S([O-])(=O)=O', 'Sodium sulfate', 'Salt'),
    ('[Mg+2].[O-]S([O-])(=O)=O', 'Magnesium sulfate', 'Salt'),
    ('[Ca+2].[O-]S([O-])(=O)=O', 'Calcium sulfate', 'Salt'),
    ('[Na+].[Na+].OP([O-])([O-])=O', 'Disodium hydrogen phosphate', 'Phosphate salt'),
    ('[Ca+2].OP([O-])([O-])=O', 'Calcium hydrogen phosphate', 'Phosphate salt'),
    ('[NH4+].[Cl-]', 'Ammonium chloride', 'Salt'),
    ('[K+].[I-]', 'Potassium iodide', 'Halide salt'),
    ('[Na+].[Br-]', 'Sodium bromide', 'Halide salt'),
    ('[Na+].[F-]', 'Sodium fluoride', 'Halide salt'),
    ('[Na+].[N-]=[N+]=[N-]', 'Sodium azide', 'Azide'),
    ('OO', 'Hydrogen peroxide', 'Oxidant'),
    ('CC(C)(C)OO', 'tert-Butyl hydroperoxide', 'Oxidant'),
    ('CC(C)(OO)c1ccccc1', 'Cumene hydroperoxide', 'Oxidant'),
    ('CC(=O)OO', 'Peracetic acid', 'Oxidant'),
    ('OOC(=O)c1cccc(Cl)c1', 'm-Chloroperbenzoic acid', 'Oxidant'),
    ('[Na+].[O-]Cl', 'Sodium hypochlorite', 'Oxidant'),
    ('[K+].[O-][Mn](=O)(=O)=O', 'Potassium permanganate', 'Oxidant'),
    ('O=O', 'Molecular oxygen', 'Oxidant'),
    ('[O-][O+]=O', 'Ozone', 'Oxidant'),
    ('[Na+].[O-]N=O', 'Sodium nitrite', 'Nitrosating agent'),
    ('[Na+].[O-][N+]([O-])=O', 'Sodium nitrate', 'Nitrate salt'),
    ('[Na+].[BH4-]', 'Sodium borohydride', 'Hydride reductant'),
    ('[Na+].OS([O-])=O', 'Sodium bisulfite', 'Reductant'),
    ('[Na+].[Na+].[O-]S([O-])=O', 'Sodium sulfite', 'Reductant'),
    ('[Na+].[Na+].[O-]S(=O)OS([O-])(=O)=O', 'Sodium metabisulfite', 'Reductant'),
]

_NAME_BY_CANON: Dict[str, str] = {}


def _build_names() -> None:
    for smi, name, _cat in LIB_ROWS:
        m = parse_smiles(smi)
        if m is None:
            continue
        c = canon(m)
        prev = _NAME_BY_CANON.get(c)
        if prev is None:
            _NAME_BY_CANON[c] = name
        elif name not in prev.split(" / "):
            _NAME_BY_CANON[c] = prev + " / " + name


_build_names()


def library_name(m: Chem.Mol) -> Optional[str]:
    n = _NAME_BY_CANON.get(canon(m))
    if not n:
        return None
    return n.split(" / ")[0] + " (or a stereoisomer)" if " / " in n else n


# --------------------------------------------------------------------------------------
# Reactivity profiles (one per functional-group family)
# --------------------------------------------------------------------------------------
PROFILES = {
    'acid': dict(name='Carboxylic Acid', cat='Carboxylic Acid', frag='-C(=O)OH', site='Carboxyl proton and carbonyl carbon',
        acid=('Low', 'Stays un-ionised; only Fischer esterification with alcohols under strong acid / heat.'),
        base=('High', 'Immediate deprotonation to the carboxylate salt (pKa ~4-5); no covalent breakdown.'),
        hyd=('Resistant', 'Hydrolytically inert terminus.'),
        photo=('Low', 'Weak UV absorber unless conjugated; alpha-aryl acids can photodecarboxylate.'),
        therm=('Low', 'Simple acids resist decarboxylation; beta-keto, malonic, alpha-EWG and electron-rich aryl acids lose CO2 on heating.'),
        ox=('Low', 'Stable to O2; radical decarboxylation only with strong oxidants / photocatalysis.'),
        sec=('High', 'Bases, amines, alcohols, metal oxides', 'Acid-base salt formation with bases/amines; esterification with alcohols and amidation with amines on heating.'),
        partners=['base', 'base_strong', 'base_weak', 'N_nuc', 'O_nuc', 'metal', 'alkylating'], provides=['acid_weak']),
    'carboxylate': dict(name='Carboxylate Salt', cat='Carboxylic Acid', frag='-C(=O)O-', site='Carboxylate anion / ion pair',
        acid=('High', 'Protonated to the free acid by stronger acids (may precipitate).'),
        base=('Resistant', 'Stable as the anion.'),
        hyd=('Resistant', 'Hydrolytically inert.'),
        photo=('Low', 'Photostable unless conjugated.'),
        therm=('Low', 'Stable to moderate heat; decarboxylates only at high temperature.'),
        ox=('Low', 'Resistant to auto-oxidation.'),
        sec=('Moderate', 'Multivalent metal ions, acids', 'Forms insoluble metal carboxylates (soaps) with Ca2+/Mg2+/Zn2+/Al3+; protonated by acids.'),
        partners=['acid_weak', 'acid_strong', 'metal', 'alkylating'], provides=['base']),
    'ester': dict(name='Carboxylic Ester', cat='Carbonyl', frag='-C(=O)O-', site='Ester carbonyl carbon and acyl-oxygen bond',
        acid=('High', 'Acid-catalysed acyl-oxygen cleavage (A_AC2) via protonated carbonyl to acid + alcohol; reversible.'),
        base=('Critical', 'Saponification (B_AC2): hydroxide attacks the carbonyl; irreversible via carboxylate formation.'),
        hyd=('Low', 'Slow neutral hydrolysis, accelerated by humidity, heat and autocatalysis by liberated acid.'),
        photo=('Low', 'Weak n→pi* absorber; simple alkyl esters are photostable.'),
        therm=('Moderate', 'Thermal transesterification; beta-H syn-elimination (ester pyrolysis) only at high temperature.'),
        ox=('Low', 'Alkoxy alpha C-H oxidation only under radical conditions.'),
        sec=('High', 'Amines, alcohols, water', 'Aminolysis to amides and transesterification with nucleophilic co-components.'),
        partners=['N_nuc', 'O_nuc', 'base', 'base_strong', 'water', 'acid_strong'], provides=['ester']),
    'ester_aryl': dict(name='Aryl / Vinyl Ester (Activated)', cat='Carbonyl', frag='Ar-O-C(=O)-', site='Ester carbonyl carbon (good phenoxide/enolate leaving group)',
        acid=('High', 'Acid-catalysed hydrolysis to acid + phenol/enol.'),
        base=('Critical', 'Rapid saponification because the phenoxide/enolate is an excellent leaving group.'),
        hyd=('High', 'Activated ester: measurable neutral hydrolysis and buffer catalysis under humidity.'),
        photo=('High', 'Photo-Fries rearrangement (aryl esters) or acyl-O homolysis on UV excitation.'),
        therm=('Moderate', 'Thermal acyl transfer to nucleophiles; Fries-type rearrangement with Lewis acids.'),
        ox=('Low', 'Resistant to ambient oxidation.'),
        sec=('Critical', 'Amines, alcohols, bases', 'Facile acyl transfer (aminolysis / transesterification) and base-catalysed cleavage.'),
        partners=['N_nuc', 'O_nuc', 'S_nuc', 'base', 'base_strong', 'water'], provides=['ester', 'acylating']),
    'lactone': dict(name='Cyclic Lactone', cat='Heterocycle / Ester', frag='C1OC(=O)CC1', site='Cyclic ester carbonyl',
        acid=('Moderate', 'Acid-catalysed ring opening to hydroxy-acid; equilibrium often favours the closed 5-/6-ring lactone.'),
        base=('Critical', 'Alkaline ring opening to the hydroxy-carboxylate salt.'),
        hyd=('Moderate', 'Ring-chain equilibrium in water; strained (4-ring) lactones hydrolyse rapidly.'),
        photo=('Moderate', 'Photodecarboxylation / photorearrangement possible for conjugated lactones.'),
        therm=('Moderate', 'Thermal ring opening/polymerisation (ROP) and dehydration.'),
        ox=('Low', 'Resistant to direct oxidation.'),
        sec=('High', 'Amines, alcohols, bases', 'Aminolysis/alcoholysis opens the ring to hydroxy-amide / hydroxy-ester.'),
        partners=['N_nuc', 'O_nuc', 'base', 'base_strong', 'water'], provides=['ester', 'acylating']),
    'amide': dict(name='Amide Bond', cat='Carbonyl / Nitrogen', frag='-C(=O)N<', site='Amide carbonyl carbon and C-N bond',
        acid=('Moderate', 'Acid-catalysed hydrolysis (A_AC2) to carboxylic acid + ammonium; needs strong acid / heat.'),
        base=('Moderate', 'Base-promoted acyl substitution; slowed by amide resonance.'),
        hyd=('Low', 'Very slow neutral hydrolysis (half-lives of years); faster at extreme pH.'),
        photo=('Moderate', 'UV-induced C-N cleavage; anilides can undergo photo-Fries rearrangement.'),
        therm=('Low', 'Thermally robust; deamidation or cyclisation only at high temperature.'),
        ox=('Low', 'Radical alpha-C-H abstraction at N-alkyl groups (N-dealkylation) under harsh conditions.'),
        sec=('Low', 'Polar co-components', 'Hydrogen-bond donor/acceptor interactions; no covalent reactivity under mild conditions.'),
        partners=['acid_strong', 'base_strong', 'water'], provides=[]),
    'lactam': dict(name='Lactam (Cyclic Amide)', cat='Heterocycle / Amide', frag='C1NC(=O)CC1', site='Cyclic amide carbonyl',
        acid=('Moderate', 'Ring opening to amino acid under strong acid / heat (6- and 7-ring faster than 5-ring).'),
        base=('Moderate', 'Base-catalysed ring opening to amino-carboxylate; slow for 5-ring.'),
        hyd=('Low', 'Stable under ambient humidity.'),
        photo=('Low', 'Photostable unless conjugated.'),
        therm=('Moderate', 'Ring-opening polymerisation of medium-ring lactams at elevated temperature.'),
        ox=('Low', 'Resistant; N-alkyl alpha C-H oxidation under radical stress.'),
        sec=('Low', 'Strong nucleophiles', 'Ring opening only with strong nucleophiles or catalysis.'),
        partners=['base_strong', 'acid_strong'], provides=[]),
    'beta_lactam': dict(name='Beta-Lactam Core', cat='Heterocycle / Strained Amide', frag='N1C(=O)CC1', site='Strained four-membered lactam carbonyl',
        acid=('Critical', 'Protonation then rapid nucleophilic ring opening (ring strain ~26 kcal/mol).'),
        base=('Critical', 'Direct hydroxide attack on the strained carbonyl; irreversible ring cleavage.'),
        hyd=('Critical', 'Spontaneous neutral solvolysis driven by ring strain.'),
        photo=('Moderate', 'Photolytic fragmentation of the four-membered ring.'),
        therm=('High', 'Thermally accelerated ring rupture.'),
        ox=('Moderate', 'Oxidation at adjacent heteroatoms / activated C-H.'),
        sec=('Critical', 'Amines, alcohols, thiols', 'Rapid aminolysis / alcoholysis opening the strained lactam.'),
        partners=['N_nuc', 'O_nuc', 'S_nuc', 'base', 'base_strong', 'water', 'acid_strong'], provides=['acylating']),
    'imide': dict(name='Imide', cat='Carbonyl / Nitrogen', frag='-C(=O)NC(=O)-', site='Imide carbonyls (electrophilic, acidic N-H)',
        acid=('Moderate', 'Acid-catalysed ring / C-N cleavage to amic acid, then diacid.'),
        base=('Critical', 'Hydroxide opens the imide to the amic acid within minutes-hours; acidic N-H forms salts.'),
        hyd=('Moderate', 'Neutral hydrolysis to amic acid on prolonged exposure to moisture.'),
        photo=('Moderate', 'Norrish-type chemistry for aryl imides.'),
        therm=('Low', 'Thermally stable.'),
        ox=('Low', 'Resistant.'),
        sec=('High', 'Amines, bases', 'Aminolysis opens the imide ring; N-H deprotonation by bases.'),
        partners=['N_nuc', 'base', 'base_strong', 'water'], provides=['acylating']),
    'urea': dict(name='Urea', cat='Carbonyl / Nitrogen', frag='-NC(=O)N-', site='Urea carbonyl carbon',
        acid=('Moderate', 'Slow acid hydrolysis to amine + CO2 (via carbamic acid).'),
        base=('Low', 'Slow base hydrolysis; N-H urea can eliminate to isocyanate at high pH/heat.'),
        hyd=('Low', 'Hydrolytically stable at ambient conditions.'),
        photo=('Low', 'Photostable unless aryl-conjugated.'),
        therm=('Moderate', 'Thermal dissociation to isocyanate + amine above ~130 C.'),
        ox=('Low', 'Resistant.'),
        sec=('Low', 'Aldehydes', 'Condensation with formaldehyde / aldehydes (methylol ureas).'),
        partners=['carbonyl'], provides=[]),
    'carbamate': dict(name='Carbamate', cat='Carbonyl / Nitrogen', frag='-NC(=O)O-', site='Carbamate carbonyl carbon',
        acid=('Moderate', 'Acid hydrolysis to amine + CO2 + alcohol; tert-alkyl carbamates cleave rapidly (A_AL1).'),
        base=('Moderate', 'Base hydrolysis, faster for aryl carbamates via E1cB isocyanate.'),
        hyd=('Low', 'Stable at neutral pH.'),
        photo=('Moderate', 'Photocleavage of aryl carbamates.'),
        therm=('Moderate', 'Thermolysis to isocyanate + alcohol (above 150 C); tert-butyl types lose isobutene.'),
        ox=('Low', 'Resistant.'),
        sec=('Moderate', 'Amines, alcohols', 'Transcarbamoylation with strong nucleophiles.'),
        partners=['N_nuc', 'O_nuc', 'base_strong'], provides=[]),
    'carbonate': dict(name='Carbonate Ester', cat='Carbonyl', frag='-OC(=O)O-', site='Carbonate carbonyl carbon',
        acid=('Moderate', 'Acid hydrolysis to alcohols + CO2.'),
        base=('Critical', 'Rapid base hydrolysis to alcohols + carbonate.'),
        hyd=('Moderate', 'Slow neutral hydrolysis.'),
        photo=('Low', 'Photostable (aliphatic).'),
        therm=('Moderate', 'Thermal decarboxylation / transesterification.'),
        ox=('Low', 'Resistant.'),
        sec=('High', 'Amines, alcohols', 'Aminolysis to carbamates; transesterification.'),
        partners=['N_nuc', 'O_nuc', 'base', 'base_strong', 'water'], provides=['acylating']),
    'anhydride': dict(name='Acid Anhydride', cat='Carbonyl', frag='-C(=O)OC(=O)-', site='Anhydride carbonyl carbons',
        acid=('Critical', 'Acid-catalysed hydrolysis to two carboxylic acids.'),
        base=('Critical', 'Instant saponification.'),
        hyd=('Critical', 'Spontaneous hydrolysis by ambient moisture.'),
        photo=('Low', 'Photostable (aliphatic).'),
        therm=('Moderate', 'Thermal disproportionation / decarboxylative fragmentation.'),
        ox=('Low', 'Resistant.'),
        sec=('Critical', 'Amines, alcohols, water', 'Vigorous acylation of amines (amides) and alcohols (esters).'),
        partners=['N_nuc', 'O_nuc', 'S_nuc', 'water', 'base'], provides=['acylating']),
    'acyl_halide': dict(name='Acyl Halide', cat='Carbonyl', frag='-C(=O)X', site='Acyl carbon (strong electrophile)',
        acid=('Critical', 'Immediate hydrolysis to the acid + HX.'),
        base=('Critical', 'Instant saponification.'),
        hyd=('Critical', 'Reacts violently with atmospheric moisture.'),
        photo=('Moderate', 'Photolytic C-X cleavage.'),
        therm=('High', 'Thermally labile; ketene formation with alpha-H.'),
        ox=('Low', 'Resistant.'),
        sec=('Critical', 'Amines, alcohols, water', 'Fast acylation of any nucleophile.'),
        partners=['N_nuc', 'O_nuc', 'S_nuc', 'water', 'base'], provides=['acylating', 'acid_strong']),
    'thioester': dict(name='Thioester', cat='Carbonyl / Sulfur', frag='-C(=O)S-', site='Thioester carbonyl carbon and sulfur',
        acid=('Moderate', 'Acid hydrolysis to acid + thiol (slower than oxoesters in acid).'),
        base=('High', 'Base hydrolysis / thiolate exchange.'),
        hyd=('Moderate', 'Slow neutral hydrolysis.'),
        photo=('Moderate', 'Photo-cleavage of the C-S bond.'),
        therm=('Low', 'Thermally stable.'),
        ox=('High', 'Sulfur oxidation to sulfoxide-type species.'),
        sec=('High', 'Amines, thiols', 'Aminolysis and thiol-thioester exchange.'),
        partners=['N_nuc', 'S_nuc', 'base'], provides=['acylating']),
    'aldehyde': dict(name='Aldehyde', cat='Carbonyl', frag='-CHO', site='Formyl carbon (oxidation / condensation centre)',
        acid=('High', 'Acid-catalysed hydration and acetal formation with alcohols.'),
        base=('High', 'Aldol self-condensation (alpha-H present) or Cannizzaro disproportionation (no alpha-H).'),
        hyd=('Moderate', 'Reversible hydration to the gem-diol in water.'),
        photo=('High', 'Norrish type I/II cleavage and photo-decarbonylation.'),
        therm=('Moderate', 'Thermal oligomerisation / decarbonylation.'),
        ox=('Critical', 'Radical auto-oxidation to the carboxylic acid via peracid.'),
        sec=('Critical', 'Amines, alcohols', 'Imine (Schiff base) formation with amines; hemiacetal/acetal with alcohols.'),
        partners=['N_nuc', 'N_nuc_prim', 'O_nuc', 'base', 'base_strong', 'oxidant'], provides=['carbonyl']),
    'ketone': dict(name='Ketone', cat='Carbonyl', frag='-C(=O)-', site='Carbonyl carbon and enolisable alpha-carbons',
        acid=('Moderate', 'Acid-catalysed enolisation, aldol condensation and ketal formation.'),
        base=('Moderate', 'Enolate formation; aldol condensation; alpha-epimerisation.'),
        hyd=('Resistant', 'Resistant to neutral hydrolysis.'),
        photo=('High', 'n→pi* absorption (280-320 nm) initiates Norrish type I/II photocleavage.'),
        therm=('Low', 'Thermally stable.'),
        ox=('Moderate', 'Baeyer-Villiger oxidation by peroxides/peracids to esters.'),
        sec=('Moderate', 'Primary amines, alcohols', 'Ketimine/enamine formation; ketalisation.'),
        partners=['N_nuc_prim', 'N_nuc', 'O_nuc', 'oxidant'], provides=['carbonyl']),
    'michael': dict(name='Alpha,Beta-Unsaturated Carbonyl (Michael Acceptor)', cat='Carbonyl / Unsaturation', frag='C=C-C=O', site='Electrophilic beta-carbon and carbonyl',
        acid=('Moderate', 'Acid-catalysed hydration / conjugate addition of water.'),
        base=('High', 'Hydroxide conjugate addition; base-catalysed polymerisation.'),
        hyd=('Low', 'Slow conjugate hydration in water.'),
        photo=('High', '[2+2] photocycloaddition and E/Z photo-isomerisation.'),
        therm=('Moderate', 'Thermal radical polymerisation / Diels-Alder dimerisation.'),
        ox=('High', 'Nucleophilic epoxidation (Weitz-Scheffer) with peroxides.'),
        sec=('Critical', 'Amines, thiols', 'Aza-/thia-Michael addition to the beta-carbon.'),
        partners=['N_nuc', 'S_nuc', 'base', 'oxidant'], provides=['michael']),
    'quinone': dict(name='Quinone / Quinone Imine', cat='Carbonyl / Redox-active', frag='O=C1C=CC(=O)C=C1', site='Electrophilic ring carbons and redox-active carbonyls',
        acid=('Low', 'Stable in mild acid; acid-catalysed conjugate addition of water.'),
        base=('High', 'Hydroxide addition / base-promoted decomposition and polymerisation.'),
        hyd=('Low', 'Slow conjugate hydration.'),
        photo=('High', 'Photoreduction and photo-addition to alkenes.'),
        therm=('Moderate', 'Thermal dimerisation / polymerisation.'),
        ox=('Low', 'Already in the oxidised state.'),
        sec=('Critical', 'Thiols, amines, reductants', 'Michael addition of thiols/amines; reduced to hydroquinone by reductants.'),
        partners=['S_nuc', 'N_nuc', 'reductant'], provides=['oxidant', 'michael']),
    'thiocarbonyl': dict(name='Thiocarbonyl (C=S)', cat='Sulfur / Carbonyl', frag='C=S', site='Thiocarbonyl carbon and sulfur',
        acid=('Moderate', 'Acid hydrolysis to the carbonyl analogue with H2S loss.'),
        base=('Moderate', 'Base hydrolysis to the carbonyl analogue.'),
        hyd=('Low', 'Slow hydrolysis.'),
        photo=('High', 'Photo-oxidation / desulfurisation.'),
        therm=('Moderate', 'Thermal rearrangement (thione-thiol).'),
        ox=('Critical', 'S-oxidation and oxidative desulfurisation to C=O.'),
        sec=('Moderate', 'Oxidants, metals', 'Oxidative desulfurisation; strong metal coordination via sulfur.'),
        partners=['oxidant', 'metal'], provides=['reductant']),
    'amine': dict(name='Aliphatic Amine', cat='Amine', frag='-NH2 / -NHR', site='Basic, nucleophilic nitrogen lone pair',
        acid=('High', 'Protonation to the ammonium salt (reversible); no covalent breakdown.'),
        base=('Low', 'Stays as the free base (nucleophilic form).'),
        hyd=('Resistant', 'Hydrolytically inert.'),
        photo=('Moderate', 'Photosensitised alpha-C-H abstraction / deamination.'),
        therm=('Moderate', 'Thermal deamination, condensation or oxidative discolouration.'),
        ox=('High', 'Air / peroxide oxidation to hydroxylamine, nitrone/imine or N-dealkylation.'),
        sec=('Critical', 'Carbonyls, reducing sugars, esters, acylating agents, nitrite', 'Schiff-base / Maillard chemistry with carbonyls; acylation by esters/anhydrides; N-alkylation; N-nitrosation of secondary amines.'),
        partners=['carbonyl', 'reducing_sugar', 'acylating', 'ester', 'alkylating', 'michael', 'nitrosating', 'acid_weak', 'acid_strong', 'oxidant'], provides=['N_nuc', 'N_nuc_prim', 'base']),
    'amine_tert': dict(name='Tertiary Amine', cat='Amine', frag='-NR3', site='Tertiary nitrogen lone pair and N-alkyl alpha C-H',
        acid=('High', 'Protonation to the tertiary ammonium salt.'),
        base=('Resistant', 'Not ionisable by bases.'),
        hyd=('Resistant', 'Hydrolytically stable.'),
        photo=('Low', 'Weak direct UV absorption (unless aryl-conjugated).'),
        therm=('Moderate', 'Hofmann-type / Cope elimination pathways at elevated temperature.'),
        ox=('Critical', 'N-oxidation to the N-oxide by peroxides / O2; oxidative N-dealkylation.'),
        sec=('Moderate', 'Peroxides, alkyl halides, nitrite', 'Peroxide N-oxidation; quaternisation by alkylating agents; nitrosative dealkylation to nitrosamines.'),
        partners=['oxidant', 'peroxide_former', 'alkylating', 'nitrosating', 'acid_weak', 'acid_strong'], provides=['base']),
    'arylamine': dict(name='Aromatic Amine (Aniline)', cat='Amine', frag='Ar-NH2 / Ar-NHR', site='Aniline nitrogen and activated ortho/para ring carbons',
        acid=('Moderate', 'Protonation to anilinium (pKa ~4-5); diazotisation with nitrite/HONO.'),
        base=('Low', 'Stable as the free base.'),
        hyd=('Resistant', 'Hydrolytically stable.'),
        photo=('Critical', 'Photo-oxidation giving coloured azo / quinone-imine products.'),
        therm=('Moderate', 'Thermal oxidative coupling / dimerisation.'),
        ox=('Critical', 'Multi-electron oxidation to hydroxylamine, nitroso, nitro, azo and quinone-imine polymers.'),
        sec=('High', 'Aldehydes, reducing sugars, nitrite, acylating agents', 'Schiff-base condensation with carbonyls; diazotisation / N-nitrosation with nitrite; acylation.'),
        partners=['carbonyl', 'reducing_sugar', 'acylating', 'nitrosating', 'oxidant', 'alkylating'], provides=['N_nuc', 'N_nuc_prim', 'base']),
    'quat': dict(name='Quaternary Ammonium', cat='Amine', frag='-N+R4', site='Quaternary nitrogen and beta C-H (Hofmann)',
        acid=('Resistant', 'Stable in acid.'),
        base=('Moderate', 'Hofmann elimination / hydroxide displacement at elevated temperature.'),
        hyd=('Resistant', 'Hydrolytically stable.'),
        photo=('Low', 'Photostable.'),
        therm=('Moderate', 'Hofmann elimination and dealkylation on heating.'),
        ox=('Low', 'Resistant to oxidation.'),
        sec=('Low', 'Anionic species', 'Ion-pairing with anionic compounds (surfactant-type complexation).'),
        partners=['acid_weak'], provides=[]),
    'ammonium': dict(name='Ammonium / Protonated Amine', cat='Amine salt', frag='-NH3+', site='Acidic N-H+',
        acid=('Resistant', 'Stable in acid.'),
        base=('Critical', 'Deprotonated to the free (nucleophilic) amine by bases.'),
        hyd=('Resistant', 'Hydrolytically inert.'),
        photo=('Low', 'Photostable.'),
        therm=('Moderate', 'Loss of HX / amine on heating (dissociation of the salt).'),
        ox=('Low', 'Resistant until free-based.'),
        sec=('Moderate', 'Bases', 'Liberation of the free amine, which then reacts as an amine.'),
        partners=['base', 'base_strong'], provides=['acid_weak']),
    'amidine': dict(name='Amidine / Guanidine', cat='Nitrogen', frag='-C(=N)N<', site='Amidinium carbon and basic nitrogen',
        acid=('High', 'Protonation to a highly stabilised amidinium/guanidinium; slow hydrolysis to amide/urea.'),
        base=('Moderate', 'Base hydrolysis to amide / urea.'),
        hyd=('Moderate', 'Slow hydrolysis to amide / urea.'),
        photo=('Low', 'Photostable.'),
        therm=('Low', 'Thermally stable.'),
        ox=('Low', 'Resistant.'),
        sec=('Moderate', 'Acids, carbonyls', 'Strong base; salt formation; condensation with 1,3-dicarbonyls.'),
        partners=['acid_weak', 'acid_strong', 'carbonyl'], provides=['base']),
    'imine': dict(name='Imine (Schiff Base)', cat='Nitrogen / C=N', frag='C=N', site='Electrophilic imine carbon',
        acid=('Critical', 'Hydrolysis to carbonyl + amine (protonated iminium).'),
        base=('Moderate', 'Base-catalysed hydrolysis / tautomerisation to enamine.'),
        hyd=('High', 'Reversible hydrolysis in water.'),
        photo=('Moderate', 'E/Z photo-isomerisation and photo-cleavage.'),
        therm=('Moderate', 'Thermal tautomerisation / oligomerisation.'),
        ox=('Moderate', 'Oxidation to amide or oxaziridine.'),
        sec=('High', 'Amines, carbonyls, reductants', 'Transimination with amines; reduction to amines.'),
        partners=['N_nuc_prim', 'N_nuc', 'reductant', 'water'], provides=[]),
    'oxime': dict(name='Oxime / Nitrone', cat='Nitrogen / C=N-O', frag='C=N-OH', site='Oxime carbon, N-O bond',
        acid=('Moderate', 'Hydrolysis to carbonyl + hydroxylamine; Beckmann rearrangement with strong acid.'),
        base=('Low', 'Stable; oximate salt formation.'),
        hyd=('Low', 'Slow hydrolysis.'),
        photo=('Moderate', 'Photo-isomerisation (E/Z) and N-O cleavage.'),
        therm=('Moderate', 'Thermal Beckmann / dehydration to nitrile.'),
        ox=('Moderate', 'Oxidative cleavage to nitro / carbonyl.'),
        sec=('Low', 'Acids', 'Acid-promoted rearrangement.'),
        partners=['acid_strong'], provides=[]),
    'hydrazone': dict(name='Hydrazone', cat='Nitrogen / C=N-N', frag='C=N-N', site='Hydrazone carbon, N-N unit',
        acid=('High', 'Hydrolysis to carbonyl + hydrazine.'),
        base=('Low', 'Stable; Wolff-Kishner-type reduction requires strong base/heat.'),
        hyd=('Moderate', 'Slow hydrolysis.'),
        photo=('Moderate', 'Photo-isomerisation / N-N cleavage.'),
        therm=('Moderate', 'Thermal decomposition.'),
        ox=('High', 'Oxidation to azo/diazo compounds.'),
        sec=('Moderate', 'Carbonyls', 'Exchange with carbonyl compounds.'),
        partners=['carbonyl'], provides=[]),
    'hydrazine': dict(name='Hydrazine (N-N)', cat='Nitrogen', frag='N-N', site='Nucleophilic N-N nitrogens',
        acid=('Moderate', 'Protonation to hydrazinium salts.'),
        base=('Low', 'Stable.'),
        hyd=('Resistant', 'Hydrolytically stable.'),
        photo=('Moderate', 'Photo-oxidation / N-N cleavage.'),
        therm=('Moderate', 'Thermal decomposition to N2/NH3.'),
        ox=('Critical', 'Oxidation to diazene / azo compounds with N2 loss.'),
        sec=('Critical', 'Carbonyls, oxidants', 'Condensation with carbonyls to hydrazones; reduces oxidants and metal ions.'),
        partners=['carbonyl', 'reducing_sugar', 'oxidant', 'acylating'], provides=['N_nuc', 'N_nuc_prim', 'reductant']),
    'hydroxylamine': dict(name='Hydroxylamine (N-OH)', cat='Nitrogen / Oxygen', frag='N-OH', site='Nucleophilic N and N-O bond',
        acid=('Moderate', 'Protonation to hydroxylammonium salts.'),
        base=('Moderate', 'Free-base is unstable; disproportionation.'),
        hyd=('Resistant', 'Hydrolytically stable.'),
        photo=('Moderate', 'Photolytic N-O cleavage.'),
        therm=('High', 'Thermal disproportionation (can be energetic).'),
        ox=('Critical', 'Oxidation to nitroso / nitrone / nitroxide.'),
        sec=('Critical', 'Carbonyls, oxidants', 'Oxime formation with carbonyls; reducing agent for oxidants / metal ions.'),
        partners=['carbonyl', 'reducing_sugar', 'oxidant', 'acylating'], provides=['N_nuc', 'reductant']),
    'nitrile': dict(name='Nitrile', cat='Nitrogen', frag='-C#N', site='Electrophilic nitrile carbon',
        acid=('Moderate', 'Hydration to amide then carboxylic acid (needs strong acid / heat).'),
        base=('Moderate', 'Base hydrolysis to amide, then carboxylate + NH3.'),
        hyd=('Resistant', 'Stable at neutral pH.'),
        photo=('Low', 'Photostable (aliphatic).'),
        therm=('Low', 'Thermally stable.'),
        ox=('Low', 'Resistant.'),
        sec=('Low', 'Alcohols + acid', 'Pinner reaction with alcohols under HCl.'),
        partners=['O_nuc', 'acid_strong'], provides=[]),
    'nitro': dict(name='Nitro Group', cat='Nitrogen / Oxygen', frag='-NO2', site='Nitro nitrogen (strong electron-withdrawing group)',
        acid=('Resistant', 'Stable in acid.'),
        base=('Low', 'Stable; aliphatic nitro compounds form nitronate salts and aromatic nitro activates SNAr / Meisenheimer chemistry.'),
        hyd=('Resistant', 'Hydrolytically stable.'),
        photo=('High', 'Photoreduction, nitro→nitrite rearrangement and ortho-nitrobenzyl photochemistry.'),
        therm=('Moderate', 'Exothermic thermal decomposition at high temperature.'),
        ox=('Resistant', 'Already fully oxidised.'),
        sec=('Moderate', 'Reductants, nucleophiles', 'Reduced to nitroso / hydroxylamine / amine; activates SNAr displacement of ortho/para leaving groups.'),
        partners=['reductant', 'N_nuc', 'base_strong'], provides=['oxidant']),
    'nitrosamine': dict(name='N-Nitroso (Nitrosamine)', cat='Nitrogen / Nitroso', frag='N-N=O', site='N-nitroso group (N-N=O)',
        acid=('High', 'Acid-mediated denitrosation / transnitrosation.'),
        base=('Low', 'Stable in base (alpha-deprotonation only with strong base).'),
        hyd=('Low', 'Hydrolytically stable.'),
        photo=('Critical', 'UV (230-350 nm) cleaves the N-N bond giving aminyl radical + NO.'),
        therm=('Moderate', 'Thermally stable at moderate temperature.'),
        ox=('Moderate', 'Alpha-hydroxylation (P450-like / radical) leading to alkyl-diazonium ions.'),
        sec=('Moderate', 'Thiols, acids', 'Transnitrosation to thiols/amines; denitrosation with strong acid.'),
        partners=['S_nuc', 'N_nuc', 'acid_strong'], provides=['nitrosating']),
    'nitroso': dict(name='C-Nitroso', cat='Nitrogen / Nitroso', frag='C-N=O', site='Nitroso nitrogen and oxygen',
        acid=('Moderate', 'Acid-catalysed rearrangement / condensation.'),
        base=('Moderate', 'Base-promoted condensation.'),
        hyd=('Low', 'Hydrolytically stable.'),
        photo=('High', 'Photo-dimerisation / photolysis to radicals.'),
        therm=('High', 'Thermal dimerisation-monomerisation equilibrium and decomposition.'),
        ox=('High', 'Oxidation to nitro.'),
        sec=('Critical', 'Amines, thiols', 'Condenses with amines (azo) and thiols; ene reactions.'),
        partners=['N_nuc', 'S_nuc', 'reductant'], provides=[]),
    'azo': dict(name='Azo (N=N)', cat='Nitrogen', frag='-N=N-', site='Azo nitrogens (chromophore)',
        acid=('Low', 'Stable (azo-hydrazone tautomerism for hydroxy-azo dyes).'),
        base=('Low', 'Stable.'),
        hyd=('Resistant', 'Hydrolytically stable.'),
        photo=('Critical', 'E/Z photo-isomerisation and photo-cleavage with N2 loss.'),
        therm=('Moderate', 'Thermal cis-trans relaxation; aliphatic azo compounds lose N2 (radical source).'),
        ox=('Moderate', 'Oxidative cleavage by strong oxidants.'),
        sec=('Moderate', 'Reductants', 'Reduced to hydrazo compounds and then amines.'),
        partners=['reductant'], provides=[]),
    'azide': dict(name='Azide', cat='Nitrogen', frag='-N3', site='Azide nitrogens (energetic)',
        acid=('Moderate', 'Protonation gives volatile toxic hydrazoic acid (HN3).'),
        base=('Low', 'Stable.'),
        hyd=('Resistant', 'Hydrolytically stable.'),
        photo=('Critical', 'Photolysis extrudes N2 to give a nitrene.'),
        therm=('High', 'Thermal N2 loss (Curtius-like); explosive hazard.'),
        ox=('Low', 'Resistant.'),
        sec=('High', 'Phosphines, alkynes, reductants', 'Staudinger reduction, azide-alkyne cycloaddition, hydrogenation to amines.'),
        partners=['reductant', 'N_nuc'], provides=[]),
    'isocyanate': dict(name='Isocyanate', cat='Nitrogen / Carbonyl', frag='-N=C=O', site='Cumulated N=C=O carbon',
        acid=('Critical', 'Hydrolysis to amine + CO2.'),
        base=('Critical', 'Base-catalysed hydrolysis / trimerisation.'),
        hyd=('Critical', 'Reacts with moisture to give amine + CO2 (then urea).'),
        photo=('Moderate', 'Photolysis to nitrene / CO.'),
        therm=('Moderate', 'Dimerisation / trimerisation on heating.'),
        ox=('Low', 'Resistant.'),
        sec=('Critical', 'Amines, alcohols', 'Addition to amines gives ureas; to alcohols gives carbamates.'),
        partners=['N_nuc', 'O_nuc', 'S_nuc', 'water'], provides=['acylating']),
    'isothiocyanate': dict(name='Isothiocyanate', cat='Nitrogen / Sulfur', frag='-N=C=S', site='Cumulated N=C=S carbon',
        acid=('Moderate', 'Hydrolysis to amine + COS.'),
        base=('High', 'Base-promoted hydrolysis.'),
        hyd=('Moderate', 'Slow hydrolysis in water.'),
        photo=('Moderate', 'Photolysis.'),
        therm=('Moderate', 'Thermal rearrangement.'),
        ox=('Moderate', 'Oxidative desulfurisation.'),
        sec=('Critical', 'Amines, thiols', 'Addition of amines gives thioureas; thiols give dithiocarbamates.'),
        partners=['N_nuc', 'S_nuc'], provides=[]),
    'n_oxide': dict(name='N-Oxide', cat='Nitrogen / Oxygen', frag='N+-O-', site='N-oxide oxygen',
        acid=('Moderate', 'Protonation to N-hydroxy ammonium.'),
        base=('Low', 'Stable.'),
        hyd=('Resistant', 'Hydrolytically stable.'),
        photo=('High', 'Photo-deoxygenation / rearrangement.'),
        therm=('High', 'Cope / Meisenheimer / Polonovski rearrangements on heating.'),
        ox=('Low', 'Already oxidised.'),
        sec=('Moderate', 'Reductants, acylating agents', 'Deoxygenation by reductants; Polonovski activation by anhydrides.'),
        partners=['reductant', 'acylating'], provides=['oxidant']),
    'alcohol': dict(name='Aliphatic Alcohol', cat='Hydroxyl', frag='-CH(OH)-', site='Carbinol carbon (dehydration / oxidation centre)',
        acid=('Moderate', 'Acid-promoted dehydration, etherification or esterification.'),
        base=('Resistant', 'Resistant to base (forms alkoxide only with strong base).'),
        hyd=('Resistant', 'Hydrolytically stable.'),
        photo=('Low', 'Transparent in the solar UV.'),
        therm=('Moderate', 'Thermal dehydration / dehydrogenation.'),
        ox=('High', 'Oxidation to aldehyde/ketone (primary/secondary carbinols).'),
        sec=('Moderate', 'Acids, esters, carbonyls', 'Esterification / transesterification; hemiacetal and acetal formation.'),
        partners=['acid_weak', 'acid_strong', 'ester', 'acylating', 'carbonyl', 'alkylating'], provides=['O_nuc']),
    'alcohol_tert': dict(name='Tertiary Alcohol', cat='Hydroxyl', frag='-C(OH)<', site='Tertiary carbinol carbon',
        acid=('High', 'E1 dehydration to alkene via tertiary carbocation.'),
        base=('Resistant', 'Resistant.'),
        hyd=('Resistant', 'Hydrolytically stable.'),
        photo=('Low', 'Transparent in the solar UV.'),
        therm=('Moderate', 'Thermal / acid-catalysed dehydration.'),
        ox=('Resistant', 'Resists oxidation without C-C cleavage.'),
        sec=('Low', 'Acylating agents', 'Sterically hindered esterification.'),
        partners=['acylating', 'acid_strong'], provides=['O_nuc']),
    'phenol': dict(name='Phenol (Ar-OH)', cat='Hydroxyl', frag='Ar-OH', site='Phenolic oxygen and activated ortho/para ring carbons',
        acid=('Resistant', 'Resistant to acid solvolysis.'),
        base=('High', 'Deprotonation to phenolate, which oxidises much faster.'),
        hyd=('Resistant', 'Hydrolytically stable.'),
        photo=('High', 'UV excitation gives phenoxyl radicals; photo-coupling to dimers.'),
        therm=('Moderate', 'Thermally accelerated oxidative coupling.'),
        ox=('Critical', 'Single-electron oxidation to phenoxyl radical, then C-C/C-O coupling or quinone formation.'),
        sec=('Moderate', 'Bases, oxidants, carbonyls, metal ions', 'Salt formation with bases; peroxide oxidation; Fe/Cu-catalysed oxidation; Mannich/aldehyde condensation.'),
        partners=['base', 'base_strong', 'oxidant', 'redox_metal', 'carbonyl', 'acylating'], provides=['O_nuc', 'acid_weak', 'reductant']),
    'enol': dict(name='Enol', cat='Hydroxyl / Unsaturation', frag='C=C-OH', site='Enolic carbon-oxygen system',
        acid=('Moderate', 'Tautomerises to the carbonyl form.'),
        base=('High', 'Enolate formation.'),
        hyd=('Low', 'Tautomerisation in water.'),
        photo=('Moderate', 'Photo-tautomerisation.'),
        therm=('Moderate', 'Thermal tautomerisation / oxidation.'),
        ox=('High', 'Air oxidation (e.g. ene-diols such as ascorbate).'),
        sec=('Moderate', 'Oxidants, metals', 'Oxidised readily; chelates metal ions.'),
        partners=['oxidant', 'redox_metal', 'metal'], provides=['reductant']),
    'ether': dict(name='Ether', cat='Ether', frag='C-O-C', site='Ether oxygen and alpha C-H bonds',
        acid=('Low', 'Cleavage only with HBr/HI or strong Lewis acids.'),
        base=('Resistant', 'Stable to base.'),
        hyd=('Resistant', 'Hydrolytically inert.'),
        photo=('Low', 'Weakly UV active.'),
        therm=('Low', 'Thermally stable.'),
        ox=('High', 'Autoxidation at alpha C-H to hydroperoxides (peroxide formation).'),
        sec=('Moderate', 'Oxidisable co-components', 'Peroxides formed on ageing oxidise co-existing S, N and phenolic groups.'),
        partners=['oxidant'], provides=['peroxide_former']),
    'aryl_ether': dict(name='Aryl Ether', cat='Ether', frag='Ar-O-C', site='Ether oxygen and O-alkyl alpha C-H',
        acid=('Low', 'Stable except with HBr/HI or BBr3.'),
        base=('Resistant', 'Stable to base.'),
        hyd=('Resistant', 'Hydrolytically inert.'),
        photo=('Moderate', 'UV-induced C-O homolysis for some aryl ethers.'),
        therm=('Low', 'Thermally stable.'),
        ox=('Moderate', 'O-dealkylation via alpha C-H abstraction (radical / enzymatic-like).'),
        sec=('Low', 'Strong Lewis acids', 'Demethylation with BBr3/HBr.'),
        partners=['acid_strong', 'lewis_acid'], provides=[]),
    'epoxide': dict(name='Epoxide (Oxirane)', cat='Ether / Strained ring', frag='C1OC1', site='Strained C-O bonds (ring-opening centre)',
        acid=('Critical', 'Acid-catalysed ring opening to the 1,2-diol (or halohydrin with HX).'),
        base=('High', 'SN2 ring opening by hydroxide to the diol.'),
        hyd=('Moderate', 'Slow neutral hydrolysis to diol.'),
        photo=('Low', 'Photostable (aliphatic).'),
        therm=('Moderate', 'Thermal rearrangement to carbonyl (Meinwald) / polymerisation.'),
        ox=('Low', 'Resistant.'),
        sec=('Critical', 'Amines, alcohols, acids, thiols', 'Ring opening by nucleophiles giving beta-substituted alcohols.'),
        partners=['N_nuc', 'O_nuc', 'S_nuc', 'acid_weak', 'acid_strong', 'base_strong', 'water'], provides=['alkylating']),
    'acetal': dict(name='Acetal / Ketal', cat='Ether / Carbonyl derivative', frag='C(OR)2', site='Acetal carbon (oxocarbenium precursor)',
        acid=('Critical', 'Rapid acid hydrolysis via oxocarbenium to carbonyl + alcohols (glycosides likewise).'),
        base=('Resistant', 'Stable to base.'),
        hyd=('Low', 'Stable at neutral pH.'),
        photo=('Low', 'Photostable.'),
        therm=('Low', 'Thermally stable in the absence of acid.'),
        ox=('Moderate', 'Oxidation of the acetal C-H to esters.'),
        sec=('Moderate', 'Alcohols + acid', 'Transacetalisation with alcohols under acid catalysis.'),
        partners=['acid_weak', 'acid_strong', 'O_nuc', 'water'], provides=[]),
    'hemiacetal': dict(name='Hemiacetal / Reducing End (Masked Aldehyde)', cat='Carbonyl derivative', frag='C(OH)(OR)', site='Anomeric carbon (ring-chain tautomerism)',
        acid=('Moderate', 'Mutarotation and glycoside formation under acid catalysis.'),
        base=('High', 'Ring opening, enolisation, epimerisation and alkaline degradation (Lobry de Bruyn).'),
        hyd=('Moderate', 'Ring-chain tautomerism (mutarotation) in water.'),
        photo=('Low', 'Photostable.'),
        therm=('High', 'Dehydration / caramelisation on heating.'),
        ox=('High', 'Oxidation of the aldehyde form to aldonic acid.'),
        sec=('Critical', 'Amines, amino acids', 'Reducing end condenses with amines (glycosylamine / Schiff base, Maillard browning).'),
        partners=['N_nuc', 'N_nuc_prim', 'oxidant', 'base'], provides=['reducing_sugar', 'carbonyl', 'reductant']),
    'peroxide': dict(name='Peroxide / Hydroperoxide', cat='Peroxide', frag='-O-O-', site='Weak O-O bond (oxidant / radical source)',
        acid=('Moderate', 'Acid-catalysed rearrangement (Hock/Criegee) or heterolysis.'),
        base=('Moderate', 'Base-promoted decomposition / perhydrolysis.'),
        hyd=('Low', 'Stable in neutral water (aqueous H2O2 disproportionates slowly).'),
        photo=('Critical', 'UV homolyses O-O to alkoxy / hydroxyl radicals.'),
        therm=('Critical', 'Thermal homolysis of O-O initiates radical chains.'),
        ox=('Resistant', 'Already an oxidant.'),
        sec=('Critical', 'Sulfides, amines, phenols, alkenes, metals', 'Oxidises S, N, phenol and C=C groups; transition metals catalyse radical (Fenton) decomposition.'),
        partners=['reductant', 'N_nuc', 'S_nuc', 'redox_metal', 'ene'], provides=['oxidant']),
    'dioxygen': dict(name='Molecular Oxygen / Ozone', cat='Oxidant', frag='O=O', site='Diradical / electrophilic oxygen',
        acid=('Resistant', 'Inert to acid.'),
        base=('Resistant', 'Inert to base.'),
        hyd=('Resistant', 'Does not hydrolyse.'),
        photo=('Moderate', 'Photosensitised singlet oxygen formation.'),
        therm=('Low', 'Stable.'),
        ox=('Resistant', 'Oxidant, not oxidised.'),
        sec=('High', 'Autoxidisable groups', 'Radical autoxidation of aldehydes, thiols, ethers, enols, phenols, amines, benzylic/allylic C-H.'),
        partners=['ene', 'S_nuc', 'N_nuc', 'reductant', 'peroxide_former'], provides=['oxidant']),
    'strong_base': dict(name='Hydroxide / Alkoxide / Metal Oxide (Strong Base)', cat='Base', frag='OH- / RO- / O2-', site='Basic oxygen anion',
        acid=('Critical', 'Neutralised by acids.'),
        base=('Resistant', 'Already basic.'),
        hyd=('Resistant', 'Does not hydrolyse.'),
        photo=('Low', 'Photostable.'),
        therm=('Low', 'Thermally stable.'),
        ox=('Resistant', 'Not oxidisable.'),
        sec=('Critical', 'Acids, esters, amides, halides', 'Raises microenvironmental pH: catalyses saponification, aldol/Cannizzaro chemistry, E2/SN2; neutralises acids.'),
        partners=['acid_weak', 'acid_strong', 'ester', 'acylating', 'alkylating', 'carbonyl'], provides=['base_strong', 'base']),
    'water': dict(name='Water', cat='Solvent / Reagent', frag='H2O', site='Nucleophile / proton donor-acceptor',
        acid=('Resistant', 'Amphoteric spectator.'),
        base=('Resistant', 'Amphoteric spectator.'),
        hyd=('Resistant', 'Reagent for hydrolysis of susceptible groups.'),
        photo=('Low', 'Transparent in the solar UV.'),
        therm=('Resistant', 'Stable.'),
        ox=('Resistant', 'Stable.'),
        sec=('Moderate', 'Electrophiles', 'Hydrolysis reagent and hydrogen-bonding medium; mobilises ions (moisture-mediated solid-state reactions).'),
        partners=['acylating', 'ester', 'alkylating'], provides=['water']),
    'thiol': dict(name='Thiol (R-SH)', cat='Sulfur', frag='-SH', site='Thiol S-H (nucleophile / redox centre)',
        acid=('Low', 'Stable in acid.'),
        base=('High', 'Thiolate formation (pKa 8-10) makes air oxidation very fast.'),
        hyd=('Resistant', 'Hydrolytically inert.'),
        photo=('Moderate', 'S-H homolysis giving thiyl radicals.'),
        therm=('Low', 'Thermally stable.'),
        ox=('Critical', 'Oxidation to disulfide, then sulfinic / sulfonic acids.'),
        sec=('Critical', 'Oxidants, electrophiles, metals', 'Disulfide formation; Michael addition; thiolate-metal binding.'),
        partners=['oxidant', 'michael', 'alkylating', 'acylating', 'redox_metal', 'metal'], provides=['S_nuc', 'reductant']),
    'thioether': dict(name='Thioether (Sulfide)', cat='Sulfur', frag='-C-S-C-', site='Divalent sulfur lone pair (S-oxidation centre)',
        acid=('Low', 'Resistant to acid cleavage.'),
        base=('Low', 'Resistant to base.'),
        hyd=('Resistant', 'Hydrolytically inert.'),
        photo=('Moderate', 'Singlet-oxygen sensitised photo-oxidation; C-S homolysis.'),
        therm=('Moderate', 'Thermal C-S bond homolysis at high temperature.'),
        ox=('Critical', 'Rapid oxidation by air / peroxides to sulfoxide and then sulfone.'),
        sec=('High', 'Peroxides, alkylating agents', 'Oxidised by peroxides (or peroxide-forming ethers); S-alkylation to sulfonium salts.'),
        partners=['oxidant', 'peroxide_former', 'alkylating'], provides=['S_nuc']),
    'disulfide': dict(name='Disulfide (S-S)', cat='Sulfur', frag='-S-S-', site='Weak S-S bond',
        acid=('Low', 'Stable in acid.'),
        base=('Moderate', 'Thiolate-disulfide exchange and beta-elimination in base.'),
        hyd=('Resistant', 'Hydrolytically stable.'),
        photo=('High', 'S-S homolysis by UV to thiyl radicals.'),
        therm=('Moderate', 'Thermal scrambling of S-S bonds.'),
        ox=('Moderate', 'Oxidation to thiosulfinate, then sulfonic acids.'),
        sec=('High', 'Thiols, reductants', 'Thiol-disulfide exchange; reduction to thiols.'),
        partners=['S_nuc', 'reductant'], provides=[]),
    'sulfoxide': dict(name='Sulfoxide', cat='Sulfur', frag='-S(=O)-', site='Sulfinyl sulfur / oxygen',
        acid=('Moderate', 'Pummerer rearrangement under acidic activation.'),
        base=('Low', 'Stable in base.'),
        hyd=('Resistant', 'Hydrolytically stable.'),
        photo=('Moderate', 'Photo-deoxygenation to sulfide or photo-oxidation.'),
        therm=('High', 'Thermal syn-elimination (sulfoxide pyrolysis) giving alkene + sulfenic acid.'),
        ox=('High', 'Further oxidation to the sulfone.'),
        sec=('Moderate', 'Reductants, anhydrides', 'Deoxygenation to sulfide; Pummerer activation by anhydrides.'),
        partners=['reductant', 'acylating', 'oxidant'], provides=[]),
    'sulfone': dict(name='Sulfone', cat='Sulfur', frag='-SO2-', site='Sulfonyl sulfur (very stable)',
        acid=('Resistant', 'Resistant.'),
        base=('Low', 'Alpha-deprotonation only with strong bases.'),
        hyd=('Resistant', 'Hydrolytically inert.'),
        photo=('Low', 'Photostable.'),
        therm=('Low', 'Thermally stable (SO2 extrusion only at very high T).'),
        ox=('Resistant', 'Fully oxidised.'),
        sec=('Resistant', 'None', 'Chemically inert under formulation conditions.'),
        partners=[], provides=[]),
    'sulfonamide': dict(name='Sulfonamide', cat='Sulfur / Nitrogen', frag='-SO2NH-', site='Sulfonamide N and S-N bond',
        acid=('Low', 'Cleaved only by boiling strong acid.'),
        base=('Moderate', 'N-H deprotonation (pKa ~6-10) gives soluble salts.'),
        hyd=('Resistant', 'Hydrolytically stable at ambient conditions.'),
        photo=('High', 'Strong UV chromophore; S-N cleavage and SO2 extrusion on photolysis.'),
        therm=('Low', 'High thermal stability.'),
        ox=('Low', 'Resistant to oxidation.'),
        sec=('Moderate', 'Bases, metal cations', 'Salt formation with bases; chelation to divalent metals.'),
        partners=['base', 'base_strong', 'metal'], provides=['acid_weak']),
    'sulfonic': dict(name='Sulfonic Acid / Sulfonate', cat='Sulfur', frag='-SO3H / -SO3-', site='Sulfonate group (strong acid)',
        acid=('Resistant', 'Stable.'),
        base=('Moderate', 'Fully ionised; forms salts.'),
        hyd=('Resistant', 'Hydrolytically stable (aryl sulfonic acids desulfonate only in hot aqueous acid).'),
        photo=('Low', 'Photostable (aliphatic).'),
        therm=('Low', 'Thermally stable.'),
        ox=('Resistant', 'Resistant.'),
        sec=('High', 'Bases, amines', 'Salt formation with amines/metal ions.'),
        partners=['base', 'N_nuc', 'metal', 'alkylating'], provides=['acid_strong']),
    'sulfonate_ester': dict(name='Sulfonate / Sulfate Ester (Alkylating)', cat='Sulfur / Ester', frag='-SO2-O-R', site='Electrophilic alkyl carbon and sulfonyl sulfur',
        acid=('Moderate', 'Acid hydrolysis to sulfonic acid + alcohol.'),
        base=('Critical', 'Hydroxide attacks carbon (SN2) or sulfur, releasing sulfonate.'),
        hyd=('Moderate', 'Solvolysis in water / alcohols (alkyl sulfonates are alkylating agents).'),
        photo=('Low', 'Photostable.'),
        therm=('Moderate', 'Thermal alkylation / elimination.'),
        ox=('Resistant', 'Resistant.'),
        sec=('Critical', 'Amines, thiols, water', 'Alkylates nucleophiles (N-, S-, O-alkylation).'),
        partners=['N_nuc', 'S_nuc', 'O_nuc', 'water', 'base_strong', 'acid_weak'], provides=['alkylating']),
    'sulfonyl_halide': dict(name='Sulfonyl Halide', cat='Sulfur', frag='-SO2X', site='Electrophilic sulfonyl sulfur',
        acid=('High', 'Hydrolysis to sulfonic acid + HX.'),
        base=('Critical', 'Rapid hydrolysis.'),
        hyd=('High', 'Moisture-sensitive; hydrolyses in water.'),
        photo=('Moderate', 'Photolytic S-X cleavage.'),
        therm=('Moderate', 'Thermal SO2 extrusion.'),
        ox=('Low', 'Resistant.'),
        sec=('Critical', 'Amines, alcohols', 'Sulfonylation of amines (sulfonamides) and alcohols (sulfonate esters).'),
        partners=['N_nuc', 'O_nuc', 'water', 'base'], provides=['acylating', 'acid_strong']),
    'sulfite': dict(name='Sulfite / Bisulfite (Reductant)', cat='Sulfur / Inorganic', frag='SO3(2-) / HSO3-', site='Sulfite sulfur (nucleophile / reductant)',
        acid=('High', 'Releases SO2 in acid.'),
        base=('Low', 'Stable as sulfite.'),
        hyd=('Resistant', 'Stable.'),
        photo=('Moderate', 'Photo-oxidation.'),
        therm=('Moderate', 'Loses SO2 on heating (bisulfite/metabisulfite).'),
        ox=('Critical', 'Rapid oxidation to sulfate by O2 / peroxides.'),
        sec=('High', 'Carbonyls, oxidants, quinones', 'Adds to aldehydes/ketones (bisulfite adducts); reduces peroxides, quinones and diazonium compounds.'),
        partners=['carbonyl', 'oxidant', 'michael'], provides=['reductant', 'S_nuc']),
    'sulfate': dict(name='Sulfate / Sulfate Ester', cat='Sulfur / Inorganic', frag='-OSO3-', site='Sulfate oxygen / S-O-C linkage',
        acid=('Moderate', 'Alkyl sulfate esters hydrolyse in acid to alcohol + hydrogen sulfate; inorganic sulfate is inert.'),
        base=('Low', 'Stable.'),
        hyd=('Low', 'Alkyl sulfates hydrolyse slowly; inorganic sulfate is inert.'),
        photo=('Low', 'Photostable.'),
        therm=('Low', 'Thermally stable.'),
        ox=('Resistant', 'Fully oxidised.'),
        sec=('Low', 'Metal ions', 'Ion pairing / precipitation with Ca2+, Ba2+, Pb2+.'),
        partners=['metal'], provides=[]),
    'phosphate': dict(name='Phosphate / Phosphonate', cat='Phosphorus', frag='-P(=O)(O)O-', site='Phosphoryl phosphorus and P-O-C ester bonds',
        acid=('Moderate', 'Acid-catalysed P-O-C hydrolysis (esters) - slow.'),
        base=('Moderate', 'Base-catalysed ester hydrolysis (triesters faster).'),
        hyd=('Low', 'Slow neutral hydrolysis.'),
        photo=('Low', 'Photostable.'),
        therm=('Low', 'Thermally stable; condensation to pyrophosphates on dehydration.'),
        ox=('Resistant', 'Fully oxidised.'),
        sec=('Moderate', 'Divalent metal ions', 'Chelation / precipitation with Ca2+, Mg2+, Al3+, Fe3+.'),
        partners=['metal', 'base_strong'], provides=[]),
    'phosphine': dict(name='Phosphine', cat='Phosphorus', frag='-PR3', site='Phosphorus lone pair',
        acid=('Moderate', 'Protonation to phosphonium in strong acid.'),
        base=('Low', 'Stable.'),
        hyd=('Resistant', 'Hydrolytically stable.'),
        photo=('Moderate', 'Photo-oxidation.'),
        therm=('Low', 'Thermally stable.'),
        ox=('Critical', 'Rapid oxidation to phosphine oxide.'),
        sec=('Critical', 'Azides, oxidants, alkyl halides', 'Staudinger reaction with azides; quaternisation with alkyl halides.'),
        partners=['oxidant', 'alkylating', 'reductant'], provides=['reductant']),
    'alkyl_halide': dict(name='Alkyl Halide', cat='Halide', frag='-C-X', site='Electrophilic carbon and C-X bond',
        acid=('Low', 'Stable; tertiary/benzylic halides solvolyse (SN1).'),
        base=('High', 'SN2 substitution or E2 dehydrohalogenation by hydroxide/bases.'),
        hyd=('Moderate', 'Solvolysis in water, fast for tertiary/benzylic/allylic and alpha-halo ethers.'),
        photo=('Moderate', 'C-X homolysis under UV (C-I > C-Br > C-Cl).'),
        therm=('Moderate', 'Thermal dehydrohalogenation.'),
        ox=('Low', 'Resistant.'),
        sec=('High', 'Amines, thiols, carboxylates', 'SN2 alkylation of nucleophiles (N-, S-, O-alkylation, quaternisation).'),
        partners=['N_nuc', 'S_nuc', 'O_nuc', 'base_strong', 'base', 'water'], provides=['alkylating']),
    'aryl_halide': dict(name='Aryl / Vinyl Halide', cat='Halide', frag='Ar-X / C=C-X', site='sp2 C-X bond',
        acid=('Resistant', 'Resistant.'),
        base=('Low', 'Stable unless activated by ortho/para electron-withdrawing groups (SNAr).'),
        hyd=('Resistant', 'Hydrolytically stable.'),
        photo=('Moderate', 'Photodehalogenation (C-I > C-Br > C-Cl much greater than C-F).'),
        therm=('Low', 'Thermally stable.'),
        ox=('Low', 'Resistant.'),
        sec=('Low', 'Strong nucleophiles, metals', 'SNAr with activated rings; metal-catalysed coupling.'),
        partners=['N_nuc', 'S_nuc', 'O_nuc', 'base_strong', 'metal'], provides=[]),
    'cf3': dict(name='Trifluoromethyl / Perfluoroalkyl', cat='Halide', frag='-CF3', site='C-F bonds (very strong)',
        acid=('Resistant', 'Resistant.'),
        base=('Low', 'Stable (haloform-type loss only when alpha to carbonyl).'),
        hyd=('Resistant', 'Stable.'),
        photo=('Low', 'Photostable.'),
        therm=('Resistant', 'Stable.'),
        ox=('Resistant', 'Resistant.'),
        sec=('Resistant', 'None', 'Chemically inert.'),
        partners=[], provides=[]),
    'halide_ion': dict(name='Halide Anion', cat='Inorganic', frag='X-', site='Halide ion',
        acid=('Resistant', 'Spectator.'),
        base=('Resistant', 'Spectator.'),
        hyd=('Resistant', 'Spectator.'),
        photo=('Low', 'Photostable (I- can photo-oxidise).'),
        therm=('Resistant', 'Stable.'),
        ox=('Moderate', 'Oxidised to X2 by strong oxidants (I- and Br- readily).'),
        sec=('Moderate', 'Oxidants, electrophiles', 'Nucleophilic catalysis; oxidised by peroxides/metal oxidants.'),
        partners=['oxidant', 'alkylating'], provides=['halide']),
    'hydrogen_halide': dict(name='Hydrogen Halide (Strong Acid)', cat='Inorganic acid', frag='H-X', site='Proton (strong acid)',
        acid=('Resistant', 'Acid itself.'),
        base=('Critical', 'Neutralised by bases.'),
        hyd=('Resistant', 'Ionises in water.'),
        photo=('Low', 'Photostable.'),
        therm=('Low', 'Stable.'),
        ox=('Moderate', 'Oxidised by strong oxidants to X2.'),
        sec=('Critical', 'Bases, acid-labile groups', 'Protonates amines/bases; catalyses ester/acetal/amide hydrolysis.'),
        partners=['base', 'base_strong', 'N_nuc', 'ester', 'acetal'], provides=['acid_strong', 'halide']),
    'hypohalite': dict(name='Hypohalite / Active Halogen (Oxidant)', cat='Oxidant', frag='-O-X / N-X', site='Electrophilic halogen (oxidant)',
        acid=('Critical', 'Releases X2 / HOX in acid.'),
        base=('Moderate', 'Disproportionates.'),
        hyd=('Moderate', 'Decomposes in water.'),
        photo=('Critical', 'Photolysis to radicals.'),
        therm=('High', 'Thermal decomposition.'),
        ox=('Resistant', 'Oxidant.'),
        sec=('Critical', 'Amines, sulfides, alkenes', 'Oxidises/chlorinates S, N, alkenes and aromatic rings.'),
        partners=['N_nuc', 'S_nuc', 'ene', 'reductant'], provides=['oxidant']),
    'alkene': dict(name='Alkene (Olefin)', cat='Unsaturation', frag='-C=C-', site='Pi bond and allylic C-H',
        acid=('Moderate', 'Acid-catalysed Markovnikov hydration or carbocation oligomerisation.'),
        base=('Low', 'Resistant unless conjugated with electron-withdrawing groups.'),
        hyd=('Resistant', 'Resistant.'),
        photo=('High', 'cis-trans photo-isomerisation and [2+2] photodimerisation.'),
        therm=('Moderate', 'Thermal isomerisation / Diels-Alder dimerisation (dienes).'),
        ox=('Critical', 'Epoxidation by peroxides; allylic autoxidation to hydroperoxides / enones.'),
        sec=('Moderate', 'Peroxides, radicals, electrophiles', 'Epoxidation by peroxides; radical addition / crosslinking.'),
        partners=['oxidant', 'peroxide_former', 'ene'], provides=['ene']),
    'alkyne': dict(name='Alkyne', cat='Unsaturation', frag='-C#C-', site='Pi system (and terminal C-H)',
        acid=('Moderate', 'Acid/metal-catalysed hydration to ketone.'),
        base=('Moderate', 'Terminal alkynes deprotonate to acetylides.'),
        hyd=('Resistant', 'Resistant.'),
        photo=('Moderate', 'Photo-cycloaddition / polymerisation.'),
        therm=('Moderate', 'Thermal oligomerisation.'),
        ox=('Moderate', 'Oxidative cleavage to acids / diketones.'),
        sec=('Moderate', 'Cu+/Ag+ ions, azides', 'Metal acetylide formation; azide-alkyne cycloaddition.'),
        partners=['metal', 'N_nuc'], provides=[]),
    'aromatic': dict(name='Aromatic Ring', cat='Aromatic', frag='c1ccccc1', site='Aromatic pi-system and ring C-H',
        acid=('Resistant', 'Resistant to acid solvolysis.'),
        base=('Resistant', 'Resistant to base.'),
        hyd=('Resistant', 'Hydrolytically inert.'),
        photo=('High', 'Strong UV absorption (254-280 nm) leading to triplet-state chemistry.'),
        therm=('Resistant', 'Aromatic stabilisation.'),
        ox=('Moderate', 'Hydroxyl-radical attack forms hydroxylated (phenolic) products; electron-rich rings oxidise faster.'),
        sec=('Moderate', 'Electrophiles, radicals', 'Electrophilic substitution (nitration, halogenation) with activated rings; pi-stacking / charge-transfer complexes.'),
        partners=['oxidant', 'acid_strong'], provides=[]),
    'hetero_rich': dict(name='Heteroaromatic Ring (Electron-Rich: Furan / Thiophene / Pyrrole / Azole)', cat='Aromatic heterocycle', frag='c1ccoc1', site='Electron-rich ring C-H and heteroatom',
        acid=('High', 'Acid-catalysed ring opening / polymerisation (furan, pyrrole).'),
        base=('Low', 'Stable; N-H azoles deprotonate with strong base.'),
        hyd=('Low', 'Stable.'),
        photo=('High', 'Photo-oxidation with singlet oxygen ([4+2]) and photoisomerisation.'),
        therm=('Low', 'Thermally stable.'),
        ox=('High', 'Oxidative ring opening / S-oxidation / polymerisation.'),
        sec=('Moderate', 'Oxidants, electrophiles', 'Oxidised/halogenated readily; azole N coordinates metal ions.'),
        partners=['oxidant', 'metal', 'alkylating'], provides=[]),
    'hetero_azine': dict(name='Heteroaromatic Ring (Pi-Deficient: Pyridine / Azine / Azole N)', cat='Aromatic heterocycle', frag='c1ccncc1', site='Ring nitrogen lone pair',
        acid=('Moderate', 'Protonation to the azinium salt (pyridine pKa ~5).'),
        base=('Resistant', 'Resistant.'),
        hyd=('Resistant', 'Hydrolytically stable.'),
        photo=('Moderate', 'Photo-isomerisation and photo-oxidation.'),
        therm=('Low', 'Thermally stable.'),
        ox=('Moderate', 'N-oxidation by peroxides/peracids.'),
        sec=('Moderate', 'Alkylating agents, metal ions, acids', 'N-alkylation (quaternisation); metal coordination; salt formation.'),
        partners=['alkylating', 'metal', 'oxidant', 'acid_weak', 'acid_strong'], provides=['base']),
    'benzylic_ch': dict(name='Benzylic / Allylic C-H', cat='Hydrocarbon', frag='Ar-CH< / C=C-CH<', site='Weak benzylic / allylic C-H bonds (BDE ~85-90 kcal/mol)',
        acid=('Resistant', 'Resistant.'),
        base=('Resistant', 'Resistant (unless acidified by conjugation).'),
        hyd=('Resistant', 'Resistant.'),
        photo=('Moderate', 'Photo-initiated hydrogen abstraction.'),
        therm=('Moderate', 'Thermal autoxidation initiated by trace radicals.'),
        ox=('High', 'Autoxidation to hydroperoxide, then alcohol / ketone / acid.'),
        sec=('Moderate', 'Radical initiators, peroxides, metals', 'Radical chain oxidation accelerated by peroxides and redox metals.'),
        partners=['oxidant', 'peroxide_former', 'redox_metal'], provides=[]),
    'metal': dict(name='Metal Cation (Main-Group)', cat='Inorganic ion', frag='M(n+)', site='Lewis-acidic cation',
        acid=('Resistant', 'Spectator.'),
        base=('Moderate', 'Precipitates as hydroxide/oxide at high pH.'),
        hyd=('Resistant', 'Spectator.'),
        photo=('Low', 'Photostable.'),
        therm=('Resistant', 'Stable.'),
        ox=('Resistant', 'Not oxidisable (fixed oxidation state).'),
        sec=('High', 'Carboxylates, phenolates, phosphates, carbonyls', 'Chelation / salt formation with anionic donors; Lewis-acid activation of carbonyls and esters (Mg2+, Ca2+, Zn2+, Al3+).'),
        partners=['acid_weak', 'acid_strong', 'N_nuc', 'carbonyl', 'ester'], provides=['metal', 'lewis_acid']),
    'metal_transition': dict(name='Transition-Metal Ion (Redox-Active)', cat='Inorganic ion', frag='M(n+)', site='Redox-active metal centre',
        acid=('Resistant', 'Spectator.'),
        base=('Moderate', 'Precipitates as hydroxide / oxide at high pH.'),
        hyd=('Resistant', 'Spectator.'),
        photo=('Moderate', 'Ligand-to-metal charge-transfer photochemistry.'),
        therm=('Low', 'Stable.'),
        ox=('High', 'Cycles between oxidation states (Fenton-type radical generation).'),
        sec=('Critical', 'Peroxides, thiols, phenols, chelators', 'Catalyses autoxidation and peroxide decomposition (radicals); coordinates N/O/S donors.'),
        partners=['oxidant', 'peroxide_former', 'S_nuc', 'N_nuc', 'acid_weak', 'enol'], provides=['metal', 'redox_metal', 'lewis_acid']),
    'metal_oxide': dict(name='Inorganic Oxide / Oxo-Metal Species', cat='Inorganic', frag='M=O', site='Oxide surface / oxo-metal centre',
        acid=('Moderate', 'Basic oxides dissolve in acid; SiO2/TiO2 are inert.'),
        base=('Moderate', 'Amphoteric / basic surface sites.'),
        hyd=('Low', 'Insoluble; basic oxides hydrate to hydroxides.'),
        photo=('Moderate', 'TiO2 / ZnO are photocatalytic (generate radicals under UV); SiO2 is inert.'),
        therm=('Resistant', 'Stable.'),
        ox=('Moderate', 'High-valent oxo-metals (MnO4-, CrO4 2-) are strong oxidants.'),
        sec=('Moderate', 'Acids, adsorbable species', 'Surface adsorption; basic oxides neutralise acids; photocatalytic oxidation (TiO2/ZnO).'),
        partners=['acid_weak', 'acid_strong'], provides=[]),
    'hydride': dict(name='Hydride Donor (Borohydride)', cat='Reductant', frag='BH4-', site='Hydridic B-H',
        acid=('Critical', 'Evolves H2 in acid.'),
        base=('Low', 'Stable in base.'),
        hyd=('High', 'Slowly hydrolyses to H2 + borate in water.'),
        photo=('Low', 'Photostable.'),
        therm=('Moderate', 'Decomposes at high temperature.'),
        ox=('Critical', 'Strong reductant (readily oxidised).'),
        sec=('Critical', 'Carbonyls, nitro, disulfides', 'Reduces aldehydes/ketones, imines, disulfides (and nitro with catalysts).'),
        partners=['carbonyl', 'acylating'], provides=['reductant']),
    'carbonate_ion': dict(name='Carbonate / Bicarbonate', cat='Inorganic base', frag='CO3(2-) / HCO3-', site='Basic carbonate oxygen',
        acid=('Critical', 'Neutralised by acid with CO2 evolution.'),
        base=('Resistant', 'Stable.'),
        hyd=('Low', 'Buffers water (pH ~8-11).'),
        photo=('Low', 'Photostable.'),
        therm=('Moderate', 'Bicarbonate decomposes to carbonate + CO2 + H2O on heating.'),
        ox=('Resistant', 'Not oxidisable.'),
        sec=('High', 'Acids, esters', 'Neutralises acids (salt formation, CO2); raises pH and catalyses base hydrolysis.'),
        partners=['acid_weak', 'acid_strong', 'ester', 'acylating'], provides=['base', 'base_strong']),
    'nitrite': dict(name='Nitrite / Nitrosating Agent', cat='Inorganic', frag='NO2-', site='Nitrite nitrogen / N=O',
        acid=('Critical', 'Forms HONO / N2O3 / NO+ (nitrosating species) in acid.'),
        base=('Low', 'Stable as nitrite.'),
        hyd=('Low', 'Stable.'),
        photo=('High', 'UV releases NO / NO2.'),
        therm=('Moderate', 'Disproportionates on heating.'),
        ox=('Moderate', 'Oxidised to nitrate.'),
        sec=('Critical', 'Secondary/tertiary amines, aromatic amines', 'N-nitrosation of amines (nitrosamines); diazotisation of aromatic amines.'),
        partners=['N_nuc', 'N_nuc_sec', 'S_nuc', 'acid_weak', 'acid_strong'], provides=['nitrosating', 'oxidant']),
    'nitrate': dict(name='Nitrate', cat='Inorganic', frag='NO3-', site='Nitrate oxygen',
        acid=('Resistant', 'Spectator (HNO3 is a strong oxidising acid).'),
        base=('Resistant', 'Spectator.'),
        hyd=('Resistant', 'Spectator.'),
        photo=('Moderate', 'Photolysis to nitrite / NO2.'),
        therm=('Low', 'Stable.'),
        ox=('Resistant', 'Fully oxidised.'),
        sec=('Low', 'Reductants', 'Weak oxidant in neutral media; strong oxidant as HNO3.'),
        partners=['reductant'], provides=['oxidant']),
    'hydrocarbon': dict(name='Aliphatic / Polar Heteroatom Framework', cat='Other', frag='C-C / Heteroatom', site='Aliphatic C-H and polar heteroatom centres',
        acid=('Low', 'No acid-labile functional group.'),
        base=('Low', 'No base-labile functional group.'),
        hyd=('Low', 'No hydrolysable group.'),
        photo=('Low', 'UV absorption depends on chromophores.'),
        therm=('Low', 'Thermal bond homolysis at high temperature.'),
        ox=('Moderate', 'Radical autoxidation across C-H positions (tertiary > secondary > primary).'),
        sec=('Low', 'Co-reactants', 'Non-covalent physical interactions only.'),
        partners=[], provides=[]),
}

# --------------------------------------------------------------------------------------
# Detection rules: (family, rule key, SMARTS, profile key, override, options)
#   * key atom = atom mapped :1 (first atom if none); one rule per (family, key atom) wins,
#     so list specific variants BEFORE the generic rule of the same family.
#   * override: dict cond -> (vulnerability, text) applied on top of the profile.
#   * options: check=callable(mol, match)->bool ; not_in=[families whose atoms must not overlap] ;
#     claim2=True: the atom mapped :2 is claimed too (a symmetric pair such as an imide / anhydride
#     is reported once).
# --------------------------------------------------------------------------------------
_EWG = "[$([N+](=O)[O-]),$(C=O),$(C#N),$(S(=O)=O)]"
_NEXC = "!$(N-C=[O,S,N]);!$(N-[SX4,PX4]=O);!$(N-[N,O])"
_ETH_C = "[CX4;!$(C(O)O)]"


def _o(**kw: Tuple[str, str]) -> Dict[str, Tuple[str, str]]:
    return dict(kw)


RULES: List[Tuple[str, str, str, str, Optional[dict], Optional[dict]]] = []


def _r(fam: str, key: str, smarts: str, pkey: str, over: Optional[dict] = None, **opts: Any) -> None:
    RULES.append((fam, key, smarts, pkey, over, opts or None))


# ---- quinones (must come first so ketone / michael / alkene rules can exclude them)
_r("quinone", "quinone_para", "[#6:1]1(=[O,N])[#6]=[#6][#6](=[O,N])[#6]=[#6]1", "quinone")
_r("quinone", "quinone_ortho", "[#6:1]1(=O)[#6](=[O,N])[#6]=[#6][#6]=[#6]1", "quinone")

# ---- acyl derivatives (key = carbonyl carbon)
_r("acyl", "acyl_halide", "[CX3:1](=O)[F,Cl,Br,I]", "acyl_halide")
_r("acyl", "anhydride", "[CX3:1](=O)[OX2][CX3:2]=O", "anhydride", claim2=True)
_r("acyl", "carbonate_ion", "[CX3:1](=O)([OX2H1,OX1-])[OX2H1,OX1-]", "carbonate_ion")
_r("acyl", "carbonate", "[CX3:1](=O)([OX2][#6])[OX2][#6]", "carbonate")
_r("acyl", "carbamate_tert", "[CX3:1](=O)([NX3])[OX2][CX4;H0]([#6])([#6])[#6]", "carbamate",
   _o(acid=("Critical", "tert-Butyl-type carbamate (Boc-like): A_AL1 / E1 cleavage loses the alkene, then CO2, giving the free amine.")))
_r("acyl", "carbamate", "[CX3:1](=O)([NX3])[OX2][#6]", "carbamate")
_r("acyl", "urea", "[CX3:1](=O)([NX3])[NX3]", "urea")
_r("acyl", "imide", "[CX3:1](=O)[NX3][CX3:2]=O", "imide", claim2=True)
_r("acyl", "beta_lactam", "[CX3;r4:1](=O)@[NX3;r4]", "beta_lactam")
_r("acyl", "lactam", "[CX3;R:1](=O)@[NX3]", "lactam")
_r("acyl", "amide_arom", "[c:1](=O)[nX3]", "amide",
   _o(acid=("Low", "Aromatic (pyridone-type) amide: resists acid hydrolysis."), base=("Low", "Aromatic (pyridone-type) amide: resists base hydrolysis."),
      hyd=("Resistant", "Aromatic (pyridone-type) amide is hydrolytically inert.")))
_r("acyl", "amide_anilide", "[CX3:1](=O)[NX3][c]", "amide",
   _o(photo=("High", "Anilides undergo photo-Fries rearrangement to amino-aryl ketones.")))
_r("acyl", "amide_formamide", "[CX3H1:1](=O)[NX3]", "amide",
   _o(hyd=("Moderate", "Formamides hydrolyse faster than other amides."), acid=("High", "Acid hydrolysis of the formamide is comparatively easy.")))
_r("acyl", "amide_hydrazide", "[CX3:1](=O)[NX3][NX3]", "amide",
   _o(acid=("High", "Hydrazides hydrolyse more readily than amides.")))
_r("acyl", "amide_tert", "[CX3:1](=O)[NX3;H0]", "amide",
   _o(base=("Low", "Tertiary amides resist base hydrolysis (no N-H, hindered).")))
_r("acyl", "amide", "[CX3:1](=O)[NX3]", "amide")
_r("acyl", "lactone_4", "[CX3;r4:1](=O)@[OX2]", "lactone",
   _o(acid=("Critical", "Strained beta-lactone: acid opens the ring rapidly."), hyd=("Critical", "Strained beta-lactone hydrolyses rapidly with ring-strain relief (~23 kcal/mol)."),
      therm=("High", "Thermal decarboxylation to alkene / polymerisation.")))
_r("acyl", "lactone_enol", "[CX3;R:1](=O)(@[OX2])[CX3]([OX2H1])=[CX3]", "lactone",
   _o(base=("Low", "Vinylogous acid (ene-diol lactone): deprotonated to the enolate, which resists hydroxide attack."), acid=("Low", "Vinylogous acid lactone resists acid opening."),
      hyd=("Low", "Vinylogous acid lactone resists neutral hydrolysis.")))
_r("acyl", "lactone_arom", "[c:1](=O)[o]", "lactone",
   _o(base=("High", "Aromatic (coumarin-type) lactone opens in base to the hydroxycinnamate."), acid=("Low", "Aromatic lactone resists acid opening."),
      photo=("High", "Coumarins undergo [2+2] photodimerisation and photo-oxidation.")))
_r("acyl", "lactone", "[CX3;R:1](=O)@[OX2][#6]", "lactone")
_r("acyl", "ester_tert", "[CX3:1](=O)!@[OX2][CX4;H0]([#6])([#6])[#6]", "ester",
   _o(acid=("Critical", "Acid-labile tert-alkyl ester: A_AL1 / E1 cleavage releases the acid + alkene (e.g. isobutene)."), base=("Low", "Sterically hindered; resists saponification."),
      therm=("High", "Thermal E1 / syn-elimination expels the alkene (isobutene) to give the acid.")))
_r("acyl", "ester_formate", "[CX3H1:1](=O)!@[OX2][#6]", "ester",
   _o(hyd=("Moderate", "Formate esters hydrolyse faster than other alkyl esters."), base=("Critical", "Formate esters saponify very rapidly.")))
_r("acyl", "ester_alphaEWG", "[CX3:1](=O)(!@[OX2][#6])[CX4;$(C(F)F),$(C(Cl)Cl),$(C[N+](=O)[O-])]", "ester",
   _o(base=("Critical", "Alpha electron-withdrawing group accelerates saponification."), hyd=("Moderate", "Alpha-EWG activation raises neutral hydrolysis rate.")))
_r("acyl", "ester_aryl", "[CX3:1](=O)!@[OX2;!$(O[OX2,#7,#16,#15])][c,$([CX3]=[CX3])]", "ester_aryl")
_r("acyl", "ester", "[CX3:1](=O)!@[OX2;!$(O[OX2,#7,#16,#15])][#6]", "ester")
_r("acyl", "thioester", "[CX3:1](=O)[SX2][#6]", "thioester")
_r("acyl", "acid_malonic", "[CX3:1](=O)([OX2H1])[CX4][CX3](=O)[OX2H1]", "acid",
   _o(therm=("High", "Malonic-type acid: thermal decarboxylation through a cyclic transition state releases CO2.")))
_r("acyl", "acid_betaketo", "[CX3:1](=O)([OX2H1])[CX4][CX3]=O", "acid",
   _o(therm=("High", "Beta-keto acid: cyclic six-membered transition state releases CO2 on gentle heating.")))
_r("acyl", "acid_alphaEWG", "[CX3:1](=O)([OX2H1])[CX4;$(C(F)F),$(C(Cl)Cl),$(C[N+](=O)[O-])]", "acid",
   _o(therm=("High", "Electron-withdrawing alpha-substituent stabilises the carbanion-like transition state: decarboxylation on heating.")))
_r("acyl", "acid_aryl_o", "[CX3:1](=O)([OX2H1])c:c[OX2H1,NX3;H1,H2]", "acid",
   _o(therm=("High", "Electron-rich aryl acid (ortho OH/NH2): protodecarboxylation on heating.")))
_r("acyl", "acid_aryl_p", "[CX3:1](=O)([OX2H1])c1ccc([OX2H1,NX3;H1,H2])cc1", "acid",
   _o(therm=("High", "Electron-rich aryl acid (para OH/NH2): protodecarboxylation on heating.")))
_r("acyl", "acid_picolinic", "[CX3:1](=O)([OX2H1])c:n", "acid",
   _o(therm=("High", "2-Pyridyl acid: zwitterionic (Hammick-type) decarboxylation on heating.")))
_r("acyl", "acid", "[CX3:1](=O)[OX2H1]", "acid")
_r("acyl", "carboxylate", "[CX3:1](=O)[OX1-]", "carboxylate")
_r("acyl", "aldehyde", "[CX3;$([CH1](=O)[#6]),$([CH2]=O):1]=O", "aldehyde")
_r("acyl", "ketone", "[#6][CX3:1](=O)[#6]", "ketone", not_in=["quinone"])
_r("michael", "michael", "[CX3;!a:1](=O)[CX3;!a]=[CX3;!a]", "michael", not_in=["quinone"])
_r("cs", "thiocarbonyl", "[CX3:1]=[SX1]", "thiocarbonyl")
_r("cs", "isocyanate", "[CX2:1](=[NX2])=O", "isocyanate")
_r("cs", "isothiocyanate", "[CX2:1](=[NX2])=S", "isothiocyanate")

# ---- nitrogen groups
_r("nitrile", "nitrile", "[CX2:1]#[NX1]", "nitrile")
_r("nitro", "nitrate", "[N+:1](=O)([O-])[O-,OX2H1]", "nitrate")
_r("nitro", "nitro_aliph", "[CX4][N+:1](=O)[O-]", "nitro",
   _o(base=("High", "Alpha-H nitro compounds ionise to nitronate salts (pKa ~10); nitronates are prone to Nef-type chemistry.")))
_r("nitro", "nitro", "[N+:1](=O)[O-]", "nitro")
_r("nitro", "nitro_hyper", "[NX3:1](=O)=O", "nitro")
_r("nitro", "nitrosamine", "[NX3][NX2:1]=O", "nitrosamine")
_r("nitro", "nitroso", "[#6][NX2:1]=O", "nitroso")
_r("nitro", "nitrite", "[NX2:1](=O)[O-,OX2]", "nitrite")
_r("azo", "azo", "[#6][NX2:1]=[NX2][#6]", "azo")
_r("azo", "azide", "[N:1]=[N+]=[N-]", "azide")
_r("azo", "azide2", "[N-:1]=[N+]=N", "azide")
_r("noxide", "n_oxide", "[#7+;!$([N+](=O)[O-]);!$([N+]=O):1][O-]", "n_oxide")
_r("amine", "amine_tert", f"[NX3;H0;{_NEXC};!$(N-a);!$(N=*):1]([#6])([#6])[#6]", "amine_tert")
_r("amine", "amine", f"[NX3;H2,H1;{_NEXC};!$(N-a):1]", "amine")
_r("amine", "ammonia", "[NX3;H3;D0:1]", "amine")
_r("amine", "arylamine", f"[NX3;{_NEXC};!$(N=*):1]-[c;$(c1ccccc1)]", "arylamine")
_r("amine", "arylamine_hetero", f"[NX3;{_NEXC};!$(N=*):1]-[a;!$(c1ccccc1)]", "arylamine",
   _o(ox=("Moderate", "Heteroaryl amine (amidine-like, delocalised lone pair): much less easily oxidised than an aniline."),
      photo=("Moderate", "Heteroaryl amines are moderately photolabile."), acid=("Low", "Weakly basic; protonation occurs on the ring nitrogen.")))
_r("amine", "quat", "[NX4+;H0:1]", "quat")
_r("amine", "ammonium", "[NX4+;!H0:1]", "ammonium")
_r("cn", "amidine", "[CX3:1](=[NX2;!$(N-[N,O]);!a])-[NX3,OX2]", "amidine")
_r("cn", "oxime", "[CX3:1]=[NX2,NX3+][OX2,OX1-]", "oxime")
_r("cn", "hydrazone", "[CX3:1]=[NX2][NX3]", "hydrazone")
_r("cn", "imine", "[CX3;!$(C-[N,O,S]):1]=[NX2;!$(N-[N,O]);!a]", "imine")
_r("nn", "hydrazine", "[NX3;!$(N-C=[O,S,N]);!$(N-[SX4,PX4]=O);!$(N=*):1]-[NX3;!$(N-C=[O,S,N]);!$(N=*)]", "hydrazine")
_r("nn", "hydroxylamine", "[NX3;!$(N-C=[O,S,N]);!$(N-[SX4,PX4]=O);!$(N=*):1]-[OX2]", "hydroxylamine")

# ---- oxygen groups (key = oxygen unless stated)
_r("oxy", "peroxide", "[OX2:1][OX2]", "peroxide")
_r("oxy", "dioxygen", "[O:1]=[O]", "dioxygen")
_r("oxy", "ozone", "[O-][O+:1]=O", "dioxygen")
_r("oxy", "hypohalite", "[OX1-:1][Cl,Br,I]", "hypohalite")
_r("oxy", "water", "[OX2H2:1]", "water")
_r("oxy", "hydroxide", "[OX1H-:1]", "strong_base")
_r("oxy", "oxide", "[O-2:1]", "strong_base")
_r("oxy", "alkoxide", "[OX1-:1][CX4]", "strong_base")
_r("oxy", "phenolate", "[OX1-:1]c", "phenol")
_r("oxy", "phenol_hindered", "[OX2H:1]c:c[CX4;H0]([CH3])([CH3])[CH3]", "phenol",
   _o(ox=("High", "Hindered phenol scavenges radicals (antioxidant) and is converted to phenoxyl / quinone-methide products.")))
_r("oxy", "phenol_ewg_o", f"[OX2H:1]c:c{_EWG}", "phenol",
   _o(ox=("Moderate", "Electron-withdrawing substituents raise the oxidation potential."), base=("Critical", "Acidified phenol (pKa below 8) is fully ionised in mild base.")))
_r("oxy", "phenol_ewg_p", f"[OX2H:1]c1ccc({_EWG})cc1", "phenol",
   _o(ox=("Moderate", "Electron-withdrawing substituents raise the oxidation potential."), base=("Critical", "Acidified phenol (pKa below 8) is fully ionised in mild base.")))
_r("oxy", "phenol", "[OX2H:1]c", "phenol")
_r("oxy", "enol", "[OX2H:1][CX3]=[CX3]", "enol")
_r("oxy", "alcohol_tert", f"[OX2H:1][CX4;H0;!$(C([OX2])[OX2])]([#6])([#6])[#6]", "alcohol_tert")
_r("oxy", "alcohol_benzylic", "[OX2H:1][CX4;!$(C([OX2])[OX2])][c,$([CX3]=[CX3])]", "alcohol",
   _o(ox=("Critical", "Benzylic / allylic carbinols oxidise very readily to carbonyls.")))
_r("oxy", "alcohol", "[OX2H:1][CX4;!$(C([OX2])[OX2])]", "alcohol")
_r("oxy", "epoxide", "[OX2;r3:1]", "epoxide")
_r("oxy", "ether_cyclic", f"[OX2;R;!r3:1]({_ETH_C})[CX4;R;!$(C(O)O)]", "ether",
   _o(ox=("Critical", "Cyclic ethers (THF, dioxane) autoxidise rapidly to hydroperoxides.")))
_r("oxy", "aryl_ether", f"[OX2;!r3:1](c){_ETH_C}", "aryl_ether")
_r("oxy", "ether", f"[OX2;!r3:1]({_ETH_C}){_ETH_C}", "ether")
_r("acetal", "acetal", "[CX4:1]([OX2][#6])[OX2][#6]", "acetal")
_r("acetal", "hemiacetal", "[CX4:1]([OX2H1])[OX2][#6]", "hemiacetal")

# ---- sulfur / phosphorus
_r("sulfur", "thiol", "[SX2H1:1]", "thiol")
_r("sulfur", "thiolate", "[SX1-:1]", "thiol")
_r("sulfur", "h2s", "[SX2H2:1]", "thiol")
_r("sulfur", "disulfide", "[SX2:1][SX2]", "disulfide")
_r("sulfur", "thioether", "[SX2:1]([#6])[#6]", "thioether")
_r("sulfur", "sulfoxide", "[SX3:1](=O)([#6])[#6]", "sulfoxide")
_r("sulfur", "sulfite", "[SX3;!$(S[#6]):1](=O)[OX2H1,OX1-]", "sulfite")
_r("sulfur", "sulfonyl_halide", "[SX4:1](=O)(=O)[F,Cl,Br,I]", "sulfonyl_halide")
_r("sulfur", "sulfonamide", "[SX4:1](=O)(=O)[NX3]", "sulfonamide")
_r("sulfur", "sulfonate_ester", "[SX4:1](=O)(=O)([#6])[OX2][#6]", "sulfonate_ester")
_r("sulfur", "dialkyl_sulfate", "[SX4:1](=O)(=O)([OX2][#6])[OX2][#6]", "sulfonate_ester")
_r("sulfur", "sulfonic", "[SX4:1](=O)(=O)([#6])[OX2H1,OX1-]", "sulfonic")
_r("sulfur", "sulfate_ester", "[SX4:1](=O)(=O)([OX2][#6])[OX2H1,OX1-]", "sulfate")
_r("sulfur", "sulfate", "[SX4:1](=O)(=O)([OX2H1,OX1-])[OX2H1,OX1-]", "sulfate")
_r("sulfur", "sulfone", "[SX4:1](=O)(=O)([#6])[#6]", "sulfone")
_r("phos", "phosphate", "[PX4:1](=O)([OX2,OX1-])[OX2,OX1-]", "phosphate")
_r("phos", "phosphine_oxide", "[PX4:1](=O)([#6])([#6])[#6]", "phosphine",
   _o(ox=("Resistant", "Phosphine oxide is already fully oxidised."), acid=("Low", "Phosphine oxide is stable in acid.")))
_r("phos", "phosphine", "[PX3:1]", "phosphine")

# ---- halogens (key = halogen atom, or carbon for CF3)
_r("cf3", "cf3", "[CX4:1](F)(F)F", "cf3")
_r("halogen", "alkyl_F", "[F:1][CX4;!$(C(F)F)]", "alkyl_halide",
   _o(acid=("Resistant", "C-F is extremely stable."), base=("Low", "C-F resists substitution."), hyd=("Resistant", "C-F is hydrolytically stable."),
      photo=("Low", "C-F is photostable."), therm=("Low", "C-F is thermally stable."), ox=("Resistant", "C-F resists oxidation.")))
_r("halogen", "alkyl_tert", "[Cl,Br,I:1][CX4;H0]([#6])([#6])[#6]", "alkyl_halide",
   _o(hyd=("High", "Tertiary halides solvolyse via SN1 in water."), acid=("Moderate", "Acid-assisted SN1 solvolysis.")))
_r("halogen", "alkyl_benzylic", "[Cl,Br,I:1][CX4][c,$([CX3]=[CX3])]", "alkyl_halide",
   _o(hyd=("High", "Benzylic / allylic halides solvolyse via SN1 in water."), acid=("Moderate", "Acid-assisted SN1 solvolysis.")))
_r("halogen", "alkyl_gem", "[Cl,Br,I:1][CX4]([Cl,Br,I])", "alkyl_halide",
   _o(base=("High", "Gem-polyhalides undergo base-induced alpha-elimination (carbene) / haloform-type cleavage.")))
_r("halogen", "alkyl_I", "[I:1][CX4]", "alkyl_halide", _o(photo=("Critical", "C-I bond homolysis under near-UV / visible light.")))
_r("halogen", "alkyl_Br", "[Br:1][CX4]", "alkyl_halide", _o(photo=("High", "C-Br homolysis under UV.")))
_r("halogen", "alkyl_halide", "[Cl:1][CX4]", "alkyl_halide",
   _o(hyd=("Low", "Primary/secondary chloride: slow neutral solvolysis; SN2 with hydroxide dominates.")))
_snar_o = ("High", "Ortho/para electron-withdrawing group activates SNAr displacement of the halide by hydroxide/amines.")
_r("halogen", "aryl_snar_p", f"[F,Cl,Br,I:1]c1ccc({_EWG})cc1", "aryl_halide", _o(base=_snar_o))
_r("halogen", "aryl_snar_o", f"[F,Cl,Br,I:1]c:c{_EWG}", "aryl_halide", _o(base=_snar_o))
_r("halogen", "aryl_snar_n", "[F,Cl,Br,I:1]c:n", "aryl_halide", _o(base=_snar_o))
_r("halogen", "aryl_I", "[I:1][c,$([CX3]=[CX3])]", "aryl_halide", _o(photo=("High", "C-I homolysis / photodehalogenation.")))
_r("halogen", "aryl_Br", "[Br:1][c,$([CX3]=[CX3])]", "aryl_halide", _o(photo=("Moderate", "C-Br photodehalogenation under UV.")))
_r("halogen", "aryl_F", "[F:1][c,$([CX3]=[CX3])]", "aryl_halide", _o(photo=("Low", "C-F is photostable.")))
_r("halogen", "aryl_halide", "[Cl:1][c,$([CX3]=[CX3])]", "aryl_halide")
_r("halogen", "halide_ion", "[F-,Cl-,Br-,I-:1]", "halide_ion")
_r("halogen", "hydrogen_halide", "[F,Cl,Br,I;H1;D0:1]", "hydrogen_halide")

# ---- carbon frameworks
_r("cc", "alkyne", "[CX2:1]#[CX2]", "alkyne")
_r("cc", "alkene", "[CX3;!a:1]=[CX3;!a]", "alkene", not_in=["quinone", "michael", "oxy"])
_r("benz", "benzylic_ch", "[CX4;H1,H2,H3;!$(C-[!#6;!#1]);$(C-[c,$([CX3]=[CX3])]):1]", "benzylic_ch")

# ---- ionic / inorganic
_r("metal", "metal_oxide", "[Ti,V,Cr,Mn,Mo,W,Ce,Si:1]=O", "metal_oxide")
_r("metal", "metal_transition", "[Ti,V,Cr,Mn,Fe,Co,Ni,Cu,Ag,Au,Pd,Pt,Mo,W,Ce:1]", "metal_transition")
_r("metal", "metal", "[Li,Na,K,Rb,Cs,Mg,Ca,Sr,Ba,Al,Ga,Zn,Cd,Hg,Sn,Pb,Bi,Sb,Zr:1]", "metal")
_r("hydride", "hydride", "[BH4-:1]", "hydride")

_COMPILED: List[Tuple[str, str, Chem.Mol, str, Optional[dict], Optional[dict], int, int]] = []
for _fam, _key, _sma, _pk, _ov, _op in RULES:
    _q = Chem.MolFromSmarts(_sma)
    if _q is None:  # fail loudly at import time rather than silently missing a group
        raise ValueError(f"invalid SMARTS for rule {_key}: {_sma}")
    _ka = next((a.GetIdx() for a in _q.GetAtoms() if a.GetAtomMapNum() == 1), 0)
    _k2 = next((a.GetIdx() for a in _q.GetAtoms() if a.GetAtomMapNum() == 2), -1)
    _COMPILED.append((_fam, _key, _q, _pk, _ov, _op, _ka, _k2))


# --------------------------------------------------------------------------------------
# Analysis
# --------------------------------------------------------------------------------------
class Site:
    __slots__ = ("rule", "pkey", "fam", "atoms", "key", "eff")

    def __init__(self, rule: str, pkey: str, fam: str, atoms: Tuple[int, ...], key: int, eff: dict):
        self.rule, self.pkey, self.fam, self.atoms, self.key, self.eff = rule, pkey, fam, atoms, key, eff


class Analysis:
    def __init__(self, mol: Chem.Mol):
        self.mol = mol
        self.sites: List[Site] = []
        self.entries: List[dict] = []
        self.classes: Set[str] = set()
        self.present: Set[str] = set()
        self._run()

    # ---- helpers on the RDKit mol
    def _h(self, i: int) -> int:
        return self.mol.GetAtomWithIdx(i).GetTotalNumHs()

    def _is_ewg(self, a: Chem.Atom) -> bool:
        m = a.GetOwningMol()
        s = a.GetSymbol()
        if s == "N":
            return a.GetFormalCharge() == 1 and any(n.GetSymbol() == "O" for n in a.GetNeighbors())
        if s == "C":
            for b in a.GetBonds():
                o = b.GetOtherAtom(a)
                if b.GetBondTypeAsDouble() == 2 and o.GetSymbol() in ("O", "S", "N"):
                    return True
                if b.GetBondTypeAsDouble() == 3:
                    return True
            return a.GetDegree() == 4 and sum(1 for n in a.GetNeighbors() if n.GetSymbol() == "F") >= 3
        if s == "S":
            return any(b.GetBondTypeAsDouble() == 2 and b.GetOtherAtom(a).GetSymbol() == "O" for b in a.GetBonds())
        del m
        return False

    def _is_donor(self, a: Chem.Atom) -> bool:
        s = a.GetSymbol()
        if s not in ("O", "N") or a.GetIsAromatic() or a.GetFormalCharge() != 0:
            return False
        if any(b.GetBondTypeAsDouble() != 1 for b in a.GetBonds()):
            return False
        for n in a.GetNeighbors():
            if n.GetSymbol() == "C" and any(b.GetBondTypeAsDouble() == 2 and b.GetOtherAtom(n).GetSymbol() in ("O", "S")
                                            for b in n.GetBonds()):
                return False
            if n.GetSymbol() in ("S", "P"):
                return False
        return True

    def _ring_sites(self) -> None:
        m = self.mol
        ri = m.GetRingInfo()
        seen: Set[Tuple[int, ...]] = set()
        for ring in ri.AtomRings():
            if len(ring) not in (5, 6):
                continue
            atoms = [m.GetAtomWithIdx(i) for i in ring]
            if not all(a.GetIsAromatic() for a in atoms):
                continue
            key = tuple(sorted(ring))
            if key in seen:
                continue
            seen.add(key)
            het = [a for a in atoms if a.GetSymbol() != "C"]
            subs = [n for a in atoms for n in a.GetNeighbors() if n.GetIdx() not in ring]
            if not het:
                donors = sum(1 for s in subs if self._is_donor(s))
                ewgs = sum(1 for s in subs if self._is_ewg(s))
                rk, over = "aromatic", None
                if donors > 0 and ewgs == 0:
                    rk, over = "aromatic_activated", _o(ox=("High", "Electron-rich (donor-substituted) ring is readily hydroxylated / oxidised by radicals."))
                elif ewgs > 0 and donors == 0:
                    rk, over = "aromatic_deactivated", _o(ox=("Low", "Electron-poor ring resists radical hydroxylation."))
                self._add_site(rk, "aromatic", "ring", tuple(ring), ring[0], over)
                continue
            if len(ring) == 6:
                if any(a.GetSymbol() == "N" for a in het):
                    self._add_site("azine", "hetero_azine", "ring", tuple(ring), ring[0], None)
                continue
            pyridine_n = any(a.GetSymbol() == "N" and a.GetTotalNumHs() == 0 and a.GetDegree() == 2 for a in het)
            pyrrole_n = any(a.GetSymbol() == "N" and (a.GetTotalNumHs() > 0 or a.GetDegree() == 3) for a in het)
            fused = any(len(r2) == 6 and r2 != ring and all(m.GetAtomWithIdx(i).GetIsAromatic() for i in r2) and set(r2) & set(ring)
                        for r2 in ri.AtomRings())
            kind, over = "azole", _o(acid=("Moderate", "Protonation of the pyridine-type N (imidazolium-type salts)."), ox=("Moderate", "Azoles are moderately oxidised."))
            if any(a.GetSymbol() == "O" for a in het):
                kind, over = "furan", _o(acid=("Critical", "Furans ring-open under aqueous acid to 1,4-dicarbonyls."), ox=("Critical", "Furans are oxidised (endoperoxide / ring opening) very readily."))
            elif any(a.GetSymbol() == "S" for a in het):
                kind, over = "thiophene", _o(acid=("Low", "Thiophene is stable in acid."), ox=("Moderate", "S-oxidation to thiophene S-oxide needs strong oxidants; the ring is fairly stable."))
            elif pyrrole_n and not pyridine_n:
                if fused:
                    kind, over = "indole", _o(acid=("High", "Indoles dimerise / oligomerise under acid (C3 protonation)."), ox=("Critical", "Electron-rich indole C2=C3 is oxidised to oxindole / cleavage products."))
                else:
                    kind, over = "pyrrole", _o(acid=("High", "Pyrroles polymerise in acid (pyrrole red)."), ox=("Critical", "Pyrroles darken by air / peroxide oxidation."))
            self._add_site(kind, "hetero_rich", "ring", tuple(ring), ring[0], over)
            if pyridine_n:
                self._add_site("azole_n", "hetero_azine", "ring", tuple(ring), ring[0], None)

    def _add_site(self, rule: str, pkey: str, fam: str, atoms: Tuple[int, ...], key: int, over: Optional[dict]) -> None:
        eff = dict(PROFILES[pkey])
        if over:
            eff.update(over)
        self.sites.append(Site(rule, pkey, fam, atoms, key, eff))
        self.present.add(rule)

    def _run(self) -> None:
        m = self.mol
        claimed: Set[Tuple[str, int]] = set()
        fam_atoms: Dict[str, Set[int]] = {}
        for fam, key, q, pk, ov, op, ka, k2i in _COMPILED:
            try:
                matches = m.GetSubstructMatches(q, uniquify=True, maxMatches=2000)
            except Exception:
                continue
            for mt in matches:
                k = mt[ka]
                k2 = mt[k2i] if (op and op.get("claim2") and k2i >= 0) else -1
                if (fam, k) in claimed or (k2 >= 0 and (fam, k2) in claimed):
                    continue
                if op:
                    ni = op.get("not_in")
                    if ni and any(set(mt) & fam_atoms.get(f, set()) for f in ni):
                        continue
                    ck = op.get("check")
                    if ck and not ck(m, mt):
                        continue
                claimed.add((fam, k))
                if k2 >= 0:
                    claimed.add((fam, k2))
                fam_atoms.setdefault(fam, set()).update(mt)
                self._add_site(key, pk, fam, tuple(mt), k, ov)
        self._ring_sites()
        if not self.sites:
            self._add_site("hydrocarbon", "hydrocarbon", "fallback", (0,), 0, None)
        self._build_entries()
        self.classes = self._derive_classes()

    # ---- merge sites of the same display name (worst case per condition)
    def _build_entries(self) -> None:
        by: Dict[str, dict] = {}
        for s in self.sites:
            e = s.eff
            cur = by.get(e["name"])
            if cur is None:
                by[e["name"]] = dict(e, n=1, rules={s.rule}, pkeys={s.pkey})
                continue
            cur["n"] += 1
            cur["rules"].add(s.rule)
            cur["pkeys"].add(s.pkey)
            for c in CKEYS:
                if vrank(e[c][0]) > vrank(cur[c][0]):
                    cur[c] = e[c]
            if vrank(e["sec"][0]) > vrank(cur["sec"][0]):
                cur["sec"] = e["sec"]
        self.entries = list(by.values())

    def _derive_classes(self) -> Set[str]:
        c: Set[str] = set()
        h = self._h
        for s in self.sites:
            p, k = s.pkey, s.key
            if p == "acid":
                c.add("acid_weak")
            elif p == "carboxylate":
                c.add("base_weak")
            elif p in ("ester", "lactone"):
                c.add("ester")
            elif p == "ester_aryl":
                c |= {"ester", "acylating"}
            elif p in ("beta_lactam", "anhydride", "acyl_halide", "isocyanate", "isothiocyanate", "carbonate", "sulfonyl_halide", "thioester", "imide"):
                c.add("acylating")
            elif p in ("aldehyde", "ketone"):
                c.add("carbonyl")
            elif p == "hemiacetal":
                c |= {"reducing_sugar", "carbonyl", "reductant"}
            elif p == "michael":
                c.add("michael")
            elif p == "quinone":
                c |= {"michael", "oxidant"}
            elif p in ("amine", "arylamine"):
                c |= {"N_nuc"}
                if p == "amine":
                    c.add("base")
                if h(k) >= 2:
                    c.add("N_nuc_prim")
                if h(k) == 1:
                    c.add("N_nuc_sec")
            elif p in ("amine_tert", "amidine", "hetero_azine"):
                c.add("base")
            elif p == "ammonium":
                c.add("acid_weak")
            elif p in ("hydrazine", "hydroxylamine"):
                c |= {"N_nuc", "N_nuc_prim", "reductant"}
            elif p in ("alcohol", "alcohol_tert"):
                c.add("O_nuc")
            elif p == "phenol":
                c |= {"O_nuc", "acid_weak"}
            elif p == "enol":
                c |= {"O_nuc", "reductant"}
            elif p == "thiol":
                c |= {"S_nuc", "reductant", "acid_weak"}
            elif p in ("epoxide", "alkyl_halide", "sulfonate_ester"):
                c.add("alkylating")
            elif p == "peroxide":
                c |= {"oxidant", "peroxide"}
            elif p in ("dioxygen", "hypohalite", "nitrate", "n_oxide"):
                c.add("oxidant")
            elif p == "ether":
                if any(h(n.GetIdx()) > 0 for n in self.mol.GetAtomWithIdx(k).GetNeighbors()):
                    c.add("peroxide_former")
            elif p == "benzylic_ch":
                c.add("peroxide_former")
            elif p in ("strong_base", "carbonate_ion"):
                c |= {"base", "base_strong"}
            elif p == "water":
                c.add("water")
            elif p == "metal":
                c |= {"metal", "lewis_acid"}
            elif p == "metal_transition":
                c |= {"metal", "redox_metal", "lewis_acid"}
                a = self.mol.GetAtomWithIdx(k)
                if a.GetFormalCharge() >= 3 or (a.GetSymbol() == "Cu" and a.GetFormalCharge() >= 2):
                    c.add("oxidant")
            elif p == "metal_oxide":
                if self.mol.GetAtomWithIdx(k).GetSymbol() in ("Mn", "Cr", "Mo", "V", "W", "Ce"):
                    c.add("oxidant")
            elif p == "nitrite":
                c.add("nitrosating")
            elif p == "hydride":
                c.add("reductant")
            elif p == "sulfite":
                c |= {"reductant", "S_nuc"}
            elif p == "halide_ion":
                c.add("halide")
            elif p == "hydrogen_halide":
                c |= {"acid_strong", "halide"}
            elif p in ("sulfonic", "sulfate"):
                if any(n.GetSymbol() == "O" and n.GetTotalNumHs() > 0 for n in self.mol.GetAtomWithIdx(k).GetNeighbors()):
                    c.add("acid_strong")
            elif p == "phosphate":
                if any(n.GetSymbol() == "O" and n.GetTotalNumHs() > 0 for n in self.mol.GetAtomWithIdx(k).GetNeighbors()):
                    c.add("acid_weak")
            elif p in ("alkene", "alkyne"):
                c.add("ene")
            elif p == "nitrosamine":
                c.add("nitrosating")
        return c

    def vuln(self, rules: Iterable[str], ck: str) -> Optional[str]:
        """worst-case vulnerability under condition ck among the present rules (None if none present)"""
        best: Optional[str] = None
        for s in self.sites:
            if s.rule in rules:
                v = s.eff[ck][0]
                if best is None or vrank(v) > vrank(best):
                    best = v
        return best

    def sec_vuln(self, rules: Iterable[str]) -> Optional[str]:
        best: Optional[str] = None
        for s in self.sites:
            if s.rule in rules:
                v = s.eff["sec"][0]
                if best is None or vrank(v) > vrank(best):
                    best = v
        return best


CLASS_TEXT = {
    "N_nuc": "nucleophilic amine / hydrazine", "N_nuc_prim": "primary amine", "N_nuc_sec": "secondary amine",
    "O_nuc": "hydroxyl nucleophile (alcohol/phenol)", "S_nuc": "thiol / sulfur nucleophile", "carbonyl": "aldehyde / ketone",
    "reducing_sugar": "reducing sugar (masked aldehyde)", "acylating": "acylating agent", "ester": "ester",
    "acid_weak": "Bronsted acid (carboxylic acid / phenol)", "acid_strong": "strong acid", "base": "base",
    "base_strong": "strong base / alkaline microenvironment", "alkylating": "alkylating agent", "michael": "Michael acceptor",
    "oxidant": "oxidant", "peroxide_former": "peroxide-forming / autoxidisable group", "reductant": "reductant",
    "metal": "metal ion", "redox_metal": "redox-active metal", "lewis_acid": "Lewis acid", "nitrosating": "nitrosating agent",
    "halide": "halide", "water": "water", "ene": "C=C / C#C unsaturation", "peroxide": "peroxide",
    "base_weak": "weakly basic carboxylate salt",
}


def refine_secondary(an: Analysis, co_name: str, co_classes: Set[str]) -> List[dict]:
    """Re-grade the cross-interaction column against the ACTUAL co-reactant(s)."""
    out = []
    for e in an.entries:
        e2 = dict(e)
        hit = [c for c in e["partners"] if c in co_classes]
        if hit:
            txt = ", ".join(dict.fromkeys(CLASS_TEXT.get(h, h) for h in hit[:3]))
            e2["sec"] = (e["sec"][0], f"{co_name}: {txt}", e["sec"][2])
        else:
            e2["sec"] = (step_vuln(e["sec"][0], -2), f"{co_name} (no complementary reactive group)",
                         "No matching reactive partner class in the co-reactant; only non-covalent (hydrogen-bond, ionic, dipolar) association expected.")
        out.append(e2)
    return out


def group_dicts(entries: List[dict]) -> List[dict]:
    """Entries -> the dictionary format used by the Streamlit app."""
    out = []
    for e in entries:
        out.append({
            "name": clean_text(e["name"]), "category": clean_text(e["cat"]), "fragment": e["frag"],
            "reactive_site": clean_text(e["site"] + (f" (x{e['n']})" if e.get("n", 1) > 1 else "")),
            "acidic": (e["acid"][0], clean_text(e["acid"][1])), "basic": (e["base"][0], clean_text(e["base"][1])),
            "hydrolysis": (e["hyd"][0], clean_text(e["hyd"][1])), "photolytic": (e["photo"][0], clean_text(e["photo"][1])),
            "thermal": (e["therm"][0], clean_text(e["therm"][1])), "oxidative": (e["ox"][0], clean_text(e["ox"][1])),
            "cross_reaction": (e["sec"][0], clean_text(e["sec"][2] if not e["sec"][1] else f"{e['sec'][1]}: {e['sec'][2]}")),
        })
    return out


def identify_functional_groups(smiles: str) -> List[dict]:
    """Functional groups of a SMILES in the app's dictionary format (never guesses: [] if unparseable)."""
    m = parse_smiles(smiles)
    if m is None:
        return []
    return group_dicts(Analysis(m).entries)


# --------------------------------------------------------------------------------------
# Reaction runner
# --------------------------------------------------------------------------------------
_HAL = {"F", "Cl", "Br", "I"}


class Template:
    __slots__ = ("rules", "ck", "key", "label", "rxn", "mech", "dG", "mult", "vuln", "loses", "kek", "maxn", "post", "same", "smarts")

    def __init__(self, rules, ck, key, label, smarts, mech, dG, mult=1.0, vuln=None, loses="", kek=False, maxn=4, post=None, same=False):
        self.rules, self.ck, self.key, self.label, self.mech, self.dG = list(rules), ck, key, label, clean_text(mech), dG
        self.mult, self.vuln, self.loses, self.kek, self.maxn, self.post, self.same, self.smarts = mult, vuln, loses, kek, maxn, post, same, smarts
        self.rxn = AllChem.ReactionFromSmarts(smarts)
        if self.rxn is None or self.rxn.GetNumReactantTemplates() < 1:
            raise ValueError(f"bad reaction SMARTS for {key}: {smarts}")


STRESS: List[Template] = []


def T(rules, ck, key, label, smarts, mech, dG, **kw) -> None:
    STRESS.append(Template(rules, ck, key, label, smarts, mech, dG, **kw))


def H3(rules, key, label, smarts, tx, dG, mult=(1.0, 1.0, 1.0), **kw) -> None:
    """the same transformation under acid / base / neutral hydrolysis (own text, energy, factor)"""
    for ck, t, g, m in zip(("acid", "base", "hyd"), tx, dG, mult):
        T(rules, ck, key, label, smarts, t, g, mult=m, **kw)


def react(rxn, reactants: Sequence[Chem.Mol], loses: str = "", maxn: int = 4, kek: bool = False,
          post: Optional[Callable[[Chem.Mol], Optional[Chem.Mol]]] = None) -> List[Chem.Mol]:
    """run a reaction SMARTS on fragment(s); keep sanitised, element-conserving, distinct products"""
    rs = []
    for r in reactants:
        r2 = Chem.Mol(r)
        if kek:
            try:
                Chem.Kekulize(r2, clearAromaticFlags=True)
            except Exception:
                return []
        rs.append(r2)
    rc: Dict[str, int] = {}
    for r in reactants:
        for k, v in element_counts(r).items():
            rc[k] = rc.get(k, 0) + v
    try:
        sets = rxn.RunReactants(tuple(rs), 300)
    except Exception:
        return []
    out: Dict[str, Chem.Mol] = {}
    for ps in sets:
        mols = []
        for p in ps:
            try:
                Chem.SanitizeMol(p)
            except Exception:
                mols = []
                break
            mols.append(p)
        if not mols:
            continue
        prod = mols[0]
        for p in mols[1:]:
            prod = Chem.CombineMols(prod, p)
        if post is not None:
            prod = post(prod)
            if prod is None:
                continue
        pc = element_counts(prod)
        ok = True
        for e in set(pc) | set(rc):
            if e in ("H", "O"):
                continue
            p_, r_ = pc.get(e, 0), rc.get(e, 0)
            lost = (e in loses) or ("X" in loses and e in _HAL)
            if p_ > r_ or (p_ < r_ and not lost):
                ok = False
                break
        if not ok:
            continue
        try:
            k = canon(prod)
        except Exception:
            continue
        out.setdefault(k, prod)
    return [out[k] for k in sorted(out)[:maxn]]


def _enol_to_keto(prod: Chem.Mol) -> Optional[Chem.Mol]:
    """vinyl-alcohol fragments released by hydrolysis tautomerise to the carbonyl compound"""
    rx = AllChem.ReactionFromSmarts("[CX3;!a:1]=[CX3;!a:2][OX2H1:3]>>[C:1][C:2]=[O:3]")
    cur = prod
    for _ in range(3):
        ps = rx.RunReactants((cur,), 5)
        if not ps:
            break
        p = ps[0][0]
        try:
            Chem.SanitizeMol(p)
        except Exception:
            break
        cur = p
    return cur


ACID_RULES = ["acid", "acid_malonic", "acid_betaketo", "acid_alphaEWG", "acid_aryl_o", "acid_aryl_p", "acid_picolinic"]
ESTER_RULES = ["ester", "ester_tert", "ester_formate", "ester_alphaEWG"]
LACTONE_RULES = ["lactone", "lactone_4", "lactone_enol"]
AMIDE_RULES = ["amide", "amide_anilide", "amide_formamide", "amide_hydrazide", "amide_tert"]
HALIDE_ALKYL = ["alkyl_halide", "alkyl_tert", "alkyl_benzylic", "alkyl_gem", "alkyl_I", "alkyl_Br"]
SNAR_RULES = ["aryl_snar_p", "aryl_snar_o", "aryl_snar_n"]

# ======================================================================================
# HYDROLYSIS (acid / base / neutral)
# ======================================================================================
H3(ESTER_RULES, "acyl-O-cleavage", "ester hydrolysis",
   "[CX3:1](=[O:2])!@[OX2:3][CX4:4]>>([C:1](=[O:2])O.[O:3][#6:4])",
   ("Acid-catalysed acyl-oxygen cleavage (A_AC2): the protonated carbonyl is attacked by water and the ester splits into carboxylic acid + alcohol.",
    "Saponification (B_AC2): hydroxide adds to the carbonyl and expels the alkoxide; irreversible once the carboxylate forms.",
    "Water-mediated hydrolysis of the ester under humidity / neutral pH (slow; accelerated by heat and by liberated acid)."),
   (-2.0, -6.0, -1.5))
H3(["ester_aryl"], "acyl-O-cleavage-aryl", "aryl/vinyl ester hydrolysis",
   "[CX3:1](=[O:2])!@[OX2:3][c,$([CX3]=[CX3]):4]>>([C:1](=[O:2])O.[O:3][*:4])",
   ("Acid-catalysed acyl-oxygen cleavage (A_AC2) to the carboxylic acid + phenol / enol (which tautomerises to the carbonyl compound).",
    "Saponification: the phenoxide / enolate is an excellent leaving group, so hydroxide cleaves the activated ester very fast.",
    "Activated ester: measurable neutral hydrolysis and buffer catalysis under humidity."),
   (-3.0, -8.0, -3.5), post=_enol_to_keto)
H3(LACTONE_RULES, "lactone-opening", "lactone hydrolysis",
   "[CX3:1](=[O:2])@[OX2:3][#6:4]>>([C:1](=[O:2])O.[O:3][#6:4])",
   ("Acid-catalysed ring opening to the hydroxy-acid (equilibrium often favours the closed 5-/6-ring lactone).",
    "Alkaline ring opening to the hydroxy-carboxylate (irreversible after ionisation).",
    "Ring-chain equilibrium in water; strained (4-ring) lactones hydrolyse rapidly."),
   (0.5, -4.0, 0.2))
H3(["lactone_arom"], "lactone-opening-arom", "coumarin-type lactone ring opening",
   "[#6:1](=[O:2])@[#8:3][#6:4]>>([#6:1](=[O:2])O.[#8:3][#6:4])",
   ("Acid-catalysed ring opening of the aromatic lactone.", "Hydroxide opens the coumarin-type lactone to the hydroxycinnamate salt.", "Slow ring opening of the aromatic lactone in water."),
   (-3.0, -7.0, -2.0), kek=True)
H3(AMIDE_RULES, "amide-cleavage", "amide hydrolysis",
   "[CX3;!$(C(=O)([!#6])[!#6]):1](=[O:2])!@[NX3;!$(N(C=O)C=O):3]>>([C:1](=[O:2])O.[N:3])",
   ("Acid-catalysed hydrolysis (A_AC2): O-protonation, water attack and C-N cleavage give the carboxylic acid + ammonium.",
    "Hydroxide attacks the amide carbonyl; C-N cleavage gives the carboxylate + amine.",
    "Neutral hydrolysis of the amide is very slow (half-lives of years) at ambient pH."),
   (-2.5, -3.0, -0.5))
H3(["lactam"], "lactam-opening", "lactam ring hydrolysis",
   "[CX3;R;!r4:1](=[O:2])@[NX3;!$(N(C=O)C=O):3]>>([C:1](=[O:2])O.[N:3])",
   ("Acid-catalysed ring opening to the amino-acid (6- and 7-membered lactams open faster than 5-membered ones).",
    "Hydroxide attacks the lactam carbonyl; ring opening to the amino-carboxylate (slow for 5-membered lactams).",
    "Neutral hydrolysis of the lactam ring is very slow at ambient pH."),
   (-2.0, -3.0, -0.8))
H3(["beta_lactam"], "beta-lactam-opening", "strained beta-lactam ring hydrolysis",
   "[CX3;r4:1](=[O:2])@[NX3;r4:3]>>([C:1](=[O:2])O.[N:3])",
   ("Protonation then rapid nucleophilic ring opening; ring strain (~26 kcal/mol) makes the 4-membered beta-lactam open fast.",
    "Hydroxide attacks the strained lactam carbonyl; irreversible ring cleavage.", "Spontaneous neutral solvolysis driven by ring strain."),
   (-8.0, -9.0, -7.5))
H3(["imide"], "imide-opening", "imide hydrolysis (amic acid formation)",
   "[CX3:1](=[O:2])[NX3:3][CX3:4]=[O:5]>>([C:1](=[O:2])O.[N:3][C:4]=[O:5])",
   ("Acid-catalysed C-N cleavage of the imide to the amic acid.", "Hydroxide opens the imide (strongly electrophilic carbonyls) to the amic acid within minutes to hours.",
    "Slow neutral ring opening of the imide to the amic acid on moisture exposure."),
   (-2.0, -6.0, -2.5))
H3(["urea"], "urea-cleavage", "urea hydrolysis (via carbamic acid, CO2 loss)",
   "[NX3:3][CX3:1](=[O:2])[NX3:4]>>([N:3].[N:4])",
   ("Acid hydrolysis of the urea carbonyl gives a carbamic acid that loses CO2 to two amines.", "Base hydrolysis (slow) gives the carbamate, which decarboxylates to the amines.",
    "Very slow neutral hydrolysis of the urea to two amines + CO2."),
   (-3.0, -2.0, -1.0), loses="C")
H3(["carbamate", "carbamate_tert"], "carbamate-cleavage", "carbamate hydrolysis (via carbamic acid, CO2 loss)",
   "[NX3:3][CX3:1](=[O:2])[OX2:4][#6:5]>>([N:3].[O:4][#6:5])",
   ("Acid hydrolysis of the carbamate gives the alcohol + carbamic acid, which loses CO2 to the free amine.",
    "Base hydrolysis (BAc2; E1cB via isocyanate for N-H aryl carbamates) releases the alcohol and, after CO2 loss, the amine.",
    "Slow neutral hydrolysis of the carbamate to amine + alcohol + CO2."),
   (-3.0, -3.5, -1.5), loses="C")
T(["carbamate_tert"], "acid", "boc-cleavage", "tert-butyl carbamate (Boc-type) cleavage",
  "[NX3:3][CX3:1](=[O:2])[OX2:4][CX4;H0:5]([CH3:6])([#6:7])[#6:8]>>([N:3].[C:5](=[CH2:6])([#6:7])[#6:8])",
  "A_AL1 / E1 acidolysis: loss of the tert-alkyl cation as an alkene and CO2 releases the free amine.", -4.5, loses="C")
T(["carbamate_tert"], "therm", "boc-cleavage", "tert-butyl carbamate thermolysis",
  "[NX3:3][CX3:1](=[O:2])[OX2:4][CX4;H0:5]([CH3:6])([#6:7])[#6:8]>>([N:3].[C:5](=[CH2:6])([#6:7])[#6:8])",
  "Thermolysis (above 150 C) expels the alkene and CO2 to give the free amine.", -4.0, loses="C")
H3(["carbonate"], "carbonate-cleavage", "carbonate hydrolysis (CO2 loss)",
   "[#6:5][OX2:3][CX3:1](=[O:2])[OX2:4][#6:6]>>([O:3][#6:5].[O:4][#6:6])",
   ("Acid-catalysed hydrolysis to two alcohols + CO2.", "Hydroxide attacks the carbonate carbonyl; two alkoxides + carbonate result.", "Slow neutral hydrolysis to two alcohols + CO2."),
   (-5.5, -6.5, -4.0), loses="C")
H3(["anhydride"], "anhydride-hydrolysis", "anhydride hydrolysis",
   "[CX3:1](=[O:2])[OX2:3][CX3:4]=[O:5]>>([C:1](=[O:2])O.[O:3][C:4]=[O:5])",
   ("Acid-catalysed hydrolysis of the anhydride to two carboxylic acids.", "Instant saponification to two carboxylates.", "Spontaneous hydrolysis by ambient moisture to two carboxylic acids."),
   (-8.5, -10.0, -9.0))
H3(["acyl_halide"], "acyl-halide-hydrolysis", "acyl halide hydrolysis",
   "[CX3:1](=[O:2])[F,Cl,Br,I:3]>>([C:1](=[O:2])O.[*:3])",
   ("Immediate hydrolysis to the carboxylic acid + HX.", "Instant saponification.", "Reacts violently with atmospheric moisture to the carboxylic acid + HX."),
   (-10.0, -11.0, -10.0))
H3(["thioester"], "thioester-hydrolysis", "thioester hydrolysis",
   "[CX3:1](=[O:2])[SX2:3][#6:4]>>([C:1](=[O:2])O.[S:3][#6:4])",
   ("Acid hydrolysis of the thioester to carboxylic acid + thiol.", "Base hydrolysis / thiolate exchange releases the carboxylate + thiol.", "Slow neutral hydrolysis to carboxylic acid + thiol."),
   (-2.0, -4.5, -2.5))
H3(["sulfonate_ester", "dialkyl_sulfate", "sulfate_ester"], "sulfonate-ester-hydrolysis", "sulfonate/sulfate ester hydrolysis (alkylating agent destroyed)",
   "[SX4:1](=[O:2])(=[O:3])[OX2:4][CX4:5]>>([S:1](=[O:2])(=[O:3])[O:4].[C:5]O)",
   ("Acid hydrolysis to the sulfonic acid + alcohol.", "Hydroxide attacks carbon (SN2) or sulfur, releasing the sulfonate + alcohol.",
    "Solvolysis in water / alcohol: the alkyl sulfonate alkylates water to give the alcohol + sulfonic acid."),
   (-5.0, -7.0, -5.0))
H3(["sulfonate_ester"], "aryl-sulfonate-hydrolysis", "aryl sulfonate ester hydrolysis",
   "[SX4:1](=[O:2])(=[O:3])[OX2:4][c:5]>>([S:1](=[O:2])(=[O:3])O.[O:4][c:5])",
   ("Acid hydrolysis of the aryl sulfonate to the sulfonic acid + phenol.", "Hydroxide attacks sulfur; the phenoxide leaves.", "Slow neutral hydrolysis of the aryl sulfonate."),
   (-5.0, -7.0, -5.0))
H3(["sulfonyl_halide"], "sulfonyl-halide-hydrolysis", "sulfonyl halide hydrolysis",
   "[SX4:1](=[O:2])(=[O:3])[F,Cl,Br,I:4]>>([S:1](=[O:2])(=[O:3])O.[*:4])",
   ("Acid-mediated hydrolysis to the sulfonic acid + HX.", "Rapid base hydrolysis to the sulfonate.", "Moisture hydrolyses the sulfonyl halide to the sulfonic acid + HX."),
   (-6.0, -8.0, -7.0))
T(["acetal", "hemiacetal"], "acid", "glycoside-hydrolysis", "glycoside / cyclic acetal cleavage (aglycone released)",
  "[CX4:1](@[OX2:2])!@[OX2:3][#6:4]>>([C:1]([O:2])O.[O:3][#6:4])",
  "Protonation of the exocyclic acetal oxygen, loss of the aglycone alcohol to the ring oxocarbenium ion and water capture (A1) give the hemiacetal + alcohol.", -4.0)
T(["acetal"], "acid", "acetal-hydrolysis", "acetal / ketal hydrolysis to the carbonyl compound + alcohols",
  "[CX4:1](!@[OX2:2][#6:5])!@[OX2:3][#6:4]>>([C:1]=O.[O:2][#6:5].[O:3][#6:4])",
  "Protonation of an acetal oxygen, loss of alcohol to the oxocarbenium ion, then water capture (A1 mechanism); very fast even at mild acidity.", -1.5)
T(["acetal"], "acid", "dioxolane-hydrolysis", "cyclic acetal (dioxolane / dioxane) hydrolysis to the carbonyl compound + diol",
  "[CX4:1](@[OX2:2][#6:5])@[OX2:3][#6:4]>>([C:1]=O.[O:2][#6:5].[O:3][#6:4])",
  "Acid-catalysed ring opening of the cyclic acetal via the oxocarbenium ion gives the carbonyl compound + diol.", -1.5)
T(["acetal"], "hyd", "acetal-hydrolysis", "acetal / ketal hydrolysis to the carbonyl compound + alcohols",
  "[CX4:1](!@[OX2:2][#6:5])!@[OX2:3][#6:4]>>([C:1]=O.[O:2][#6:5].[O:3][#6:4])",
  "Neutral hydrolysis is slow, but trace acid autocatalysis in the presence of moisture can cleave the acetal.", -1.5, mult=0.5)
H3(["imine", "oxime", "hydrazone"], "c=n-hydrolysis", "C=N hydrolysis to the carbonyl compound",
   "[CX3;!a:1]=[NX2;!a:2]>>([C:1]=O.[N:2])",
   ("Protonation gives an iminium that is attacked by water; the carbinolamine collapses to the carbonyl compound + amine / hydroxylamine / hydrazine.",
    "Base-catalysed hydrolysis of the C=N bond to the carbonyl compound + amine.", "Reversible hydrolysis of the C=N bond in water (equilibrium shifts with dilution)."),
   (-2.0, -1.0, -1.0))
H3(["amidine"], "amidine-hydrolysis", "amidine / guanidine hydrolysis to the amide / urea + amine",
   "[CX3:1](=[NX2;!$(N-[N,O]);!a:2])-[NX3:3]>>([C:1](=O)[N:3].[N:2])",
   ("Protonated amidinium is attacked by water; the intermediate expels ammonia / amine to give the amide (urea for guanidines).",
    "Hydroxide adds to the amidine carbon and expels the amine to give the amide / urea (needs heat).", "Slow neutral hydrolysis of the amidine / guanidine to the amide / urea + amine."),
   (-2.0, -2.5, -1.5))
T(["nitrile"], "acid", "nitrile-hydration", "nitrile hydration to the primary amide", "[CX2:1]#[NX1:2]>>[C:1](=O)[N:2]",
  "Acid-catalysed hydration of the C#N to the primary amide; further hydrolysis to the acid + NH4+ needs forcing conditions.", -4.0)
T(["nitrile"], "base", "nitrile-hydration", "nitrile hydration to the primary amide", "[CX2:1]#[NX1:2]>>[C:1](=O)[N:2]",
  "Hydroxide adds to the nitrile carbon; tautomerisation gives the primary amide (acid + NH3 on forcing).", -4.5)
H3(["epoxide"], "epoxide-opening", "epoxide ring opening to the 1,2-diol",
   "[C:1]1[O:2][C:3]1>>[C:1]([O:2])[C:3]O",
   ("Protonated epoxide is opened by water at the more substituted carbon (SN1-like) to the vicinal diol.", "Hydroxide opens the strained ring by SN2 at the less hindered carbon to the diol.",
    "Slow neutral hydrolysis of the strained oxirane to the diol."),
   (-9.0, -8.0, -7.0))
_HSUB = ("SN2 (primary/secondary) or SN1 (tertiary/benzylic/allylic) displacement of the halide by hydroxide gives the alcohol.",)
T(HALIDE_ALKYL, "base", "halide-substitution", "alkyl halide substitution by hydroxide", "[CX4:2][F,Cl,Br,I]>>[C:2]O",
  _HSUB[0], -4.0, loses="X", mult=0.85)
T(HALIDE_ALKYL, "hyd", "halide-substitution", "alkyl halide solvolysis", "[CX4:2][F,Cl,Br,I]>>[C:2]O",
  "Solvolysis: ionisation (SN1) of tertiary / benzylic / allylic halides, or slow SN2 by water, gives the alcohol + HX.", -4.0, loses="X", mult=0.8)
T(["alkyl_tert", "alkyl_benzylic"], "acid", "halide-substitution", "alkyl halide solvolysis", "[CX4:2][F,Cl,Br,I]>>[C:2]O",
  "Acid-assisted SN1 solvolysis of the labile halide to the alcohol + HX.", -4.0, loses="X", mult=0.8)
T(HALIDE_ALKYL, "base", "dehydrohalogenation", "dehydrohalogenation (E2/E1) to the alkene", "[CX4:1]([Cl,Br,I])[CX4;!H0:2]>>[C:1]=[C:2]",
  "Base-induced elimination of HX (Zaitsev alkene); dominant for tertiary halides, minor for primary.", -3.0, loses="X", mult=0.6, maxn=2)
T(SNAR_RULES, "base", "snar-hydroxide", "SNAr displacement of the halide by hydroxide (phenol formation)", "[F,Cl,Br,I][c:1]>>[c:1]O",
  "Hydroxide adds to the ring carbon activated by an ortho/para electron-withdrawing group (Meisenheimer complex) and expels halide.", -6.0, loses="X")
T(SNAR_RULES, "hyd", "snar-hydroxide", "SNAr displacement of the halide by water", "[F,Cl,Br,I][c:1]>>[c:1]O",
  "Activated aryl halides are slowly hydrolysed to the phenol by water / weak nucleophiles.", -5.0, loses="X", mult=0.5)
T(["alkene"], "acid", "alkene-hydration", "Markovnikov hydration of the C=C bond",
  "[CX3;!$(C-C=O);!$(C-C#N);!$(C-[N+](=O)[O-]);!$(C-[O,N]);H0,H1:1]=[CX3;!$(C-C=O);!$(C-C#N);!$(C-[O,N]);H1,H2:2]>>[C:1](O)[C:2]",
  "Protonation forms the more stable carbocation, trapped by water (Markovnikov hydration to the alcohol).", -1.8, maxn=2)
T(["alkene"], "acid", "enol-ether-hydrolysis", "enol ether / enamine / enamide hydrolysis to the carbonyl compound",
  "[CX3;!a:1](=[CX3;!a:2])!@[OX2,NX3;!$([OX2H1]):3]>>([C:1](=O)[C:2].[*:3])",
  "C-protonation of the electron-rich alkene gives an oxocarbenium / iminium ion that water converts to the carbonyl compound + alcohol / amine.", -2.5, vuln="High")
T(["alkene"], "hyd", "enol-ether-hydrolysis", "enol ether / enamine hydrolysis",
  "[CX3;!a:1](=[CX3;!a:2])!@[OX2,NX3;!$([OX2H1]):3]>>([C:1](=O)[C:2].[*:3])",
  "Moisture-driven hydrolysis of the heteroatom-substituted alkene to the carbonyl compound.", -1.0, vuln="Moderate")
H3(["michael"], "conjugate-hydration", "conjugate (oxa-Michael) hydration to the beta-hydroxy carbonyl",
   "[CX3;!a:1](=[O:2])[CX3;!a:3]=[CX3;!a:4]>>[C:1](=[O:2])[C:3][C:4]O",
   ("Acid-catalysed conjugate addition of water to the beta-carbon (reversible).", "Hydroxide adds to the beta-carbon (oxa-Michael); reversible retro-Michael / retro-aldol follows.",
    "Slow reversible conjugate hydration in water."),
   (0.8, 0.5, 1.2), mult=(0.5, 0.5, 0.5))
H3(["isocyanate"], "isocyanate-hydrolysis", "isocyanate hydrolysis to the amine (CO2 loss)", "[NX2:1]=[CX2:2]=[O:3]>>[N:1]",
   ("Water adds to N=C=O to give the carbamic acid, which loses CO2 to the amine.", "Hydroxide adds to the isocyanate carbon; the carbamate decarboxylates to the amine.",
    "Reacts with moisture: carbamic acid, then amine + CO2 (the amine can add to remaining isocyanate to form ureas)."),
   (-8.0, -9.0, -8.0), loses="C")
T(["quat"], "base", "hofmann", "Hofmann elimination to the alkene + tertiary amine", "[NX4+:1][CX4:2][CX4;!H0:3]>>[N+0:1].[C:2]=[C:3]",
  "Hydroxide removes a beta-hydrogen (E2, Hofmann rule) and expels the neutral tertiary amine; needs heat.", -3.0, mult=0.5, maxn=2)
T(["quat"], "base", "quat-dealkylation", "SN2 dealkylation of the quaternary ammonium by hydroxide", "[NX4+:1][CH3:2]>>[N+0:1].[C:2]O",
  "Hydroxide displaces the neutral tertiary amine from a methyl group (SN2; high temperature).", -3.0, mult=0.4, maxn=1)
H3(["phosphate"], "phosphoanhydride-hydrolysis", "phosphoanhydride (P-O-P) hydrolysis",
   "[PX4:1](=[O:2])[OX2:3][PX4:4]>>([P:1](=[O:2])O.[O:3][P:4])",
   ("Protonation of the bridging oxygen and water attack at phosphorus cleave the P-O-P bond to two phosphates.", "Hydroxide attack at phosphorus (slowed by charge repulsion of the polyanion).",
    "Water attacks phosphorus of the pyrophosphate-type linkage (metal-ion and pH dependent)."),
   (-7.0, -7.0, -7.0), mult=(1.0, 0.8, 1.0))
H3(["phosphate"], "phosphate-ester-hydrolysis", "phosphate ester hydrolysis (P-O cleavage)",
   "[PX4:1](=[O:2])[OX2:3][#6:4]>>([P:1](=[O:2])O.[O:3][#6:4])",
   ("Acid-catalysed hydrolysis of the P-O-C ester to the phosphate diester / monoester + alcohol (slow).",
    "Hydroxide attacks phosphorus (SN2@P) and expels the alkoxide (triesters fast, diesters very slow).", "Slow neutral hydrolysis of the phosphate ester."),
   (-3.0, -3.5, -2.5))
T(["alkyne"], "acid", "alkyne-hydration", "alkyne hydration to the ketone / aldehyde (Markovnikov)", "[CX2;H0:1]#[CX2;!H0:2]>>[C:1](=O)[C:2]",
  "Acid (usually Hg2+/Au+ assisted) hydration via the enol gives the methyl ketone (aldehyde for ethyne).", -5.0, mult=0.6, maxn=2)
T(["furan"], "acid", "furan-ring-opening", "acid-catalysed furan ring opening to the 1,4-dicarbonyl",
  "[#8:1]1[#6:2]=[#6:3][#6:4]=[#6:5]1>>[O:1]=[C:2][C:3][C:4][C:5]=O",
  "C2-protonation of the furan gives an oxocarbenium ion; water addition and ring opening give the saturated 1,4-dicarbonyl (reverse Paal-Knorr).", -3.0, kek=True, maxn=2)

# ======================================================================================
# OXIDATION
# ======================================================================================
T(["thioether"], "ox", "s-oxidation", "sulfide to sulfoxide (S-oxidation)", "[SX2:1]([#6:2])[#6:3]>>[S:1](=O)([#6:2])[#6:3]",
  "Electrophilic oxygen transfer from peroxide / O2 to the nucleophilic sulfur lone pair gives the sulfoxide (-S(=O)-).", -7.5)
T(["thioether"], "ox", "s-oxidation-sulfone", "sulfide to sulfone (over-oxidation)", "[SX2:1]([#6:2])[#6:3]>>[S:1](=O)(=O)([#6:2])[#6:3]",
  "Excess oxidant converts the sulfoxide onward to the sulfone (-SO2-).", -6.5, mult=0.5)
T(["sulfoxide"], "ox", "sulfoxide-to-sulfone", "sulfoxide to sulfone", "[SX3:1](=[O:4])([#6:2])[#6:3]>>[S:1](=[O:4])(=O)([#6:2])[#6:3]",
  "Further oxygen transfer to the sulfoxide sulfur gives the sulfone.", -6.0)
T(["thiol", "thiolate"], "ox", "disulfide-formation", "oxidative dimerisation to the disulfide", "[SX2H1:1][#6:2].[SX2H1:3][#6:4]>>[#6:2][S:1][S:3][#6:4]",
  "Two thiols are oxidised (thiolate / thiyl radical coupling; O2, peroxides, trace metals) to the S-S disulfide.", -5.0, same=True, maxn=1)
T(["thiol"], "ox", "sulfonic-acid", "over-oxidation to the sulfonic acid", "[SX2H1:1][#6:2]>>[S:1](=O)(=O)(O)[#6:2]",
  "Strong oxidants take the thiol through sulfenic and sulfinic acids to the sulfonic acid.", -4.5, mult=0.45)
T(["amine_tert", "arylamine"], "ox", "n-oxide", "tertiary amine N-oxidation",
  "[NX3;H0;!$(N-C=[O,S,N]);!$(N-[SX4,PX4]=O);!$(N-[N,O]);!$(N=*);!$(N(-a)-a):1]>>[N+:1][O-]",
  "Electrophilic oxygen transfer from peroxides / O2 to the nucleophilic nitrogen lone pair gives the polar N-oxide (R3N+-O-).", -4.5, maxn=2)
T(["amine_tert", "amine", "arylamine"], "ox", "n-dealkylation", "oxidative N-dealkylation (carbinolamine collapse)",
  "[NX3;H0,H1;!$(N-C=[O,S,N]);!$(N-[SX4,PX4]=O);!$(N-[N,O]);!$(N=*):1][CH3:2]>>([N:1].[C:2]=O)",
  "Hydrogen abstraction at the N-alpha C-H, oxygen rebound and collapse of the carbinolamine release the dealkylated amine + formaldehyde.", -5.0, mult=0.9, maxn=1)
T(["amine_tert", "amine", "arylamine"], "ox", "n-dealkylation", "oxidative N-dealkylation (carbinolamine collapse)",
  "[NX3;H0,H1;!$(N-C=[O,S,N]);!$(N-[SX4,PX4]=O);!$(N-[N,O]);!$(N=*):1][CX4;H2:2]>>([N:1].[C:2]=O)",
  "Hydrogen abstraction at the N-alpha C-H, oxygen rebound and collapse of the carbinolamine release the dealkylated amine + aldehyde.", -5.0, mult=0.7, maxn=2)
T(["amine"], "ox", "oxidative-deamination", "oxidative deamination to the carbonyl compound + NH3", "[NX3;H2;!$(N-a):1][CX4;!H0:2]>>([N:1].[C:2]=O)",
  "Hydrogen abstraction at the alpha C-H, oxygen rebound and collapse of the carbinolamine release the carbonyl compound + ammonia.", -5.0, mult=0.5, maxn=1)
T(["amine", "arylamine"], "ox", "n-hydroxylation", "N-hydroxylation to the hydroxylamine",
  "[NX3;H2,H1;!$(N-C=[O,S,N]);!$(N-[SX4,PX4]=O);!$(N-[N,O]):1]>>[N:1]O",
  "Two-electron N-oxidation of the amine N-H gives the hydroxylamine (further oxidation gives nitroso / nitrone species).", -2.5, maxn=2)
T(["arylamine"], "ox", "nitroso", "aromatic amine oxidation to the nitroso compound", "[NX3;H2:1][c;$(c1ccccc1):2]>>[N:1](=O)[c:2]",
  "Stepwise oxidation (hydroxylamine, then nitroso) of the aniline nitrogen gives Ar-N=O.", -3.0, mult=0.8, maxn=1)
T(["arylamine"], "ox", "azo-dimer", "oxidative coupling to the azo dimer (Ar-N=N-Ar)", "[NX3;H2:1][c;$(c1ccccc1):2].[NX3;H2:3][c;$(c1ccccc1):4]>>[c:2][N:1]=[N:3][c:4]",
  "One-electron oxidation gives anilino radicals that couple N-N to hydrazo and then azo compounds (coloured impurities).", -6.0, mult=0.5, same=True, maxn=1)
T(["azine", "azole_n"], "ox", "aza-n-oxide", "heteroaromatic N-oxidation", "[nX2;r6:1]>>[n+:1][O-]",
  "Peroxides / peracids oxidise the pyridine-type ring nitrogen to the N-oxide.", -3.0, mult=0.9, maxn=1)
T(["hemiacetal"], "ox", "aldose-oxidation", "oxidation of the reducing end to the lactone / aldonic acid", "[CX4;H1:1]([OX2H1:2])[OX2:3][#6:4]>>[C:1](=[O:2])[O:3][#6:4]",
  "The masked aldehyde (hemiacetal) is oxidised at the anomeric carbon to the lactone (aldonolactone) / aldonic acid.", -8.0, maxn=1)
T(["aldehyde"], "ox", "aldehyde-autoxidation", "aldehyde autoxidation to the carboxylic acid", "[CX3;H1,H2:1]=[O:2]>>[C:1](=[O:2])O",
  "Radical chain autoxidation: the acyl radical adds O2 to give a peracid, which oxidises a second aldehyde; net formation of the carboxylic acid.", -8.5)
T(["alcohol", "alcohol_benzylic"], "ox", "alcohol-oxidation", "alcohol oxidation to the aldehyde", "[CX4;H2,H3;!$(C(O)O):1][OX2H:2]>>[C:1]=[O:2]",
  "Hydride / hydrogen-atom abstraction from the carbinol C-H (peroxide, O2, metal catalysis) gives the carbonyl compound.", -3.5, maxn=2)
T(["alcohol", "alcohol_benzylic"], "ox", "alcohol-oxidation", "alcohol oxidation to the ketone", "[CX4;H1;!$(C(O)O):1][OX2H:2]>>[C:1]=[O:2]",
  "Hydride / hydrogen-atom abstraction from the carbinol C-H (peroxide, O2, metal catalysis) gives the carbonyl compound.", -4.0, maxn=2)
T(["phenol", "phenol_ewg_o", "phenol_ewg_p"], "ox", "quinone", "oxidation to the para-quinone / quinone-imine",
  "[OX2H:1][c:2]1[c:3][c:4][c:5]([OX2H,$([NX3;H1,H2;!$(N-S=O)]):6])[c:7][c:8]1>>[O:1]=[C:2]1[C:3]=[C:4][C:5](=[*:6])[C:7]=[C:8]1",
  "Two-electron, two-proton oxidation (or SET then disproportionation) of the hydroquinone / aminophenol-type ring gives the quinone / quinone-imine.", -2.5, maxn=1)
T(["phenol", "phenol_ewg_o", "phenol_ewg_p"], "ox", "quinone", "oxidation to the ortho-quinone / quinone-imine",
  "[OX2H:1][c:2]1[c:3]([OX2H,$([NX3;H1,H2;!$(N-S=O)]):4])[c:5][c:6][c:7][c:8]1>>[O:1]=[C:2]1[C:3](=[*:4])[C:5]=[C:6][C:7]=[C:8]1",
  "Oxidation of the catechol / aminophenol ring to the ortho-quinone (reactive electrophile).", -2.5, maxn=1)
T(["phenol", "phenol_ewg_o", "phenol_ewg_p"], "ox", "phenol-coupling", "oxidative ortho-ortho C-C coupling to the biaryl-diol",
  "[OX2H:1][c:2]:[cH:3].[OX2H:4][c:5]:[cH:6]>>[OH:1][c:2]:[c:3]-[c:6]:[c:5][OH:4]",
  "Phenoxyl radicals (from SET / H-abstraction) couple at the ortho carbons; tautomerisation restores aromaticity to give the 2,2'-biaryl diol.", -3.0, mult=0.5, same=True, maxn=1)
T(["enol"], "ox", "enediol-oxidation", "ene-diol oxidation to the 1,2-dicarbonyl (dehydro form)", "[OX2H:1][CX3:2]=[CX3:3][OX2H:4]>>[O:1]=[C:2][C:3]=[O:4]",
  "Two-electron, two-proton oxidation of the ene-diol (via the radical anion; metal ions catalyse) gives the 1,2-dicarbonyl.", -4.0, maxn=1)
T(["alkene"], "ox", "epoxidation", "alkene epoxidation", "[CX3;!a:1]=[CX3;!a:2]>>[C:1]1O[C:2]1",
  "Electrophilic oxygen transfer from peroxide / peracid (or radical addition-cyclisation) converts the C=C bond to the epoxide.", -4.5, maxn=2)
T(["michael"], "ox", "epoxidation", "nucleophilic epoxidation of the enone", "[CX3;!a:3](=[O:5])[CX3;!a:1]=[CX3;!a:2]>>[C:3](=[O:5])[C:1]1O[C:2]1",
  "Nucleophilic (Weitz-Scheffer) epoxidation of the electron-poor C=C by hydroperoxide anions.", -4.5, mult=0.5, maxn=1)
T(["ether", "ether_cyclic"], "ox", "ether-hydroperoxide", "ether autoxidation to the alpha-hydroperoxide",
  "[CX4;!H0;!$(C(O)O):1][OX2:2][CX4:3]>>[C:1](OO)[O:2][C:3]",
  "Radical autoxidation: H-abstraction at the alpha C-H next to oxygen, O2 capture and chain transfer give the alpha-hydroperoxide (peroxide formation on storage).", -3.0, maxn=2)
T(["aryl_ether"], "ox", "o-dealkylation", "oxidative O-dealkylation to the phenol + aldehyde", "[c:1][OX2:2][CX4;!H0:3]>>([c:1][O:2].[C:3]=O)",
  "Alpha C-H hydroxylation of the O-alkyl group gives a hemiacetal that collapses to the phenol + aldehyde.", -4.0, mult=0.7, maxn=2)
T(["benzylic_ch"], "ox", "benzylic-hydroperoxide", "benzylic/allylic autoxidation to the hydroperoxide",
  "[CX4;H1;!$(C-[!#6;!#1]):1]-[c,$([CX3]=[CX3]):2]>>[C:1](OO)-[*:2]",
  "Radical chain: H-abstraction at the weak benzylic / allylic C-H and O2 capture give the hydroperoxide (e.g. cumene hydroperoxide).", -3.0, maxn=2)
T(["benzylic_ch"], "ox", "benzylic-hydroxylation", "benzylic/allylic hydroxylation to the alcohol",
  "[CX4;H2,H3;!$(C-[!#6;!#1]):1]-[c,$([CX3]=[CX3]):2]>>[C:1](O)-[*:2]",
  "Autoxidation via the hydroperoxide followed by reduction / decomposition gives the benzylic / allylic alcohol.", -4.0, maxn=2)
T(["benzylic_ch"], "ox", "benzylic-oxidation", "benzylic oxidation to the ketone",
  "[CX4;H2;!$(C-[!#6;!#1]):1](-[#6:3])-[c,$([CX3]=[CX3]):2]>>[C:1](=O)([#6:3])-[*:2]",
  "Further oxidation of the secondary hydroperoxide / alcohol gives the ketone.", -4.5, mult=0.6, maxn=2)
T(["hydrazine"], "ox", "hydrazine-oxidation", "hydrazine oxidation to the diazene", "[NX3;!H0:1][NX3;!H0:2]>>[N:1]=[N:2]",
  "Two-electron oxidation of the N-N unit gives the diazene (N=N), which can lose N2.", -5.0, maxn=1)
T(["hydroxylamine"], "ox", "hydroxylamine-oxidation", "hydroxylamine oxidation to the nitroso compound", "[NX3;!H0:1][OX2H1:2]>>[N:1]=[O:2]",
  "Oxidation of the N-hydroxy-amine removes two hydrogens to give the C-nitroso compound.", -3.0, maxn=1)
T(["nitroso"], "ox", "nitroso-to-nitro", "nitroso oxidation to the nitro compound", "[#6:2][NX2:1]=[O:3]>>[#6:2][N+:1](=[O:3])[O-]",
  "Oxidation of the nitroso nitrogen (peroxide / O2) gives the nitro group.", -8.0)
T(["thiocarbonyl"], "ox", "oxidative-desulfurisation", "oxidative desulfurisation (C=S to C=O)", "[CX3:1]=[SX1]>>[C:1]=O",
  "S-oxidation to the sulfine / S-oxide is followed by loss of sulfur, converting the thiocarbonyl to the carbonyl.", -5.0, loses="S")
T(["phosphine"], "ox", "p-oxidation", "phosphine oxidation to the phosphine oxide", "[PX3:1]>>[P:1]=O",
  "Air / peroxide oxidises the P(III) lone pair to the strong P=O bond.", -9.0)
T(["ketone"], "ox", "baeyer-villiger", "Baeyer-Villiger oxidation to the ester / lactone",
  "[#6:2][CX3:1](=[O:3])[#6;$([CX4;H0,H1]),$(c):4]>>[#6:2][C:1](=[O:3])O[#6:4]",
  "Peroxide / peracid adds to the carbonyl (Criegee intermediate); the higher-aptitude group migrates to oxygen, inserting O between carbonyl and that group.", -7.0, mult=0.6, maxn=2)

T(["ketone"], "ox", "baeyer-villiger", "Baeyer-Villiger oxidation to the ester / lactone",
  "[#6:2][CX3:1](=[O:3])[#6;!$([CX4;H0,H1]);!$(c):4]>>[#6:2][C:1](=[O:3])O[#6:4]",
  "Peroxide / peracid adds to the carbonyl (Criegee intermediate); the group migrates to oxygen, inserting O between carbonyl and that group.", -7.0, mult=0.4, maxn=2)
T(["azole_n"], "ox", "aza-n-oxide", "heteroaromatic N-oxidation", "[nX2;r5:1]>>[n+:1][O-]",
  "Peroxides / peracids oxidise the pyridine-type ring nitrogen to the N-oxide.", -3.0, mult=0.5, maxn=1)
T(["alkene"], "acid", "alkene-hydration", "hydration of the symmetric C=C bond", "[CX3;H2:1]=[CX3;H2:2]>>[C:1](O)[C:2]",
  "Protonation forms the carbocation, trapped by water (hydration to the alcohol).", -1.8, maxn=1)
T(["alkyne"], "acid", "alkyne-hydration", "alkyne hydration to the aldehyde / ketone", "[CX2;H1:1]#[CX2;H1:2]>>[C:1](=O)[C:2]",
  "Acid (usually Hg2+/Au+ assisted) hydration via the enol gives the carbonyl compound.", -5.0, mult=0.6, maxn=1)

# ======================================================================================
# PHOTOLYSIS
# ======================================================================================
_ARYL_FRIES = ("UV homolysis of the aryl ester C(acyl)-O bond gives a radical pair in the solvent cage; recombination at the {} ring position gives the hydroxyaryl ketone.")
T(["ester_aryl"], "photo", "photo-fries-ortho", "photo-Fries rearrangement (ortho-hydroxyaryl ketone)", "[CX3:1](=[O:2])!@[OX2:3][c:4]:[cH:5]>>[C:1](=[O:2])[c:5]:[c:4][O:3]",
  _ARYL_FRIES.format("ortho"), 1.5, maxn=1)
T(["ester_aryl"], "photo", "photo-fries-para", "photo-Fries rearrangement (para-hydroxyaryl ketone)", "[CX3:1](=[O:2])!@[OX2:3][c:4]1[c:5][c:6][cH:7][c:8][c:9]1>>[C:1](=[O:2])[c:7]1[c:6][c:5][c:4]([O:3])[c:9][c:8]1",
  _ARYL_FRIES.format("para"), 1.5, maxn=1)
T(["amide_anilide"], "photo", "photo-fries-anilide", "photo-Fries rearrangement of the anilide (ortho-aminoaryl ketone)", "[CX3:1](=[O:2])!@[NX3:3][c:4]:[cH:5]>>[C:1](=[O:2])[c:5]:[c:4][N:3]",
  "UV homolysis of the anilide C(acyl)-N bond and radical-cage recombination at the ortho ring position gives the amino-aryl ketone.", 2.0, mult=0.8, maxn=1)
T(["amide_anilide"], "photo", "photo-fries-anilide-p", "photo-Fries rearrangement of the anilide (para-aminoaryl ketone)", "[CX3:1](=[O:2])!@[NX3:3][c:4]1[c:5][c:6][cH:7][c:8][c:9]1>>[C:1](=[O:2])[c:7]1[c:6][c:5][c:4]([N:3])[c:9][c:8]1",
  "UV homolysis of the anilide C(acyl)-N bond and radical-cage recombination at the para ring position gives the amino-aryl ketone.", 2.0, mult=0.8, maxn=1)
T(["ketone"], "photo", "norrish-II", "Norrish type II photocleavage (methyl ketone + alkene)", "[CX3:1](=[O:2])[CX4:3][CX4:4][CX4;!H0:5]>>([C:1](=[O:2])[C:3].[C:4]=[C:5])",
  "n to pi* excitation, intramolecular gamma-hydrogen abstraction (1,4-biradical) and beta-scission give a shorter ketone (via the enol) + an alkene.", 2.0, maxn=2)
T(["aldehyde"], "photo", "photodecarbonylation", "photodecarbonylation (Norrish type I) to the alkane / arene", "[#6:1][CX3;H1:2]=[O:3]>>[#6:1]",
  "n to pi* excitation cleaves the acyl C-C bond (Norrish I); loss of CO and H-atom transfer gives the decarbonylated hydrocarbon.", 1.0, mult=0.6, loses="C", maxn=1)
T(["acid"], "photo", "photodecarboxylation", "photodecarboxylation of the alpha-aryl acid", "[CX4:1]([c:5])[CX3](=O)[OX2H1]>>[C:1][c:5]",
  "UV excitation of the aryl chromophore drives electron transfer / decarboxylation, releasing CO2 and a benzylic radical that abstracts H.", 0.5, mult=0.7, loses="C", maxn=1)
T(["aryl_I", "aryl_Br", "aryl_F", "aryl_halide", "aryl_snar_p", "aryl_snar_o", "aryl_snar_n"], "photo", "photodehalogenation", "photodehalogenation (Ar-X to Ar-H)", "[c:1][F,Cl,Br,I]>>[c:1]",
  "UV homolysis of the C-X bond (weakest for I, then Br, Cl) gives an aryl radical that abstracts hydrogen from the medium.", 2.5, loses="X", maxn=1)
T(HALIDE_ALKYL, "photo", "photo-hydrodehalogenation", "photo-hydrodehalogenation (R-X to R-H)", "[CX4:1][Cl,Br,I]>>[C:1]",
  "UV homolysis of the C-X bond (weakest for I, then Br, Cl) gives an alkyl radical that abstracts hydrogen from the medium.", 2.0, loses="X", maxn=1)
T(["nitro"], "photo", "nitro-photoreduction", "nitroarene photoreduction to the nitroso compound", "[c:1][N+:2](=[O:3])[O-]>>[c:1][N+0:2]=[O:3]",
  "Excited nitroarene (n,pi* triplet) abstracts hydrogen / transfers an O atom, giving the nitroso compound (further reduction gives hydroxylamine / amine).", 3.5, maxn=1)
T(["alkene"], "photo", "[2+2]-dimer", "[2+2] photodimerisation to the cyclobutane", "[CX3;!a:1](=[CX3;!a:2])[c,$(C=O):5].[CX3;!a:3]=[CX3;!a:4]>>[C:1]1([*:5])[C:2][C:4][C:3]1",
  "Excited-state alkene adds to a ground-state alkene ([2+2] cycloaddition) giving the cyclobutane dimer.", 2.5, mult=0.7, same=True, maxn=2)
T(["michael"], "photo", "[2+2]-dimer", "[2+2] photodimerisation of the enone to the cyclobutane", "[CX3;!a:1](=[O:5])[CX3;!a:2]=[CX3;!a:3].[CX3;!a:4]=[CX3;!a:6]>>[C:1](=[O:5])[C:2]1[C:3][C:6][C:4]1",
  "Triplet enone adds to a ground-state alkene ([2+2] photocycloaddition) forming the cyclobutane dimer.", 2.5, mult=0.7, same=True, maxn=2)
T(["sulfonamide"], "photo", "so2-extrusion", "photo-desulfonylation (SO2 extrusion, Ar-N bond formation)", "[c:1][SX4](=O)(=O)[NX3:3]>>[c:1][N:3]",
  "UV cleavage of the S-N / S-C bonds with SO2 extrusion (Smiles-type rearrangement) joins the aryl carbon to nitrogen.", 2.0, mult=0.8, loses="S", maxn=1)
T(["azide", "azide2"], "photo", "azide-photolysis", "azide photolysis (nitrene) to the amine", "[#6:1][N:2]=[N+]=[N-]>>[#6:1][N:2]",
  "UV photolysis extrudes N2 to the nitrene, which abstracts hydrogen to give the primary amine (or inserts / rearranges).", -1.0, loses="N", maxn=1)
T(["nitrosamine"], "photo", "denitrosation", "photolytic N-N cleavage (denitrosation) to the secondary amine", "[NX3:1][NX2]=O>>[N:1]",
  "UV (230-350 nm) cleaves the N-N bond giving an aminyl radical + NO; H-abstraction gives the parent amine.", 2.0, loses="N", maxn=1)
T(["nitrosamine"], "acid", "denitrosation", "acid-mediated denitrosation to the secondary amine", "[NX3:1][NX2]=O>>[N:1]",
  "Protonated nitrosamine undergoes nucleophile-assisted (Br-, SCN-, thiourea) transnitrosation, releasing the amine + NO+.", 1.0, mult=0.7, loses="N", maxn=1)
T(["disulfide"], "photo", "disulfide-homolysis", "disulfide S-S photolysis to the thiols", "[SX2:1][SX2:2]>>([S:1].[S:2])",
  "UV homolysis of the S-S bond gives thiyl radicals that abstract hydrogen to give two thiols (or scramble / recombine).", 3.0, maxn=1)
T(["n_oxide"], "photo", "n-oxide-deoxygenation", "N-oxide photo-deoxygenation to the amine", "[#7+;!$([N+](=O)[O-]):1][O-]>>[#7+0:1]",
  "Photolysis of the N-O bond releases atomic oxygen and regenerates the tertiary amine / pyridine (competes with rearrangement).", 2.0, maxn=1)
T(["sulfoxide"], "photo", "sulfoxide-deoxygenation", "sulfoxide photo-deoxygenation to the sulfide", "[SX3:1](=O)([#6:2])[#6:3]>>[S:1]([#6:2])[#6:3]",
  "UV excitation cleaves the S=O bond (photo-deoxygenation) to give the sulfide.", 2.5, mult=0.7, maxn=1)
_PEROX_P = "UV homolysis of the weak O-O bond gives alkoxy / hydroxyl radicals that abstract hydrogen to give the alcohol (acid for peracids)."
T(["peroxide"], "photo", "peroxide-homolysis", "peroxide O-O homolysis (reduction to alcohol/acid)", "[OX2H1][OX2:2][#6:3]>>[O:2][#6:3]", _PEROX_P, -2.0, loses="", maxn=1)
T(["peroxide"], "photo", "peroxide-homolysis", "peroxide O-O homolysis (reduction to alcohols)", "[#6:3][OX2:1][OX2:2][#6:4]>>([#6:3][O:1].[#6:4][O:2])", _PEROX_P, -2.0, maxn=1)
T(["peroxide"], "therm", "peroxide-homolysis", "peroxide thermolysis (O-O homolysis)", "[OX2H1][OX2:2][#6:3]>>[O:2][#6:3]",
  "Thermal O-O homolysis initiates radical chains and gives the alcohol / acid after H-abstraction.", -2.0, maxn=1)
T(["peroxide"], "therm", "peroxide-homolysis", "peroxide thermolysis (O-O homolysis)", "[#6:3][OX2:1][OX2:2][#6:4]>>([#6:3][O:1].[#6:4][O:2])",
  "Thermal O-O homolysis initiates radical chains and gives the alcohols after H-abstraction.", -2.0, maxn=1)

# ======================================================================================
# THERMAL
# ======================================================================================
T(["acid_betaketo"], "therm", "decarboxylation", "decarboxylation (beta-keto acid)", "[CX3](=O)([OX2H1])[CX4:4][CX3:5]=O>>[C:4][C:5]=O",
  "Cyclic six-membered transition state (beta-keto acid) releases CO2 and gives the enol, which tautomerises to the ketone.", -4.0, loses="C", maxn=1)
T(["acid_malonic"], "therm", "decarboxylation", "decarboxylation (malonic-type acid)", "[CX3;$(C(=O)[OX2H1])]([OX2H1])(=O)[CX4:4][CX3:5](=[O:6])[OX2H1:7]>>[C:4][C:5](=[O:6])[O:7]",
  "Malonic-type acids lose CO2 through a cyclic transition state to give the mono-acid.", -3.5, loses="C", maxn=1)
T(["acid_alphaEWG"], "therm", "decarboxylation", "decarboxylation (alpha-EWG acid)", "[CX3](=O)([OX2H1])[CX4:4]>>[C:4]",
  "The electron-withdrawing alpha-substituent stabilises the carbanion-like transition state, releasing CO2.", -1.0, loses="C", maxn=1)
T(["acid_aryl_o", "acid_aryl_p"], "therm", "decarboxylation", "protodecarboxylation of the electron-rich aryl acid", "[CX3](=O)([OX2H1])[c:4]>>[c:4]",
  "Ipso-protonation of the electron-rich aryl acid (protodecarboxylation) releases CO2.", -1.0, loses="C", maxn=1)
T(["acid_picolinic"], "therm", "decarboxylation", "decarboxylation of the 2-pyridyl acid (Hammick-type)", "[CX3](=O)([OX2H1])[c:4]:n>>[c:4]:n",
  "The zwitterionic 2-pyridyl acid loses CO2 to give the ylide / pyridine.", -1.0, loses="C", maxn=1)
T(["ester_tert"], "acid", "tbu-cleavage", "tert-alkyl ester cleavage", "[CX3:1](=[O:2])!@[OX2:3][CX4;H0:4]([CH3:5])([#6:6])[#6:7]>>([C:1](=[O:2])[O:3].[C:4](=[CH2:5])([#6:6])[#6:7])",
  "A_AL1 / E1: protonation and loss of the tertiary carbocation as an alkene (e.g. isobutene) leaves the free carboxylic acid; no water is consumed.", -4.0)
T(["ester_tert"], "therm", "tbu-cleavage", "tert-alkyl ester thermolysis", "[CX3:1](=[O:2])!@[OX2:3][CX4;H0:4]([CH3:5])([#6:6])[#6:7]>>([C:1](=[O:2])[O:3].[C:4](=[CH2:5])([#6:6])[#6:7])",
  "Thermal E1 / syn-elimination expels the alkene and leaves the carboxylic acid.", -4.0)
T(["alcohol", "alcohol_tert", "alcohol_benzylic"], "therm", "aldol-dehydration", "dehydration of the beta-hydroxy carbonyl to the enone (E1cB)",
  "[OX2H:1][CX4:2][CX4;!H0:3][CX3:4]=[O:5]>>[C:2]=[C:3][C:4]=[O:5]",
  "Enolisation and E1cB loss of water from the aldol-type beta-hydroxy carbonyl gives the conjugated alpha,beta-unsaturated carbonyl.", -1.0, mult=0.9, maxn=1)
T(["alcohol_tert"], "therm", "dehydration", "acid- / heat-catalysed dehydration to the alkene", "[OX2H:1][CX4;H0:2][CX4;!H0:3]>>[C:2]=[C:3]",
  "E1 dehydration through the tertiary carbocation (trace acid, heat) gives the alkene + water.", 0.8, maxn=2)
T(["alcohol_benzylic"], "therm", "dehydration", "acid- / heat-catalysed dehydration to the alkene", "[OX2H:1][CX4;H1:2]([c:4])[CX4;!H0:3]>>[C:2]([c:4])=[C:3]",
  "E1 dehydration through the benzylic carbocation (trace acid, heat) gives the conjugated alkene + water.", 0.8, mult=0.7, maxn=2)
T(["sulfoxide"], "therm", "sulfoxide-syn-elimination", "sulfoxide pyrolysis (syn-elimination) to alkene + sulfenic acid", "[SX3:1](=[O:2])[CX4:3][CX4;!H0:4]>>([S:1][O:2].[C:3]=[C:4])",
  "Concerted five-membered cyclic (Ei) syn-elimination expels the sulfenic acid and forms the alkene.", -0.5, maxn=1)
T(["carbamate", "carbamate_tert"], "therm", "carbamate-to-isocyanate", "carbamate thermolysis to isocyanate", "[NX3;!H0:1][CX3:2](=[O:3])[OX2:4][#6:5]>>([N:1]=[C:2]=[O:3].[O:4][#6:5])",
  "Reversible thermal dissociation of an N-H carbamate into isocyanate + alcohol (above ~150 C).", 4.0, maxn=1)
T(["urea"], "therm", "urea-dissociation", "thermal dissociation of the urea to isocyanate + amine", "[NX3;!H0:1][CX3:2](=[O:3])[NX3:4]>>([N:1]=[C:2]=[O:3].[N:4])",
  "N-H ureas dissociate reversibly (above 130 C) into an isocyanate and an amine.", 3.5, maxn=2)


def _cyc(ring: int, nuc: str) -> None:
    chain = "[#6;!$([#6]=O):{}]"
    n = ring - 2                       # chain carbons between nucleophile and carbonyl
    idx = list(range(5, 5 + n))
    left = "".join(chain.format(i) for i in idx)
    closing = "".join(f"[C:{i}]" for i in idx)
    isN = nuc == "N"
    nu_s = "[NX3;H2,H1;!$(NC=[O,S,N]);!$(N-a):4]" if isN else "[OX2H1:4]"
    prod_nu = "[N:4]" if isN else "[O:4]"
    smarts = f"{nu_s}{left}[CX3:1](=[O:2])[OX2:3]>>({prod_nu}1{closing}[C:1]1=[O:2].[O:3])"
    T(ACID_RULES + ESTER_RULES + ["ester_aryl"], "therm", f"cyclisation-{nuc}{ring}",
      "intramolecular lactamisation" if isN else "intramolecular lactonisation", smarts,
      f"Intramolecular nucleophilic acyl substitution: the {'amine' if isN else 'hydroxyl'} attacks the carboxylic acid / ester carbonyl through a {ring}-membered transition state, closing the {ring}-membered {'lactam' if isN else 'lactone'} and expelling water or the alcohol (heat / acid catalysis).",
      -2.5 if isN else 0.5, vuln="High" if isN else "Moderate", maxn=1)


for _ring in (5, 6):
    _cyc(_ring, "N")
    _cyc(_ring, "O")


# --------------------------------------------------------------------------------------
# Special graph edits (regiochemistry decided in Python)
# --------------------------------------------------------------------------------------
def _aryl_hydroxylation(an: Analysis) -> List[Tuple[Chem.Mol, str]]:
    """hydroxyl-radical / peroxide attack at the most electron-rich free C-H of a benzenoid ring"""
    m = an.mol
    best = None
    for s in an.sites:
        if s.pkey != "aromatic":
            continue
        ring = list(s.atoms)
        if len(ring) != 6:
            continue
        for ia, a in sorted(enumerate(ring), key=lambda t: t[1]):
            at = m.GetAtomWithIdx(a)
            if at.GetTotalNumHs() < 1:
                continue
            sc = 0.0
            for ik, k in enumerate(ring):
                d = min(abs(ia - ik), 6 - abs(ia - ik))
                if d == 0:
                    continue
                for nb in m.GetAtomWithIdx(k).GetNeighbors():
                    if nb.GetIdx() in ring:
                        continue
                    o = d in (1, 3)
                    if an._is_donor(nb):
                        sc += 2 if o else 0.3
                    elif nb.GetSymbol() == "C" and not nb.GetIsAromatic() and all(b.GetBondTypeAsDouble() == 1 for b in nb.GetBonds()) and not an._is_ewg(nb):
                        sc += 1 if o else 0.2
                    elif nb.GetSymbol() in _HAL:
                        sc += 0.4 if o else -0.2
                    elif an._is_ewg(nb):
                        sc += -2 if o else 0.3
            if best is None or sc > best[0]:
                best = (sc, a, s)
    if best is None:
        return []
    _, a, s = best
    rw = Chem.RWMol(m)
    o = rw.AddAtom(Chem.Atom(8))
    rw.AddBond(a, o, Chem.BondType.SINGLE)
    try:
        Chem.SanitizeMol(rw)
    except Exception:
        return []
    return [(rw.GetMol(), "activated" if s.rule == "aromatic_activated" else "deactivated" if s.rule == "aromatic_deactivated" else "plain")]


def _hydrocarbon_oxidation(an: Analysis) -> List[Chem.Mol]:
    m = an.mol
    cs = [a for a in m.GetAtoms() if a.GetSymbol() == "C" and not a.GetIsAromatic() and a.GetTotalNumHs() >= 1
          and all(b.GetBondTypeAsDouble() == 1 for b in a.GetBonds())]
    if not cs:
        return []
    cs.sort(key=lambda a: (a.GetTotalNumHs(), a.GetIdx()))
    rw = Chem.RWMol(m)
    o = rw.AddAtom(Chem.Atom(8))
    rw.AddBond(cs[0].GetIdx(), o, Chem.BondType.SINGLE)
    try:
        Chem.SanitizeMol(rw)
    except Exception:
        return []
    return [rw.GetMol()]


def gen_stress(an: Analysis) -> List[dict]:
    """all intrinsic stress pathways of the analysed compound"""
    out: List[dict] = []
    frs = frag_mols(an.mol)

    def emit(ck: str, key: str, label: str, mol: Chem.Mol, mech: str, dG: float, vuln: str, mult: float = 1.0) -> None:
        out.append(dict(cond=COND_LABEL[ck], ck=ck, key=key, label=label, mol=mol, mech=clean_text(mech), dG=dG, vuln=vuln, mult=mult,
                        source=SRC_STRESS, co=None, physical=False, merged=False))

    for t in STRESS:
        v = t.vuln or an.vuln(t.rules, t.ck)
        if v is None:
            continue
        if t.same:
            groups = [(f, f) for f in frs]
        else:
            groups = [(f,) for f in frs]
        for reac in groups:
            for p in react(t.rxn, reac, t.loses, t.maxn, t.kek, t.post):
                emit(t.ck, t.key, t.label, p, t.mech, t.dG, v, t.mult)
    for p, kind in _aryl_hydroxylation(an):
        emit("ox", "aromatic-hydroxylation", "aromatic hydroxylation by hydroxyl radical / peroxide (phenol formation)", p,
             "Electrophilic HO. / peroxide attack at the most electron-rich free ring C-H (ortho/para to donors, away from electron-withdrawing groups) gives the phenol.",
             -2.0, an.vuln(["aromatic", "aromatic_activated", "aromatic_deactivated"], "ox") or "Moderate", 0.9 if kind == "activated" else 0.3 if kind == "deactivated" else 0.6)
    if len(an.sites) == 1 and an.sites[0].rule == "hydrocarbon":
        for p in _hydrocarbon_oxidation(an):
            emit("ox", "alkane-autoxidation", "C-H autoxidation to the alcohol (tertiary above secondary above primary)", p,
                 "Radical chain autoxidation abstracts the weakest C-H; O2 capture and decomposition of the hydroperoxide give the alcohol.", -3.0, "Moderate", 0.7)
    return out


# --------------------------------------------------------------------------------------
# Co-reactant (secondary compound) interaction templates
# --------------------------------------------------------------------------------------
class XTemplate:
    __slots__ = ("e_rules", "n_rules", "key", "label", "rxn", "mech", "dG", "mult", "vuln", "cond", "physical", "loses", "maxn", "acid_dep", "kek")

    def __init__(self, e_rules, n_rules, key, label, smarts, mech, dG, cond="Thermal Degradation", mult=1.0, vuln=None,
                 physical=False, loses="", maxn=2, acid_dep=False, kek=False):
        self.e_rules, self.n_rules, self.key, self.label = list(e_rules), list(n_rules), key, label
        self.mech, self.dG, self.mult, self.vuln, self.cond = clean_text(mech), dG, mult, vuln, cond
        self.physical, self.loses, self.maxn, self.acid_dep, self.kek = physical, loses, maxn, acid_dep, kek
        self.rxn = AllChem.ReactionFromSmarts(smarts)
        if self.rxn is None or self.rxn.GetNumReactantTemplates() != 2:
            raise ValueError(f"bad cross reaction SMARTS for {key}: {smarts}")


XT: List[XTemplate] = []


def X(*a, **k) -> None:
    XT.append(XTemplate(*a, **k))


NUC_N = ["amine", "ammonia", "arylamine", "arylamine_hetero", "hydrazine", "hydroxylamine"]
NUC_O = ["alcohol", "alcohol_tert", "alcohol_benzylic", "phenol", "phenol_hindered", "phenol_ewg_o", "phenol_ewg_p"]
NUC_S = ["thiol"]
_N = "[NX3;!H0;!$(N-C=[O,S,N]);!$(N-[SX4,PX4]=O):5]"            # amine / hydrazine / hydroxylamine with a removable H
_N1 = "[NX3;H2;!$(N-C=[O,S,N]);!$(N-[SX4,PX4]=O):5]"           # primary
_O = "[OX2H1:5][#6;!$(C=O)]"
_S = "[SX2H1:5]"
_NUC_ANY = "[N,S;!H0;!$(N-C=[O,S,N]);!$(N-[SX4,PX4]=O):5]"
_BASE_N = "[NX3;!$(N-C=[O,S,N]);!$(N-[SX4,PX4]=O);!$(N-[N,O]);!$(N-a);!$(N=*):5]"

# ---- acyl transfer: aminolysis / transesterification / thioester formation
_ACYL = [
    (ESTER_RULES, "ester", "[CX3:1](=[O:2])!@[OX2:3][CX4:4]", "([C:1](=[O:2]){nu}.[O:3][#6:4])", -3.5, 0.4, 0.9),
    (["ester_aryl"], "aryl/vinyl ester", "[CX3:1](=[O:2])!@[OX2:3][c,$([CX3]=[CX3]):4]", "([C:1](=[O:2]){nu}.[O:3][*:4])", -5.0, -2.0, 1.0),
    (LACTONE_RULES, "lactone", "[CX3:1](=[O:2])@[OX2:3][#6:4]", "([C:1](=[O:2]){nu}.[O:3][#6:4])", -3.0, 0.4, 0.9),
    (["anhydride"], "anhydride", "[CX3:1](=[O:2])[OX2:3][CX3:4]=[O:6]", "([C:1](=[O:2]){nu}.[O:3][C:4]=[O:6])", -9.0, -8.0, 1.0),
    (["acyl_halide"], "acyl halide", "[CX3:1](=[O:2])[F,Cl,Br,I:3]", "([C:1](=[O:2]){nu})", -9.0, -9.0, 1.0),
    (["thioester"], "thioester", "[CX3:1](=[O:2])[SX2:3][#6:4]", "([C:1](=[O:2]){nu}.[S:3][#6:4])", -4.0, -2.0, 0.9),
    (["carbonate"], "carbonate", "[#6:6][OX2:7][CX3:1](=[O:2])[OX2:3][#6:4]", "([C:1](=[O:2])([O:3][#6:4]){nu}.[O:7][#6:6])", -6.0, -3.0, 0.9),
    (["beta_lactam"], "beta-lactam", "[CX3:1](=[O:2])@[NX3:3]", "([C:1](=[O:2]){nu}.[N:3])", -9.0, -8.0, 1.0),
]
_O9 = "[OX2H1:5][#6;!$(C=O):9]"
for _er, _nm, _ep, _pp, _dgN, _dgO, _w in _ACYL:
    _loses = "X" if _nm == "acyl halide" else ""
    X(_er, NUC_N, "aminolysis", f"aminolysis: {_nm} + amine to amide", f"{_ep}.{_N}>>{_pp.format(nu='[N:5]')}",
      "Nucleophilic addition-elimination: the amine nitrogen of {N} attacks the electrophilic " + _nm + " carbonyl of {E}, expels the leaving group and forms a covalent amide adduct.",
      _dgN, mult=_w, loses=_loses)
    X(_er, NUC_O, "transesterification", f"acyl transfer: {_nm} + alcohol/phenol (transesterification)", f"{_ep}.{_O9}>>{_pp.format(nu='[O:5][#6:9]')}",
      "Intermolecular acyl transfer: the hydroxyl oxygen of {N} attacks the " + _nm + " carbonyl of {E} (accelerated by heat, moisture and trace acid / base), giving a new ester and the released alcohol.",
      _dgO, mult=_w * 0.6, loses=_loses)
    if _nm in ("ester", "aryl/vinyl ester", "anhydride", "acyl halide", "carbonate"):
        X(_er, NUC_S, "thioester-formation", f"acyl transfer: {_nm} + thiol (thioester formation)", f"{_ep}.{_S}>>{_pp.format(nu='[S:5]')}",
          "Thiolate / thiol of {N} attacks the " + _nm + " carbonyl of {E} giving a thioester.", -1.5, mult=_w * 0.7, loses=_loses)
X(["isocyanate"], NUC_N + NUC_O + NUC_S, "isocyanate-addition", "isocyanate + nucleophile (urea / carbamate / thiocarbamate)",
  "[NX2:1]=[CX2:2]=[O:3]." + "[N,O,S;!H0;!$(N-C=[O,S,N]):5]>>[N:1][C:2](=[O:3])[*:5]",
  "Nucleophilic addition of the nucleophile of {N} to the cumulated N=C=O carbon of {E} gives the urea / carbamate / thiocarbamate (very fast).", -9.0)

# ---- acid + alcohol / amine (esterification / amidation)
X(ACID_RULES, NUC_O, "esterification", "Fischer esterification (acid + alcohol)",
  "[CX3:1](=[O:2])[OX2H1:3].[OX2H1:5][#6;!$(C=O):9]>>[C:1](=[O:2])[O:5][#6:9]",
  "Acid-catalysed condensation of the carboxylic acid of {E} with the hydroxyl of {N} gives an ester + water (equilibrium; driven by heat / water loss).", 0.8, mult=0.6, vuln="Moderate")
X(ACID_RULES, NUC_N, "amidation", "thermal amidation (acid + amine)", "[CX3:1](=[O:2])[OX2H1:3]." + _N + ">>[C:1](=[O:2])[N:5]",
  "Heating the carboxylic acid of {E} with the amine of {N} (via the ammonium carboxylate salt) dehydrates to the amide.", -1.0, mult=0.6, vuln="Moderate")

# ---- carbonyl + amine (Schiff base) and reducing sugar + amine (Maillard)
X(["aldehyde"], NUC_N, "condensation", "aldehyde + primary amine condensation to the imine",
  "[CX3;H1,H2:1]=[O:2]." + _N1 + ">>[C:1]=[N:5]",
  "Nucleophilic addition of the primary amine of {N} to the aldehyde carbonyl of {E} and dehydration gives the imine (Schiff base; hydrazone / oxime for hydrazines / hydroxylamines); reversible, water removal drives it.", 1.5)
X(["ketone"], NUC_N, "condensation", "ketone + primary amine condensation to the ketimine",
  "[#6][CX3;H0:1](=[O:2])[#6]." + _N1 + ">>[#6][C:1]([#6])=[N:5]",
  "Nucleophilic addition of the primary amine of {N} to the ketone of {E} and dehydration gives the ketimine (reversible).", 3.0, vuln="Moderate")
X(["hemiacetal"], NUC_N, "maillard", "Maillard condensation: reducing sugar + amine to the N-glycosylamine",
  "[CX4:1]([OX2H1:2])[OX2:3][#6:4]." + _N + ">>[C:1]([O:3][#6:4])[N:5]",
  "The masked aldehyde (anomeric hemiacetal) of {E} condenses with the unprotonated amine of {N} to the N-glycosylamine / Schiff base; the irreversible Amadori rearrangement and browning (melanoidins) that follow on heating / moisture pull the equilibrium forward.", -3.0)

# ---- epoxide ring opening / alkylation / Michael addition
for _nr, _ns, _nm, _w in ((NUC_N, _N, "amine", 1.0), (NUC_O, _O, "alcohol", 0.5), (NUC_S, _S, "thiol", 1.0)):
    _tail = "[N:5]" if _nm == "amine" else "[O:5][#6:9]" if _nm == "alcohol" else "[S:5]"
    _ns2 = _ns.replace("[#6;!$(C=O)]", "[#6;!$(C=O):9]")
    X(["epoxide"], _nr, "epoxide-opening", f"epoxide ring opening by {_nm}", f"[C:1]1[O:2][C:3]1.{_ns2}>>[C:1]([O:2])[C:3]{_tail}",
      "SN2 opening of the strained epoxide of {E} by the nucleophile of {N} gives the beta-substituted alcohol.", -8.0, mult=_w)
    X(HALIDE_ALKYL, _nr, "alkylation", f"alkylation by the alkyl halide ({_nm})", f"[CX4:1][Cl,Br,I:2].{_ns2}>>[C:1]{_tail}",
      "SN2 displacement of halide from {E} by the nucleophile of {N} gives the alkylated product + HX (alkyl halides are alkylating, genotoxic-type reagents).", -6.0,
      mult=0.4 if _nm == "alcohol" else 1.0, loses="X")
    X(["sulfonate_ester", "dialkyl_sulfate"], _nr, "alkylation-sulfonate", f"alkylation by the sulfonate / sulfate ester ({_nm})",
      f"[CX4:1][OX2:2][SX4:3](=[O:4])(=[O:6])[*:7].{_ns2}>>[C:1]{_tail}.[O:2][S:3](=[O:4])(=[O:6])[*:7]",
      "{E} is an alkylating agent: SN2 attack by the nucleophile of {N} on the alkyl carbon releases the sulfonic acid and gives the alkylated adduct.", -7.0,
      mult=0.4 if _nm == "alcohol" else 1.0)
    if _nm != "alcohol":
        X(["michael"], _nr, "michael-addition", f"{'aza' if _nm == 'amine' else 'thia'}-Michael addition to the alpha,beta-unsaturated carbonyl",
          f"[CX3;!a:1](=[O:2])[CX3;!a:3]=[CX3;!a:4].{_ns}>>[C:1](=[O:2])[C:3][C:4]{_tail}",
          "Conjugate addition of the nucleophile of {N} to the electrophilic beta-carbon of {E} gives the beta-substituted carbonyl adduct.", -5.0)
X(HALIDE_ALKYL, ["amine_tert"], "quaternisation", "quaternisation of the tertiary amine (Menshutkin)",
  "[CX4:1][Cl,Br,I:2]." + _BASE_N.replace(";!$(N=*)", ";H0;!$(N=*)") + ">>[C:1][N+:5]",
  "Menshutkin reaction: the tertiary amine of {N} displaces halide at {E}, giving a quaternary ammonium salt.", -6.0, mult=0.8, loses="X")

# ---- N-nitrosation / diazotisation by a nitrosating co-reactant (nitrite)
_SEC_N = "[NX3;H1;!$(N-C=[O,S,N]);!$(N-[SX4,PX4]=O);!$(N-[N,O]):1]"
X(["nitrite"], ["amine", "arylamine", "arylamine_hetero"], "n-nitrosation", "N-nitrosation of the secondary amine (nitrosamine formation)",
  "[NX2:2](=[O:3])[O-:4]." + _SEC_N + ">>[N:1][N:2]=[O:3]",
  "Under acidic conditions nitrite ({E}) forms HONO / N2O3 / NO+, which nitrosates the secondary amine of {N} to the N-nitrosamine (a mutagenic impurity class).", -8.5,
  cond="Acidic Hydrolysis", vuln="Critical", acid_dep=True)
X(["nitrite"], ["amine_tert", "arylamine"], "nitrosative-dealkylation", "nitrosative dealkylation of the tertiary amine to the N-nitrosamine",
  "[NX2:2](=[O:3])[O-:4].[NX3;H0;!$(N-C=[O,S,N]);!$(N-[SX4,PX4]=O);!$(N-[N,O]);!$(N=*):1][CX4;!H0:6]>>([N:1][N:2]=[O:3].[C:6]=O)",
  "Nitrosation of the tertiary amine of {N} by {E} forms an unstable N-nitrosammonium species that loses an aldehyde (dealkylation), leaving the secondary N-nitrosamine.", -5.0,
  cond="Acidic Hydrolysis", mult=0.4, vuln="Moderate", acid_dep=True)
X(["nitrite"], ["arylamine"], "diazotisation", "diazotisation of the primary aromatic amine (aryl diazonium ion)",
  "[NX2:2](=[O:3])[O-:4].[NX3;H2:1][c;$(c1ccccc1):6]>>[N+:1](#[N:2])[c:6]",
  "In acid, nitrite ({E}) diazotises the primary aromatic amine of {N} to the aryl diazonium ion (which couples, hydrolyses to phenols, or loses N2).", -6.0,
  cond="Acidic Hydrolysis", vuln="High", acid_dep=True, loses="")

# ---- proton transfer salts (physical / ionic).  cond = (E is primary, E is co-reactant)
_ACC = "[NX3;!$(N-C=[O,S,N]);!$(N-[SX4,PX4]=O);!$(N-[N,O]);!$(N-a);!$(N=*):5]"
_SALT_TXT = "Proton transfer from the {what} of {E} to the basic nitrogen of {N} gives an ammonium-anion salt (changes solubility, hygroscopicity and the pH of the microenvironment; may precipitate)."
X(ACID_RULES, ["amine", "amine_tert", "amidine", "ammonia"], "proton-transfer-salt", "acid-base salt formation (ionic pair, proton transfer)",
  "[CX3:1](=[O:2])[OX2H1:3]." + _ACC + ">>[C:1](=[O:2])[O-:3].[N+:5]", _SALT_TXT.replace("{what}", "carboxylic acid"), -4.5,
  cond=("Basic Hydrolysis", "Acidic Hydrolysis"), vuln="High", physical=True)
X(["sulfonic", "sulfate", "hydrogen_halide"], ["amine", "amine_tert", "amidine", "ammonia", "arylamine", "arylamine_hetero", "azine"], "proton-transfer-salt-strong",
  "strong acid - base salt formation (ionic pair)", "[Cl,Br,I,F;H1;D0:1]." + "[NX3;!$(N-C=[O,S,N]);!$(N-[SX4,PX4]=O);!$(N-[N,O]);!$(N=*):5]>>[*-:1].[N+:5]",
  _SALT_TXT.replace("{what}", "strong acid"), -6.5, cond=("Basic Hydrolysis", "Acidic Hydrolysis"), vuln="Critical", physical=True)
X(["sulfonic"], ["amine", "amine_tert", "amidine", "ammonia", "arylamine", "arylamine_hetero", "azine"], "proton-transfer-salt-strong",
  "sulfonic acid - base salt formation (ionic pair)", "[SX4:1](=O)(=O)([#6:2])[OX2H1:3]." + "[NX3;!$(N-C=[O,S,N]);!$(N-[SX4,PX4]=O);!$(N-[N,O]);!$(N=*):5]>>[S:1](=O)(=O)([#6:2])[O-:3].[N+:5]",
  _SALT_TXT.replace("{what}", "sulfonic acid"), -6.5, cond=("Basic Hydrolysis", "Acidic Hydrolysis"), vuln="Critical", physical=True)
X(["azine"], ["hydrogen_halide"], "proton-transfer-salt-strong", "hydrohalide salt of the heteroaromatic base", "[Cl,Br,I,F;H1;D0:1].[nX2:5]>>[*-:1].[n+:5]", _SALT_TXT.replace("{what}", "strong acid"), -5.0,
  cond=("Acidic Hydrolysis", "Basic Hydrolysis"), vuln="High", physical=True)
X(["acid_alphaEWG", "acid_aryl_o", "acid_aryl_p", "acid_picolinic"], ["carboxylate"], "acid-base-exchange", "acid-base exchange (stronger acid displaces the weaker carboxylate)",
  "[CX3:1](=[O:2])([OX2H1:3])[c:9].[CX3:6](=[O:7])([OX1-:8])[CX4:10]>>[C:1](=[O:2])([O-:3])[c:9].[C:6](=[O:7])([O:8])[C:10]",
  "{N} is the salt of a weak acid (a mild base); the stronger aromatic / activated carboxylic acid of {E} protonates the carboxylate, liberating the free weaker acid and forming the salt of {E} with the metal cation (alters solubility and raises the local pH).",
  -3.0, cond=("Basic Hydrolysis", "Acidic Hydrolysis"), vuln="High", physical=True)
X(["acid"], ["carboxylate"], "acid-base-exchange", "acid-base exchange (aromatic acid displaces an aliphatic carboxylate)",
  "[CX3:1](=[O:2])([OX2H1:3])[c:9].[CX3:6](=[O:7])([OX1-:8])[CX4:10]>>[C:1](=[O:2])([O-:3])[c:9].[C:6](=[O:7])([O:8])[C:10]",
  "{N} is the salt of a weak acid (a mild base); the stronger aromatic carboxylic acid of {E} protonates the carboxylate, liberating the free weaker acid and forming the salt of {E} with the metal cation (alters solubility and raises the local pH).",
  -3.0, cond=("Basic Hydrolysis", "Acidic Hydrolysis"), vuln="High", physical=True)


def _min_v(a: str, b: str) -> str:
    return a if vrank(a) <= vrank(b) else b


# ---- reductions: co-reactant supplies the reducing power (gated on rules present in the co-reactant)
RED: List[Tuple[List[str], List[str], str, str, str, str, float, float, str, Any]] = []


def _red(a_rules, gate, key, label, smarts, mech, dG, mult, vuln, loses="") -> None:
    RED.append((a_rules, gate, key, label, AllChem.ReactionFromSmarts(smarts), clean_text(mech), dG, mult, vuln, loses))


_HYD = ["hydride"]
_RED_ANY = ["hydride", "hydrazine", "sulfite", "thiol", "hydroxylamine", "enol"]
_red(["aldehyde", "ketone"], _HYD, "hydride-reduction", "reduction of the carbonyl to the alcohol by {B}", "[CX3;!$(C(=O)[!#6]):1]=[O:2]>>[C:1][O:2]",
     "Hydride transfer from {B} to the carbonyl carbon of {A}, then protonation, gives the alcohol (redox interaction).", -8.0, 1.0, "Critical")
_red(["quinone_para"], _RED_ANY, "quinone-reduction", "reduction of the quinone to the hydroquinone by {B}",
     "[#6:1]1(=[O:7])[#6:2]=[#6:3][#6:4](=[O:8])[#6:5]=[#6:6]1>>[c:1]1([OH:7])[c:2][c:3][c:4]([OH:8])[c:5][c:6]1",
     "Two-electron, two-proton reduction of the quinone of {A} by {B} gives the aromatic hydroquinone.", -6.0, 0.9, "High")
_red(["quinone_para"], _RED_ANY, "quinone-imine-reduction", "reduction of the quinone-imine to the aminophenol by {B}",
     "[#6:1]1(=[O:7])[#6:2]=[#6:3][#6:4](=[N:8])[#6:5]=[#6:6]1>>[c:1]1([OH:7])[c:2][c:3][c:4]([NH:8])[c:5][c:6]1",
     "Two-electron, two-proton reduction of the quinone-imine of {A} by {B} gives the aromatic aminophenol.", -6.0, 0.9, "High")
_red(["disulfide"], _RED_ANY, "disulfide-reduction", "disulfide reduction to thiols by {B}", "[SX2:1][SX2:2]>>([S:1].[S:2])",
     "Thiol-disulfide exchange / hydride reduction by {B} cleaves the S-S bond of {A} to two thiols.", -4.0, 1.0, "High")
_red(["nitro"], ["hydride", "hydrazine", "sulfite"], "nitro-reduction", "nitro group reduction to the amine by {B}", "[c:1][N+:2](=O)[O-]>>[c:1][N+0:2]",
     "{B} reduces the aromatic nitro group of {A} through nitroso and hydroxylamine to the aniline (usually needs a catalyst / forcing conditions).", -9.0, 0.6, "Moderate")
_red(["peroxide"], ["hydride", "hydrazine", "sulfite", "thiol"], "peroxide-reduction", "hydroperoxide reduction to the alcohol by {B}", "[OX2H1][OX2:2][#6:3]>>[O:2][#6:3]",
     "{B} reduces the O-O bond of the hydroperoxide of {A} to the alcohol.", -12.0, 1.0, "Critical")
_red(["azide", "azide2"], ["hydride", "hydrazine", "sulfite"], "azide-reduction", "azide reduction to the amine by {B}", "[#6:1][N:2]=[N+]=[N-]>>[#6:1][N:2]",
     "{B} reduces the azide of {A} to the primary amine with N2 loss.", -10.0, 0.8, "High", "N")


def metal_salt(A: Analysis, B: Analysis, nameA: str, nameB: str) -> List[dict]:
    """carboxylic acid + basic metal species (oxide / hydroxide / carbonate) -> metal carboxylate (soap) salt"""
    if "base_strong" not in B.classes:
        return []
    ms = [s for s in B.sites if s.pkey in ("metal", "metal_transition")]
    if not ms or not any(s.pkey in ("strong_base", "carbonate_ion") for s in B.sites):
        return []
    at = B.mol.GetAtomWithIdx(ms[0].key)
    z = max(1, abs(at.GetFormalCharge()))
    acids = [s for s in A.sites if s.pkey == "acid"]
    if not acids:
        return []
    fr = None
    for f in frag_mols(A.mol):
        if any(a.GetSymbol() == "O" and a.GetTotalNumHs() == 1 for a in f.GetAtoms()):
            fr = f
            break
    if fr is None:
        return []
    q = Chem.MolFromSmarts("[CX3](=O)[OX2H1:1]")
    m0 = fr.GetSubstructMatch(q)
    if not m0:
        return []
    o_idx = m0[2]
    pieces = []
    for _ in range(min(z, 3)):
        rw = Chem.RWMol(fr)
        a = rw.GetAtomWithIdx(o_idx)
        a.SetFormalCharge(-1)
        a.SetNumExplicitHs(0)
        a.SetNoImplicit(True)
        pieces.append(rw.GetMol())
    prod = pieces[0]
    for p in pieces[1:]:
        prod = Chem.CombineMols(prod, p)
    metal = Chem.RWMol()
    ma = Chem.Atom(at.GetAtomicNum())
    ma.SetFormalCharge(z)
    metal.AddAtom(ma)
    prod = Chem.CombineMols(prod, metal.GetMol())
    try:
        Chem.SanitizeMol(prod)
    except Exception:
        return []
    return [dict(cond="Basic Hydrolysis", ck="base", key="metal-carboxylate", label=f"metal carboxylate (soap) formation with {at.GetSymbol()}{z}+", mol=prod,
                 mech=clean_text(f"The basic {nameB} (oxide / hydroxide / carbonate) neutralises the carboxylic acid of {nameA} stoichiometrically, forming the {at.GetSymbol()}({z}+) carboxylate salt (ionic; Lewis-acid cation coordination; often insoluble)."),
                 dG=-6.0, vuln="Critical", mult=1.0, source=SRC_CROSS, co=nameB, physical=True)]


def environment_boost(stress: List[dict], B: Analysis, nameB: str) -> List[dict]:
    """copies of the primary compound's own pathways, re-graded because the co-reactant changes the microenvironment"""
    out: List[dict] = []
    cl = B.classes
    strong_acid = "acid_strong" in cl
    photocat = any(s.pkey == "metal_oxide" and B.mol.GetAtomWithIdx(s.key).GetSymbol() == "Ti" for s in B.sites)
    zno = any(s.pkey == "metal" and B.mol.GetAtomWithIdx(s.key).GetSymbol() == "Zn" for s in B.sites) and any(s.pkey == "strong_base" for s in B.sites)

    def push(r: dict, why: str) -> None:
        r2 = dict(r)
        r2.update(vuln=step_vuln(r["vuln"], 1), source=SRC_CROSS, co=nameB, mech=clean_text(why + " " + r["mech"]), physical=False)
        out.append(r2)

    for r in stress:
        c = r["cond"]
        if c == "Basic Hydrolysis" and "base_strong" in cl:
            push(r, f"Alkaline microenvironment created by {nameB} (local pH much greater than 9) accelerates this pathway.")
        elif c == "Basic Hydrolysis" and "base_weak" in cl and "metal" in cl:
            push(r, f"{nameB} is the metal salt of a weak acid: dissolving in the moisture film it raises the local pH (mildly alkaline microenvironment) and its cation acts as a Lewis acid, accelerating base-catalysed cleavage.")
        elif c == "Acidic Hydrolysis" and strong_acid:
            push(r, f"Acidic microenvironment created by {nameB} accelerates this pathway.")
        elif c == "Hydrolysis" and "water" in cl:
            push(r, f"Water / moisture supplied by {nameB} drives this hydrolysis.")
        elif c == "Oxidation" and "oxidant" in cl:
            push(r, f"{nameB} is an oxidant / contains peroxide and supplies the oxidising equivalents.")
        elif c == "Oxidation" and "peroxide_former" in cl:
            push(r, f"Autoxidation of {nameB} generates trace hydroperoxides / radicals that oxidise susceptible groups.")
        elif c == "Oxidation" and "redox_metal" in cl:
            push(r, f"Redox-active metal ions in {nameB} catalyse radical (Fenton-type) oxidation.")
        elif c == "Photodegradation" and (photocat or zno):
            push(r, f"{nameB} (TiO2 / ZnO-type semiconductor) is a photocatalyst that generates reactive oxygen species under UV.")
    return out


def gen_cross(A: Analysis, B: Analysis, nameA: str, nameB: str, stress: List[dict], ambient_acid: bool = False) -> List[dict]:
    """all interaction pathways between the primary compound A and one co-reactant B"""
    out: List[dict] = []
    acid_present = ambient_acid or bool({"acid_weak", "acid_strong"} & (A.classes | B.classes))
    fa, fb = frag_mols(A.mol), frag_mols(B.mol)
    for (E, N, nE, nN, e_is_a, fes, fns) in ((A, B, nameA, nameB, True, fa, fb), (B, A, nameB, nameA, False, fb, fa)):
        for t in XT:
            ve = E.sec_vuln(t.e_rules)
            vn = N.sec_vuln(t.n_rules) if t.n_rules else ve
            if ve is None or vn is None:
                continue
            v = t.vuln or _min_v(ve, vn)
            cond = t.cond if isinstance(t.cond, str) else (t.cond[0] if e_is_a else t.cond[1])
            mult = t.mult * (1.0 if (not t.acid_dep or acid_present) else 0.55)
            txt = t.mech.replace("{E}", nE).replace("{N}", nN)
            for fe in fes:
                for fn in fns:
                    for p in react(t.rxn, (fe, fn), t.loses, t.maxn, t.kek):
                        out.append(dict(cond=cond, ck="acid", key=t.key, label=t.label, mol=p, mech=txt, dG=t.dG, vuln=v, mult=mult,
                                        source=SRC_CROSS, co=nameB, physical=t.physical))
    for a_rules, gate, key, label, rxn, mech, dG, mult, vuln, loses in RED:
        if not any(g in B.present or g in {s.pkey for s in B.sites} for g in gate):
            continue
        if A.vuln(a_rules, "ox") is None:
            continue
        for f in fa:
            for p in react(rxn, (f,), loses, 2):
                out.append(dict(cond="Oxidation", ck="ox", key=key, label=label.replace("{B}", nameB), mol=p,
                                mech=mech.replace("{B}", nameB).replace("{A}", nameA), dG=dG, vuln=vuln, mult=mult, source=SRC_CROSS, co=nameB, physical=False))
    out += metal_salt(A, B, nameA, nameB)
    out += environment_boost(stress, B, nameB)
    for r in out:
        r["co"] = nameB
    return out


# --------------------------------------------------------------------------------------
# Candidate finalisation, naming, ranking
# --------------------------------------------------------------------------------------
_TRIVIAL = {"O", "O=C=O", "[C-]#[O+]", "N#N", "O=S=O", "Cl", "Br", "I", "F", "[Cl-]", "[Br-]", "[I-]", "[F-]", "[H][H]"}


def _prune_trivial(mol: Chem.Mol) -> Tuple[Optional[Chem.Mol], List[str]]:
    """drop trivial co-products (H2O, CO2, CO, N2, SO2, HX, halide ions) and return the rest as one molecule"""
    keep, dropped = [], []
    for f in frag_mols(mol):
        c = canon(f)
        if c in _TRIVIAL:
            dropped.append({"O": "H2O", "O=C=O": "CO2", "[C-]#[O+]": "CO", "N#N": "N2", "O=S=O": "SO2"}.get(c, "H" + c.strip("[]-") if c in ("Cl", "Br", "I", "F") else c))
        else:
            keep.append(f)
    if not keep:
        return None, dropped
    out = keep[0]
    for f in keep[1:]:
        out = Chem.CombineMols(out, f)
    return out, dropped


def _fragment_name(f: Chem.Mol) -> str:
    n = library_name(f)
    if n:
        return n
    try:
        an = Analysis(f)
        names = [e["name"] for e in an.entries if not re.match(r"^(Aromatic Ring|Benzylic|Aliphatic / Polar|Aliphatic Alcohol)", e["name"])]
        main = names[0] if names else (an.entries[0]["name"] if an.entries else "species")
    except Exception:
        main = "species"
    main = re.sub(r"\s*\(.*?\)", "", main.split(" / ")[0]).lower()
    return f"{formula_of(f)} {main}"


def _cap(s: str) -> str:
    return s[:1].upper() + s[1:]


def _finalize(rx: dict, parent_canon: Set[str], ref_mass: float) -> Optional[dict]:
    pr, dropped = _prune_trivial(rx["mol"])
    if pr is None:
        return None
    try:
        Chem.SanitizeMol(pr)
    except Exception:
        return None
    frs = sorted(frag_mols(pr), key=lambda f: (-f.GetNumHeavyAtoms(), canon(f)))
    strs = [canon(f) for f in frs]
    if all(s in parent_canon for s in strs):
        return None
    smiles = ".".join(strs)
    main = frs[0]
    mass = mono_mass(main)
    like = max(0.01, min(0.99, VSCORE[rx["vuln"]] * rx.get("mult", 1.0)))
    return dict(rx, smiles=smiles, mol=pr, key_smiles=".".join(sorted(strs)), names=" + ".join(_fragment_name(f) for f in frs),
                dropped=dropped, main_formula=formula_of(main), main_mass=mass, delta_mass=mass - ref_mass, like=like, also=[])


def _fmt_mass(x: float) -> str:
    return ("+" if x >= 0 else "") + f"{x:.4f}"


def _interaction_type(cross_pool: List[dict], has_co: bool) -> str:
    if not has_co:
        return "None"
    chem = [c for c in cross_pool if not c["physical"]]
    if chem:
        return "Chemical"
    return "Physical" if cross_pool else "None"


# --------------------------------------------------------------------------------------
# Public API
# --------------------------------------------------------------------------------------
def _display_name(m: Chem.Mol, fallback: str) -> str:
    n = library_name(m)
    return n.replace(" (or a stereoisomer)", "") if n else fallback


def predict(primary_smiles: str, secondary_smiles_list: Sequence[str], method: str = "Both") -> Dict[str, Any]:
    """
    Degradation / interaction prediction.
    Returns dict(ok, reason, functional_groups, impurities, chain_of_thought, interaction_type, mechanism)
    with impurities already ranked (top 5) in the field names the Streamlit app uses.
    """
    pm = parse_smiles(primary_smiles)
    if pm is None:
        return {"ok": False, "reason": "The primary SMILES could not be parsed as a valid molecule.", "functional_groups": [], "impurities": [],
                "chain_of_thought": "", "interaction_type": "None", "mechanism": ""}
    if pm.GetNumHeavyAtoms() > 400:
        return {"ok": False, "reason": "The structure has more than 400 heavy atoms, which is beyond the size the reaction engine is designed for (polymers / biomacromolecules).",
                "functional_groups": [], "impurities": [], "chain_of_thought": "", "interaction_type": "None", "mechanism": ""}
    an = Analysis(pm)
    pname = _display_name(pm, "the primary compound")
    porigin = pname if pname != "the primary compound" else "Primary compound"
    p_frags = frag_mols(pm)
    parent_canon = {canon(f) for f in p_frags}
    big = max(p_frags, key=lambda f: f.GetNumHeavyAtoms())
    ref_mass = mono_mass(big)
    p_formula = formula_of(big)

    co: List[Tuple[str, str, Chem.Mol, Analysis]] = []
    unresolved: List[str] = []
    for i, s in enumerate(secondary_smiles_list or []):
        if not s or not s.strip():
            continue
        m = parse_smiles(s)
        if m is None:
            unresolved.append(s.strip())
            continue
        co.append((s.strip(), _display_name(m, f"Co-reactant {len(co) + 1}"), m, Analysis(m)))

    cands: List[dict] = []

    def add(rxs: Iterable[dict], extra_parents: Set[str]) -> None:
        for r in rxs:
            try:
                c = _finalize(r, parent_canon | extra_parents, ref_mass)
            except Exception:
                c = None
            if c:
                cands.append(c)

    stress = gen_stress(an)[:800]
    add(stress, set())
    ambient_acid = any({"acid_strong", "acid_weak"} & a.classes for _, _, _, a in co) or bool({"acid_strong", "acid_weak"} & an.classes)
    for smi, nm, m, a in co:
        co_parents = {canon(f) for f in frag_mols(m)}
        add(gen_cross(an, a, pname, nm, stress, ambient_acid), co_parents)

    # ---- de-duplicate identical products (keep the most likely origin, remember other conditions)
    by: Dict[str, dict] = {}
    for c in cands:
        k = c["key_smiles"] + ("|ionic" if c["physical"] else "")
        cur = by.get(k)
        if cur is None:
            by[k] = c
            continue
        cross_c, cross_cur = c["source"] == SRC_CROSS, cur["source"] == SRC_CROSS
        if c["like"] > cur["like"] + 1e-9 or (abs(c["like"] - cur["like"]) <= 1e-9 and cross_c and not cross_cur):
            c["also"] = [x for x in dict.fromkeys(cur["also"] + [cur["cond"]]) if x != c["cond"]]
            by[k] = c
        elif c["cond"] != cur["cond"] and c["cond"] not in cur["also"]:
            cur["also"].append(c["cond"])
    unique = list(by.values())

    # ---- Boltzmann distribution (T = 298.15 K), numerically stabilised: P_i = exp(-dG_i/RT) / sum_j exp(-dG_j/RT)
    R = 0.0019872  # kcal/(mol*K)
    T_K = 298.15
    RT = R * T_K
    raw = [-c["dG"] / RT for c in unique]
    mx = max(raw) if raw else 0.0
    exps = [math.exp(max(-500.0, min(500.0, e - mx))) for e in raw]
    tot = sum(exps) or 1.0

    ranked: List[dict] = []
    for i, c in enumerate(unique):
        p_b = round(min(0.99, max(0.01, exps[i] / tot)), 4)
        p_h = round(min(0.99, max(0.01, c["like"])), 4)
        if method == "Boltzmann":
            prob = p_b
        elif method == "Heuristic":
            prob = p_h
        else:
            prob = round((p_b + p_h) / 2.0, 4)
        cross = c["source"] == SRC_CROSS
        sub = f" Fragments: {c['names']}." if "." in c["smiles"] else ""
        trivial = f" Co-product(s) not listed: {', '.join(dict.fromkeys(c['dropped']))}." if c["dropped"] else ""
        desc = (f"Ionic / non-covalent association product ({c['main_formula']}; no covalent bond change).{sub}" if c["physical"] else
                f"{_cap(c['label'])}. Main species {c['main_formula']}, monoisotopic mass {c['main_mass']:.4f} Da ({_fmt_mass(c['delta_mass'])} Da vs {p_formula} parent).{sub}{trivial}")
        also = f" Also formed under: {', '.join(c['also'])}." if c["also"] else ""
        ranked.append({
            "iupacName": clean_text(f"{_cap(c['label'])}: {c['names']}"), "smiles": c["smiles"], "structureDescription": clean_text(desc),
            "origin": f"{porigin} + {c['co']}" if cross else porigin, "condition": c["cond"], "source": c["source"],
            "mechanismExplanation": clean_text(c["mech"] + also), "deltaG": round(c["dG"], 2), "kineticLikelihood": p_h,
            "probability": prob, "probabilityBoltzmann": p_b, "probabilityHeuristic": p_h, "_cross": cross, "_h": c["like"], "_physical": c["physical"], "_label": c["label"],
            "_mech": c["mech"],
        })
    ranked.sort(key=lambda x: x["probability"], reverse=True)
    top = ranked[:5]
    pool = [x for x in ranked if x["_cross"] and x["_h"] >= 0.5]
    if co and pool and not any(t in pool for t in top):
        top = ranked[:4] + [pool[0]]

    # interaction classification uses every co-reactant pathway before de-duplication
    cross_all = [c for c in cands if c["source"] == SRC_CROSS and c["like"] >= 0.5]
    cross_all.sort(key=lambda c: -c["like"])
    itype = _interaction_type(cross_all, bool(co))
    best = next((c for c in cross_all if not c["physical"]), cross_all[0] if cross_all else None)

    # ---- functional-group table (cross column evaluated against the actual co-reactants)
    fg_entries = an.entries
    if co:
        union: Set[str] = set()
        for _, _, _, a in co:
            union |= a.classes
        fg_entries = refine_secondary(an, " + ".join(nm for _, nm, _, _ in co), union)

    # ---- narrative
    def worst(e: dict, k: str) -> str:
        return e[k][0]

    lines = []
    for k in CKEYS:
        hits = sorted([e for e in an.entries if vrank(worst(e, k)) >= 2], key=lambda e: -vrank(worst(e, k)))[:3]
        lines.append(f"   - {COND_TITLE[k]}: " + ("; ".join(f"{e['name']} [{worst(e, k)}]" for e in hits) if hits else "no group above Moderate susceptibility"))
    if co:
        blocks = []
        for smi, nm, m, a in co:
            blocks.append(f"   - Co-reactant evaluated: {nm} (SMILES: {smi})\n   - Co-reactant functional groups: {', '.join(e['name'] for e in a.entries) or 'none detected'}\n"
                          f"   - Reactive classes detected: {'; '.join(CLASS_TEXT.get(c, c) for c in sorted(a.classes)) or 'none'}")
        co_block = ("CROSS-FUNCTIONAL INTERACTIONS WITH CO-REACTANTS:\n" + "\n".join(blocks) +
                    "\n   - Templates evaluated: acyl transfer (aminolysis / transesterification), Schiff-base / Maillard condensation, esterification / amidation, epoxide and alkyl-halide alkylation, "
                    "Michael addition, isocyanate addition, N-nitrosation / diazotisation, proton-transfer salts, metal carboxylates, redox (oxidant / reductant) and microenvironment (pH, moisture, peroxide, photocatalyst) effects."
                    f"\n   - Result: {best['label'] + (' (ionic / physical)' if best['physical'] else ' (covalent)') if best else 'no complementary reactive pair with meaningful likelihood; only non-covalent association expected'}.")
    else:
        co_block = "MULTI-COMPONENT MATRIX: Single-compound intrinsic forced degradation analysis."
    if unresolved:
        co_block += "\n   - Unparseable co-reactant SMILES ignored (no structure assumed): " + ", ".join(f'"{u}"' for u in unresolved) + "."
    groups_txt = ", ".join(e["name"] for e in an.entries[:4])
    cot = clean_text(f"""[Systematic Functional Group Reactivity & Computational Degradation Assessment]

PRIMARY MOLECULAR FUNCTIONAL GROUP INVENTORY:
   - Target Structure: {pname} (SMILES: {primary_smiles}; {p_formula}, monoisotopic {ref_mass:.4f} Da)
   - Identified Functional Groups: {', '.join(f"{e['name']} [{e['cat']}]" for e in an.entries)}
   - Identified Reactive Centers: {'; '.join(e['site'] + (f" (x{e['n']})" if e.get('n', 1) > 1 else '') for e in an.entries)}

REACTION SUSCEPTIBILITY BY STRESS CONDITION (most vulnerable groups):
{chr(10).join(lines)}

{co_block}

THERMODYNAMIC BOLTZMANN PARTITION & KINETIC PROBABILITIES (T = 298.15 K):
   - {len(cands)} candidate pathways were generated from reaction templates; {len(unique)} distinct products remained after de-duplication (identical products from several conditions are merged).
   - Relative formation free energies (dG, kcal/mol) are reaction-class estimates on a compressed scale (ring strain, bond-energy and ionisation terms), not quantum-chemical values - confirm key pathways with DFT if needed.
   - Boltzmann distribution derived via P_i = exp(-dG_i / RT) / sum exp(-dG_j / RT); heuristic likelihood derives from the vulnerability grade of the reacting group under the same condition (Critical 0.93, High 0.80, Moderate 0.58, Low 0.30, Resistant 0.10) scaled by pathway-specific factors. Method used: {method}.
   - Ranked top {len(top)} dominant degradation impurities and reaction adducts.""")
    if not unique:
        mech = f"{pname} shows no functional-group liability under the evaluated stress conditions ({groups_txt}); it is chemically inert or ionic under acid, base, neutral hydrolysis, light, heat and oxidative stress."
    elif co:
        mech = (f"Cross-functional interaction ({best['label']}): {best['mech']}" if best else
                f"No complementary reactive functional-group pair was found between {pname} and {' + '.join(nm for _, nm, _, _ in co)}; interaction is limited to non-covalent association. Intrinsic stress degradation is governed by {groups_txt}.")
    else:
        mech = f"Intrinsic stress degradation governed by hydrolytic, oxidative, photolytic and thermal reactivity of the functional groups present in {pname} ({groups_txt})."
    for t in top:
        for k in [k for k in t if k.startswith("_")]:
            del t[k]
    return {"ok": True, "reason": "", "functional_groups": group_dicts(fg_entries), "impurities": top, "chain_of_thought": cot,
            "interaction_type": itype, "mechanism": clean_text(mech), "n_candidates": len(cands), "n_unique": len(unique)}
