import numpy as np
import pandas as pd
import json
import warnings
warnings.filterwarnings('ignore')

from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors, AllChem
from sklearn.linear_model import LinearRegression, Ridge
from sklearn.ensemble import GradientBoostingRegressor
from sklearn.metrics import r2_score, mean_squared_error
from sklearn.model_selection import LeaveOneOut, cross_val_predict
from sklearn.preprocessing import StandardScaler
import plotly.graph_objects as go
import plotly.express as px
import plotly.io as pio
from plotly.subplots import make_subplots

# ── Dataset from Brusseau 2019 (Water Research) Table 1 ──────────────────────
raw = [
    # name,       SMILES,                                                                         Ki_cm,   class,  n_CF, branched
    ("PFAA",    "OC(=O)C(F)(F)F",                                                                1.52e-7, "PFCA", 2,  False),
    ("PFPrA",   "OC(=O)C(F)(F)C(F)(F)F",                                                        7.14e-7, "PFCA", 3,  False),
    ("PFBA",    "OC(=O)C(F)(F)C(F)(F)C(F)(F)F",                                                 1.72e-6, "PFCA", 4,  False),
    ("PFPeA",   "OC(=O)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F",                                         5.84e-6, "PFCA", 5,  False),
    ("PFHxA",   "OC(=O)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F",                                  2.20e-5, "PFCA", 6,  False),
    ("PFHpA",   "OC(=O)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F",                          5.76e-5, "PFCA", 7,  False),
    ("PFOA",    "OC(=O)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F",                   2.32e-4, "PFCA", 8,  False),
    ("PFNA",    "OC(=O)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F",            9.34e-4, "PFCA", 9,  False),
    ("PFDA",    "OC(=O)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F",    3.72e-3, "PFCA", 10, False),
    ("PFUnA",   "OC(=O)" + "C(F)(F)"*11 + "F",                                                  1.90e-2, "PFCA", 11, False),
    ("Iso-PFOA","OC(=O)C(F)(C(F)(F)F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F",                   2.39e-4, "PFCA", 8,  True),
    ("Iso-PFDA","OC(=O)C(F)(C(F)(F)F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F",    6.94e-3, "PFCA", 10, True),
    ("PFBS",    "OS(=O)(=O)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F",                                     1.78e-5, "PFSA", 4,  False),
    ("PFHxS",   "OS(=O)(=O)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F",                      9.72e-5, "PFSA", 6,  False),
    ("PFHpS",   "OS(=O)(=O)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F",              5.14e-4, "PFSA", 7,  False),
    ("PFOS",    "OS(=O)(=O)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F",       2.30e-3, "PFSA", 8,  False),
    ("PFNS",    "OS(=O)(=O)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F",3.70e-2, "PFSA", 9,  False),
]

cols = ["name","smiles","Ki_cm","head_class","n_CF","branched"]
df = pd.DataFrame(raw, columns=cols)
df["log_Ki"] = np.log10(df["Ki_cm"])

# ── Compute RDKit descriptors ─────────────────────────────────────────────────
def get_descriptors(smi):
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        return {}
    return {
        "MW":        Descriptors.MolWt(mol),
        "LogP":      Descriptors.MolLogP(mol),
        "TPSA":      Descriptors.TPSA(mol),        # key: head-group polarity proxy
        "MolMR":     Descriptors.MolMR(mol),       # molar refractivity ~ Vm proxy
        "HeavyAtoms":mol.GetNumHeavyAtoms(),
        "n_F":       sum(1 for a in mol.GetAtoms() if a.GetAtomicNum() == 9),
        "n_O":       sum(1 for a in mol.GetAtoms() if a.GetAtomicNum() == 8),
        "RotBonds":  rdMolDescriptors.CalcNumRotatableBonds(mol),
    }

desc = df["smiles"].apply(get_descriptors).apply(pd.Series)
df = pd.concat([df, desc], axis=1)

# Physics-informed manual descriptors
df["head_type"]  = (df["head_class"] == "PFSA").astype(int)   # 0=PFCA, 1=PFSA
df["is_branched"]= df["branched"].astype(int)
df["F_per_C"]    = df["n_F"] / df["n_CF"]                    # fluorination density of tail

# Approximate Vm = MW / density (density: PFCA~1.65, PFSA~1.8 g/cm³)
df["density_est"] = np.where(df["head_class"]=="PFCA", 1.65, 1.80)
df["Vm_est"] = df["MW"] / df["density_est"]

print("Dataset shape:", df.shape)
print("\nKey descriptors:")
print(df[["name","head_class","n_CF","log_Ki","Vm_est","TPSA","LogP","MolMR"]].to_string())

# ── MODELS ────────────────────────────────────────────────────────────────────
# Drop rows with NaN
df_clean = df.dropna(subset=["Vm_est","TPSA","LogP","n_CF"]).copy()
y = df_clean["log_Ki"].values

# Model 1: Baseline — molar volume only (Brusseau 2019 benchmark)
X1 = df_clean[["Vm_est"]].values
m1 = LinearRegression().fit(X1, y)
loo = LeaveOneOut()
y_pred_m1 = cross_val_predict(m1, X1, y, cv=loo)
r2_m1 = r2_score(y, y_pred_m1)
rmse_m1 = np.sqrt(mean_squared_error(y, y_pred_m1))

# Model 2: Physics-informed linear — Vm + TPSA + head_type + n_CF + branching
X2 = df_clean[["Vm_est","TPSA","head_type","n_CF","is_branched"]].values
scaler2 = StandardScaler()
X2s = scaler2.fit_transform(X2)
m2 = Ridge(alpha=1.0).fit(X2s, y)
y_pred_m2 = cross_val_predict(m2, X2s, y, cv=loo)
r2_m2 = r2_score(y, y_pred_m2)
rmse_m2 = np.sqrt(mean_squared_error(y, y_pred_m2))

# Model 3: Full physics-informed — all descriptors + GBR
feat_cols = ["Vm_est","TPSA","LogP","MolMR","n_CF","head_type","is_branched","F_per_C","n_O","RotBonds"]
X3 = df_clean[feat_cols].values
scaler3 = StandardScaler()
X3s = scaler3.fit_transform(X3)
m3 = GradientBoostingRegressor(n_estimators=100, max_depth=2, learning_rate=0.1, random_state=42)
y_pred_m3 = cross_val_predict(m3, X3s, y, cv=loo)
r2_m3 = r2_score(y, y_pred_m3)
rmse_m3 = np.sqrt(mean_squared_error(y, y_pred_m3))

print(f"\n{'Model':<35} {'LOO R²':>8} {'LOO RMSE':>10}")
print("-"*55)
print(f"{'Baseline (Vm only)':35} {r2_m1:8.4f} {rmse_m1:10.4f}")
print(f"{'Physics-informed linear':35} {r2_m2:8.4f} {rmse_m2:10.4f}")
print(f"{'Full physics-informed (GBR)':35} {r2_m3:8.4f} {rmse_m3:10.4f}")

# ── HEAD-GROUP DISCRIMINATION TEST ────────────────────────────────────────────
# Paired PFCA/PFSA at same chain length
pairs = [(4,"PFBA","PFBS"),(6,"PFHxA","PFHxS"),(7,"PFHpA","PFHpS"),(8,"PFOA","PFOS"),(9,"PFNA","PFNS")]
print("\nHead-Group Discrimination Analysis:")
print(f"{'Pair':20} {'True Ki PFCA':>14} {'True Ki PFSA':>14} {'PFSA>PFCA':>10} | Vm Correct | TPSA Correct")
print("-"*90)

for n_cf, ca_name, sa_name in pairs:
    row_ca = df_clean[df_clean["name"]==ca_name]
    row_sa = df_clean[df_clean["name"]==sa_name]
    if len(row_ca)==0 or len(row_sa)==0:
        continue
    logKi_ca = row_ca["log_Ki"].values[0]
    logKi_sa = row_sa["log_Ki"].values[0]
    true_order = logKi_sa > logKi_ca  # PFSA should be higher
    # Vm comparison
    vm_ca = row_ca["Vm_est"].values[0]; vm_sa = row_sa["Vm_est"].values[0]
    vm_pred_ca = m1.predict([[vm_ca]])[0]; vm_pred_sa = m1.predict([[vm_sa]])[0]
    vm_correct = (vm_pred_sa > vm_pred_ca) == true_order
    # TPSA comparison
    tpsa_ca = row_ca["TPSA"].values[0]; tpsa_sa = row_sa["TPSA"].values[0]
    tpsa_correct = (tpsa_sa > tpsa_ca) == true_order
    print(f"{ca_name+'/'+sa_name:20} {logKi_ca:14.3f} {logKi_sa:14.3f} {str(true_order):>10} | {'✓' if vm_correct else '✗':^10} | {'✓' if tpsa_correct else '✗':^12}")

# ── BRANCHING DISCRIMINATION TEST ────────────────────────────────────────────
print("\nBranching Discrimination Analysis (same n_CF, different branching):")
branch_pairs = [(8,"PFOA","Iso-PFOA"),(10,"PFDA","Iso-PFDA")]
for n_cf, lin_name, iso_name in branch_pairs:
    row_l = df_clean[df_clean["name"]==lin_name]
    row_i = df_clean[df_clean["name"]==iso_name]
    if len(row_l)==0 or len(row_i)==0: continue
    logKi_l = row_l["log_Ki"].values[0]; logKi_i = row_i["log_Ki"].values[0]
    vm_l = row_l["Vm_est"].values[0]; vm_i = row_i["Vm_est"].values[0]
    vm_pred_l = m1.predict([[vm_l]])[0]; vm_pred_i = m1.predict([[vm_i]])[0]
    true_higher = iso_name if logKi_i > logKi_l else lin_name
    vm_higher = iso_name if vm_pred_i > vm_pred_l else lin_name
    print(f"  {lin_name} log Ki={logKi_l:.3f} | {iso_name} log Ki={logKi_i:.3f}")
    print(f"  True higher: {true_higher} | Vm predicts higher: {vm_higher} | Match: {true_higher==vm_higher}")

df_clean.to_csv("pfas_ki_dataset.csv", index=False)
print("\nDataset saved to pfas_ki_dataset.csv")